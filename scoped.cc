
// Daniel Pitzl, DESY, Sep 2017, Apr 2018
// telescope analysis with eudaq and ROC4Sens edge-on

// make scoped
// needs runs.dat
// needs align_31972.dat from tele

// scoped 31972

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TFile.h>
#include <TH1.h> // counting
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <set> // for hotset
#include <cmath>

using namespace std;
using namespace eudaq;

class pixel {
 public:
  int col;
  int row;
  double ph;
  double q;
  int ord;
  bool big;
};

class rowsrt {
public:
  bool operator() ( const pixel p1,  const pixel p2 )
  {
    return p1.row < p2.row;
  }
};

class colsrt {
public:
  bool operator() ( const pixel p1,  const pixel p2 )
  {
    return p1.col < p2.col;
  }
};

struct cluster {
  vector <pixel> vpix; // Armin Burgmeier: list
  int size;
  int ncol, nrow;
  double col, row;
  double charge;
  bool big;
};

struct triplet {
  double xm;
  double ym;
  double zm;
  double sx;
  double sy;
  bool lk;
  double ttdmin;
  vector <double> vx;
  vector <double> vy;
};

//------------------------------------------------------------------------------
vector < cluster > getClus( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with local coordinates
  // next-neighbour topological clustering (allows fCluCut-1 empty pixels)

  vector < cluster > v;
  if( pb.size() == 0 ) return v;

  int * gone = new int[pb.size()]{0};

  unsigned seed = 0;

  while( seed < pb.size() ) {

    // start a new cluster:

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do {
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( (   dr>=-fCluCut) && (dr<=fCluCut) 
		&& (dc>=-fCluCut) && (dc<=fCluCut) ) {
              c.vpix.push_back(pb[i]);
	      gone[i] = 1;
              growing = 1;
              break; // important!
            }
          } // loop over vpix
        } // not gone
      } // loop over all pix
    }
    while( growing );

    // added all I could. determine position and append it to the list of clusters:

    c.size = c.vpix.size();
    c.col = 0;
    c.row = 0;
    double sumQ = 0;
    c.big = 0;
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      double Qpix = p->q; // calibrated [Vcal]
      if( Qpix < 0 ) Qpix = 1; // DP 1.7.2012
      //Qpix = 3; // [ke] for shallow angle
      sumQ += Qpix;
      c.col += (*p).col*Qpix;
      c.row += (*p).row*Qpix;
      if( p->big ) c.big = 1;
      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( sumQ > 0 ) {
      c.col /= sumQ;
      c.row /= sumQ;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with non-positive charge" << endl;
    }

    c.charge = sumQ;
    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;

    v.push_back(c); // add cluster to vector

    // look for a new seed = unused pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  // nothing left, return clusters

  delete gone;
  return v;
}

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc == 1 ) {
    cout << "give run number" << endl;
    return 1;
  }

  // run number = last arg

  string runnum( argv[argc-1] );
  int run = atoi( argv[argc-1] );

  cout << "run " << run << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // further arguments:

  int lev = 900200100; // last event

  double thr = 0; // offline pixel threshold [ADC]

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

    if( !strcmp( argv[i], "-t" ) )
      thr = atof( argv[++i] ); // [ke]

  } // argc

  cout << "apply offline pixel threshold at " << thr << " ke" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  string geoFileName( "geo.dat" );
  double DUTtilt0 = 19.3;
  double pbeam = 4.8;
  int chip0 = 110;
  string gainFileName( "gain.dat" );

  ifstream runsFile( "runs.dat" );

  if( runsFile.bad() || ! runsFile.is_open() ) {
    cout << "Error opening runs.dat" << endl;
    return 1;
  }
  // can there be instructions between if and else ? no

  else {

    cout << "read runs from runs.dat" << endl;

    string hash( "#" );
    string RUN( "run" );
    string GEO( "geo" );
    string GeV( "GeV" );
    string CHIP( "chip" );
    string GAIN( "gain" );
    string TILT( "tilt" );
    bool found = 0;

    while( ! runsFile.eof() ) {

      string line;
      getline( runsFile, line );

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == RUN )  {
	int ival;
	tokenizer >> ival;
	if( ival == run ) {
	  found = 1;
	  break; // end file reading
	}
      }

      if( tag == TILT ) {
	tokenizer >> DUTtilt0;
	continue;
      }

      if( tag == GAIN ) {
	tokenizer >> gainFileName;
	continue;
      }

      if( tag == GEO ) {
	tokenizer >> geoFileName;
	continue;
      }

      if( tag == GeV ) {
	tokenizer >> pbeam;
	continue;
      }

      if( tag == CHIP ) {
	tokenizer >> chip0;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    if( found )
      cout 
	<< "settings for run " << run << ":" << endl
	<< "  beam " << pbeam << " GeV" << endl
	<< "  geo file " << geoFileName << endl
	<< "  nominal DUT tilt " << DUTtilt0 << " deg" << endl
	<< "  DUT chip " << chip0 << endl
	<< "  DUT gain file " << gainFileName << endl
	<< endl;
    else {
      cout << "run " << run << " not found in runs.dat" << endl;
      return 1;
    }

  } // runsFile

  runsFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // geometry:

  int nx[9]; // x-pixels per plane
  int ny[9]; // y-pixels per plane
  double sizex[9]; // x size per plane
  double sizey[9]; // y size per plane
  double ptchx[9]; // x-pixel size
  double ptchy[9]; // y-pixel size
  double midx[9]; // x mid
  double midy[9]; // y mid

  double zz[9];

  for( int ipl = 0; ipl < 9; ++ipl )
    nx[ipl] = 0; // missing plane flag

  ifstream geoFile( geoFileName );

  cout << endl;

  if( geoFile.bad() || ! geoFile.is_open() ) {
    cout << "Error opening " << geoFileName << endl;
    return 1;
  }

  cout << "read geometry from " << geoFileName << endl;

  { // open local scope

    string hash( "#" );
    string plane( "plane" );
    string type( "type" );
    string sizexs( "sizex" );
    string sizeys( "sizey" );
    string npixelx( "npixelx" );
    string npixely( "npixely" );
    string zpos( "zpos" );

    int ipl = 0;
    string chiptype;

    while( ! geoFile.eof() ) {

      string line;
      getline( geoFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane ) {
	tokenizer >> ipl;
	continue;
      }

      if( ipl < 0 || ipl >= 9 ) {
	cout << "wrong plane number " << ipl << endl;
	continue;
      }

      if( tag == type ) {
	tokenizer >> chiptype;
	continue;
      }

      if( tag == sizexs ) {
	double val;
	tokenizer >> val;
	sizex[ipl] = val;
	continue;
      }

      if( tag == sizeys ) {
	double val;
	tokenizer >> val;
	sizey[ipl] = val;
	continue;
      }

      if( tag == npixelx ) {
	int val;
	tokenizer >> val;
	nx[ipl] = val;
	continue;
      }

      if( tag == npixely ) {
	int val;
	tokenizer >> val;
	ny[ipl] = val;
	continue;
      }

      if( tag == zpos ) {
	double val;
	tokenizer >> val;
	zz[ipl] = val;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    for( int ipl = 0; ipl < 9; ++ipl ) {
      if( nx[ipl] == 0 ) continue; // missing plane flag
      ptchx[ipl] = sizex[ipl] / nx[ipl]; // pixel size
      ptchy[ipl] = sizey[ipl] / ny[ipl];
      midx[ipl] = 0.5 * sizex[ipl]; // mid plane
      midy[ipl] = 0.5 * sizey[ipl]; // mid plane
    }

  } // geo scope

  geoFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignments:

  int aligniteration = 0;
  double alignx[9];
  double aligny[9];
  double rotx[9];
  double roty[9];

  ostringstream alignFileName; // output string stream

  alignFileName << "align_" << run << ".dat";

  ifstream ialignFile( alignFileName.str() );

  cout << endl;

  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    cout << "Error opening " << alignFileName.str() << endl
	 << "  please do: tele -g " << geoFileName << " " << run << endl
	 << endl;
    return 1;
  }
  else {

    cout << "read alignment from " << alignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string plane( "plane" );
    string shiftx( "shiftx" );
    string shifty( "shifty" );
    string rotxvsy( "rotxvsy" );
    string rotyvsx( "rotyvsx" );

    int ipl = 0;

    while( ! ialignFile.eof() ) {

      string line;
      getline( ialignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) 
	tokenizer >> aligniteration;

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 0 || ipl >= 9 ) {
	cout << "wrong plane number " << ipl << endl;
	continue;
      }

      double val;
      tokenizer >> val;
      if(      tag == shiftx )
	alignx[ipl] = val;
      else if( tag == shifty )
	aligny[ipl] = val;
      else if( tag == rotxvsy )
	rotx[ipl] = val;
      else if( tag == rotyvsx )
	roty[ipl] = val;

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  ialignFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // hot pixels:

  ostringstream hotFileName; // output string stream

  hotFileName << "hot_" << run << ".dat";

  ifstream ihotFile( hotFileName.str() );

  set <int> hotset[9];

  if( ihotFile.bad() || ! ihotFile.is_open() ) {
    cout << "no " << hotFileName.str() << " (created by tele)" << endl;
  }
  else {

    cout << "read hot pixel list from " << hotFileName.str() << endl;

    string hash( "#" );
    string plane( "plane" );
    string pix( "pix" );

    int ipl = 0;

    while( ! ihotFile.eof() ) {

      string line;
      getline( ihotFile, line );
      //cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 0 || ipl >= 6 ) {
	//cout << "wrong plane number " << ipl << endl;
	continue;
      }

      if( tag == pix ) {
	int ix, iy;
	tokenizer >> ix;
	tokenizer >> iy;
	int ipx = ix*ny[ipl]+iy;
	hotset[ipl].insert(ipx);
      }

    } // while getline

  } // hotFile

  ihotFile.close();

  for( int ipl = 0; ipl < 6; ++ipl )
    cout << ipl << ": hot " << hotset[ipl].size() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT:

  const double log10 = log(10);
  const double wt = atan(1.0) / 45.0; // pi/180 deg

  const double qwid = 1.2; // [ke] for Moyal

  bool rot90 = 0; // straight
  if( chip0 == 106 ) rot90 = 1;
  if( chip0 == 107 ) rot90 = 1;
  if( chip0 == 108 ) rot90 = 1;
  if( chip0 == 109 ) rot90 = 1;
  if( chip0 == 110 ) rot90 = 1;
  if( chip0 == 111 ) rot90 = 1;
  if( chip0 == 112 ) rot90 = 1;
  if( chip0 == 113 ) rot90 = 1;
  if( chip0 == 114 ) rot90 = 1;
  if( chip0 == 115 ) rot90 = 1;
  if( chip0 == 116 ) rot90 = 1;
  if( chip0 == 117 ) rot90 = 1;
  if( chip0 == 118 ) rot90 = 1;

  if( chip0 == 120 ) rot90 = 1; /// irrad
  if( chip0 == 124 ) rot90 = 1; /// irrad
  if( chip0 == 129 ) rot90 = 1; /// irrad
  if( chip0 == 130 ) rot90 = 1; /// irrad
  if( chip0 == 136 ) rot90 = 1; /// irrad
  if( chip0 == 137 ) rot90 = 1; /// irrad

  if( chip0 == 139 ) rot90 = 1;
  if( chip0 == 142 ) rot90 = 1;
  if( chip0 == 143 ) rot90 = 1;
  if( chip0 == 144 ) rot90 = 1;
  if( chip0 == 147 ) rot90 = 1;
  if( chip0 == 155 ) rot90 = 1;

  if( !rot90 ) {
    cout << "only for DUT 90 degree rotated" << endl;
    return 1;
  }

  bool fifty = 0;
  if( chip0 == 102 ) fifty = 1;
  if( chip0 == 106 ) fifty = 1;
  if( chip0 == 111 ) fifty = 1;
  if( chip0 == 117 ) fifty = 1;
  if( chip0 == 118 ) fifty = 1;
  if( chip0 == 122 ) fifty = 1; // irr
  if( chip0 == 123 ) fifty = 1; // irr
  if( chip0 == 133 ) fifty = 1; // irr
  if( chip0 == 139 ) fifty = 1;
  if( chip0 == 142 ) fifty = 1;
  if( chip0 == 143 ) fifty = 1;
  if( chip0 == 144 ) fifty = 1;
  if( chip0 == 147 ) fifty = 1;
  if( chip0 == 155 ) fifty = 1;

  int iDUT = 7;

  int DUTaligniteration = 0;
  double DUTalignx = 0.0;
  double DUTaligny = 0.0;
  double DUTrot = 0.0;
  double DUTturn = 0;
  double DUTtilt = DUTtilt0; // [deg]
  double DUTz = 0.5*( zz[2] + zz[3] );
  double DUTthickness = 0.150; // [mm]
  if( chip0 == 155 )
    DUTthickness = 0.180; // deep diffused

  ostringstream DUTalignFileName; // output string stream

  DUTalignFileName << "alignedg_" << run << ".dat";

  ifstream iDUTalignFile( DUTalignFileName.str() );

  cout << endl;

  if( iDUTalignFile.bad() || ! iDUTalignFile.is_open() ) {
    cout << "no " << DUTalignFileName.str() << ", use runs.dat" << endl;
  }
  else {

    cout << "read DUTalignment from " << DUTalignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string alignx( "alignx" );
    string aligny( "aligny" );
    string rot( "rot" );
    string tilt( "tilt" );
    string turn( "turn" );
    string dz( "dz" );

    while( ! iDUTalignFile.eof() ) {

      string line;
      getline( iDUTalignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) 
	tokenizer >> DUTaligniteration;

      double val;
      tokenizer >> val;
      if(      tag == alignx )
	DUTalignx = val;
      else if( tag == aligny )
	DUTaligny = val;
      else if( tag == rot )
	DUTrot = val;
      else if( tag == tilt )
	DUTtilt = val;
      else if( tag == turn )
	DUTturn = val;
      else if( tag == dz )
	DUTz = val + zz[2];

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  iDUTalignFile.close();

  double co = cos( DUTturn ); // rad
  double so = sin( DUTturn );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT gain:

  double p0[155][160]; // Fermi
  double p1[155][160];
  double p2[155][160];
  double p3[155][160];

  ifstream gainFile( gainFileName );

  if( ! gainFile ) {
    cout << "gain file " << gainFileName << " not found" << endl;
    return 1;
  }
  else {
    cout << endl << "using DUT gain file " << gainFileName << endl;

    while( ! gainFile.eof() ) {

      int icol;
      int irow;
      gainFile >> icol;
      gainFile >> irow;
      gainFile >> p0[icol][irow];
      gainFile >> p1[icol][irow];
      gainFile >> p2[icol][irow];
      gainFile >> p3[icol][irow];

    } // while

  } // gainFile

  double ke = 0.0367; // default

  double dphcut = 12; // 31619 dy8cq 4.77

  if( run >= 31635 && run <= 32266 ) dphcut = 24; // gain_2 2018

  //if( chip0 == 109 ) dphcut = 50; // like edges2

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "scoped" << run << ".root";

  TFile* histoFile = new TFile( rootFileName.str(  ).c_str(  ), "RECREATE" );

  // book histos:

  double f = 4.8/pbeam;

  TH1I hdttlu( "dttlu", "TLU time between events;TLU time between events log_{10}(#Deltat [s]);events",
	       60, -4, 2 );
  TH1I hdtdtb( "dtdtb", "DTB time between events;DTB time between events log_{10}(#Deltat [s]);events",
	       60, -4, 2 );

  TH1I hddt( "ddt", "#Deltadt TLU - DTB;TLU - DTB #Deltadt [ms];events", 200, -1, 1 );
  TProfile ddtvsev1( "ddtvsev1", "#Deltadt TLU - DTB;event;<#Deltadt TLU - DTB> [ms]",
		    100, 0, 50000, -100, 100 );
  TProfile ddtvsev2( "ddtvsev2", "#Deltadt TLU - DTB;event;<#Deltadt TLU - DTB> [ms]",
		    1000, 0, 1000*1000, -100, 100 );

  TH1I t1Histo( "t1", "event time;event time [s];events", 100, 0, 1 );
  TH1I t2Histo( "t2", "event time;event time [s];events", 300, 0, 300 );
  TH1I t3Histo( "t3", "event time;event time [s];events", 150, 0, 1500 );
  TH1I t4Histo( "t4", "event time;event time [s];events", 600, 0, 6000 );
  TH1I t5Histo( "t5", "event time;event time [s];events", 600, 0, 60000 );
  TH1I t6Histo( "t6", "event time;event time [h];events", 1000, 0, 50 );

  TH1I * hcol[9];
  TH1I * hrow[9];
  TH1I * hnpx[9];
  TH2I * hmap[9];

  TH1I * hncl[9];
  TH1I * hsiz[9];
  TH1I * hncol[9];
  TH1I * hnrow[9];

  for( int ipl = 0; ipl < 9; ++ipl ) {

    hcol[ipl] = new TH1I( Form( "col%i", ipl ),
			  Form( "%i col;col;plane %i pixels", ipl, ipl ), 
			  max( 155, nx[ipl]/4 ), 0, nx[ipl] );
    hrow[ipl] = new TH1I( Form( "row%i", ipl ),
		      Form( "%i row;row;plane %i pixels", ipl, ipl ),
		      max( 160, ny[ipl]/2 ), 0, ny[ipl] );
    hmap[ipl] = new TH2I( Form( "map%i", ipl ),
			  Form( "%i map;col;row;plane %i pixels", ipl, ipl ),
			  max( 155, nx[ipl]/4 ), 0, nx[ipl], max( 160, ny[ipl]/2 ), 0, ny[ipl] );

    hnpx[ipl] = new TH1I( Form( "npx%i", ipl ),
		      Form( "%i pixel per event;pixels;plane %i events", ipl, ipl ),
		      200, 0, 200 );

    hncl[ipl] = new TH1I( Form( "ncl%i", ipl ),
		      Form( "plane %i cluster per event;cluster;plane %i events", ipl, ipl ),
		      51, -0.5, 50.5 );
    hsiz[ipl] = new TH1I( Form( "clsz%i", ipl ),
		      Form( "%i cluster size;pixels/cluster;plane %i clusters", ipl, ipl ),
		      51, -0.5, 50.5 );
    hncol[ipl] = new TH1I( Form( "ncol%i", ipl ), 
		       Form( "%i cluster size x;columns/cluster;plane %i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );
    hnrow[ipl] = new TH1I( Form( "nrow%i", ipl ),
		       Form( "%i cluster size y;rows/cluster;plane %i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );

  } // planes

  TProfile phvsprev( "phvsprev", "Tsunami;previous PH [ADC];<PH> [ADC]", 80, 0, 800, -999, 1999 );
  TH1I dutphHisto( "dutph", "DUT PH;ADC-PED [ADC];pixels", 500, -100, 900 );
  TH1I dutdphHisto( "dutdph", "DUT #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );
  TH1I dutdph0Histo( "dutdph0", "DUT #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );

  TProfile dutvvsdph( "dutvvsdph", "gain;#DeltaPH [ADC];Fermi Vcal [mv]", 500, -100, 900, -999, 1999 );
  TH1I dutgHisto( "dutg", "DUT gain;1/gain [mv/ADC];pixels", 200, 0, 2 );
  TProfile linvvsdph( "linvvsdph", "linear gain;#DeltaPH [ADC];linear Vcal [mv]", 500, -100, 900, -999, 1999 );

  TProfile dutnpxvst2( "dutnpxvst2",
	      "DUT pixels vs time;time [s];DUT pixels per event",
	      150, 0, 1500, -0.5, 99.5 );
  TProfile dutyldvst2( "dutyldvst2",
	      "DUT yield vs time;time [s];DUT events with pixels",
	      150, 0, 1500, -0.5, 1.5 );
  TProfile dutyldvst6( "dutyldvst6",
	      "DUT yield vs time;time [h];DUT events with pixels",
	      1000, 0, 50, -0.5, 1.5 );

  // DUT:

  TH1I dutadcHisto( "dutadc",
		    "DUT pixel ADC;pixel pulse height [ADC];pixels",
		    400, 0, 800 );
  TH1I dutcolHisto( "dutcol",
		    "DUT pixel column;pixel column;pixels",
		    nx[iDUT], -0.5, nx[iDUT]-0.5 );
  TH1I dutrowHisto( "dutrow",
		    "DUT pixel row;pixel row;pixels",
		    ny[iDUT], -0.5, ny[iDUT]-0.5 );

  TH1I dutq0Histo( "dutq0",
		   "normal fiducial cluster charge;normal fiducial cluster charge [ke];fiducial clusters",
		   160, 0, 80 );

  TH1I dutnpxHisto( "dutnpx",
		     "DUT cluster size;cluster size [pixels];clusters",
		     200, 0.5, 200.5 );
  TH1I dutncolHisto( "dutncol",
		     "DUT cluster size;cluster size [columns];clusters",
		     80, 0.5, 80.5 );
  TH1I dutnrowHisto( "dutnrow",
		     "DUT cluster size;cluster size [rows];clusters",
		     80, 0.5, 80.5 );
  TH1I dutcolminHisto( "dutcolmin",
		       "DUT first cluster column;first cluster column;clusters",
		       155, -0.5, 154.5 );
  TH1I dutcolmaxHisto( "dutcolmax",
		       "DUT last cluster column;last cluster column;clusters",
		       155, -0.5, 154.5 );
  TH1I dutcol0qHisto( "dutcol0q",
		      "DUT first column charge;first column charge [ke];clusters",
		      100, 0, 50 );
  TH1I dutcol0oddqHisto( "dutcol0oddq",
			 "DUT odd first column charge;odd first column charge [ke];clusters",
			 100, 0, 50 );
  TH1I dutcol0eveqHisto( "dutcol0eveq",
			 "DUT eve first column charge;eve first column charge [ke];clusters",
			 100, 0, 50 );
  TH1I dutcol9qHisto( "dutcol9q",
		      "DUT last column charge;last column charge [ke];clusters",
		      100, 0, 50 );

  // triplets:

  TH1I hdx02( "dx02", "0-2 dx;0-2 dx [mm];cluster pairs", 100, -f, f );
  TH1I hdy02( "dy02", "0-2 dy;0-2 dy [mm];cluster pairs", 100, -f, f );

  TH1I htridx( "tridx", "triplet dx;triplet dx [mm];triplets", 100, -0.1, 0.1 );
  TH1I htridy( "tridy", "triplet dy;triplet dy [mm];triplets", 100, -0.1, 0.1 );

  TH1I htridxc( "tridxc", "triplet dx;triplet dx [mm];triplets", 100, -0.05, 0.05 );
  TH1I htridyc( "tridyc", "triplet dy;triplet dy [mm];triplets", 100, -0.05, 0.05 );

  TH1I htridxc1( "tridxc1", "triplet dx 1-col;1-col triplet dx [mm];1-col triplets",
		 100, -0.05, 0.05 );
  TH1I htridxc2( "tridxc2", "triplet dx 2-col;2-col triplet dx [mm];2-col triplets",
			100, -0.05, 0.05 );
  TH1I htridxc3( "tridxc3", "triplet dx 3-col;3-col triplet dx [mm];3-col triplets",
			100, -0.05, 0.05 );
  TH1I htridxc4( "tridxc4", "triplet dx 4-col;4-col triplet dx [mm];4-col triplets",
			100, -0.05, 0.05 );
  TH1I htridxc5( "tridxc5", "triplet dx 5-col;5-col triplet dx [mm];5-col triplets",
			100, -0.05, 0.05 );

  TH1I htridxs1( "tridxs1", "triplet dx 1-px;1-px triplet dx [mm];1-px triplets",
			100, -0.05, 0.05 );
  TH1I htridxs2( "tridxs2", "triplet dx 2-px;2-px triplet dx [mm];2-px triplets",
			100, -0.05, 0.05 );
  TH1I htridxs3( "tridxs3", "triplet dx 3-px;3-px triplet dx [mm];3-px triplets",
			100, -0.05, 0.05 );
  TH1I htridxs4( "tridxs4", "triplet dx 4-px;4-px triplet dx [mm];4-px triplets",
			100, -0.05, 0.05 );
  TH1I htridxs5( "tridxs5", "triplet dx 5-px;5-px triplet dx [mm];5-px triplets",
			100, -0.05, 0.05 );
  TProfile tridxvsy( "tridxvsy",
		     "triplet dx vs y;triplet yB [mm];<triplet #Deltax> [mm]",
		     110, -5.5, 5.5, -0.05, 0.05 );
  TProfile tridxvstx( "tridxvstx",
		      "triplet dx vs slope x;triplet slope x [rad];<triplet #Deltax> [mm]",
		      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile tridxvst3( "tridxvst3",
		      "triplet dx vs time;time [s];<triplet #Deltax> [mm]",
		      300, 0, 6000, -0.05, 0.05 );
  TProfile tridxvst6( "tridxvst6",
		      "triplet dx vs time;time [h];<triplet #Deltax> [mm]",
		      1000, 0, 50, -0.05, 0.05 );

  TProfile tridyvsx( "tridyvsx",
		     "triplet dy vs x;triplet xB [mm];<triplet #Deltay> [mm]",
		     110, -11, 11, -0.05, 0.05 );
  TProfile tridyvsty( "tridyvsty",
		      "triplet dy vs slope y;triplet slope y [rad];<triplet #Deltay> [mm]",
		      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile tridyvst3( "tridyvst3",
		      "triplet dy vs time;time [s];<triplet #Deltay> [mm]",
		      300, 0, 6000, -0.05, 0.05 );
  TProfile tridyvst6( "tridyvst6",
		      "triplet dy vs time;time [h];<triplet #Deltay> [mm]",
		      1000, 0, 50, -0.05, 0.05 );

  TH1I trixHisto( "trix", "triplets x;x [mm];triplets",
			  240, -12, 12 );
  TH1I triyHisto( "triy", "triplets y;y [mm];triplets",
			  120, -6, 6 );
  TH2I * trixyHisto = new
    TH2I( "trixy", "triplets x-y;x [mm];y [mm];triplets",
	  240, -12, 12, 120, -6, 6 );
  TH1I tritxHisto( "tritx", "triplet slope x;slope x [rad];triplets",
			    100, -0.005, 0.005 );
  TH1I trityHisto( "trity", "triplet slope y;slope y [rad];triplets",
			    100, -0.005, 0.005 );

  TH1I ntriHisto( "ntri", "triplets;triplets;events", 51, -0.5, 50.5 );

  TH1I ttdxHisto( "ttdx", "telescope triplets;triplet #Deltax [mm];triplet pairs",
			 100, -5, 5 );
  TH1I ttdx1Histo( "ttdx1", "telescope triplets;triplet #Deltax [mm];triplet pairs",
			  100, -0.5, 0.5 );
  TH1I ttdmin1Histo( "ttdmin1",
			    "telescope triplets isolation;triplet min #Delta_{xy} [mm];triplet pairs",
			    100, 0, 1 );
  TH1I ttdmin2Histo( "ttdmin2",
			    "telescope triplets isolation;triplet min #Delta_{xy} [mm];triplet pairs",
			    150, 0, 15 );

  // DUT pixel vs triplets:

  TH1I triycHisto( "triyc", "triplets y at DUT;y [mm];triplets at DUT",
		   120, -6, 6 );
  TH1I triyc0Histo( "triyc0", "triplets y at DUT 0 pix;y [mm];triplets at DUT 0 pix",
		    120, -6, 6 );
  TH1I triyc1Histo( "triyc1", "triplets y at DUT 1 pix;y [mm];triplets at DUT 1 pix",
		    120, -6, 6 );
  TH1I triyc2Histo( "triyc2", "triplets y at DUT 2 pix;y [mm];triplets at DUT 2 pix",
		    120, -6, 6 );
  TH1I triyc3Histo( "triyc3", "triplets y at DUT 3 pix;y [mm];triplets at DUT 3 pix",
		    120, -6, 6 );
  TH1I triyc4Histo( "triyc4", "triplets y at DUT 4 pix;y [mm];triplets at DUT 4 pix",
		    120, -6, 6 );
  TH1I triyc6Histo( "triyc6", "triplets y at DUT 6 pix;y [mm];triplets at DUT 6 pix",
		    120, -6, 6 );
  TH1I triyc8Histo( "triyc8", "triplets y at DUT 8 pix;y [mm];triplets at DUT 8 pix",
		    120, -6, 6 );
  TH1I triycnHisto( "triycn", "triplets y at DUT n pix;y [mm];triplets at DUT n pix",
		    120, -6, 6 );

  TH1I cmssxaHisto( "cmssxa",
		    "DUT + Telescope x;pixel + triplet #Sigmax [mm];pixels",
		    440, -11, 11 );
  TH1I cmsdxaHisto( "cmsdxa",
		    "DUT - Telescope x;pixel - triplet #Deltax [mm];pixels",
		    440, -11, 11 );

  TH1I cmssyaHisto( "cmssya",
		    "DUT + Telescope y;pixel + triplet #Sigmay [mm];pixels",
		    198, -99, 99 ); // shallow needs wide range

  TH1I cmsdxHisto( "cmsdx",
		   "DUT - Telescope x;pixel - triplet #Deltax [mm];pixels",
		   200, -1, 1 );

  TProfile cmsdxvsx( "cmsdxvsx",
		     "#Deltax vs x;x track [mm];<pixel - triplet #Deltax> [mm]",
		     50, -3.75, 3.75, -0.5, 0.5 );
  TProfile cmsdxvsz( "cmsdxvsz",
		     "#Deltax vs z;z [mm];<pixel - triplet #Deltax> [mm]",
		     80, -4, 4, -0.5, 0.5 );
  TProfile cmsdxvstx( "cmsdxvstx",
		      "#Deltax vs #theta_{x};x track slope [mrad];<pixel - triplet #Deltax> [#mum]",
		      40, -2, 2, -200, 200 );

  TH1I triyclkHisto( "triyclk", "linked triplets y at DUT;y [mm];linked triplets at DUT",
		     120, -6, 6 );

  TH1I cmsdxcHisto( "cmsdxc",
		    "DUT - Telescope x;pixel - triplet #Deltax [mm];pixels",
		    200, -0.5, 0.5 );

  TH1I cmsdyHisto( "cmsdy",
		    "DUT - Telescope y;pixel - triplet #Deltay [mm];pixels",
		    200, -2, 2 );
  TH1I cmsdycHisto( "cmsdyc",
		    "DUT - Telescope y;pixel - triplet #Deltay [mm];pixels",
		    200, -0.5, 0.5 );
  TProfile cmsdyvsx( "cmsdyvsx",
		     "#Deltay vs x;x track [mm];<pixel - triplet #Deltay> [mm]",
		     50, -3.75, 3.75, -0.2, 0.2 );
  TProfile cmsdyvsz( "cmsdyvsz",
		     "#Deltay vs z;z [mm];<pixel - triplet #Deltay> [mm]",
		     80, -4, 4, -0.2, 0.2 );

  TH1I cmspxpHisto( "cmspxp",
		    "DUT pixel PH linked;pixel PH [ADC];linked pixels",
		    500, 0, 1000 );
  TH1I cmspxqHisto( "cmspxq",
		    "DUT pixel charge linked;pixel charge [ke];linked pixels",
		    250, 0, 25 );

  TH1I nlnkHisto( "nlnk",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );
  TProfile2D nlnkvsxy( "nlnkvsxy",
		       "link map;dx track [mm];dy [mm];<links>",
		       100, -5, 5, 100, -0.5, 0.5, -1, 999 );
  TProfile2D nlnkvsdy( "nlnkvsdy",
		       "link map;dy0 track [mm];dy9 [mm];<links>",
		       100, -1, 1, 100, -1, 1, -1, 999 );
  TProfile nlnkvsdy0( "nlnkvsdy0",
		      "links vs dy0;dy0 [mm];<links>",
		      100, -1, 1, -1, 999 );
  TProfile nlnkvsdy9( "nlnkvsdy9",
		      "links vs dy9;dy9 [mm];<links>",
		      100, -1, 1, -1, 999 );
  TH1I nlnk0Histo( "nlnk0",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );
  TH1I nlnk9Histo( "nlnk9",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );
  TH1I nlnk1Histo( "nlnk1",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );
  TH1I nlnk2Histo( "nlnk2",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );
  TH1I nlnk3Histo( "nlnk3",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );
  TH1I nlnk4Histo( "nlnk4",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );

  TH1I zminHisto( "zmin",
		  "first link;zmin [mm];tracks",
		  80, -4, 4 );
  TH1I zmaxHisto( "zmax",
		  "first link;zmax [mm];tracks",
		  80, -4, 4 );
  TH1I zlngHisto( "zlng",
		  "length;z length [mm];tracks",
		  80, 0, 8 );
  TH1I zlngcHisto( "zlngc",
		  "edge-on length;z length [mm];tracks",
		  80, 0, 8 );
  TProfile zlngvsy( "zlngvsy",
		    "edge-on track length;y [mm];length [mm]",
		    40, -0.1, 0.1, -8, 8 );

  TProfile ymaxvsz( "ymaxvsz",
		    "entry point;z [mm];y [mm]",
		    80, -4, 4, -5, 5 );
  TProfile yminvsz( "yminvsz",
		    "exit point;z [mm];y [mm]",
		    80, -4, 4, -5, 5 );

  TProfile dmaxvsz( "dmaxvsz",
		    "entry point;z [mm];y [mm]",
		    80, -4, 4, -0.5, 0.5 );
  TProfile dminvsz( "dminvsz",
		    "exit point;z [mm];y [mm]",
		    80, -4, 4, -0.5, 0.5 );

  TProfile zlngvsz( "zlngvsz",
		    "track length;z [mm];length [mm]",
		    80, -4, 4, -8, 8 );

  TH1I cmsp0Histo( "cmsp0",
		   "trackcharge;track charge [ADC/mm];tracks",
		   100, 0, 8000 );
  TH1I cmsq0Histo( "cmsq0",
		   "track charge;track charge [ke/mm];tracks",
		   100, 0, 400 );
  TProfile q0vszlng( "q0vszlng",
		     "track charge;track length [mm];<charge/length> [ke/mm]",
		     80, 0, 8, 0, 800 );

  TH1I colpHisto( "colp",
		    "DUT linked column PH;coumn PH [ADC];linked columns",
		    500, 0, 1000 );
  TH1I colqHisto( "colq",
		    "DUT linked column charge;column charge [ke];linked columns",
		    400, 0, 40 );
  TProfile colqvsy( "colqvsy",
		    "column charge vs depth;y [mm];<column charge> [ke]",
		    40, -0.1, 0.1, -1, 80 );
  TH1I colqyHisto( "colqy",
		    "DUT linked column charge;y [mm];sum column charge [ke]",
		    40, -0.1, 0.1 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  cout << endl;

  FileReader * reader;
  if(      run <    100 )
    reader = new FileReader( runnum.c_str(), "data/run0000$2R$X" );
  else if( run <   1000 )
    reader = new FileReader( runnum.c_str(), "data/run000$3R$X" );
  else if( run <  10000 )
    reader = new FileReader( runnum.c_str(), "data/run00$4R$X" );
  else if( run < 100000 )
    reader = new FileReader( runnum.c_str(), "data/run0$5R$X" );
  else
    reader = new FileReader( runnum.c_str(), "data/run$6R$X" );

  // DUT R4S:

  string evFileName = Form( "roi%06i.txt", run );
  cout << "try to open  " << evFileName;
  ifstream evFile( evFileName.c_str() );
  if( !evFile ) {
    cout << " : failed " << endl;
    return 2;
  }
  cout << " : succeed " << endl;

  string START {"START"};
  string hd;
  while( hd != START ) {
    getline( evFile, hd ); // read one line into string
    cout << "  " << hd << endl;
  }
  bool readnext = 1;
  string DUTev;

  string F {"F"}; // filled flag
  string E {"E"}; // empty  flag
  string A {"A"}; // added  flag

  int iev = 0;
  int nresync = 0;

  uint64_t tlutime0 = 0;
  const double fTLU = 384E6; // 384 MHz TLU clock
  uint64_t prevtlutime = 0;

  const double fDTB = 39.917E6; // 40 MHz DTB clock
  uint64_t prevdtbtime = 0;

  bool ldbt = 0;

  do {
    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() ) {
      eudaq::PluginManager::Initialize(evt);
      reader->NextEvent();
      evt = reader->GetDetectorEvent();
    }

    if( evt.IsEORE() ) {
      cout << "EORE" << endl;
      break;
    }

    bool ldb = 0;

    if( ldb ) std::cout << "debug ev " << iev << endl << flush;

    if( iev <  0 )
      ldb = 1;

    if( lev < 10 )
      ldb = 1;

    uint64_t tlutime = evt.GetTimestamp(); // 384 MHz = 2.6 ns
    if( iev < 2  )
      tlutime0 = tlutime;
    double evsec = (tlutime - tlutime0) / fTLU;
    t1Histo.Fill( evsec );
    t2Histo.Fill( evsec );
    t3Histo.Fill( evsec );
    t4Histo.Fill( evsec );
    t5Histo.Fill( evsec );
    t6Histo.Fill( evsec/3600 );

    double tludt = (tlutime - prevtlutime) / fTLU;
    if( tludt > 1e-6 )
      hdttlu.Fill( log(tludt)/log10 );
    prevtlutime = tlutime;

    if( iev < 10 )
      cout << "scoped processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "scoped processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "scoped processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 )
      cout << "scoped processing  " << run << "." << iev << "  taken " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    vector < cluster > cl[9];

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldb ) std::cout << "PLANE " << plane.ID();

      // /home/pitzl/eudaq/main/include/eudaq/CMSPixelHelper.hh

      int ipl = plane.ID();

      if( run > 28000 && ipl > 0 && ipl < 7 ) // 2017, eudaq 1.6: Mimosa 1..6, DUT 7, REF 8, QAD 9
	ipl -= 1; // 0..5
      if( ipl > 8 ) ipl = 6; // QUAD

      if( ldb ) cout << " = ipl " << ipl << ", size " << pxl.size() << flush;

      vector <pixel> pb; // for clustering

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldb )
	  std::cout << ", " << plane.GetX(ipix)
		    << " " << plane.GetY(ipix)
		    << " " << plane.GetPixel(ipix) << flush;

	int ix = plane.GetX(ipix); // global column 0..415
	int iy = plane.GetY(ipix); // global row 0..159
	int ph = plane.GetPixel(ipix); // ADC 0..255

	// skip hot pixels:

	int ipx = ix*ny[ipl] + iy;
	if( hotset[ipl].count(ipx) ) {
	  if( ldb ) cout << " hot" << flush;
	  continue;
	}

	double q = ph;

	hcol[ipl]->Fill( ix+0.5 );
	hrow[ipl]->Fill( iy+0.5 );
	hmap[ipl]->Fill( ix+0.5, iy+0.5 );

	// fill pixel block

	pixel px;
	px.col = ix; // col
	px.row = iy; // row
	px.ph = ph;
	px.q = q;
	px.ord = pb.size(); // readout order
	px.big = 0;
	pb.push_back(px);

	if( pb.size() > 990 ) {
	  cout << "pixel buffer overflow in plane " << ipl
	       << ", event " << iev
	       << endl;
	  break;
	}

      } // pix

      hnpx[ipl]->Fill( pb.size() );

      if( ldb ) std::cout << std::endl;

      // clustering:

      cl[ipl] = getClus(pb);

      if( ldb ) cout << "clusters " << cl[ipl].size() << endl;

      hncl[ipl]->Fill( cl[ipl].size() );

      for( vector<cluster>::iterator c = cl[ipl].begin(); c != cl[ipl].end(); ++c ) {

	hsiz[ipl]->Fill( c->size );
	hncol[ipl]->Fill( c->ncol );
	hnrow[ipl]->Fill( c->nrow );

      } // cl

    } // eudaq planes

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // read DUT stream:

    if( readnext )
      getline( evFile, DUTev );

    if( evFile.eof() ) {
      cout << evFileName << " EOF" << endl;
      break; // event loop
    }

    if( ldb ) cout << "  DUT ev " << DUTev << endl << flush;

    istringstream iss( DUTev ); // tokenize string

    int dut_ev;
    iss >> dut_ev; // DUT event

    string filled;
    iss >> filled;

    int iblk; // event block number: 100, 200, 300...
    iss >> iblk;

    unsigned long dtbtime = prevdtbtime;

    if( filled == F )
      iss >> dtbtime; // from run 456 = 31093

    else if( filled == E )
      iss >> dtbtime; // from run 456 = 31093

    double dtbdt = ( dtbtime - prevdtbtime ) / fDTB;
    if( dtbdt > 1e-6 )
      hdtdtb.Fill( log(dtbdt)/log10 );

    hddt.Fill( (tludt - dtbdt)*1E3 ); // [ms]
    ddtvsev1.Fill( iev, (tludt - dtbdt)*1E3 ); // [ms]
    ddtvsev2.Fill( iev, (tludt - dtbdt)*1E3 ); // [ms]

    if( iev > 65100 )
      ldbt = 0;

    if( ldbt && tludt > 0.5 ) cout << endl;
    if( ldbt )
      cout << "\t" << iev << " TLU " << tludt*1E3
	   << ", DTB " << dtbdt*1e3
	   << endl; // [ms]

    // large time gap = DTB readout of one data block

    while( iev > 88 && tludt > 0.5 && dtbdt < 0.2*tludt ) {

      prevdtbtime = dtbtime;
      ++nresync;
      //if( ldbt )
	cout << "  resync " << nresync << endl;
      if( filled == F )
	getline( evFile, DUTev ); // hits
      getline( evFile, DUTev ); // next event
      istringstream iss( DUTev ); // tokenize string
      iss >> dut_ev;
      iss >> filled;
      iss >> iblk;
      if( filled == F || filled == E ) {
	iss >> dtbtime; // from run 456 = 31093
	dtbdt = ( dtbtime - prevdtbtime ) / fDTB;
	prevdtbtime = dtbtime;
	if( ldbt )
	  cout << "\t" << iev << " TLU " << tludt*1E3
	       << ", DTB " << dtbdt*1e3
	       << endl; // [ms]
      }
      else { // added event
	if( ldbt )
	  cout << "\t DTB added" << endl;
	if( chip0 == 109 && run == 31095 && iev == 712020 )
	  continue;
	else if( chip0 == 109 && run == 31095 && iev == 722668 )
	  continue;
	else
	  break;
      }

    } // tludt

    vector <pixel> pb; // for clustering

    if( readnext && filled == F ) {

      string roi;
      getline( evFile, roi );
      istringstream roiss( roi ); // tokenize string

      int ipx = 0;
      vector <pixel> vpx;
      vpx.reserve(35);

      if( ldb ) cout << "  px";

      while( ! roiss.eof() ) {

	int col;
	int row;
	double ph;
	roiss >> col;
	roiss >> row;
	roiss >> ph;
	if( ldb ) cout << " " << col << " " << row << " " << ph;
	hcol[iDUT]->Fill( col+0.5 );
	hrow[iDUT]->Fill( row+0.5 );
	dutphHisto.Fill( ph );

	pixel px { col, row, ph, ph, ipx, 0 };
	vpx.push_back(px);
	++ipx;

      } // roi px

      // columns-wise common mode correction:

      set <pixel,rowsrt> compx[155]; // per column, sorted along row

      for( unsigned ipx = 0; ipx < vpx.size(); ++ipx ) {
	int col = vpx[ipx].col;
	pixel px { col, vpx[ipx].row, vpx[ipx].ph, vpx[ipx].q, vpx[ipx].ord, 0 };
	compx[col].insert(px); // sorted along row
      }

      for( unsigned col = 0; col < 155; ++col ) {

	if( compx[col].size() < 2 ) continue;

	auto px1 = compx[col].begin();
	auto px7 = compx[col].end(); --px7; // last
				       
	int row1 = px1->row;
	int row7 = px7->row;
	double ph1 = px1->ph;
	double ph7 = px7->ph;

	auto px4 = px1;
	double phprev = px4->ph;
	++px4;
	for( ; px4 != px7; ++px4 ) { // between 1 and 7, exclusively

	  int col4 = px4->col;
	  int row4 = px4->row;
	  double ph4 = px4->ph; // pedestal already subtracted online

	  phvsprev.Fill( phprev, ph4 );

	  if( run >= 31635 && run <= 32266 ) { // 2018 gain_2
	    if( chip0 == 109 )
	      ph4 -= 0.11*phprev; // Tsunami
	    if( chip0 == 118 )
	      ph4 -= 0.21*phprev; // Tsunami
	  }

	  double dph;
	  if( row4 - row1 < row7 - row4 )
	    dph = ph4 - ph1;
	  else
	    dph = ph4 - ph7;
	  phprev = px4->ph;

	  dutdphHisto.Fill( dph );
	  if( row4 > 10 && row4 < 150 ) dutdph0Histo.Fill( dph ); // skip noisy regions

	  // r4scal.C

	  double U = ( dph - p3[col4][row4] ) / p2[col4][row4];

	  if( U >= 1 )
	    U = 0.9999999; // avoid overflow

	  double vcal = p0[col4][row4] - p1[col4][row4] * log( (1-U)/U ); // inverse Fermi

	  dutvvsdph.Fill( dph, vcal );

	  // linear gain from Vcal:
	  double v = 300; // fresh
	  if( chip0 >= 119 && chip0 <= 138 )
	    v = 200; // irrad
	  if( fifty )
	    v = 100;
	  double t = ( v - p0[col4][row4] ) / p1[col4][row4]; // from Vcal
	  double a = p3[col4][row4] + p2[col4][row4] / ( 1 + exp(-t) ); // [ADC]
	  double g = v / a; // mv/ADC
	  dutgHisto.Fill( g );
	  linvvsdph.Fill( dph, g*dph );

	  vcal = g*dph; // overwrite !!

	  pixel px;

	  if( fifty ) {
	    px.col = col4; // 50x50
	    px.row = row4;
	  }
	  else {
	    px.col = (col4+1)/2; // 100 um
	    if( col4%2 ) 
	      px.row = 2*row4 + 0;
	    else
	      px.row = 2*row4 + 1;
	  }
	  px.ph = dph;
	  px.q = ke*vcal; // calibrated
	  px.ord = pb.size(); // readout order
	  px.big = 0;

	  if( dph > dphcut ) {

	    hmap[iDUT]->Fill( px.col+0.5, px.row+0.5 );
	    pb.push_back(px);

	  } // dph cut

	} // p4

      } // cols

      if( ldb ) cout << " npx " << pb.size() << endl << flush;

    } // filled

    readnext = 1;
    if( iev > 88 && dtbdt > 0.5 && tludt < 0.2*dtbdt ) {
      readnext = 0;
      if( ldbt ) cout << "repeat DTB event " << DUTev << endl;
    }
    if( readnext )
      prevdtbtime = dtbtime;

    hnpx[iDUT]->Fill( pb.size() );

    dutnpxvst2.Fill( evsec, pb.size() );

    int DUTyld = 0;
    if( pb.size() ) DUTyld = 1; // no double counting: events with at least one px
    dutyldvst2.Fill( evsec, DUTyld );
    dutyldvst6.Fill( evsec/3600, DUTyld );

    if( ldb ) cout << "  DUT px " << pb.size() << endl << flush;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT:

    int colmin = 999;
    int colmax = -1;
    int rowmin = 999;
    int rowmax = -1;

    double qsum = 0;
    double sumcol = 0;
    double sumrow = 0;

    double qcol[nx[iDUT]];
    for( int icol = 0; icol < nx[iDUT]; ++icol ) qcol[icol] = 0;

    double qrow[ny[iDUT]];
    for( int irow = 0; irow < ny[iDUT]; ++irow ) qrow[irow] = 0;

    for( vector<pixel>::iterator px = pb.begin(); px != pb.end(); ++px ) {

      int icol = px->col;
      int irow = px->row;

      dutadcHisto.Fill( px->ph );
      dutcolHisto.Fill( icol );
      dutrowHisto.Fill( irow );

      if( icol < colmin ) colmin = icol;
      if( icol > colmax ) colmax = icol;
      if( irow < rowmin ) rowmin = irow;
      if( irow > rowmax ) rowmax = irow;

      double q = px->q;
      if( q < 0 ) continue;

      qsum += q;
      qcol[icol] += q;
      qrow[irow] += q;

      sumcol += icol*q;
      sumrow += irow*q;

    } // pix

    int ncol = colmax - colmin + 1;
    int nrow = rowmax - rowmin + 1;

    if( colmin > 0 && colmax < nx[iDUT]-2 && rowmin > 0 && rowmax < ny[iDUT]-1 ) {

      dutq0Histo.Fill( qsum / ncol );
      dutncolHisto.Fill( ncol );
      dutnrowHisto.Fill( nrow );

      dutcolminHisto.Fill( colmin ); // strong even peaks at large turn 14848
      dutcolmaxHisto.Fill( colmax ); // weaker even peaks

      dutcol0qHisto.Fill(qcol[colmin]);
      if( colmin%2 )
	dutcol0oddqHisto.Fill(qcol[colmin]);
      else
	dutcol0eveqHisto.Fill(qcol[colmin]);
      dutcol9qHisto.Fill(qcol[colmax]);

    } // DUT clusters

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make triplets 2+0-1:

    vector <triplet> triplets;

    double triCut = 0.1; // [mm]

    for( vector<cluster>::iterator cA = cl[0].begin(); cA != cl[0].end(); ++cA ) {

      double xA = cA->col*ptchx[0] - alignx[0];
      double yA = cA->row*ptchy[0] - aligny[0];
      double xmid = xA - midx[0];
      double ymid = yA - midy[0];
      xA = xmid - ymid*rotx[0];
      yA = ymid + xmid*roty[0];

      for( vector<cluster>::iterator cC = cl[2].begin(); cC != cl[2].end(); ++cC ) {

	double xC = cC->col*ptchx[2] - alignx[2];
	double yC = cC->row*ptchy[2] - aligny[2];
	double xmid = xC - midx[2];
	double ymid = yC - midy[2];
	xC = xmid - ymid*rotx[2];
	yC = ymid + xmid*roty[2];

	double dx2 = xC - xA;
	double dy2 = yC - yA;
	double dz02 = zz[2] - zz[0]; // from 0 to 2 in z
	hdx02.Fill( dx2 );
	hdy02.Fill( dy2 );

	if( fabs( dx2 ) > 0.005 * dz02 ) continue; // angle cut ?
	if( fabs( dy2 ) > 0.005 * dz02 ) continue; // angle cut

	double avx = 0.5 * ( xA + xC ); // mid
	double avy = 0.5 * ( yA + yC );
	double avz = 0.5 * ( zz[0] + zz[2] ); // mid z
 
	double slpx = ( xC - xA ) / dz02; // slope x
	double slpy = ( yC - yA ) / dz02; // slope y

	// middle plane B = 1:

	for( vector<cluster>::iterator cB = cl[1].begin(); cB != cl[1].end(); ++cB ) {

	  double xB = cB->col*ptchx[1] - alignx[1];
	  double yB = cB->row*ptchy[1] - aligny[1];
	  double xmid = xB - midx[1];
	  double ymid = yB - midy[1];
	  xB = xmid - ymid*rotx[1];
	  yB = ymid + xmid*roty[1];

	  // interpolate track to B:

	  double dz = zz[1] - avz;
	  double xk = avx + slpx * dz; // triplet at k
	  double yk = avy + slpy * dz;

	  double dx3 = xB - xk;
	  double dy3 = yB - yk;
	  htridx.Fill( dx3 );
	  htridy.Fill( dy3 );

	  if( fabs( dy3 ) < 0.05 ) {

	    htridxc.Fill( dx3 );
	    tridxvsy.Fill( yB, dx3 );
	    tridxvstx.Fill( slpx, dx3 );
	    tridxvst3.Fill( evsec, dx3 );
	    tridxvst6.Fill( evsec/3600, dx3 );

	    if(      cB->size == 1 )
	      htridxs1.Fill( dx3 ); // 4.2 um
	    else if( cB->size == 2 )
	      htridxs2.Fill( dx3 ); // 4.0 um
	    else if( cB->size == 3 )
	      htridxs3.Fill( dx3 ); // 3.8 um
	    else if( cB->size == 4 )
	      htridxs4.Fill( dx3 ); // 4.3 um
	    else
	      htridxs5.Fill( dx3 ); // 3.6 um

	    if(      cB->ncol == 1 )
	      htridxc1.Fill( dx3 ); // 4.0 um
	    else if( cB->ncol == 2 )
	      htridxc2.Fill( dx3 ); // 4.1 um
	    else if( cB->ncol == 3 )
	      htridxc3.Fill( dx3 ); // 3.6 um
	    else if( cB->ncol == 4 )
	      htridxc4.Fill( dx3 ); // 3.5 um
	    else
	      htridxc5.Fill( dx3 ); // 4.1 um

	  } // dy

	  if( fabs( dx3 ) < 0.05 ) {
	    htridyc.Fill( dy3 );
	    tridyvsx.Fill( xB, dy3 );
	    tridyvsty.Fill( slpy, dy3 );
	    tridyvst3.Fill( evsec, dy3 );
	    tridyvst6.Fill( evsec/3600, dy3 );
	  }

	  // telescope triplet cuts:

	  if( fabs(dx3) > triCut ) continue;
	  if( fabs(dy3) > triCut ) continue;

	  triplet tri;
	  tri.xm = avx;
	  tri.ym = avy;
	  tri.zm = avz;
	  tri.sx = slpx;
	  tri.sy = slpy;
	  tri.lk = 0;
	  tri.ttdmin = 99.9; // isolation [mm]

	  vector <double> ux(3);
	  ux[0] = xA;
	  ux[1] = xB;
	  ux[2] = xC;
	  tri.vx = ux;

	  vector <double> uy(3);
	  uy[0] = yA;
	  uy[1] = yB;
	  uy[2] = yC;
	  tri.vy = uy;

	  triplets.push_back(tri);

	  trixHisto.Fill( avx );
	  triyHisto.Fill( avy );
	  trixyHisto->Fill( avx, avy );
	  tritxHisto.Fill( slpx );
	  trityHisto.Fill( slpy );

	} // cl B

      } // cl C

    } // cl A

    ntriHisto.Fill( triplets.size() );
    if( ldb ) cout << "  triplets " << triplets.size() << endl << flush;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets:

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

      double xmA = triplets[iA].xm;
      double ymA = triplets[iA].ym;
      double zmA = triplets[iA].zm;
      double sxA = triplets[iA].sx; // track slope = angle
      double syA = triplets[iA].sy;
      //double txy = sqrt( sxA*sxA + syA*syA ); // angle

      double zA = DUTz - zmA; // z DUT from mid of triplet
      double xA = xmA + sxA * zA; // triplet impact point on DUT
      double yA = ymA + syA * zA;

      if( ldb ) cout << "  triplet " << iA << endl << flush;

      triycHisto.Fill( yA );
      if(      pb.size() == 0 ) triyc0Histo.Fill( yA );
      else if( pb.size() == 1 ) triyc2Histo.Fill( yA );
      else if( pb.size() == 2 ) triyc2Histo.Fill( yA );
      else if( pb.size() == 3 ) triyc3Histo.Fill( yA );
      else if( pb.size() <= 5 ) triyc4Histo.Fill( yA );
      else if( pb.size() <= 7 ) triyc6Histo.Fill( yA );
      else if( pb.size() <= 9 ) triyc8Histo.Fill( yA );
      else                      triycnHisto.Fill( yA );

      // DUT pixels:

      int nlnk = 0;
      double zmin =  99;
      double zmax = -99;
      double ymin =  99;
      double ymax = -99;
      double dmin =  99;
      double dmax = -99;
      int col0 = 155;
      int col9 =   0;

      double qcol[nx[iDUT]];
      for( int icol = 0; icol < nx[iDUT]; ++icol ) qcol[icol] = 0;

      double pcol[nx[iDUT]];
      for( int icol = 0; icol < nx[iDUT]; ++icol ) pcol[icol] = 0;

      double ycol[nx[iDUT]];
      for( int icol = 0; icol < nx[iDUT]; ++icol ) {
	double zpix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // pixel center
	zpix *= -1; // invert
	double dy = yA - DUTaligny - zpix*tan( DUTtilt ); // w.r.t. mid sensor
	ycol[icol] = dy;
      }

      double sump = 0;
      double sumq = 0;

      for( vector<pixel>::iterator px = pb.begin(); px != pb.end(); ++px ) {

	int icol = px->col;
	int irow = px->row;
	double q = px->q;

	double xpix;
	double zpix;

	xpix = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // pixel center
	zpix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // pixel center
	zpix *= -1; // invert

	double xx = co*xpix - so*zpix; // trn in DUT plane around center
	double zz =+so*xpix + co*zpix;

	cmsdxaHisto.Fill( xx - xA );
	cmssxaHisto.Fill( xx + xA );

	double dx = xx - DUTalignx + xA;

	cmsdxHisto.Fill( dx ); // alignx
	cmsdxvsx.Fill( xA, dx );
	cmsdxvstx.Fill( sxA*1E3, dx*1E3 ); // slope = -dz
	cmsdxvsz.Fill( zz, dx ); // tan(trn) = slope

	if( fabs( dx ) < 0.1 )
	  triyclkHisto.Fill( yA );

	double dy = yA - DUTaligny - zpix*tan( DUTtilt ); // rad

	if( fabs( dy ) < 0.1 )
	  cmsdxcHisto.Fill( dx );

	if( fabs( dx ) < 0.1 ) {
	  cmsdyHisto.Fill( dy ); // aligny
	  cmsdycHisto.Fill( dy ); // aligny
	  cmsdyvsx.Fill( xA, dy );
	  cmsdyvsz.Fill( zz, dy ); // tilt
	}

	if( fabs( dx ) < 0.1 && fabs( dy ) < 0.1 ) {

	  cmspxpHisto.Fill( px->ph );
	  cmspxqHisto.Fill( q );

	  ++nlnk;

	  sumq += q;
	  sump += px->ph;
	  pcol[icol] += px->ph; // project cluster onto cols
	  qcol[icol] += q; // project cluster onto cols

	  if( icol < col0 ) col0 = icol;
	  if( icol > col9 ) col9 = icol;

	  if( zz > zmax ) {
	    zmax = zz;
	    ymax = yA;
	    dmax = dy;
	  }
	  if( zz < zmin ) {
	    zmin = zz;
	    ymin = yA;
	    dmin = dy;
	  }

	} // linked

      } // px

      nlnkHisto.Fill( nlnk );

      // first:

      int icol = 0;
      double zpix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // pixel center
      zpix *= -1; // invert
      double dy0 = yA - DUTaligny - zpix*tan( DUTtilt ); // w.r.t. mid sensor
      if( dy0 > 0.1 ) nlnk0Histo.Fill( nlnk ); // empty

      // mid:

      icol = nx[iDUT]/2;
      zpix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // pixel center
      zpix *= -1; // invert
      double dym = yA - DUTaligny - zpix*tan( DUTtilt ); // w.r.t. mid sensor

      // last:

      icol = nx[iDUT] - 1;
      zpix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // pixel center
      zpix *= -1; // invert
      double dy9 = yA - DUTaligny - zpix*tan( DUTtilt ); // w.r.t. mid sensor
      if( dy9 < -0.1 ) nlnk9Histo.Fill( nlnk ); // empty

      // mid:

      int irow = ny[iDUT]/2;
      double xpix = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // pixel center
      double dxm = xpix - DUTalignx + xA;

      nlnkvsxy.Fill( dxm, dym, nlnk );

      if( fabs( dxm ) > 4.0 )
	nlnk1Histo.Fill( nlnk ); // long
      else if( fabs( dxm ) > 3.8 )
	nlnk2Histo.Fill( nlnk ); // 
      else {
	nlnkvsdy0.Fill( dy0, nlnk );
	nlnkvsdy9.Fill( dy9, nlnk );
	nlnkvsdy.Fill( dy0, dy9, nlnk ); // diag
	if( dy0 < -0.1 && dy9 > 0.1 )
	  nlnk3Histo.Fill( nlnk ); // long, some out-of-time
	else
	  nlnk4Histo.Fill( nlnk ); // short

      }

      double zlng = zmax-zmin;

      if( nlnk > 4 ) { // in-time

	zminHisto.Fill( zmin );
	zmaxHisto.Fill( zmax );
	zlngHisto.Fill( zlng );

	if( zlng > 0.15 ) {
	  q0vszlng.Fill( zlng, sumq/zlng );
	  cmsq0Histo.Fill( sumq/zlng );
	  cmsp0Histo.Fill( sump/zlng );
	}

	if( zmin < -3.7 ) { // edge-on
	  zlngcHisto.Fill( zlng ); // same shape
	  zlngvsy.Fill( yA - DUTaligny + 4*tan( DUTtilt ), zlng ); // y at zz = -4 mm
	}

	ymaxvsz.Fill( zmax, ymax );
	yminvsz.Fill( zmin, ymin );
	dmaxvsz.Fill( zmax, dmax );
	dminvsz.Fill( zmin, dmin );
	zlngvsz.Fill( zmin, zlng );

	for( int icol = 0; icol < nx[iDUT]; ++icol ) {

	  if( icol > col0 && icol < col9 ) {

	    colpHisto.Fill( pcol[icol] );
	    colqHisto.Fill( qcol[icol] );

	  } // filled

	  if( fabs( dxm ) < 3.8  &&
	      dy0 < -0.1 && dy9 > 0.1 ) {
	    colqvsy.Fill( ycol[icol], qcol[icol] );
	    colqyHisto.Fill( ycol[icol], qcol[icol] );
	  }

	} // cols

      } // nlnk

    } // loop triplets iA

    if( ldb ) cout << "done ev " << iev << endl << flush;

    ++iev;

  } while( reader->NextEvent() && iev < lev );

  delete reader;

  cout << "done after " << iev << " events" << endl;
  cout << "resyncs " << nresync << endl;

  histoFile->Write();
  histoFile->Close();

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  return 0;
}

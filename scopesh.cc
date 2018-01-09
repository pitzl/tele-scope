
// Daniel Pitzl, DESY, Sep 2017
// telescope analysis with eudaq and ROC4Sens shallow

// make scopesh
// needs runs.dat
// needs align_24500.dat from tele

// scopesh 31110
// scopesh 31336  120 V
// scopesh 31342   20 V
// scopesh 31343   10 V

#include "eudaq/FileReader.hh"

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TFile.h>
#include <TH1.h> // counting
//#include <TH1D.h> // weighted counts
#include <TH2.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <set>
#include <cmath>

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  double adc;
  double q;
  int ord;
  bool big;
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

// globals:

pixel pb[999]; // global declaration: array of pixel hits
int fNHit; // global

//------------------------------------------------------------------------------
vector < cluster > getClus()
{
  // returns clusters with local coordinates
  // decodePixels should have been called before to fill pixel buffer pb 
  // simple clusterization
  // cluster search radius fCluCut ( allows fCluCut-1 empty pixels)

  const int fCluCut = 1; // clustering: 1 = no gap (15.7.2012)
  //const int fCluCut = 2;

  vector < cluster > v;
  if( fNHit == 0 ) return v;

  int* gone = new int[fNHit];

  for( int i = 0; i < fNHit; ++i )
    gone[i] = 0;

  int seed = 0;

  while( seed < fNHit ) {

    // start a new cluster

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do {
      growing = 0;
      for( int i = 0; i < fNHit; ++i ) {
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

    // look for a new seed = used pixel:

    while( (++seed < fNHit) && gone[seed] );

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
  if( chip0 == 119 ) rot90 = 1;
  if( chip0 == 139 ) rot90 = 1;
  if( chip0 == 142 ) rot90 = 1;
  if( chip0 == 143 ) rot90 = 1;
  if( chip0 == 144 ) rot90 = 1;
  if( chip0 == 147 ) rot90 = 1;
  if( chip0 == 155 ) rot90 = 1;

  bool fifty = 0;
  if( chip0 == 102 ) fifty = 1;
  if( chip0 == 106 ) fifty = 1;
  if( chip0 == 111 ) fifty = 1;
  if( chip0 == 117 ) fifty = 1;
  if( chip0 == 118 ) fifty = 1;
  if( chip0 == 139 ) fifty = 1;
  if( chip0 == 142 ) fifty = 1;
  if( chip0 == 143 ) fifty = 1;
  if( chip0 == 144 ) fifty = 1;
  if( chip0 == 147 ) fifty = 1;
  if( chip0 == 155 ) fifty = 1;

  double upsignx =  1; // w.r.t. telescope
  double upsigny =  1;

  if( rot90 ) {
    upsignx = -1;
    upsigny =  1;
  }

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

  DUTalignFileName << "alignDUT_" << run << ".dat";

  ifstream iDUTalignFile( DUTalignFileName.str() );

  cout << endl;

  if( iDUTalignFile.bad() || ! iDUTalignFile.is_open() ) {
    cout << "no " << DUTalignFileName.str() << ", will bootstrap" << endl;
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

  if( DUTaligniteration <= 1 )
    DUTtilt = DUTtilt0; // from runs.dat

  if( rot90 )
    cout << "DUT 90 degree rotated" << endl;

  // normal vector on DUT surface:
  // N = ( 0, 0, -1 ) on DUT, towards -z
  // transform into tele system:
  // tilt alpha around x
  // turn omega around y

  const double co = cos( DUTturn*wt );
  const double so = sin( DUTturn*wt );
  const double ca = cos( DUTtilt*wt );
  const double sa = sin( DUTtilt*wt );
  const double cf = cos( DUTrot );
  const double sf = sin( DUTrot );

  const double Nx =-ca*so;
  const double Ny = sa;
  const double Nz =-ca*co;

  const double norm = cos( DUTturn*wt ) * cos( DUTtilt*wt ); // length of Nz

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

  double ke = 0.039; // Landau peak at 11 ke

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "scopesh" << run << ".root";

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

  TH1I dutphHisto( "dutph", "DUT PH;ADC-PED [ADC];pixels", 500, -100, 900 );
  TH1I dutdphHisto( "dutdph", "DUT #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );

  TProfile dutnpxvst2( "dutnpxvst2",
	      "DUT pixels vs time;time [s];DUT pixels per event",
	      150, 0, 1500, -0.5, 99.5 );
  TProfile dutnclvst2( "dutnclvst2",
	      "DUT clusters vs time;time [s];DUT clusters per event with pixels",
	      150, 0, 1500, -0.5, 99.5 );
  TProfile dutyldvst2( "dutyldvst2",
	      "DUT yield vs time;time [s];DUT events with pixels",
	      150, 0, 1500, -0.5, 1.5 );
  TProfile dutyldvst6( "dutyldvst6",
	      "DUT yield vs time;time [h];DUT events with pixels",
	      1000, 0, 50, -0.5, 1.5 );

  // DUT:

  TH1I dutpxq1stHisto( "dutpxq1st",
		       "DUT pixel charge 1st;1st pixel charge [ke];1st pixels",
		       100, 0, 25 );

  TH1I dutpxq2ndHisto( "dutpxq2nd",
		       "DUT pixel charge 2nd;2nd pixel charge [ke];2nd pixels",
		       100, 0, 25 );

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
		     52, 0.5, 52.5 );
  TH1I dutncolHisto( "dutncol",
		     "DUT cluster size;cluster size [columns];clusters",
		     52, 0.5, 52.5 );
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
  TH1I dutcol1qHisto( "dutcol1q",
		      "DUT 2nd column charge;2nd column charge [ke];clusters",
		      100, 0, 50 );
  TH1I dutcol2qHisto( "dutcol2q",
		      "DUT 3rd column charge;3rd column charge [ke];clusters",
		      100, 0, 50 );
  TH1I dutcol3qHisto( "dutcol3q",
		      "DUT 4th column charge;4th column charge [ke];clusters",
		      100, 0, 50 );
  TH1I dutcol4qHisto( "dutcol4q",
		      "DUT 5th column charge;5th column charge [ke];clusters",
		      100, 0, 50 );
  TH1I dutcol5qHisto( "dutcol5q",
		      "DUT 6th column charge;6th column charge [ke];clusters",
		      100, 0, 50 );
  TH1I dutcol6qHisto( "dutcol6q",
		      "DUT 7th column charge;7th column charge [ke];clusters",
		      100, 0, 50 );
  TH1I dutcol7qHisto( "dutcol7q",
		      "DUT 8th column charge;8th column charge [ke];clusters",
		      100, 0, 50 );
  TH1I dutcol8qHisto( "dutcol8q",
		      "DUT 9th column charge;9th column charge [ke];clusters",
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

  TH1I z3Histo( "z3",
		       "z3 should be zero;z3 [mm];triplets",
		       100, -0.01, 0.01 );

  TH1I cmssxaHisto( "cmssxa",
			   "DUT + Telescope x;cluster + triplet #Sigmax [mm];clusters",
			   440, -11, 11 );
  TH1I cmsdxaHisto( "cmsdxa",
			   "DUT - Telescope x;cluster - triplet #Deltax [mm];clusters",
			   440, -11, 11 );

  TH1I cmssyaHisto( "cmssya",
			   "DUT + Telescope y;cluster + triplet #Sigmay [mm];clusters",
			   198, -99, 99 ); // shallow needs wide range
  TH1I cmsdyaHisto( "cmsdya",
			   "DUT - Telescope y;cluster - triplet #Deltay [mm];clusters",
			   160, -40, 40 );
  TH1I cmsycHisto( "cmsyc",
		   "x-linked triplet at DUT y;triplet y at DUT [mm];x-linked triplets",
		   120, -6, 6 );

  TH2I * cmsxvsx = new TH2I( "cmsxvsx",
			     "DUT vs Telescope x;track x [mm];DUT x [mm];track-cluster combinations",
			     160, -4, 4, 160, -4, 4 );
  TH2I * cmsyvsy = new TH2I( "cmsyvsy",
			     "DUT vs Telescope y;track y [mm];DUT y [mm];track-cluster combinations",
			     160, -4, 4, 160, -4, 4 );

  TH1I cmsdxHisto( "cmsdx",
			   "DUT - Telescope x;cluster - triplet #Deltax [mm];clusters",
			   200, -0.5, 0.5 );
  TH1I cmsdyHisto( "cmsdy",
			   "DUT - Telescope y;cluster - triplet #Deltay [mm];clusters",
			   200, -10, 10 );

  TH2I * cmsdxvsev1 = new TH2I( "cmsdxvsev1",
				"DUT - Telescope x;event;#Deltax [mm]",
				100, 0, 50000, 100, -5, 5 );
  TH2I * cmsdxvsev2 = new TH2I( "cmsdxvsev2",
				"DUT - Telescope x;event;#Deltax [mm]",
				1000, 0, 1000*1000, 50, -5, 5 );

  TH1I cmsdxcHisto( "cmsdxc",
			   "DUT - Telescope x, cut dy;cluster - triplet #Deltax [mm];clusters",
			   200, -0.5, 0.5 );

  TProfile cmsdxvsx( "cmsdxvsx",
		     "#Deltax vs x;x track [mm];<cluster - triplet #Deltax> [mm]",
		     50, -3.75, 3.75, -0.5, 0.5 );
  TProfile cmsdxvsy( "cmsdxvsy",
		     "#Deltax vs y;y track [mm];<cluster - triplet #Deltax> [mm]",
		     76, -3.8, 3.8, -0.5, 0.5 );
  TProfile cmsdxvstx( "cmsdxvstx",
		      "#Deltax vs #theta_{x};x track slope [rad];<cluster - triplet #Deltax> [mm]",
		      80, -0.004, 0.004, -0.5, 0.5 );

  TH1I cmsdycHisto( "cmsdyc",
		    "#Deltay cut x;cluster - triplet #Deltay [mm];fiducial clusters",
		     200, -10, 10 );

  TProfile cmsdyvsx( "cmsdyvsx",
		     "DUT #Deltay vs x;x track [mm];<cluster - triplet #Deltay> [mm]",
		     50, -3.75, 3.75, -1, 1 );
  TProfile cmsdyvsxn( "cmsdyvsxn",
		      "DUT #Deltay vs x, long;x track [mm];<cluster - triplet #Deltay> [mm]",
		      50, -3.75, 3.75, -1, 1 );
  TProfile cmsdyvsy( "cmsdyvsy",
		     "DUT #Deltay vs y;y track [mm];<cluster - triplet #Deltay> [mm]",
		     80, -8, 8, -1, 1 );
  TProfile cmsdyvsty( "cmsdyvsty",
		      "DUT #Deltay vs #theta_{y};y track slope [rad];<cluster - triplet #Deltay> [mm]",
		      80, -0.002, 0.002, -1, 1 );

  TH1I cmsqc3Histo( "cmsqc3",
		   "column charge;column charge [ke];x-linked columns",
		   100, 0, 20 );
  TH1I cmsqc4Histo( "cmsqc4",
		   "column charge;column charge [ke];x-linked columns",
		   100, 0, 20 );
  TProfile cmsqcvsd( "cmsqcvsd",
		     "DUT charge vs depth;depth [#mum];<column charge> [ke]",
		     150/ptchx[iDUT]*fabs( tan(DUTtilt*wt) ), 0, 150, 0, 50 );
  TH2I * cmsqcvsxHisto = new
    TH2I( "cmsqcvsx",
	  "column charge vs x;x track [mm];column charge [ke];x-linked columns",
	  76, -3.8, 3.8, 50, 0, 10 );
  TProfile cmsqcvsy( "cmsqcvsy",
		     "column charge vs y;y track [mm];<column charge> [ke]",
		     76, -3.8, 3.8, 0, 50 );

  TH1I cmslkxHisto( "cmslkx",
			   "linked triplet at DUT x;triplet x at DUT [mm];linked triplets",
			   220, -11, 11 );
  TH1I cmslkyHisto( "cmslky",
			   "linked triplet at DUT y;triplet y at DUT [mm];linked triplets",
			   120, -6, 6 );

  TH1I cmscolHisto( "cmscol",
		    "DUT linked columns;DUT linked cluster column;linked clusters",
		    nx[iDUT], -0.5, nx[iDUT]-0.5 );
  TH1I cmsrowHisto( "cmsrow",
		    "DUT linked rows;DUT linked cluster row;linked clusters",
		    ny[iDUT], -0.5, ny[iDUT]-0.5 );
  TH1I cmsrowcHisto( "cmsrowc",
		    "x-DUT linked long rows;DUT x- linked long cluster row;x-linked long clusters",
		    ny[iDUT], -0.5, ny[iDUT]-0.5 );

  TH1I cmscolminHisto( "cmscolmin",
		    "DUT x-linked min column;cluster 1st column;x-linked clusters",
		    nx[iDUT], -0.5, nx[iDUT]-0.5 );
  TH1I cmscolmaxHisto( "cmscolmax",
		    "DUT x-linked max column;cluster lst column;x-linked clusters",
		    nx[iDUT], -0.5, nx[iDUT]-0.5 );
  TH1I cmsrowminHisto( "cmsrowmin",
		    "DUT x-linked min row;cluster 1st row;x-linked clusters",
		    ny[iDUT], -0.5, ny[iDUT]-0.5 );
  TH1I cmsrowmaxHisto( "cmsrowmax",
		    "DUT x-linked max row;cluster lst row;x-linked clusters",
		    ny[iDUT], -0.5, ny[iDUT]-0.5 );

  TH1I cmsnpxHisto( "cmsnpx",
		     "linked DUT cluster size;cluster size [pixels];linked fiducial clusters",
		     80, 0.5, 80.5 );
  TH1I cmsncolHisto( "cmsncol",
		     "linked DUT cluster size;cluster size [columns];linked fiducial clusters",
		     52, 0.5, 52.5 );
  TH1I cmsnrowHisto( "cmsnrow",
		     "linked DUT cluster size;cluster size [rows];linked fiducial clusters",
		     80, 0.5, 80.5 );
  TH1I cmsnrowcHisto( "cmsnrowc",
		     "x-linked DUT cluster size;cluster size [rows];x-linked clusters",
		     80, 0.5, 80.5 );

  TProfile cmsncolvsxm( "cmsncolvsxm",
			"DUT cluster size vs xmod;x track mod 0.3 [mm];<cluster size> [columns]",
			100, 0, 0.2, 0, 20 );
  TProfile cmsnrowvsxm( "cmsnrowvsxm",
			"DUT cluster size vs xmod;x track mod 0.3 [mm];<cluster size> [rows]",
			100, 0, 0.2, 0, 80 );

  TProfile cmsncolvsym( "cmsncolvsym",
			"DUT cluster size vs ymod;y track mod 0.2 [mm];<cluster size> [columns]",
			100, 0, 0.2, 0, 20 );
  TProfile cmsnrowvsym( "cmsnrowvsym",
			"DUT cluster size vs ymod;y track mod 0.2 [mm];<cluster size> [rows]",
			100, 0, 0.2, 0, 80 );
  TProfile2D * cmsnpxvsxmym = new
    TProfile2D( "cmsnpxvsxmym",
	      "DUT cluster size vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];<cluster size> [pixels]",
		40, 0, 0.2, 40, 0, 0.2, 0, 20 );

  TH1I cmsq0Histo( "cmsq0",
		   "normal fiducial cluster charge;normal cluster charge [ke];linked fiducial clusters",
		   160, 0, 80 );
  TProfile cmsqxvsx( "cmsqxvsx",
		     "DUT cluster charge vs x;x track [mm];<cluster charge> [ke]",
		     50, -3.75, 3.75, 0, 0.015 ); // cutoff at 5 ke
  TProfile cmsqxvsy( "cmsqxvsy",
		     "DUT cluster charge vs y;y track [mm];<cluster charge> [ke]",
		     76, -3.8, 3.8, 0, 0.015 );
  TProfile2D * cmsqxvsxy = new
		      TProfile2D( "cmsqxvsxy",
				  "DUT cluster charge vs xy;x track [mm];y track [mm];<cluster charge> [ke]",
				  50, -3.75, 3.75, 76, -3.8, 3.8, 0, 0.015 );
  TProfile cmsqxvsxm( "cmsqxvsxm",
		      "DUT cluster charge vs xmod;x track mod 0.3 [mm];<cluster charge> [ke]",
		      100, 0, 0.2, 0, 0.015 );
  TProfile cmsqxvsym( "cmsqxvsym",
		      "DUT cluster charge vs ymod;y track mod 0.2 [mm];<cluster charge> [ke]",
		      100, 0, 0.2, 0, 0.015 );
  TProfile cmsqxvsymn( "cmsqxvsymn",
		       "DUT cluster charge vs ymod no dot;y track mod 0.2 [mm];<no dot cluster charge> [ke]",
		       100, 0, 0.2, 0, 0.015 );
  TProfile2D * cmsqxvsxmym = new
    TProfile2D( "cmsqxvsxmym",
	      "DUT cluster charge vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];<cluster charge> [ke]",
		40, 0, 0.2, 40, 0, 0.2, 0, 0.015 );

  TH1I cmspxqHisto( "cmspxq",
		    "DUT pixel charge linked;pixel charge [ke];linked pixels",
		    100, 0, 25 );
  TH1I cmspxqlHisto( "cmspxql",
		    "DUT pixel charge long linked;pixel charge [ke];pixels in long linked clusters",
		    100, 0, 25 );

  TH1I npxrdHisto = TH1I( "npxrd",
			  "pixels in track road;pixels in track road;tracks",
			  80, 0.5, 80.5 );

  TH1I npxrdqHisto = TH1I( "npxrdq",
			  "charged pixels in track road;charged pixels in track road;tracks",
			  80, 0.5, 80.5 );
  TH1I dcol0Histo = TH1I( "dcol0",
			  "cluster-road 1st pixel;1st cluster-road [pixels];roads",
			  41, -20.5, 20.5 );
  TH1I dcol9Histo = TH1I( "dcol9",
			  "cluster-road lst pixel;lst cluster-road [pixels];roads",
			  41, -20.5, 20.5 );
  TH1I ncolrdHisto = TH1I( "ncolrd",
			   "columns in track road;pixel columns in track road;roads",
			   80, 0.5, 80.5 );

  TH1I dutpxqmidHisto =
    TH1I( "dutpxqmid",
	  "DUT pixel charge middle;pixel charge [ke];middle pixels",
	  100, 0, 25 );
  TH1I dutpxqedgHisto =
    TH1I( "dutpxqedg",
	  "DUT pixel charge edge;pixel charge [ke];edge pixels",
	  100, 0, 25 );

 TH1I yrd0Histo = TH1I( "yrd0",
			"track road entry;track road entry [mm];x-linked tracks",
			200, -10, 10 );
 TH1I yrd9Histo = TH1I( "yrd9",
			"track road exit;track road exit [mm];x-linked tracks",
			200, -10, 10 );
  TH1I cmsdyexHisto =
    TH1I( "cmsdyex",
	  "DUT exit dy;track exit - last pixel #Deltay [mm];tracks",
	 160, -2, 2 );
  TH1I cmsdyinHisto =
    TH1I( "cmsdyin",
	  "DUT entry dy;track entry - last pixel #Deltay [mm];tracks",
	 160, -2, 2 );

  TH1I cmsdyrowHisto =
    TH1I( "cmsdyrow",
	  "DUT row dy;#Deltay [mm];rows",
	 160, -2, 2 );
  TH1I cmsdyrowcHisto =
    TH1I( "cmsdyrowc",
	  "DUT row dy;#Deltay [mm];rows",
	 160, -2, 2 );
  TProfile cmsqyvsy =
    TProfile( "cmsqyvsy",
	      "DUT charge vs track y;dy [mm];<charge> [ke]",
	      80, -2, 2, 0, 99 );
  TH1I cmsdHisto =
    TH1I( "cmsd",
	  "track depth;track depth [#mum];pixels in track road",
	  100, -100, 100 );
  TH1I cmsdqHisto =
    TH1I( "cmsdq",
	  "track depth charged;track depth [#mum];charged pixels in track road",
	  100, -100, 100 );

  TH1I cmsqyHisto =
    TH1I( "cmsqy",
	  "DUT y charge;y charge [ke];y bins in track road",
	  100, 0, 25 );
  TH1I cmsqy1Histo =
    TH1I( "cmsqy1",
	  "DUT 1st y charge;y charge [ke];1st y bin in track road",
	  100, 0, 25 );
  TH1I cmsqy2Histo =
    TH1I( "cmsqy2",
	  "DUT 2nd y charge;y charge [ke];2nd y bin in track road",
	  100, 0, 25 );
  TH1I cmsqy8Histo =
    TH1I( "cmsqy8",
	  "DUT last-1 y charge;y charge [ke];lst-1 y bin in track road",
	  100, 0, 25 );
  TH1I cmsqy9Histo =
    TH1I( "cmsqy9",
	  "DUT last y charge;y charge [ke];lst y bin in track road",
	  100, 0, 25 );
  TH1I cmsqymHisto =
    TH1I( "cmsqym",
	  "DUT mid y charge;y charge [ke];mid y bin in track road",
	  100, 0, 25 );

  TProfile cmsqyvsd =
    TProfile( "cmsqyvsd",
	      "DUT charge vs track depth;track depth [#mum];<pixel charge> [ke]",
	      40, -100, 100, 0, 50 );
  TProfile cmsqyvsdxp =
    TProfile( "cmsqyvsdxp",
	      "DUT charge vs track depth, x > 0;track depth [#mum];<pixel charge> [ke]",
	      40, -100, 100, 0, 50 );
  TProfile cmsqyvsdxm =
    TProfile( "cmsqyvsdxm",
	      "DUT charge vs track depth, x < 0;track depth [#mum];<pixel charge> [ke]",
	      40, -100, 100, 0, 50 );
  TProfile cmsqy1vsd =
    TProfile( "cmsqy1vsd",
	      "DUT charge > 0 vs track depth;track depth [#mum];<pixel charge> [ke]",
	      40, -100, 100, 0, 99 );

  TH2I * cmsqydHisto = new
    TH2I( "cmsqyd",
	  "DUT charge vs track depth;track depth [#mum];pixel charge [ke];entries",
	  40, -100, 100, 249, 0.2, 50 );
  TH1I cmsqyposHisto =
    TH1I( "cmsqypos",
	  "DUT charge per pixel row above mid plane;pixel charge [ke];pixels in track road above mid",
	  100, 0, 25 );
  TH1I cmsqynegHisto =
    TH1I( "cmsqyneg",
	  "DUT charge per pixel row below mid plane;pixel charge [ke];pixels in track road below mid",
	  100, 0, 25 );

  TH1I cmsqyexHisto =
    TH1I( "cmsqyex",
	  "DUT pixel charge, d > 120 #mum;pixel charge [ke];d > 120 #mum pixels",
	  100, 0.25, 25.25 );
  TH1I cmsqy90Histo =
    TH1I( "cmsqy90",
	  "DUT pixel charge, d 90-120 #mum;pixel charge [ke];d 90-120 #mum pixels",
	  100, 0.25, 25.25 );
  TH1I cmsqy60Histo =
    TH1I( "cmsqy60",
	  "DUT pixel charge, d 60-90 #mum;pixel charge [ke];d 60-90 #mum pixels",
	  100, 0.25, 25.25 );
  TH1I cmsqy30Histo =
    TH1I( "cmsqy30",
	  "DUT pixel charge, d 30-60 #mum;pixel charge [ke];d 30-60 #mum pixels",
	  100, 0.25, 25.25 );
  TH1I cmsqy00Histo =
    TH1I( "cmsqy00",
	  "DUT pixel charge, d 0-30 #mum;pixel charge [ke];d 0-30 #mum pixels",
	  100, 0.25, 25.25 );

  TProfile cmsqyvsrm =
    TProfile( "cmsqyvsrm",
	      "DUT pixel charge vs rmod;x track mod 0.2 [mm];<pixel charge> [ke]",
	      100, 0, 0.2, 0, 99 );
  TProfile cmsqy0vsrm =
    TProfile( "cmsqy0vsrm",
	      "DUT deep pixel charge vs rmod;x track mod 0.2 [mm];<pixel charge d < -50 #mum> [ke]",
	      100, 0, 0.2, 0, 99 );      

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
  double prevtludt = 0;

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

    if( chip0 == 114 && run == 31095 && iev == 145447 ) tludt = 0.001; // see cmsdxvsev2 and use ldbt
    if( chip0 == 114 && run == 31095 && iev == 146442 ) tludt = 0.001;
    if( chip0 == 114 && run == 31095 && iev == 147437 ) tludt = 0.001;
    //if( tludt > 0.5 && prevtludt > 0.5 ) tludt = 0.001;
    prevtludt = tludt;

    if( iev < 10 )
      cout << "scopesh processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "scopesh processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "scopesh processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 )
      cout << "scopesh processing  " << run << "." << iev << "  taken " << evsec << endl;

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

      int npx = 0;

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldb )
	  std::cout << ", " << plane.GetX(ipix)
		    << " " << plane.GetY(ipix)
		    << " " << plane.GetPixel(ipix) << flush;

	int ix = plane.GetX(ipix); // global column 0..415
	int iy = plane.GetY(ipix); // global row 0..159
	int adc = plane.GetPixel(ipix); // ADC 0..255

	// skip hot pixels:

	int ipx = ix*ny[ipl] + iy;
	if( hotset[ipl].count(ipx) ) {
	  if( ldb ) cout << " hot" << flush;
	  continue;
	}

	double q = adc;

	hcol[ipl]->Fill( ix+0.5 );
	hrow[ipl]->Fill( iy+0.5 );
	hmap[ipl]->Fill( ix+0.5, iy+0.5 );

	// fill pixel block for clustering:

	pb[npx].col = ix; // col
	pb[npx].row = iy; // row
	pb[npx].adc = adc;
	pb[npx].q = q;
	pb[npx].ord = npx; // readout order
	pb[npx].big = 0;
	++npx;

	if( npx > 990 ) {
	  cout << "pixel buffer overflow in plane " << ipl
	       << ", event " << iev
	       << endl;
	  break;
	}

      } // pix

      hnpx[ipl]->Fill(npx);

      if( ldb ) std::cout << std::endl;

      // clustering:

      fNHit = npx; // for cluster search

      cl[ipl] = getClus();

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

    int npx = 0;

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

      for( unsigned ipx = 0; ipx < vpx.size(); ++ipx ) {

	int col4 = vpx[ipx].col;
	int row4 = vpx[ipx].row;
	double ph4 = vpx[ipx].adc;

	int row1 = row4;
	int row7 = row4;
	double ph1 = ph4;
	double ph7 = ph4;

	for( unsigned jpx = 0; jpx < vpx.size(); ++jpx ) {

	  if( jpx == ipx ) continue;
	  if( vpx[jpx].col != col4 ) continue; // want same column

	  int jrow = vpx[jpx].row;

	  if( jrow < row1 ) {
	    row1 = jrow;
	    ph1 = vpx[jpx].adc;
	  }

	  if( jrow > row7 ) {
	    row7 = jrow;
	    ph7 = vpx[jpx].adc;
	  }

	} // jpx

	if( row4 == row1 ) continue; // Randpixel
	if( row4 == row7 ) continue;

	double dph;
	if( row4 - row1 < row7 - row4 )
	  dph = ph4 - ph1;
	else
	  dph = ph4 - ph7;

	dutdphHisto.Fill( dph );

	//if( dph > 16 ) { // cmspxq left edge 1.3 ke
	if( dph > 12 ) { // cmspxq left edge 1.0 ke, some noise

	  // r4scal.C

	  double U = ( dph - p3[col4][row4] ) / p2[col4][row4];

	  if( U >= 1 )
	    U = 0.9999999; // avoid overflow

	  double vcal = p0[col4][row4] - p1[col4][row4] * log( (1-U)/U ); // inverse Fermi

	  if( fifty ) {
	    pb[npx].col = col4; // 50x50
	    pb[npx].row = row4;
	  }
	  else {
	    pb[npx].col = (col4+1)/2; // 100 um
	    if( col4%2 ) 
	      pb[npx].row = 2*row4 + 0;
	    else
	      pb[npx].row = 2*row4 + 1;
	  }
	  pb[npx].adc = dph;
	  pb[npx].q = ke*vcal; // should become Vcal
	  pb[npx].ord = npx; // readout order
	  pb[npx].big = 0;
	  hmap[iDUT]->Fill( pb[npx].col+0.5, pb[npx].row+0.5 );
	  ++npx;

	  if( npx > 990 ) {
	    cout << "R4S pixel buffer overflow in event " << iev
		 << endl;
	    break;
	  }

	} // dph cut

      } // ipx

      if( ldb ) cout << " npx " << npx << endl << flush;

    } // filled

    readnext = 1;
    if( iev > 88 && dtbdt > 0.5 && tludt < 0.2*dtbdt ) {
      readnext = 0;
      if( ldbt ) cout << "repeat DTB event " << DUTev << endl;
    }
    if( readnext )
      prevdtbtime = dtbtime;

    // DUT clustering:

    hnpx[iDUT]->Fill( npx );

    fNHit = npx; // for cluster search

    cl[iDUT] = getClus();

    hncl[iDUT]->Fill( cl[iDUT].size() );

    for( unsigned icl = 0; icl < cl[iDUT].size(); ++ icl ) {

      hsiz[iDUT]->Fill( cl[iDUT][icl].size );
      //hclph.Fill( cl[iDUT][icl].sum );
      //hclmap->Fill( cl[iDUT][icl].col, cl[iDUT][icl].row );

    } // icl

    dutnpxvst2.Fill( evsec, npx );
    dutnclvst2.Fill( evsec, cl[iDUT].size() );

    int DUTyld = 0;
    if( npx ) DUTyld = 1; // no double counting: events with at least one px
    dutyldvst2.Fill( evsec, DUTyld );
    dutyldvst6.Fill( evsec/3600, DUTyld );

    if( ldb ) cout << "  DUT cl " << cl[iDUT].size() << endl << flush;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT:

    for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

      if( c->vpix.size() < 2 ) continue; // skip 1-pix clusters

      if( c->vpix.size() > 2 ) continue; // want 2-pix clusters

      vector<pixel>::iterator pxa = c->vpix.begin();
      vector<pixel>::iterator pxb = pxa; // lower row (sorted along 80*col + row)
      ++pxb; // higher row in clustering

      if( pxa->ord > 0 && pxb->ord > 0 ) continue; // none is 1st

      double q1 = pxa->q; // read out first
      double q2 = pxb->q;
      if( pxb->ord == 0 ) {
	q1 = pxb->q; // read out first
	q2 = pxa->q;
      }

      dutpxq1stHisto.Fill( q1 ); // read out first on ROC
      dutpxq2ndHisto.Fill( q2 );

    } // DUT cl

    // tsunami corrected clusters:

    for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

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

      for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	int icol = px->col;
	int irow = px->row;

	dutadcHisto.Fill( px->adc );
	dutcolHisto.Fill( icol );
	dutrowHisto.Fill( irow );

	if( icol < colmin ) colmin = icol;
	if( icol > colmax ) colmax = icol;
	if( irow < rowmin ) rowmin = irow;
	if( irow > rowmax ) rowmax = irow;

	double q = px->q; // corrected
	if( q < 0 ) continue;

	qsum += q;
	qcol[icol] += q; // project cluster onto cols
	qrow[irow] += q; // project cluster onto rows

	sumcol += icol*q;
	sumrow += irow*q;

      } // pix

      int ncol = colmax - colmin + 1;
      int nrow = rowmax - rowmin + 1;

      if( colmin > 0 && colmax < nx[iDUT]-2 && rowmin > 0 && rowmax < ny[iDUT]-1 ) {
	dutq0Histo.Fill( qsum * norm );
	dutnpxHisto.Fill( c->size );
	dutncolHisto.Fill( ncol );
	dutnrowHisto.Fill( nrow );
      }

      if( ncol > 2 ) {

	dutcolminHisto.Fill( colmin ); // strong even peaks at large turn 14848
	dutcolmaxHisto.Fill( colmax ); // weaker even peaks

	dutcol0qHisto.Fill(qcol[colmin]);
	if( colmin%2 )
	  dutcol0oddqHisto.Fill(qcol[colmin]);
	else
	  dutcol0eveqHisto.Fill(qcol[colmin]);
	dutcol9qHisto.Fill(qcol[colmax]);

	if( ncol > 1 ) dutcol1qHisto.Fill(qcol[colmin+1]);
	if( ncol > 2 ) dutcol2qHisto.Fill(qcol[colmin+2]);
	if( ncol > 3 ) dutcol3qHisto.Fill(qcol[colmin+3]);
	if( ncol > 4 ) dutcol4qHisto.Fill(qcol[colmin+4]);
	if( ncol > 5 ) dutcol5qHisto.Fill(qcol[colmin+5]);
	if( ncol > 6 ) dutcol6qHisto.Fill(qcol[colmin+6]);
	if( ncol > 7 ) dutcol7qHisto.Fill(qcol[colmin+7]);
	if( ncol > 8 ) dutcol8qHisto.Fill(qcol[colmin+8]);

      } // long

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

    double xcut = 0.1;
    if( run == 31336 )
      xcut = 0.25; // skewed?

    double ycut = 0.2;
    if( fabs(DUTtilt) > 60 )
      ycut = 8;

    int rowcut = 30; // run 31110
    if( run >= 31111 ) rowcut = 20;

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

      double xmA = triplets[iA].xm;
      double ymA = triplets[iA].ym;
      double zmA = triplets[iA].zm;
      double sxA = triplets[iA].sx;
      double syA = triplets[iA].sy;
      //double txy = sqrt( sxA*sxA + syA*syA ); // angle

      double zA = DUTz - zmA; // z DUT from mid of triplet
      double xA = xmA + sxA * zA; // triplet impact point on DUT
      double yA = ymA + syA * zA;

      if( ldb ) cout << "  triplet " << iA << endl << flush;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // tri vs tri: isolation at DUT

      double ttdmin = 99.9;

      for( unsigned int jj = 0; jj < triplets.size(); ++jj ) {

	if( jj == iA ) continue;

	double zj = DUTz - triplets[jj].zm;
	double xj = triplets[jj].xm + triplets[jj].sx * zj; // triplet impact point on DUT
	double yj = triplets[jj].ym + triplets[jj].sy * zj;

	double dx = xA - xj;
	double dy = yA - yj;
	double dd = sqrt( dx*dx + dy*dy );
	if( dd < ttdmin ) ttdmin = dd;

	ttdxHisto.Fill( dx );
	ttdx1Histo.Fill( dx );

      } // jj

      ttdmin1Histo.Fill( ttdmin );
      ttdmin2Histo.Fill( ttdmin );
      triplets[iA].ttdmin = ttdmin;

      bool liso = 0;
      if( ttdmin > 0.6 ) liso = 1; // isolated track

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // intersect inclined track with tilted DUT plane:

      double zc = (Nz*zA - Ny*ymA - Nx*xmA) / (Nx*sxA + Ny*syA + Nz); // from zmA
      double yc = ymA + syA * zc;
      double xc = xmA + sxA * zc;

      double dzc = zc + zmA - DUTz; // from DUT z0 [-8,8] mm

      // transform into DUT system: (passive).
      // large rotations don't commute: careful with order

      double x1 = co*xc - so*dzc; // turn o
      double y1 = yc;
      double z1 = so*xc + co*dzc;

      double zoffs = 0;
      z1 += zoffs; // R4S space above Cu block axis

      double x2 = x1;
      double y2 = ca*y1 + sa*z1; // tilt a
      double z2 =-sa*y1 + ca*z1; // should be zero (in DUT plane). is zero

      double x3 = cf*x2 + sf*y2; // rot
      double y3 =-sf*x2 + cf*y2;
      double z3 = z2 - zoffs; // should be zero (in DUT plane). is zero

      z3Histo.Fill( z3 ); // is zero

      double x4 = upsignx*x3 + DUTalignx; // shift to mid
      double y4 = upsigny*y3 + DUTaligny; // invert y, shift to mid

      // reduce to 100x100 um region:

      double xmod = fmod( 9.000 + x4, 0.2 ); // [0,0.2] mm
      double ymod = fmod( 9.000 + y4, 0.2 ); // [0,0.2] mm

      double rmod = fmod( 9.000 + y4, 0.2 ); // [0,0.2] mm
      if( rot90 ) {
     	rmod = fmod( 9.000 + x4, 0.2 ); // [0,0.2] mm
      }

      // define track road through pixels: (shallow tilt, no turn)

      double yrd9 = y4 + 0.5*DUTthickness * fabs( tan( DUTtilt*wt ) );
      double yrd0 = y4 - 0.5*DUTthickness * fabs( tan( DUTtilt*wt ) );

      double qrd[155][160]; // charge
      double drd[155][160]; // depth
      int npxrd = 0;
      int rdrowmin = 159;
      int rdrowmax =  0;
      int rdcolmin = 154;
      int rdcolmax =  0;

      for( int irow = 0; irow < 160; ++irow ) // x
	for( int icol = 0; icol < 155; ++icol ) // x
	  qrd[icol][irow] = -1; // pixel is outside track road

      if( rot90 ) {

	for( int irow = 0; irow < 160; ++irow ) { // x

	  double xpix = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT];

	  if( fabs( x4 - xpix ) > 1.1 * ptchy[iDUT] ) continue; // 2-row road
	  //if( fabs( x4 - xpix ) > 1.6 * ptchy[iDUT] ) continue; // 3-row road
	  //if( fabs( x4 - xpix ) > 4.1 * ptchy[iDUT] ) continue; // 5-row road, poor x align

	  for( int icol = 0; icol < 155; ++icol ) { // y

	    double ypix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // [mm]

	    drd[icol][irow] = ( y4 - ypix ) / fabs( tan( DUTtilt*wt ) ); // depth, 0 = mid

	    // |   |   |   |   |

	    if( ypix > yrd9 + 1.5 * ptchx[iDUT] ) continue; // with smearing
	    if( ypix < yrd0 - 1.5 * ptchx[iDUT] ) continue;
	    if( ypix > yrd9 + 0.5 * ptchx[iDUT] ) continue; // tighter
	    if( ypix < yrd0 - 0.5 * ptchx[iDUT] ) continue;

	    ++npxrd;

	    qrd[icol][irow] = 0; // pixel is inside track road

	    if( irow < rdrowmin ) rdrowmin = irow;
	    if( irow > rdrowmax ) rdrowmax = irow;
	    if( icol < rdcolmin ) rdcolmin = icol;
	    if( icol > rdcolmax ) rdcolmax = icol;

	  } // col

	} // row

      } // rot90

      else {

	for( int icol = 0; icol < 155; ++icol ) { // x

	  double xpix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT];

	  if( fabs( x4 - xpix ) > 1.1 * ptchx[iDUT] ) continue; // 2-col road

	  for( int irow = 0; irow < 160; ++irow ) { // y

	    double ypix = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // [mm]

	    drd[icol][irow] = ( y4 - ypix ) / fabs( tan( DUTtilt*wt ) ); // depth, 0 = mid

	    // |   |   |   |   |

	    if( ypix > yrd9 + 1.5 * ptchy[iDUT] ) continue; // with smearing
	    if( ypix < yrd0 - 1.5 * ptchy[iDUT] ) continue;
	    if( ypix > yrd9 + 0.5 * ptchy[iDUT] ) continue; // tighter
	    if( ypix < yrd0 - 0.5 * ptchy[iDUT] ) continue;

	    ++npxrd;

	    qrd[icol][irow] = 0; // pixel is inside track road

	    if( irow < rdrowmin ) rdrowmin = irow;
	    if( irow > rdrowmax ) rdrowmax = irow;
	    if( icol < rdcolmin ) rdcolmin = icol;
	    if( icol > rdcolmax ) rdcolmax = icol;

	  } // row

	} // col

      } // not rot90

      npxrdHisto.Fill( npxrd );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT pixel clusters:

      int pxrowmin = 159;
      int pxrowmax =  0;
      int pxcolmin = 154;
      int pxcolmax =  0;

      for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;

	double Q0 = c->charge * norm; // cluster charge normalized to vertical incidence
	double Qx = exp(-Q0/qwid);

	int colmin = 999;
	int colmax = -1;
	int rowmin = 999;
	int rowmax = -1;

	double qcol[nx[iDUT]];
	for( int icol = 0; icol < nx[iDUT]; ++icol ) qcol[icol] = 0;

	double qrow[ny[iDUT]];
	for( int irow = 0; irow < ny[iDUT]; ++irow ) qrow[irow] = 0;

	for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	  int icol = px->col;
	  if( icol <  0 ) continue;
	  if( icol >= nx[iDUT] ) continue;

	  int irow = px->row;
	  if( irow <  0 ) continue;
	  if( irow >= ny[iDUT] ) continue;

	  if( icol < colmin ) colmin = icol;
	  if( icol > colmax ) colmax = icol;
	  if( irow < rowmin ) rowmin = irow;
	  if( irow > rowmax ) rowmax = irow;

	  double q = px->q; // [ke] corrected
	  if( q < 0 ) continue;

	  if( qrd[icol][irow] > -0.5 ) qrd[icol][irow] += q; // road

	  qcol[icol] += q; // project cluster onto cols
	  qrow[irow] += q; // project cluster onto rows

	  // match to road in x:

	  if( rot90 ) {
	    double xpix = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT];
	    if( fabs( x4 - xpix ) > 2.1 * ptchy[iDUT] ) continue;
	  }
	  else {
	    double xpix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT];
	    if( fabs( x4 - xpix ) > 2.1 * ptchx[iDUT] ) continue;
	  }

	  if( irow < pxrowmin ) pxrowmin = irow;
	  if( irow > pxrowmax ) pxrowmax = irow;
	  if( icol < pxcolmin ) pxcolmin = icol;
	  if( icol > pxcolmax ) pxcolmax = icol;

	} // pix

	bool fiducial = 1;
	if( colmin < 1 ) fiducial = 0;
	if( rowmin < 1 ) fiducial = 0;
	if( colmax >= nx[iDUT]-2 ) fiducial = 0; // nx 156, col to 154
	if( rowmax >  ny[iDUT]-2 ) fiducial = 0; // ny 160, row to 159

	if( ! fiducial ) continue;

	int ncol = colmax - colmin + 1;
	int nrow = rowmax - rowmin + 1;
	int nlng = nrow;
	if( rot90 ) nlng = ncol;

	if( nlng < 11 ) continue; // !!

	if( rot90 )
	  ccol = 0.5 * ( colmin + colmax ); // mid from head and tail, overwrite!
	else
	  crow = 0.5 * ( rowmin + rowmax ); // mid from head and tail, overwrite!

	bool lq = 1; // Landau peak
	if( ( Q0 < 5 || Q0 > 18 ) ) lq = 0; // r102

	// DUT - triplet:

	double cmsx = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	double cmsy = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm

	if( rot90 ) {
	  cmsx = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	  cmsy = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	}

	cmsxvsx->Fill( x4, cmsx );
	cmsyvsy->Fill( y4, cmsy );

	// residuals for pre-alignment:

	cmssxaHisto.Fill( cmsx + x3 ); // rot, tilt and turn but no shift
	cmsdxaHisto.Fill( cmsx - x3 ); // peak

	double cmsdx = cmsx - x4; // triplet extrapol

	if( fabs(cmsdx) < xcut ) {
	  cmsycHisto.Fill( yc );
	  cmssyaHisto.Fill( cmsy + y3 );
	  cmsdyaHisto.Fill( cmsy - y3 );
	}

	double cmsdy = cmsy - y4; // 2012 version

	cmsdxHisto.Fill( cmsdx );
	cmsdyHisto.Fill( cmsdy );
	cmsdxvsev1->Fill( iev, cmsdx ); // sync stability
	cmsdxvsev2->Fill( iev, cmsdx ); // sync stability

	if( fabs(cmsdy) < ycut ) {

	  cmsdxcHisto.Fill( cmsdx ); // align: same sign
	  cmsdxvsx.Fill( x4, cmsdx ); // align: same sign
	  cmsdxvsy.Fill( y4, cmsdx ); // align: opposite sign
	  cmsdxvstx.Fill( sxA, cmsdx );

	  if( fabs( cmsdx ) < 0.02 )
	    if( ldbt )
	      cout << "\t\t dx " << cmsdx << endl;

	} // ycut

	// for dy:

	if( fabs(cmsdx) < xcut )
	  cmsnrowcHisto.Fill( nrow );

	if( fabs(cmsdx) < xcut
	    //&& nrow > rowcut
	    ) {

	  cmsrowcHisto.Fill( crow );
	  cmsrowminHisto.Fill( rowmin );
	  cmsrowmaxHisto.Fill( rowmax );
	  cmscolminHisto.Fill( colmin );
	  cmscolmaxHisto.Fill( colmax );
	  yrd0Histo.Fill( yrd0 );
	  yrd9Histo.Fill( yrd9 );

	  double cmsymax = ( rowmax + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	  if( rot90 )
	    cmsymax = ( colmax + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm

	  cmsdyexHisto.Fill( cmsymax-yrd9 );
	  cmsdyinHisto.Fill( cmsymax-yrd0 );

	  cmsdycHisto.Fill( cmsdy );
	  cmsdyvsx.Fill( x4, cmsdy );
	  if( ncol > 28 && ncol < 33 ) // 31336
	    cmsdyvsxn.Fill( x4, cmsdy );
	  cmsdyvsy.Fill( y4, cmsdy );
	  cmsdyvsty.Fill( syA, cmsdy );

	  // cluster road:

	  if( rot90 ) {

	    if( colmin > 0 && colmax < nx[iDUT]-1 && nrow <= 2 ) {

	      for( int icol = colmin; icol <= colmax; ++icol ) {
		double dyc = ( 0.5+colmax-icol ) * ptchx[iDUT]; // [mm]
		double d = dyc / fabs( tan( DUTtilt*wt ) ); // [mm] depth
		cmsqcvsd.Fill( d*1E3, qcol[icol] );
	      }

	      for( int icol = colmin+1; icol < colmax; ++icol ) {
		cmsqcvsy.Fill( y4, qcol[icol] ); // plateau only +- 3 mm
		if( fabs( y4 ) < 3 ) {
		  cmsqcvsxHisto->Fill( x4, qcol[icol] );
		  cmsqc3Histo.Fill( qcol[icol] );
		}
		else
		  cmsqc4Histo.Fill( qcol[icol] );
	      }

	    } // colmin

	  } // rot90

	} // dx

	// xy cuts:

	if( fabs(cmsdx) < xcut  && fabs(cmsdy) < ycut ) {

	  cmslkxHisto.Fill( xc ); // telescope coordinates
	  cmslkyHisto.Fill( yc );
	  cmscolHisto.Fill( ccol ); // map
	  cmsrowHisto.Fill( crow );

	  cmsnpxHisto.Fill( c->size );
	  cmsncolHisto.Fill( ncol );
	  cmsnrowHisto.Fill( nrow );

	  cmsncolvsxm.Fill( xmod, ncol );
	  cmsnrowvsxm.Fill( xmod, nrow );

	  cmsncolvsym.Fill( ymod, ncol ); // within pixel
	  cmsnrowvsym.Fill( ymod, nrow ); // within pixel

	  cmsnpxvsxmym->Fill( xmod, ymod, c->size ); // cluster size map

	  cmsq0Histo.Fill( Q0 ); // Landau

	  // cluster charge profiles, exponential weighting Qx:

	  cmsqxvsx.Fill( x4, Qx );
	  cmsqxvsy.Fill( y4, Qx );
	  cmsqxvsxy->Fill( x4, y4, Qx );
	  cmsqxvsxm.Fill( xmod, Qx ); // Q within pixel
	  cmsqxvsym.Fill( ymod, Qx ); // Q within pixel
	  cmsqxvsxmym->Fill( xmod, ymod, Qx ); // cluster charge profile

	  for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	    cmspxqHisto.Fill( px->q );

	    if( rot90 ) {
	      if( ncol > rowcut )
		cmspxqlHisto.Fill( px->q );
	    }
	    else {
	      if( nrow > rowcut )
		cmspxqlHisto.Fill( px->q );
	    }

	  } // pix

	} // cut xy

      } // loop DUT clusters

      // look at road:

      if( liso &&
 	  yrd0 > -3.9 && // fiducial contained road
 	  yrd9 < 3.9 ) {

	// telescope has out-of-time pile-up: require some match

	int npxrdq = 0;
	for( int irow = rdrowmin; irow <= rdrowmax; ++irow )
	  for( int icol = rdcolmin; icol <= rdcolmax; ++icol )
	    if( qrd[icol][irow] > 0.1 )
	      ++npxrdq;

 	if( rot90 ) { // col is y

 	  if( rdrowmin > 0 && rdrowmax < 159 ) { // fiducial road in x

 	    npxrdqHisto.Fill( npxrdq );

	    if( npxrd > 0 && npxrdq == 0 && iev < 999 )
	      cout
		<< "  ev " << iev
		<< ": empty road"
		<< " row " << rdrowmin
		<< " to " << rdrowmax
		<< ", col " << rdcolmin
		<< " to " << rdcolmax
		<< endl;
	  }

 	  if( rdrowmin > 0 &&
 	      rdrowmax < 159 && // fiducial road in x
 	      npxrdq > 11 ) {

	    if( iev < 999 ) {
	      cout
		<< "  ev " << iev
		<< ": filld road"
		<< " row " << rdrowmin
		<< " to " << rdrowmax
		<< ", col " << rdcolmin
		<< " to " << rdcolmax
		<< " active " << npxrdq
		<< ", col " << pxcolmin
		<< " to " << pxcolmax
		<< " clus";
	      for( unsigned icl = 0; icl < cl[iDUT].size(); ++icl )
		cout << " " << cl[iDUT][icl].size;
	      cout << endl;
	    }

	    dcol0Histo.Fill( rdcolmin - pxcolmin );
	    dcol9Histo.Fill( rdcolmax - pxcolmax );
	    ncolrdHisto.Fill( pxcolmax-pxcolmin+1 );

	    for( int icol = rdcolmin; icol <= rdcolmax; ++icol ) { // y

 	      double ycol = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // [mm]
 	      double dy = y4 - ycol;

 	      double Qy = 0;
 	      double Dp = drd[icol][rdrowmin]; // depth

 	      for( int irow = rdrowmin; irow <= rdrowmax; ++irow ) {

 		if( qrd[icol][irow] < -0.5 ) continue; // pixel is outside track road

 		Qy += qrd[icol][irow];
 		Dp = drd[icol][irow];

 	      } // rows in this col

 	      if( Qy > 0.1 ) {
 		cmsdyrowHisto.Fill( dy );
 		if( yrd0 > -3.9 && yrd9 < 3.9 )
 		  cmsdyrowcHisto.Fill( dy ); // not sharper
 	      }
 	      cmsdHisto.Fill( Dp*1E3 );
 	      if( Qy > 0.1 )
 		cmsdqHisto.Fill( Dp*1E3 );

 	      cmsqyHisto.Fill( Qy );

 	      if(      icol <= rdcolmin )
 		cmsqy1Histo.Fill( Qy );
 	      else if( icol <= rdcolmin+1 )
 		cmsqy2Histo.Fill( Qy );
 	      else if( icol >= rdcolmax )
 		cmsqy9Histo.Fill( Qy );
 	      else if( icol >= rdcolmax-1 )
 		cmsqy8Histo.Fill( Qy );
 	      else
 		cmsqymHisto.Fill( Qy );

 	      cmsqyvsy.Fill( dy, Qy );

 	      cmsqyvsd.Fill( Dp*1E3, Qy );

	      if( x4 > 0 )
		cmsqyvsdxp.Fill( Dp*1E3, Qy );
	      else
		cmsqyvsdxm.Fill( Dp*1E3, Qy );

 	      if( Qy > 0.1 )
 		cmsqy1vsd.Fill( Dp*1E3, Qy );

 	      cmsqydHisto->Fill( Dp*1E3, Qy ); // 2D

 	      if( Dp > 0 )
 		cmsqyposHisto.Fill( Qy );
 	      else
 		cmsqynegHisto.Fill( Qy );

 	      if(      Dp > 0.120 )
 		cmsqyexHisto.Fill( Qy );
 	      else if( Dp > 0.090 )
 		cmsqy90Histo.Fill( Qy );
 	      else if( Dp > 0.060 )
 		cmsqy60Histo.Fill( Qy );
 	      else if( Dp > 0.030 )
 		cmsqy30Histo.Fill( Qy );
 	      else if( Dp > 0.000 )
 		cmsqy00Histo.Fill( Qy );

 	      cmsqyvsrm.Fill( rmod, Qy ); // Q within pixel
 	      if( Dp < -0.050 )
 		cmsqy0vsrm.Fill( rmod, Qy ); // Q within pixel

 	    } // loop road y

 	  } // fid x

 	} // rot90

 	else { // row is y

 	  if( rdcolmin > 0 && rdcolmax < 154 ) // fiducial road in x
 	    npxrdqHisto.Fill( npxrdq );

 	  if( rdcolmin > 0 &&
 	      rdcolmax < 154 && // fiducial road in x
 	      npxrdq > 11 ) {

 	    bool ldbrd = 0;

 	    for( int irow = rdrowmin; irow <= rdrowmax; ++irow ) { // y

 	      double yrow = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // [mm]
 	      double dy = y4 - yrow;

 	      double Qy = 0;
 	      double Dp = drd[rdcolmin][irow]; // depth

 	      for( int icol = rdcolmin; icol <= rdcolmax; ++icol ) { // x

 		if( qrd[icol][irow] < -0.5 ) continue; // pixel is outside track road

 		Qy += qrd[icol][irow];
 		Dp = drd[icol][irow];

 	      } // cols in this row

 	      if( Qy > 0.1 )
 		cmsdyrowHisto.Fill( dy );
 	      cmsdHisto.Fill( Dp*1E3 );
 	      if( Qy > 0.1 )
 		cmsdqHisto.Fill( Dp*1E3 );
 	      cmsqyHisto.Fill( Qy );
 	      if(      irow <= rdrowmin )
 		cmsqy1Histo.Fill( Qy );
 	      else if( irow <= rdrowmin+1 )
 		cmsqy2Histo.Fill( Qy );
 	      else if( irow >= rdrowmax )
 		cmsqy9Histo.Fill( Qy );
 	      else if( irow >= rdrowmax-1 )
 		cmsqy8Histo.Fill( Qy );
 	      else {
 		cmsqymHisto.Fill( Qy );
 		if( Qy < 0.1 )
 		  ldbrd = 0;
 	      }

 	      cmsqyvsy.Fill( dy, Qy );

 	      cmsqyvsd.Fill( Dp*1E3, Qy );
 	      if( Qy > 0.1 )
 		cmsqy1vsd.Fill( Dp*1E3, Qy );

 	      cmsqydHisto->Fill( Dp*1E3, Qy ); // 2D

 	      if( Dp > 0 )
 		cmsqyposHisto.Fill( Qy );
 	      else
 		cmsqynegHisto.Fill( Qy );

 	      if(      Dp > 0.120 )
 		cmsqyexHisto.Fill( Qy );
 	      else if( Dp > 0.090 )
 		cmsqy90Histo.Fill( Qy );
 	      else if( Dp > 0.060 )
 		cmsqy60Histo.Fill( Qy );
 	      else if( Dp > 0.030 )
 		cmsqy30Histo.Fill( Qy );
 	      else if( Dp > 0.000 )
 		cmsqy00Histo.Fill( Qy );

 	    } // loop y

 	    if( ldbrd ) {
 	      cout << iev << endl
 		   << "  road at x " << x4
 		   << " from col " << rdcolmin
 		   << " to " << rdcolmax
 		   << ", y " << y4
 		   << " from row " << rdrowmin
 		   << " to " << rdrowmax
 		   << ", area " << npxrd
 		   << endl;
 	      for( int irow = rdrowmin; irow <= rdrowmax; ++irow ) {
 		cout << "  " << irow;
 		for( int icol = rdcolmin; icol <= rdcolmax; ++icol )
 		  cout << "  " << fixed << setprecision(1) << qrd[icol][irow];
 		cout << setprecision(6) << endl;
 	      }

 	    } // ldb

 	  } // fid x

 	} // not rot90

      } // fid road

    } // loop triplets iA

    if( ldb ) cout << "done ev " << iev << endl << flush;

    ++iev;

  } while( reader->NextEvent() && iev < lev );

  delete reader;

  cout << "done after " << iev << " events" << endl;
  cout << "resyncs " << nresync << endl;

  histoFile->Write();
  histoFile->Close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT alignment:

  double newDUTalignx = DUTalignx;
  double newDUTaligny = DUTaligny;

  if( cmsdxaHisto.GetEntries() > 999 ) {

    if( cmsdxaHisto.GetMaximum() > cmssxaHisto.GetMaximum() ) {
      cout << endl << cmsdxaHisto.GetTitle()
	   << " bin " << cmsdxaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -10, 10 );
      double xpk = cmsdxaHisto.GetBinCenter( cmsdxaHisto.GetMaximumBin() );
      fgp0->SetParameter( 0, cmsdxaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, cmsdxaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, cmsdxaHisto.GetBinContent( cmsdxaHisto.FindBin(xpk-1) ) ); // BG
      cmsdxaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newDUTalignx = fgp0->GetParameter(1);
    }
    else {
      cout << endl << cmssxaHisto.GetTitle()
	   << " bin " << cmssxaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, cmssxaHisto.GetMaximum() ); // amplitude
      double xpk = cmssxaHisto.GetBinCenter( cmssxaHisto.GetMaximumBin() );
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, cmssxaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, cmssxaHisto.GetBinContent( cmssxaHisto.FindBin(xpk-1) ) ); // BG
      cmssxaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newDUTalignx = fgp0->GetParameter(1);
    }

    if( cmsdyaHisto.GetMaximum() > cmssyaHisto.GetMaximum() ) {
      cout << endl << cmsdyaHisto.GetTitle()
	   << " bin " << cmsdyaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, cmsdyaHisto.GetMaximum() ); // amplitude
      double xpk = cmsdyaHisto.GetBinCenter( cmsdyaHisto.GetMaximumBin() );
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, cmsdyaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, cmsdyaHisto.GetBinContent( cmsdyaHisto.FindBin(xpk-1) ) ); // BG
      cmsdyaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newDUTaligny = fgp0->GetParameter(1);
    }
    else {
      cout << endl << cmssyaHisto.GetTitle()
	   << " bin " << cmssyaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, cmssyaHisto.GetMaximum() ); // amplitude
      double xpk = cmssyaHisto.GetBinCenter( cmssyaHisto.GetMaximumBin() );
      fgp0->SetParameter( 1, xpk );
      fgp0->SetParameter( 2, cmssyaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, cmssyaHisto.GetBinContent( cmssyaHisto.FindBin(xpk-1) ) ); // BG
      cmssyaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newDUTaligny = fgp0->GetParameter(1);
    }

  } // cmsdxa

  // finer alignment:

  if( DUTaligniteration > 0 && fabs( newDUTalignx - DUTalignx ) < 0.1 &&
      cmsdxHisto.GetEntries() > 999 ) {

    cout << endl << cmsdxHisto.GetTitle()
	 << " bin " << cmsdxHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -0.5, 0.5 );
    fgp0->SetParameter( 0, cmsdxHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, cmsdxHisto.GetBinCenter( cmsdxHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 8*cmsdxHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, cmsdxHisto.GetBinContent(1) ); // BG
    cmsdxHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newDUTalignx = DUTalignx + fgp0->GetParameter(1);

    // dxvsy -> rot:

    if( cmsdxvsy.GetEntries() > 999 ) {
      cmsdxvsy.Fit( "pol1", "q", "", -midx[iDUT]+0.2, midx[iDUT]-0.2 );
      TF1 * fdxvsy = cmsdxvsy.GetFunction( "pol1" );
      cout << endl << cmsdxvsy.GetTitle()
	   << ": extra rot " << -fdxvsy->GetParameter(1) << endl;
      if( rot90 )
	DUTrot -= upsigny*fdxvsy->GetParameter(1);
      else
	DUTrot += upsigny*fdxvsy->GetParameter(1);
    }

    // dxvstx -> dz:

    if( cmsdxvstx.GetEntries() > 999 ) {
      cmsdxvstx.Fit( "pol1", "q", "", -0.002, 0.002 );
      TF1 * fdxvstx = cmsdxvstx.GetFunction( "pol1" );
      cout << endl << cmsdxvstx.GetTitle()
	   << ": z shift " << -fdxvstx->GetParameter(1)
	   << " mm"
	   << endl;
      DUTz -= upsigny*fdxvstx->GetParameter(1);
    }

  }

  if( DUTaligniteration > 0 && fabs( newDUTaligny - DUTaligny ) < 0.1 &&
      cmsdyHisto.GetEntries() > 999 ) {

    cout << endl << cmsdyHisto.GetTitle()
	 << " bin " << cmsdyHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -0.5, 0.5 );
    fgp0->SetParameter( 0, cmsdyHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, cmsdyHisto.GetBinCenter( cmsdyHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 5*cmsdyHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, cmsdyHisto.GetBinContent(1) ); // BG
    cmsdyHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newDUTaligny = DUTaligny + fgp0->GetParameter(1);

    // dyvsy -> tilt:

    if( fabs(DUTtilt0) > 5 && cmsdyvsy.GetEntries() > 999 ) {
      cmsdyvsy.Fit( "pol1", "q", "", -2, 2 );
      TF1 * fdyvsy = cmsdyvsy.GetFunction( "pol1" );
      cout << endl << cmsdyvsy.GetTitle()
	   << ": slope " << fdyvsy->GetParameter(1)
	   << ", extra tilt " << fdyvsy->GetParameter(1)/wt*sa
	   << " deg"
	   << endl;
      DUTtilt += fdyvsy->GetParameter(1)/wt*sa;
    }

  } // iter

  // write new DUT alignment:
  /*
  ofstream DUTalignFile( DUTalignFileName.str() );

  DUTalignFile << "# DUT alignment for run " << run << endl;
  ++DUTaligniteration;
  DUTalignFile << "iteration " << DUTaligniteration << endl;
  DUTalignFile << "alignx " << newDUTalignx << endl;
  DUTalignFile << "aligny " << newDUTaligny << endl;
  DUTalignFile << "rot " << DUTrot << endl;
  DUTalignFile << "tilt " << DUTtilt << endl;
  DUTalignFile << "turn " << DUTturn << endl;
  DUTalignFile << "dz " << DUTz - zz[2] << endl;

  DUTalignFile.close();

  cout << endl << "wrote DUT alignment iteration " << DUTaligniteration
       << " to " << DUTalignFileName.str() << endl
       << "  alignx " << newDUTalignx << endl
       << "  aligny " << newDUTaligny << endl
       << "  rot    " << DUTrot << endl
       << "  tilt   " << DUTtilt << endl
       << "  turn   " << DUTturn << endl
       << "  dz     " << DUTz - zz[2] << endl
    ;

  cout << endl
       << "DUT efficiency " << 100*effvst3.GetMean(2) << "%"
       << endl;
  */
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  return 0;
}

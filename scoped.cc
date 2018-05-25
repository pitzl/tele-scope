
// Daniel Pitzl, DESY, Sep 2017, Apr 2018
// telescope analysis with eudaq and ROC4Sens edge-on

// make scoped
// needs runs.dat
// needs align_31972.dat from tele

// scoped 31942  irrad 136
// scoped 31972  edge-on turn 2
// scoped 32251  shallow

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
  double DUTtilt0 = 0;
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
  //const double wt = atan(1.0) / 45.0; // pi/180 deg

  //const double qwid = 1.2; // [ke] for Moyal

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
  if( chip0 == 145 ) rot90 = 1;
  if( chip0 == 146 ) rot90 = 1;
  if( chip0 == 147 ) rot90 = 1;
  if( chip0 == 148 ) rot90 = 1;
  if( chip0 == 149 ) rot90 = 1;
  if( chip0 == 150 ) rot90 = 1;
  if( chip0 == 151 ) rot90 = 1;
  if( chip0 == 153 ) rot90 = 1;
  if( chip0 == 155 ) rot90 = 1;
  if( chip0 == 156 ) rot90 = 1;
  if( chip0 == 157 ) rot90 = 1;

  if( chip0 == 163 ) rot90 = 1;
  if( chip0 == 164 ) rot90 = 1;
  if( chip0 == 165 ) rot90 = 1;

  if( chip0 == 174 ) rot90 = 1;
  if( chip0 == 179 ) rot90 = 1;
  if( chip0 == 193 ) rot90 = 1;

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
  //double DUTthickness = 0.150; // [mm]
  //if( chip0 == 157 ) DUTthickness = 0.180; // deep diffused
  //double dycut = 2*DUTthickness/3; // 100 or 120 um

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

  double cf = cos( DUTrot ); // rad
  double sf = sin( DUTrot );
  double ca = cos( DUTtilt ); // rad
  double sa = sin( DUTtilt );
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

  if( chip0 == 116 ) dphcut = 16; // 4 sigma
  if( chip0 == 157 ) dphcut = 14; // 4 sigma

  dphcut = 24; // fresh gain_1, 2 ke, below is noise
  dphcut = 16; // fresh gain_1, 
  dphcut = 12; // fresh gain_1, 

  if( chip0 >= 119 && chip0 <= 138 ) // irrad
    //dphcut = 24; // gain_1
    dphcut = 16; // gain_1 irrad, more nrow (same shape), colq tails not worse

  if( run >= 31635 && run <= 32266 ) { // Feb-Apr 2018
    dphcut = 24; // gain_2
    if( chip0 >= 119 && chip0 <= 138 )
      dphcut = 33; // gain_2 irrad
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT hot pixels:

  cout << endl;

  ostringstream DUThotFileName; // output string stream

  DUThotFileName << "hotDUT_" << run << ".dat";

  ifstream iDUThotFile( DUThotFileName.str() );

  if( iDUThotFile.bad() || ! iDUThotFile.is_open() ) {
    cout << "no " << DUThotFileName.str() << endl;
  }
  else {

    cout << "read DUT hot pixel list from " << DUThotFileName.str() << endl;

    string hash( "#" );
    string pix( "pix" );

    while( ! iDUThotFile.eof() ) {

      string line;
      getline( iDUThotFile, line );
      //cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == pix ) {
	int ix, iy;
	tokenizer >> ix;
	tokenizer >> iy;
	int ipx = ix * 160 + iy;
	hotset[iDUT].insert(ipx);
      }

    } // while getline

  } // hotFile

  iDUThotFile.close();

  cout << "DUT hot " << hotset[iDUT].size() << endl;

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

  TH1I hcol[9];
  TH1I hrow[9];
  TH1I hnpx[9];
  TH2I * hmap[9];

  TH1I hncl[9];
  TH1I hsiz[9];
  TH1I hncol[9];
  TH1I hnrow[9];

  for( int ipl = 0; ipl < 9; ++ipl ) {

    if( ipl == iDUT ) { // R4S
      hcol[ipl] = TH1I( Form( "col%i", ipl ),
			Form( "%i col;col;plane %i pixels", ipl, ipl ), 
			155, 0, 155 );
      hrow[ipl] = TH1I( Form( "row%i", ipl ),
			Form( "%i row;row;plane %i pixels", ipl, ipl ),
			160, 0, 160 );
      hmap[ipl] = new TH2I( Form( "map%i", ipl ),
			    Form( "%i map;col;row;plane %i pixels", ipl, ipl ),
			    155, 0, 155, 160, 0, 160 );
    }
    else {
      hcol[ipl] = TH1I( Form( "col%i", ipl ),
			Form( "%i col;col;plane %i pixels", ipl, ipl ), 
			max( 52, nx[ipl]/4 ), 0, nx[ipl] );
      hrow[ipl] = TH1I( Form( "row%i", ipl ),
			Form( "%i row;row;plane %i pixels", ipl, ipl ),
			max( 80, ny[ipl]/2 ), 0, ny[ipl] );
      hmap[ipl] = new TH2I( Form( "map%i", ipl ),
			    Form( "%i map;col;row;plane %i pixels", ipl, ipl ),
			    max( 52, nx[ipl]/4 ), 0, nx[ipl], max( 80, ny[ipl]/2 ), 0, ny[ipl] );
    }
    hnpx[ipl] = TH1I( Form( "npx%i", ipl ),
		      Form( "%i pixel per event;pixels;plane %i events", ipl, ipl ),
		      200, 0, 200 );

    hncl[ipl] = TH1I( Form( "ncl%i", ipl ),
		      Form( "plane %i cluster per event;cluster;plane %i events", ipl, ipl ),
		      51, -0.5, 50.5 );
    hsiz[ipl] = TH1I( Form( "clsz%i", ipl ),
		      Form( "%i cluster size;pixels/cluster;plane %i clusters", ipl, ipl ),
		      51, -0.5, 50.5 );
    hncol[ipl] = TH1I( Form( "ncol%i", ipl ), 
		       Form( "%i cluster size x;columns/cluster;plane %i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );
    hnrow[ipl] = TH1I( Form( "nrow%i", ipl ),
		       Form( "%i cluster size y;rows/cluster;plane %i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );

  } // planes

  TH2I * hDUTmap = new TH2I( "dutmap",
			     "cool DUT map;col;row;cool DUT pixels",
			     155, 0, 155, 160, 0, 160 );

  TProfile phvsprev( "phvsprev", "Tsunami;previous PH [ADC];<PH> [ADC]", 80, 0, 800, -999, 1999 );
  TH1I dutphHisto( "dutph", "DUT PH;ADC-PED [ADC];pixels", 500, -100, 900 );
  TH1I dutdphHisto( "dutdph", "DUT #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );
  TH1I dutdpiHisto( "dutdpi", "DUT -#DeltaPH;-#DeltaPH [ADC];pixels", 500, -100, 900 );
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
		   "DUT charge / columns;DUT sum charge / columns [ke / 100 #mum];events",
		   200, 0, 40 );

  TH1I dutnpxHisto( "dutnpx",
		     "DUT pixels;DUT size [pixels];events",
		     200, 0.5, 200.5 );
  TH1I dutncolHisto( "dutncol",
		     "DUT columns;DUT size [columns];events",
		     80, 0.5, 80.5 );
  TH1I dutnrowHisto( "dutnrow",
		     "DUT rows;DUT size [rows];events",
		     80, 0.5, 80.5 );
  TH1I dutcol0Histo( "dutcol0",
		     "DUT first column;first DUT column;events",
		     155, -0.5, 154.5 );
  TH1I dutcol9Histo( "dutcol9",
		     "DUT last column;last DUT column;events",
		     155, -0.5, 154.5 );
  TH1I dutcol0qHisto( "dutcol0q",
		      "DUT first column charge;first column charge [ke];events",
		      200, 0, 20 );
  TH1I dutcol9qHisto( "dutcol9q",
		      "DUT last column charge;last column charge [ke];events",
		      200, 0, 20 );
  TH1I dutcol1qHisto( "dutcol1q",
		      "DUT 2nd column charge;2nd column charge [ke];events",
		      400, 0, 80 );
  TH1I dutcol8qHisto( "dutcol8q",
		      "DUT vorletzte column charge;vorletzte column charge [ke];events",
		      400, 0, 80 );

  TH1I dutcolszHisto( "dutcolsz",
		      "DUT column size;column size [rows];columns",
		      21, -0.5, 20.5 );
  TH1I dutcolqHisto( "dutcolq",
		     "DUT column charge;column charge [ke / 100 #mum];columns",
		     200, 0, 20 );
  TProfile2D dutcolqmap( "dutcolqmap",
			 "column charge map;columns;rows;<column charge> [ke / 100 #mum]",
			 78, 0, 78, 80, 0, 320, -2, 20 );

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

  TH1I trixcHisto( "trixc", "triplets x at DUT;x [mm];triplets at DUT",
		   240, -12, 12 );
  TH1I triycHisto( "triyc", "triplets y at DUT;y [mm];triplets at DUT",
		   120, -6, 6 );

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

  TH1I ycolHisto( "ycol",
		  "track depth;track depth [mm];linked pixels",
		  160, -8, 8 );

  TH1I cmsdxmHisto( "cmsdxm",
		    "dxm;dxm [mm];pixels",
		    400, -20, 20 );

  TH1I cmsdxaHisto( "cmsdxa",
		    "DUT - Telescope x;pixel - triplet #Deltax [mm];pixels",
		    440, -11, 11 );
  TH1I cmssxaHisto( "cmssxa",
		    "DUT + Telescope x;pixel + triplet #Sigmax [mm];pixels",
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
  TProfile cmsdxvstx( "cmsdxvstx",
		      "#Deltax vs #theta_{x};x track slope [mrad];<pixel - triplet #Deltax> [#mum]",
		      40, -2, 2, -200, 200 );
  TProfile cmsdxvsz( "cmsdxvsz",
		     "#Deltax vs z;z [mm];<pixel - triplet #Deltax> [mm]",
		     80, -4, 4, -0.5, 0.5 );

  TH1I cmsdxcHisto( "cmsdxc",
		    "DUT - Telescope x;pixel - triplet #Deltax [mm];pixels",
		    200, -0.5, 0.5 );
  TProfile cmsdxvsev( "cmsdxvsev", "DUT - Telescope y vs time;time [events];<#Deltax> [mm]",
		      9000, 0, 1000*9000, -0.1, 0.1 );
  TProfile cmsdxvsxc( "cmsdxvsxc",
		     "#Deltax vs x;x track [mm];<pixel - triplet #Deltax> [mm]",
		     50, -3.75, 3.75, -0.5, 0.5 );
  TProfile cmsdxvstxc( "cmsdxvstxc",
		      "#Deltax vs #theta_{x};x track slope [mrad];<pixel - triplet #Deltax> [#mum]",
		      40, -2, 2, -200, 200 );
  TProfile cmsdxvszc( "cmsdxvszc",
		     "#Deltax vs z;z [mm];<pixel - triplet #Deltax> [mm]",
		     80, -4, 4, -0.5, 0.5 );

  TH1I cmsdyHisto( "cmsdy",
		    "DUT - Telescope y;pixel - triplet #Deltay [mm];pixels",
		    200, -2, 2 );
  TH1I cmsdycHisto( "cmsdyc",
		    "DUT - Telescope y;pixel - triplet #Deltay [mm];pixels",
		    200, -0.5, 0.5 );
  TProfile cmsdyvsev( "cmsdyvsev", "DUT - Telescope y vs time;time [events];<#Deltay> [mm]",
		      9000, 0, 1000*9000, -0.1, 0.1 );
  TProfile cmsdyvsx( "cmsdyvsx",
		     "#Deltay vs x;x track [mm];<pixel - triplet #Deltay> [mm]",
		     50, -3.75, 3.75, -0.2, 0.2 );
  TProfile cmsdyvsz( "cmsdyvsz",
		     "#Deltay vs z;z [mm];<pixel - triplet #Deltay> [mm]",
		     80, -4, 4, -0.2, 0.2 );
  TH2I * cmsdyzHisto = new
    TH2I( "cmsdyz", "pixels dy vs z;z_{DUT} [mm];#Deltay [mm];pixels",
	  80, -4, 4, 100, -0.25, 0.25 );

  TH1I cmsx6Histo( "cmsx6",
		   "track x at pixel;track x at pixel [mm];pixels",
		   100, -5, 5 );
  TH1I cmsy6Histo( "cmsy6",
		   "track y at pixel;track y at pixel [mm];pixels",
		   100, -0.25, 0.25 );
  TH1I cmsz6Histo( "cmsz6",
		   "track z at pixel;track z at pixel [mm];pixels",
		   100, -5, 5 );
  TH1I cmsdy6Histo( "cmsdy6",
		    "check   y6;y6-dy [mm];pixels",
		    100, -0.025, 0.025 );
  TH1I cmsdz6Histo( "cmsdz6",
		    "check   z6;z6-zpix [mm];pixels",
		    100, -0.05, 0.05 );

  TH1I cmsx6cHisto( "cmsx6c",
		   "linked track x at pixel;track x at pixel [mm];linked pixels",
		   100, -5, 5 );
  TH1I cmsy6cHisto( "cmsy6c",
		   "linked track y at pixel;track y at pixel [mm];linked pixels",
		   100, -0.25, 0.25 );
  TH1I cmsz6cHisto( "cmsz6c",
		   "linked track z at pixel;track z at pixel [mm];linked pixels",
		   100, -5, 5 );
  TH2I * triyzlkHisto = new
    TH2I( "triyzlk", "linked tracks y-z;z_{tele} [mm];y_{tele} [mm];linked tracks",
	  80, -4, 4, 100, -0.25, 0.25 );
  TH2I * cmsyzlkHisto = new
    TH2I( "cmsyzlk", "linked pixels y-z;z_{DUT} [mm];y_{DUT} [mm];linked pixels",
	  80, -4, 4, 100, -0.25, 0.25 );

  TH1I pxdyHisto( "pxdy",
		  "track depth;track depth [mm];linked pixels",
		  100, -0.25, 0.25 );
  TH1I cmscolHisto( "col",
		    "DUT linked columns;column;linked pixels",
		    nx[iDUT], -0.5, nx[iDUT]-0.5 );
  TH1I cmsrowHisto( "cmsrow",
		    "DUT linked rows;row;linked pixels",
		    ny[iDUT], -0.5, ny[iDUT]-0.5 );

  TH1I cmspxpHisto( "cmspxp",
		    "DUT pixel PH linked;pixel PH [ADC];linked pixels",
		    500, 0, 1000 );
  TH1I cmspxqHisto( "cmspxq",
		    "DUT pixel charge linked;pixel charge [ke];linked pixels",
		    250, 0, 25 );
  TProfile cmspxqvsxm( "cmspxqvsxm",
		       "pixel charge vs x mod 50;track x mod 50 [#mum];pixel charge [ke]",
		       50, 0, 50, 0, 20 );

  TH1I nlnkHisto( "nlnk",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );

  TH1I trixclkHisto( "trixclk", "linked triplets x at DUT;x [mm];linked triplets at DUT",
		   240, -12, 12 );
  TH1I triyclkHisto( "triyclk", "linked triplets y at DUT;y [mm];linked triplets at DUT",
		     120, -6, 6 );

  TProfile2D nlnkvsxyA( "nlnkvsxyA",
			"link map tele;xc track [mm];yc track [mm];<pixel links>",
			160, -8, 8, 80, -4, 4, -1, 999 );
  TProfile2D nlnkvsxy( "nlnkvsxy",
		       "link map DUT;x6 track [mm];y6 track [mm];<pixel links>",
		       100, -5, 5, 100, -0.25, 0.25, -1, 999 );
  TH1I nlnk1Histo( "nlnk1",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );
  TH1I nlnk2Histo( "nlnk2",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );
  TH1I nlnk3Histo( "nlnk3",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );
  TProfile nlnkvsy( "nlnkvsy",
		    "link map;track y(midz) [mm];<pixel links>",
		    80, -0.2, 0.2, -1, 999 );
  TH1I nlnk4Histo( "nlnk4",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );
  TH1I nlnk5Histo( "nlnk5",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );
  TH1I nlnk6Histo( "nlnk6",
		  "linked pixels;linked pixels;tracks",
		  200, 0.5, 200.5 );

  TH1I zminHisto( "zmin",
		  "first link;zmin [mm];tracks",
		  80, -4, 4 );
  TH1I zmaxHisto( "zmax",
		  "last link;zmax [mm];tracks",
		  80, -4, 4 );
  TH1I zlngHisto( "zlng",
		  "length;z length [mm];tracks",
		  80, 0, 8 );

  TH1I densHisto( "dens",
		  "hit density;hit density [1/mm];tracks",
		  220, 0, 44 );

  TProfile q0vszlng( "q0vszlng",
		     "track charge;track length [mm];<charge/length> [ke/mm]",
		     80, 0, 8, 0, 800 );
  TH1I cmsp0Histo( "cmsp0",
		   "track charge;track charge [ADC/mm];tracks",
		   100, 0, 8000 );
  TH1I cmsq0Histo( "cmsq0",
		   "track charge;track charge [ke/mm];tracks",
		   100, 0, 400 );

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

  TH1I coldxHisto( "coldx",
		   "DUT - Telescope x;column - triplet #Deltax [mm];columns",
		   200, -0.2, 0.2 );
  TH2I coldxvsxHisto( "coldxvsx",
		     "columns dx vs x;track x [mm];columns #Deltax [mm];column-tracks",
		     240, -12, 12, 30, -0.15, 0.15 );

  TProfile colqvsx( "colqvsx",
		    "charge vs x;x [mm];<column charge> [ke / 100 #mum]",
		    100, -5, 5, -2, 20 );
  TProfile2D colqvsxz( "colqvsxz",
		      "column charge xz;x [mm];z [mm];<column charge> [ke / 100 #mum]",
		       88, -4.4, 4.4, 80, -4, 4, -2, 20 );

  TH1I coldxcHisto( "coldxc",
		    "DUT - Telescope x;column - triplet #Deltax [mm];columns",
		    200, -0.2, 0.2 );
  TProfile coldxvsz( "coldxvsz",
		     "#Deltax vs z;z [mm];<column - triplet #Deltax> [mm]",
		     80, -4, 4, -0.1, 0.1 );
  TH1I colpHisto( "colp",
		    "DUT linked column PH;column PH [ADC];linked columns",
		    500, 0, 1000 );
  TH1I colvHisto( "colv",
		    "DUT linked column charge;column charge [Vcal / 100 #mum];linked columns",
		    500, 0, 1000 );
  TH1I colqHisto( "colq",
		    "DUT linked column charge;column charge [ke / 100 #mum];linked columns",
		    400, 0, 80 );

  TH1I nrowHisto( "nrow",
		    "DUT linked column size;column size [rows];linked columns",
		    20, 0.5, 20.5 );
  TProfile nrowvsxm( "nrowvsxm",
		     "columns size vs x mod 50;track x mod 50 [#mum];<column size> [rows]",
		     50, 0, 50, -1, 20 );
  TProfile2D nrowvsxmz( "nrowvsxmz",
			"columns size vs x mod 50 and z;track x mod 50 [#mum];column;<column size> [rows]",
			50, 0, 50, 78, -0.5, 77.5, -1, 20 );
  TProfile nrowvsxmo( "nrowvsxmo",
		      "columns size vs x mod 50;shallow track x mod 50 [#mum];<column size> [rows]",
		      50, 0, 50, -1, 20 );
  TProfile nrowvsxmu( "nrowvsxmu",
		      "columns size vs x mod 50;deep track x mod 50 [#mum];<column size> [rows]",
		      50, 0, 50, -1, 20 );

  TProfile2D nrowvsxy( "nrowvsxy",
		       "column size vs xy;x_{at DUT} [mm];y_{in DUT} [mm];<column size> [rows]",
		       88, -4.4, 4.4, 80, -0.2, 0.2, -1, 20 );
  TProfile2D colqvsxy( "colqvsxy",
		       "column charge xy;x_{at DUT} [mm];y_{in DUT} [mm];<column charge> [ke / 100 #mum]",
		       88, -4.4, 4.4, 30, -0.15, 0.15, -2, 20 );

  TH1I colyHisto( "coly",
		  "track depth;track depth [mm];columns",
		  100, -0.25, 0.25 );
  TH2I colyzHisto( "colyz",
		   "track depth vs z;z [mm];track depth [mm];columns",
		   80, -4, 4, 80, -0.2, 0.2 );
  TH2I coly6zHisto( "coly6z",
		    "y6-z;z [mm];y6 [mm];columns",
		    80, -4, 4, 80, -0.2, 0.2 );
  TH1I coldy6Histo( "coldy6",
		  "track depth - y6;track depth - y6 [mm];columns",
		  100, -0.05, 0.05 );
  TH1I coldz6Histo( "coldz6",
		    "check z6;z6-zcol [mm];columns",
		    100, -0.25, 0.25 );

  TProfile nrowvsy( "nrowvsy",
		    "column size vs depth;y [mm];<column size> [rows]",
		    80, -0.2, 0.2, -1, 20 );
  TProfile nrowvsyi( "nrowvsyi",
		     "iso column size vs depth;y [mm];iso <column size> [rows]",
		     80, -0.2, 0.2, -1, 20 );
  TProfile2D nrowvsyz( "nrowvsyz",
		       "column size vs yz;z [mm];y [mm];<column size> [rows]",
		       80, -4, 4, 80, -0.2, 0.2, -1, 20 );
  TProfile2D nrowvsyzi( "nrowvsyzi",
			"iso column size vs yz;z [mm];y [mm];iso <column size> [rows]",
			80, -4, 4, 80, -0.2, 0.2, -1, 20 );
  TProfile2D nrowvsxmy( "nrowvsxmy",
			"column size vs xmod 50 and depth;track x mod 50 [#mum];y [mm];<column size> [rows]",
			20, 0, 50, 40, -0.1, 0.1, -1, 20 );

  TProfile colpvsy( "colpvsy",
		    "column charge vs depth;y [mm];<column charge> [ADC / 100 #mum]",
		    80, -0.2, 0.2, -50, 500 );
  TProfile colqvsy( "colqvsy",
		    "column charge vs depth;y [mm];<column charge> [ke / 100 #mum]",
		    80, -0.2, 0.2, -2, 20 );
  TProfile2D colqvsyz( "colqvsyz",
		      "column charge yz;z [mm];y [mm];<column charge> [ke / 100 #mum]",
		       80, -4, 4, 80, -0.2, 0.2, -2, 20 );
  TProfile colqvsz( "colqvsz",
		    "charge vs z;z [mm];<column charge> [ke / 100 #mum]",
		    100, -5, 5, -2, 20 );
  TProfile colqvsyc( "colqvsyc",
		     "long column charge vs depth;y [mm];<column charge> [ke / 100 #mum]",
		     80, -0.2, 0.2, -2, 20 );
  TProfile colqvsyl( "colqvsyl",
		     "long column charge vs depth;y [mm];<column charge> [ke / 100 #mum]",
		     80, -0.2, 0.2, -2, 20 );
  TProfile colqvsyd( "colqvsyd",
		     "dense column charge vs depth;y [mm];<column charge> [ke / 100 #mum]",
		     80, -0.2, 0.2, -2, 20 );

  TH1I colq90mHisto( "colq90m",
		    "DUT linked column charge, -150 < d < -80 #mum;column charge [ke / 100 #mum];-150 < d < -80 #mum linked columns",
		    400, 0, 80 );
  TH1I colq80mHisto( "colq80m",
		    "DUT linked column charge, -80 < d < -70 #mum;column charge [ke / 100 #mum];-80 < d < -70 #mum linked columns",
		    400, 0, 80 );
  TH1I colq70mHisto( "colq70m",
		    "DUT linked column charge, -70 < d < -60 #mum;column charge [ke / 100 #mum];-70 < d < -60 #mum linked columns",
		    400, 0, 80 );
  TH1I colq60mHisto( "colq60m",
		    "DUT linked column charge, -60 < d < -50 #mum;column charge [ke / 100 #mum];-60 < d < -50 #mum linked columns",
		    400, 0, 80 );
  TH1I colq50mHisto( "colq50m",
		    "DUT linked column charge, -50 < d < -40 #mum;column charge [ke / 100 #mum];-50 < d < -40 #mum linked columns",
		    400, 0, 80 );
  TH1I colq40mHisto( "colq40m",
		    "DUT linked column charge, -40 < d < -30 #mum;column charge [ke / 100 #mum];-40 < d < -30 #mum linked columns",
		    400, 0, 80 );
  TH1I colq30mHisto( "colq30m",
		    "DUT linked column charge, -30 < d < -20 #mum;column charge [ke / 100 #mum];-30 < d < -20 #mum linked columns",
		    400, 0, 80 );
  TH1I colq20mHisto( "colq20m",
		    "DUT linked column charge, -20 < d < -10 #mum;column charge [ke / 100 #mum];-20 < d < -10 #mum linked columns",
		    400, 0, 80 );
  TH1I colq10mHisto( "colq10m",
		    "DUT linked column charge, -10 < d < -0 #mum;column charge [ke / 100 #mum];-10 < d < -0 #mum linked columns",
		    400, 0, 80 );

  TH1I colq10pHisto( "colq10p",
		    "DUT linked column charge, 0 < d < 10 #mum;column charge [ke / 100 #mum];0 < d < 10 #mum linked columns",
		    400, 0, 80 );
  TH1I colq20pHisto( "colq20p",
		    "DUT linked column charge, 10 < d < 20 #mum;column charge [ke / 100 #mum];10 < d < 20 #mum linked columns",
		    400, 0, 80 );
  TH1I colq30pHisto( "colq30p",
		    "DUT linked column charge, 20 < d < 30 #mum;column charge [ke / 100 #mum];20 < d < 30 #mum linked columns",
		    400, 0, 80 );
  TH1I colq40pHisto( "colq40p",
		    "DUT linked column charge, 30 < d < 40 #mum;column charge [ke / 100 #mum];30 < d < 40 #mum linked columns",
		    400, 0, 80 );
  TH1I colq50pHisto( "colq50p",
		    "DUT linked column charge, 50 < d < 50 #mum;column charge [ke / 100 #mum];40 < d < 50 #mum linked columns",
		    400, 0, 80 );
  TH1I colq60pHisto( "colq60p",
		    "DUT linked column charge, 50 < d < 60 #mum;column charge [ke / 100 #mum];50 < d < 60 #mum linked columns",
		    400, 0, 80 );
  TH1I colq70pHisto( "colq70p",
		    "DUT linked column charge, 60 < d < 70 #mum;column charge [ke / 100 #mum];60 < d < 70 #mum linked columns",
		    400, 0, 80 );
  TH1I colq80pHisto( "colq80p",
		    "DUT linked column charge, 70 < d < 80 #mum;column charge [ke / 100 #mum];70 < d < 80 #mum linked columns",
		    400, 0, 80 );
  TH1I colq90pHisto( "colq90p",
		    "DUT linked column charge, 80 < d #mum;column charge [ke / 100 #mum];80 < d #mum linked columns",
		    400, 0, 80 );

  TH1I cmsncolHisto( "cmsncol",
		     "track length columns;track length [columns];tracks",
		     80, 0.5, 80.5 );
  TProfile q0vsncol( "q0vsncol",
		     "charge vs length;track length [columns];<charge/length> [ke / 0.1 mm]",
		     78, 0.5, 78.5, 0, 80 );
  TH1I colp0Histo( "colp0",
		   "track charge;track charge [ADC / 0.,1 mm];tracks",
		   100, 0, 800 );
  TH1I colq0Histo( "colq0",
		   "track charge;track charge [ke / 0.1 mm];tracks",
		   400, 0, 80 );

  TH1I colq1Histo( "colq1",
		      "DUT 2nd column charge;2nd column charge [ke];tracks",
		      400, 0, 80 );
  TH1I colq8Histo( "colq8",
		      "DUT vorletzte column charge;vorletzte column charge [ke];tracks",
		      400, 0, 80 );

  TH1I hq[78];
  for( int lng = 1; lng < 78; ++lng )
    hq[lng] = TH1I( Form( "q%i", lng ),
		    Form( "q/l for l = %i;charge/length [ke/columns];columns", lng ),
		    //100, 0, 20 );
		    400, 0, 80 ); // full range

  TH2I * q1q0Histo = new
    TH2I( "q1q0", "next column correlation;q0 [ke];q1 [ke];column pairs",
	  80, 0, 40, 80, 0, 40 );
  TH2I * qdqsHisto = new
    TH2I( "qdqs", "next column correlation;q0 + q1 [ke];q1 - q0 [ke];column pairs",
	  80, 0, 80, 80, -40, 40 );

  TH2I * q2q0Histo = new
    TH2I( "q2q0", "2nd column correlation;q0 [ke];q2 [ke];column pairs",
	  80, 0, 40, 80, 0, 40 );

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

    // Get next eudaq event:

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

      if( run > 28000 && run <= 32516 &&
	  ipl > 0 && ipl < 7 ) // 2017, eudaq 1.6: Mimosa 1..6, DUT 7, REF 8, QAD 9
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

	hcol[ipl].Fill( ix+0.5 );
	hrow[ipl].Fill( iy+0.5 );
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

      hnpx[ipl].Fill( pb.size() );

      if( ldb ) std::cout << std::endl;

      // clustering:

      cl[ipl] = getClus(pb);

      if( ldb ) cout << "clusters " << cl[ipl].size() << endl;

      hncl[ipl].Fill( cl[ipl].size() );

      for( vector<cluster>::iterator c = cl[ipl].begin(); c != cl[ipl].end(); ++c ) {

	hsiz[ipl].Fill( c->size );
	hncol[ipl].Fill( c->ncol );
	hnrow[ipl].Fill( c->nrow );

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
	hcol[iDUT].Fill( col+0.5 );
	hrow[iDUT].Fill( row+0.5 );
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

	  bool cool = 1;
	  int hpx = col4 * 160 + row4;
	  if( hotset[iDUT].count(hpx) ) {
	    if( ldb ) cout << " hot" << flush;
	    cool = 0;
	  }

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

	  // r4scal.C

	  double U = ( dph - p3[col4][row4] ) / p2[col4][row4];

	  if( U >= 1 )
	    U = 0.9999999; // avoid overflow

	  double vcal = p0[col4][row4] - p1[col4][row4] * log( (1-U)/U ); // inverse Fermi

	  // linear gain from Vcal at Landau peak:

	  double v = 160; // fresh gain_1
	  if( chip0 >= 119 && chip0 <= 138 )
	    v = 130; // irrad gain_1 low IA
	  if( fifty )
	    v = 100;
	  double t = ( v - p0[col4][row4] ) / p1[col4][row4]; // from Vcal
	  double a = p3[col4][row4] + p2[col4][row4] / ( 1 + exp(-t) ); // [ADC]
	  double g = v / a; // mv/ADC

	  if( cool ) {
	    phvsprev.Fill( phprev, ph4 ); // Tsunami
	    dutdphHisto.Fill( dph );
	    dutdpiHisto.Fill( -dph );
	    if( row4 > 10 && row4 < 150 ) dutdph0Histo.Fill( dph ); // skip noisy regions
	    dutvvsdph.Fill( dph, vcal );
	    dutgHisto.Fill( g );
	    linvvsdph.Fill( dph, g*dph );
	  }

	  //vcal = g*dph; // overwrite !!

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

	  phprev = px4->ph;

	  if( dph > dphcut ) {

	    hmap[iDUT]->Fill( col4+0.5, row4+0.5 ); // R4S pixels, for hot channel masking

	    if( cool ) {
	      pb.push_back(px);
	      hDUTmap->Fill( col4+0.5, row4+0.5 );
	    }

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

    hnpx[iDUT].Fill( pb.size() );

    dutnpxvst2.Fill( evsec, pb.size() );

    int DUTyld = 0;
    if( pb.size() ) DUTyld = 1; // no double counting: events with at least one px
    dutyldvst2.Fill( evsec, DUTyld );
    dutyldvst6.Fill( evsec/3600, DUTyld );

    if( ldb ) cout << "  DUT px " << pb.size() << endl << flush;

    ++iev;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT:

    int colmin = 999;
    int colmax = -1;
    int rowmin = 999;
    int rowmax = -1;

    double qsum = 0;
    double sumcol = 0;
    double sumrow = 0;

    int zcol[nx[iDUT]]; // size [rows]
    double qcol[nx[iDUT]];
    double rcol[nx[iDUT]];
    for( int icol = 0; icol < nx[iDUT]; ++icol ) {
      zcol[icol] = 0;
      qcol[icol] = 0;
      rcol[icol] = 0;
    }

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

      qsum += q;
      qcol[icol] += q;
      rcol[icol] += irow*q; // for COG
      ++zcol[icol];

      qrow[irow] += q;

      sumcol += icol*q;
      sumrow += irow*q;

    } // pix

    int ncol = colmax - colmin + 1;
    int nrow = rowmax - rowmin + 1;

    dutnpxHisto.Fill( pb.size() );

    dutq0Histo.Fill( qsum / ncol );
    dutncolHisto.Fill( ncol );
    dutnrowHisto.Fill( nrow );

    dutcol0Histo.Fill( colmin );
    dutcol9Histo.Fill( colmax );

    if( ncol > 0 ) {
      dutcol0qHisto.Fill( qcol[colmin] );
      dutcol9qHisto.Fill( qcol[colmax] );
    }
    if( ncol > 12 ) {
      dutcol1qHisto.Fill( qcol[colmin+1] );
      dutcol8qHisto.Fill( qcol[colmax-1] );
    }

    for( int icol = colmin+1; icol < colmax; ++icol ) {

      dutcolszHisto.Fill( zcol[icol] );
      dutcolqHisto.Fill( qcol[icol] );

      double crow = -1;
      if( zcol[icol] > 0 )
	crow = rcol[icol]/qcol[icol];

      dutcolqmap.Fill( icol, crow, qcol[icol] );

    } // cols

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

    //if( triplets.size() > 15 ) continue; // next event, does not help

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets:

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

      double xmA = triplets[iA].xm;
      double ymA = triplets[iA].ym;
      double zmA = triplets[iA].zm;
      double sxA = triplets[iA].sx; // track slope = angle
      double syA = triplets[iA].sy;
      //double txy = sqrt( sxA*sxA + syA*syA ); // angle

      double dzA = DUTz - zmA; // z from mid of triplet to mid DUT
      double xAc = xmA + sxA * dzA; // track at z_mid(DUT)
      double yAc = ymA + syA * dzA;

      if( ldb ) cout << "  triplet " << iA << endl << flush;

      trixcHisto.Fill( xAc );
      triycHisto.Fill( yAc );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // tri vs tri: isolation at DUT

      double ttdmin = 99.9;

      for( unsigned int jB = 0; jB < triplets.size(); ++jB ) {

	if( jB == iA ) continue;

	double zB = DUTz - triplets[jB].zm;
	double xB = triplets[jB].xm + triplets[jB].sx * zB; // triplet impact point on DUT
	double yB = triplets[jB].ym + triplets[jB].sy * zB;

	double dx = xAc - xB;
	double dy = yAc - yB;
	double dd = sqrt( dx*dx + dy*dy );
	if( dd < ttdmin ) ttdmin = dd;

	ttdxHisto.Fill( dx );
	ttdx1Histo.Fill( dx );

      } // jB

      ttdmin1Histo.Fill( ttdmin );
      ttdmin2Histo.Fill( ttdmin );
      triplets[iA].ttdmin = ttdmin;

      bool liso = 0;
      //if( ttdmin > 0.3 ) liso = 1;
      //if( ttdmin > 0.6 ) liso = 1;
      if( ttdmin > 0.9 ) liso = 1;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT pixels:

      double ycol[nx[iDUT]];
      double dxmin[nx[iDUT]];
      for( int icol = 0; icol < nx[iDUT]; ++icol ) {
	ycol[icol] = 9; // [mm] out of range
	dxmin[icol] = 99; // [mm] out of range
      }

      double dxmcut = 0.2 + 4*fabs(DUTturn); // [mm] cut on x of row at zpix = 0, but allow for turn
      //double dxmcut = 0.4 + 4*fabs(DUTturn); // does not help

      // y from track in each column: slow

      const double ypix = 0; // mid, dummy place holder for symmetry

      for( int irow = 0; irow < ny[iDUT]; ++irow ) {

	double xpix = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // pixel center

	double dxm = xpix - DUTalignx + xAc; // at z_mid(DUT)

	if( DUTaligniteration > 1 && fabs( dxm ) > dxmcut ) continue; // speedup

	for( int icol = 0; icol < nx[iDUT]; ++icol ) {

	  double zpix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // pixel center
	  zpix *= -1; // rot90: invert

	  double x1 = co*xpix - so*zpix; // turn omega
	  double y1 = ypix;
	  double z1 = so*xpix + co*zpix;

	  double x2 = x1;
	  double y2 = ca*y1 + sa*z1; // tilt alpha
	  double z2 =-sa*y1 + ca*z1;

	  double x3 = cf*x2 + sf*y2; // rot
	  double y3 =-sf*x2 + cf*y2;
	  double z3 = z2;

	  double dz = z3 + DUTz - zmA; // z of pixel from mid of triplet A
	  double xA = xmA + sxA * dz; // triplet impact point on DUT
	  double yA = ymA + syA * dz; // track A at pixel

	  double dx = xA + x3 - DUTalignx; // opposite x sign

	  if( fabs(dx) < fabs( dxmin[icol] ) ) {
	    dxmin[icol] = dx;
	    ycol[icol] = yA - y3 - DUTaligny; // dy
	  }

	} // cols

      } // rows

      for( int icol = 0; icol < nx[iDUT]; ++icol )
	ycolHisto.Fill( ycol[icol] );

      // pixels with charge:

      int nlnk = 0;
      double zmin =  99;
      double zmax = -99;
      double ymin =  99;
      double ymax = -99;
      double dmin =  99;
      double dmax = -99;
      int col0 = 155;
      int col9 =   0;

      double sump = 0;
      double sumq = 0;

      double pcol[nx[iDUT]];
      double qcol[nx[iDUT]];
      double rcol[nx[iDUT]];
      int    nrow[nx[iDUT]];
      for( int icol = 0; icol < nx[iDUT]; ++icol ) {
	pcol[icol] = 0;
	qcol[icol] = 0;
	rcol[icol] = 0;
	nrow[icol] = 0;
      }

      for( vector<pixel>::iterator px = pb.begin(); px != pb.end(); ++px ) {

	int icol = px->col;
	int irow = px->row;

	double xpix = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // pixel center

	double dxm = xpix - DUTalignx + xAc; // at z_mid(DUT)
	cmsdxmHisto.Fill( dxm );

	if( DUTaligniteration > 1 && fabs( dxm ) > dxmcut ) continue; // speedup (allow +-0.2/4 turn)

	double zpix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // pixel center
	zpix *= -1; // rot90: invert

	// transform pixel into telescope system:

	double x1 = co*xpix - so*zpix; // trn in DUT plane around center
	double y1 = ypix;
	double z1 = so*xpix + co*zpix;

	double x2 = x1;
	double y2 = ca*y1 + sa*z1; // tilt alpha
	double z2 =-sa*y1 + ca*z1;

	double x3 = cf*x2 + sf*y2; // rot
	double y3 =-sf*x2 + cf*y2;
	double z3 = z2;

	double dz = z3 + DUTz - zmA; // z of pixel from mid of triplet A
	double xA = xmA + sxA * dz; // triplet impact point on DUT
	double yA = ymA + syA * dz; // track A at pixel

	cmsdxaHisto.Fill( xA - x3 );
	cmssxaHisto.Fill( xA + x3 );

	double dx = xA + x3 - DUTalignx;

	double dy = yA - y3 - DUTaligny;

	cmsdxHisto.Fill( dx ); // alignx
	cmsdxvsx.Fill( xA, dx );
	cmsdxvstx.Fill( sxA*1E3, dx*1E3 ); // slope = -alignz
	cmsdxvsz.Fill( z3, dx ); // tan(trn) = slope

	if( fabs( dy ) < 0.07 ) {
	  cmsdxcHisto.Fill( dx );
	  cmsdxvsev.Fill( iev, dx );
	  cmsdxvsxc.Fill( xA, dx );
	  cmsdxvstxc.Fill( sxA*1E3, dx*1E3 ); // slope = -alignz
	  cmsdxvszc.Fill( z3, dx ); // tan(trn) = slope
	}

	if( fabs( dx ) < 0.100 ) {
	  cmsdyHisto.Fill( dy ); // aligny
	  cmsdycHisto.Fill( dy ); // aligny
	  cmsdyvsev.Fill( iev, dy );
	  cmsdyvsx.Fill( xA, dy );
	  cmsdyvsz.Fill( z1, dy ); // tilt
	  cmsdyzHisto->Fill( z3, dy ); // OK
	  pxdyHisto.Fill( dy );
	}

	// for xmod: transform track into DUT system:

	double xt = -xA + DUTalignx;
	double yt =  yA - DUTaligny;
	double zt = z3;

	double x4 = cf*xt - sf*yt; // rot
	double y4 = sf*xt + cf*yt;
	double z4 = zt;

	double x5 = x4;
	double y5 = ca*y4 - sa*z4; // tilt alpha
	double z5 = sa*y4 + ca*z4;

	double x6 = co*x5 + so*z5; // trn in DUT plane around center
	double y6 = y5; // track y at pixel, should be dy
	double z6 =-so*x5 + co*z5; // should be zpix

	cmsx6Histo.Fill( x6 );
	cmsy6Histo.Fill( y6 );
	cmsz6Histo.Fill( z6 );
	cmsdy6Histo.Fill( dy-y6 ); // should be zero, tails +-0.01
	cmsdz6Histo.Fill( z6-zpix ); // tails +-0.03
	if( fabs( dx ) < 0.100 )
	  cmsyzlkHisto->Fill( z3, y6 ); // OK

	double xmod5 = fmod( 9.000 + x6, 0.05 );

	// for depth profile:

	if( fabs( dx ) < 0.100 && // wide (25 um pitch)
	    //if( fabs( dx ) < 0.050 && // tight
	    fabs( dy ) < 0.150 ) { // defines plot range

	  cmsx6cHisto.Fill( x6 );
	  cmsy6cHisto.Fill( y6 );
	  cmsz6cHisto.Fill( z6 );

	  triyzlkHisto->Fill( z3, yA-DUTaligny );

	  cmscolHisto.Fill( icol );
	  cmsrowHisto.Fill( irow );

	  double q = px->q;

	  cmspxpHisto.Fill( px->ph );
	  cmspxqHisto.Fill( q );
	  cmspxqvsxm.Fill( xmod5*1E3, q );

	  ++nlnk;

	  sumq += q;
	  sump += px->ph;
	  pcol[icol] += px->ph; // project cluster onto cols
	  qcol[icol] += q; // project cluster onto cols
	  rcol[icol] += irow*q; // for COG
	  ++nrow[icol];

	  if( icol < col0 ) col0 = icol;
	  if( icol > col9 ) col9 = icol;

	  if( z3 > zmax ) {
	    zmax = z3;
	    ymax = yA;
	    dmax = dy;
	  }
	  if( z3 < zmin ) {
	    zmin = z3;
	    ymin = yA;
	    dmin = dy;
	  }

	} // linked

      } // px

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      nlnkHisto.Fill( nlnk ); // pixel links per triplet track

      if( nlnk > 9 ) {
	trixclkHisto.Fill( xAc );
	triyclkHisto.Fill( yAc );
      }

      // transform track into DUT system:

      double xt = -xAc + DUTalignx; // track at zDUT
      double yt =  yAc - DUTaligny;
      double zt = 0; // z mid

      double x4 = cf*xt - sf*yt; // rot
      double y4 = sf*xt + cf*yt;
      double z4 = zt;

      double x5 = x4;
      double y5 = ca*y4 - sa*z4; // tilt alpha
      double z5 = sa*y4 + ca*z4;

      double x6 = co*x5 + so*z5; // trn in DUT plane around center
      double y6 = y5; // track y at pixel
      //double z6 =-so*x5 + co*z5;

      nlnkvsxyA.Fill( xAc, yAc, nlnk ); // tele system
      nlnkvsxy.Fill( x6, y6, nlnk ); // DUT system

      if( fabs( x6 ) > 4.0 )
	nlnk1Histo.Fill( nlnk ); // short

      else if( fabs( x6 ) > 3.8 )
	nlnk2Histo.Fill( nlnk ); // few

      else {

	nlnk3Histo.Fill( nlnk ); // long, some out-of-time
	nlnkvsy.Fill( y6, nlnk );

	if(      fabs( y6 ) < 0.075 )
	  nlnk4Histo.Fill( nlnk );
	else if( fabs( y6 ) < 0.100 )
	  nlnk5Histo.Fill( nlnk );
	else
	  nlnk6Histo.Fill( nlnk );

      } // x6

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

      if( nlnk > 1 ) { // in-time track (or noise or random coincidence)
	//if( nlnk > 4 ) { // not better

	double zlng = ptchx[iDUT] + zmax - zmin; // [mm]
	double dens = nlnk/zlng; // hit density [1/mm]

	ymaxvsz.Fill( zmax, ymax );
	yminvsz.Fill( zmin, ymin );
	dmaxvsz.Fill( zmax, dmax );
	dminvsz.Fill( zmin, dmin );

	if( fabs( x6 ) < 3.8 && fabs( y6 ) < 0.07 ) {

	  zminHisto.Fill( zmin );
	  zmaxHisto.Fill( zmax );
	  zlngHisto.Fill( zlng );
	  zlngvsz.Fill( zmin, zlng );

	  densHisto.Fill( dens ); // peak at 12.5

	  if( zlng > 0.15 ) {
	    q0vszlng.Fill( zlng, sumq/zlng );
	    cmsp0Histo.Fill( sump/zlng );
	    cmsq0Histo.Fill( sumq/zlng );
	  }

	  if( zmin < -3.7 ) { // edge-on
	    zlngcHisto.Fill( zlng ); // same shape
	    zlngvsy.Fill( yAc - DUTaligny + 4*tan( DUTtilt ), zlng ); // y at zz = -4 mm
	  }

	} // x-y cuts

	int ncol = 0;
	int nup = 0; // outliers

	for( int icol = 0; icol < nx[iDUT]; ++icol ) { // rot90: along the track

	  //double crow = 0; // left, causes underflows
	  double crow = ny[iDUT]/2; // mid
	  if( nrow[icol] > 0 ) crow = rcol[icol] / qcol[icol]; // COG

	  double xcol = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // col center
	  double zcol = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // col center
	  zcol *= -1; // rot90: invert

	  double x1 = co*xcol - so*zcol; // trn in DUT plane around center
	  double y1 = ypix;
	  double z1 = so*xcol + co*zcol;

	  double x2 = x1;
	  double y2 = ca*y1 + sa*z1; // tilt alpha
	  double z2 =-sa*y1 + ca*z1; // track z at col

	  double x3 = cf*x2 + sf*y2; // rot
	  //double y3 =-sf*x2 + cf*y2;
	  //double z3 = z2;

	  double dz = z2 + DUTz - zmA; // z of col from mid of triplet A
	  double xA = xmA + sxA * dz; // triplet impact point on DUT
	  double yA = ymA + syA * dz; // track A at col

	  double dx = xA + x3 - DUTalignx;

	  if( nrow[icol] > 0 ) {
	    coldxHisto.Fill( dx );
	    coldxvsxHisto.Fill( xA, dx );
	  }

	  // for xmod: transform track into DUT system:

	  double xt = -xA + DUTalignx;
	  double yt =  yA - DUTaligny;
	  double zt = z2;

	  double x4 = cf*xt - sf*yt; // rot
	  double y4 = sf*xt + cf*yt;
	  double z4 = zt;

	  double x5 = x4;
	  double y5 = ca*y4 - sa*z4; // tilt alpha
	  double z5 = sa*y4 + ca*z4;

	  double x6 = co*x5 + so*z5; // trn in DUT plane around center
	  double y6 = y5; // should be ycol
	  double z6 =-so*x5 + co*z5; // should be zcol

	  coldy6Histo.Fill( y6-ycol[icol] ); // is zero
	  coldz6Histo.Fill( z6-zcol ); // tails +-0.07 mm

	  double xmod5 = fmod( 9.000 + x6, 0.05 );

	  if( icol > col0 && icol < col9 ) { // inside road with charge

	    colqvsx.Fill( x6, qcol[icol] );
	    colqvsxz.Fill( x6, zcol, qcol[icol] );

	  } // road

	  if( qcol[icol] > 0.1 ) ++ncol;

	  nrowvsxy.Fill( x6, ycol[icol], nrow[icol] );
	  colqvsxy.Fill( x6, ycol[icol], qcol[icol] );

	  //if( fabs( x6 ) < 3.2 ) { // does not help
	  if( fabs( x6 ) < 3.8 ) {

	    colyHisto.Fill( ycol[icol] );
	    colyzHisto.Fill( zcol, ycol[icol] ); // inclined ?
	    coly6zHisto.Fill( zcol, y6 );

	    nrowvsy.Fill( ycol[icol], nrow[icol] );
	    nrowvsyz.Fill( zcol, ycol[icol], nrow[icol] );
	    nrowvsxmy.Fill( xmod5*1E3, ycol[icol], nrow[icol] );
	    if( liso ) {
	      nrowvsyi.Fill( ycol[icol], nrow[icol] );
	      nrowvsyzi.Fill( zcol, ycol[icol], nrow[icol] );
	      if( ycol[icol] > 0.1 && nrow[icol] ) ++nup;
	    }

	    colpvsy.Fill( ycol[icol], pcol[icol] );
	    colqvsy.Fill( ycol[icol], qcol[icol] );
	    colqvsyz.Fill( zcol, ycol[icol], qcol[icol] );

	    if( fabs( ycol[icol] ) < 0.070 ) {

	      colqvsz.Fill( zcol, qcol[icol] );

	      if( icol > col0 && icol < col9 ) { // inside road with charge

		coldxcHisto.Fill( dx );
		coldxvsz.Fill( zcol, dx ); // for turn angle alignment: coldxvsz->Fit("pol1")

		colpHisto.Fill( pcol[icol] );
		colvHisto.Fill( qcol[icol]/ke );
		colqHisto.Fill( qcol[icol] );

		nrowHisto.Fill( nrow[icol] );
		nrowvsxm.Fill( xmod5*1E3, nrow[icol] );
		nrowvsxmz.Fill( xmod5*1E3, icol, nrow[icol] ); // alignment check
		if( ycol[icol] < 0 )
		  nrowvsxmo.Fill( xmod5*1E3, nrow[icol] );
		else
		  nrowvsxmu.Fill( xmod5*1E3, nrow[icol] );

	      } // road

	    } // y cut

	    if( x6 > -3.8 && x6 < 2 && zcol < 3.4 )
	      colqvsyc.Fill( ycol[icol], qcol[icol] );

	    if( nlnk > 10 )
	      colqvsyl.Fill( ycol[icol], qcol[icol] );

	    if( dens > 8 && dens < 19 )
	      colqvsyd.Fill( ycol[icol], qcol[icol] );

	    if(      ycol[icol] < -0.080 )
	      colq90mHisto.Fill( qcol[icol] );
	    else if( ycol[icol] < -0.070 )
	      colq80mHisto.Fill( qcol[icol] );
	    else if( ycol[icol] < -0.060 )
	      colq70mHisto.Fill( qcol[icol] );
	    else if( ycol[icol] < -0.050 )
	      colq60mHisto.Fill( qcol[icol] );
	    else if( ycol[icol] < -0.040 )
	      colq50mHisto.Fill( qcol[icol] );
	    else if( ycol[icol] < -0.030 )
	      colq40mHisto.Fill( qcol[icol] );
	    else if( ycol[icol] < -0.020 )
	      colq30mHisto.Fill( qcol[icol] );
	    else if( ycol[icol] < -0.010 )
	      colq20mHisto.Fill( qcol[icol] );
	    else if( ycol[icol] < -0.000 )
	      colq10mHisto.Fill( qcol[icol] );
	    else if( ycol[icol] <  0.010 )
	      colq10pHisto.Fill( qcol[icol] );
	    else if( ycol[icol] <  0.020 )
	      colq20pHisto.Fill( qcol[icol] );
	    else if( ycol[icol] <  0.030 )
	      colq30pHisto.Fill( qcol[icol] );
	    else if( ycol[icol] <  0.040 )
	      colq40pHisto.Fill( qcol[icol] );
	    else if( ycol[icol] <  0.050 )
	      colq50pHisto.Fill( qcol[icol] );
	    else if( ycol[icol] <  0.060 )
	      colq60pHisto.Fill( qcol[icol] );
	    else if( ycol[icol] <  0.070 )
	      colq70pHisto.Fill( qcol[icol] );
	    else if( ycol[icol] <  0.080 )
	      colq80pHisto.Fill( qcol[icol] );
	    else
	      colq90pHisto.Fill( qcol[icol] );

	  } // fiducial x

	} // loop icol

	if( nup > 99 )
	  cout << iev << ": nup " << nup
	       << ", x " << -xAc + DUTalignx
	       << ", y " << yAc - DUTaligny
	       << ", y6 " << y6
	       << ", syA " << syA*1E3
	       << ", lnk " << nlnk
	       << ", col0 " << col0
	       << ", col9 " << col9
	       << ", Q/L " << sumq / ncol // expect 10 ke/100 um
	       << ", dens " << dens
	       << endl;

	// Landau vs length:

	if( fabs( x6 ) < 3.8 ) {

	  if( fabs( y6 ) < 0.07 ) {

	    cmsncolHisto.Fill( ncol );

	    if( ncol > 0 ) {
	      q0vsncol.Fill( ncol, sumq/ncol );
	      colp0Histo.Fill( sump/ncol );
	      colq0Histo.Fill( sumq/ncol );
	    }
	    if( ncol > 3 ) {
	      colq1Histo.Fill( qcol[col0+1] );
	      colq8Histo.Fill( qcol[col9-1] );
	    }

	  } // y6

	  int mxl = colmax-colmin-1; // skip 1st and lst

	  for( int lng = 1; lng <= mxl; ++lng ) {

	    for( int col1 = colmin+1; col1 <= colmax-lng; col1 += lng ) { // step along track

	      int nc = 0;
	      double sumq = 0;

	      for( int col = col1; col < col1+lng; ++col ) // columns in lng

		if( fabs( ycol[col] ) < 0.070 ) {
		  ++nc;
		  sumq += qcol[col];
		}

	      if( nc ) hq[lng].Fill( sumq/nc ); // normalized

	    } // col1

	  } // lng

	  // correlations:

	  if( mxl > 1 ) {

	    for( int col0 = colmin+1; col0 <= colmax-2; col0 += 2 ) { // step along track

	      int col1 = col0 + 1;

	      if( fabs( ycol[col0] ) < 0.070 &&
		  fabs( ycol[col1] ) < 0.070 ) {
		q1q0Histo->Fill( qcol[col0], qcol[col1] );
		qdqsHisto->Fill( qcol[col0] + qcol[col1], qcol[col1] - qcol[col0] );
	      }
	    } // col0

	  } // mxl

	  if( mxl > 2 ) {

	    for( int col0 = colmin+1; col0 <= colmax-3; col0 += 2 ) { // step along track

	      int col2 = col0 + 2; // 2nd

	      if( fabs( ycol[col0] ) < 0.070 &&
		  fabs( ycol[col2] ) < 0.070 )
		q2q0Histo->Fill( qcol[col0], qcol[col2] );

	    } // col0

	  } // mxl

	} // x6

      } // nlnk

    } // loop triplets iA

    if( ldb ) cout << "done ev " << iev << endl << flush;

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

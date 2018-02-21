
// Daniel Pitzl, DESY, Sep 2017, Jan 2018: run-wise
// telescope analysis with eudaq and ROC4Sens and MOD

// make scopes
// needs runs.dat
// needs align_24500.dat from tele

// scopes 24500
// scopes -l 99999 24500
// scopes -t 2.2 25200
// scopes -s 24093
// scopes -l 99999 28027
// scopes -s 28037
// scopes 31210 # 3D
// scopes -l 312100 31425

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
#include <set> // for hotset
#include <list>
#include <cmath>
#include <time.h> // clock_gettime

using namespace std;
using namespace eudaq;

struct evInfo {
  uint64_t tlutime;
  uint64_t dtbtime;
  string filled;
};

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
  double xm; // mid point [mm]
  double ym;
  double zm;
  double sx; // slope [rad]
  double sy;
  bool lk; // linked
  double ttdmin; // distance top next track [mm]
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

  // nothing left

  delete gone;

  return v; // vector of clusters
}

//------------------------------------------------------------------------------
int main( int argc, char * argv[] )
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
  string modgainFileName( "/home/pitzl/psi/dtb/tst400/D4028-ia25-trim40-2016-04-gaincal.dat" );
  int weib = 3;

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
    string MODGAIN( "modgain" );
    string WEIB( "weib" );
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

      if( tag == CHIP ) {
	tokenizer >> chip0;
	continue;
      }

      if( tag == TILT ) {
	tokenizer >> DUTtilt0;
	continue;
      }

      if( tag == GeV ) {
	tokenizer >> pbeam;
	continue;
      }

      if( tag == GEO ) {
	tokenizer >> geoFileName;
	continue;
      }

      if( tag == GAIN ) {
	tokenizer >> gainFileName;
	continue;
      }

      if( tag == MODGAIN ) {
	tokenizer >> modgainFileName;
	continue;
      }

      if( tag == weib ) {
	tokenizer >> weib;
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
	<< "  MOD gain file " << modgainFileName << endl
	<< "  Weibull version " << weib << endl
	<< endl;
    else {
      cout << "run " << run << " not found in runs.dat" << endl;
      return 1;
    }

  } // runsFile

  runsFile.close();

  const double fTLU = 384E6; // 384 MHz TLU clock

  double fDTB = 39.99679E6; // 40 MHz DTB clock from ddt and ddtvsdt
  if( run >= 31592 )  // 16.12.2017 different DTB
    fDTB = 39.996943E6; // 40 MHz DTB clock from ddt and ddtvsdt: fDTB*(1-slope*1E-3)

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
  // telescope alignment:

  int MIMaligniteration = 0;
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
	tokenizer >> MIMaligniteration;

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

  double qwid = 1.1; // [ke] for Moyal in 150 um
  if( chip0 >= 300 ) qwid = 1.6; // 230 um 3D

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
  if( chip0 == 151 ) rot90 = 1;
  if( chip0 == 155 ) rot90 = 1;

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
  if( chip0 == 149 ) fifty = 1;
  if( chip0 == 151 ) fifty = 1;
  if( chip0 == 152 ) fifty = 1;
  if( chip0 == 155 ) fifty = 1;
  if( chip0 == 158 ) fifty = 1;
  if( chip0 == 159 ) fifty = 1;
  if( chip0 == 160 ) fifty = 1;
  if( chip0 >= 300 ) fifty = 1; // 3D

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

  double norm = cos( DUTturn*wt ) * cos( DUTtilt*wt ); // length of Nz

  if( rot90 )
    cout << "DUT 90 degree rotated" << endl;

  // DUT Cu window in x: from sixdtvsx

  double xminCu = -6.5;
  double xmaxCu =  6.5;

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

  double ke = 0.037; // Landau peak at 11 ke  chip 102  1002.dat
  if( run >= 31148 ) ke = 0.039; // chip 111  1004.dat
  if( run >= 31163 ) ke = 0.035; // chip 117  1006.dat
  if( run >= 31173 ) ke = 0.037; // chip 117  1008.dat
  if( run >= 31210 ) ke = 0.068; // chip 332  3D 230 um at 17.2
  if( run >= 31237 ) ke = 0.050; // chip 352  3D 230 um at 17.2
  if( run >= 31295 ) ke = 0.040; // chip 143  at 11 ke
  if( run >= 31304 ) ke = 0.0385; // chip 144  at 11 ke
  if( run >= 31309 ) ke = 0.0274; // chip 142  at 11 ke
  if( run >= 31351 ) ke = 0.025; // chip 142  at 11 ke
  if( run >= 31355 ) ke = 0.039; // chip 144  at 11 ke
  if( run >= 31390 ) ke = 0.031; // chip 102  at 11 ke
  if( run >= 31425 ) ke = 0.036; // chip 144  at 11 ke
  if( run >= 31473 ) ke = 0.0296; // chip 144  at 11 ke
  if( run >= 31488 ) ke = 0.035; // chip 160  at 11 ke
  if( run >= 31512 ) ke = 0.068; // chip 332  3D 230 um
  if( run >= 31562 ) ke = 0.050; // chip 352  3D 230 um
  if( run >= 31592 ) ke = 0.031; // chip 152  deep diff 180 um at 11 ke
  if( run >= 31610 ) ke = 0.035; // chip 160  at 11 ke
  if( run >= 31614 ) ke = 0.030; // chip 152  deep diff 180 um at 11 ke

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOD:

  int iMOD = 6;

  int MODaligniteration = 0;
  double MODalignx = 0.0;
  double MODaligny = 0.0;
  double MODrot = 0.0;
  double MODtilt = 17.2; // [deg]
  double MODturn =-27.0; // [deg]
  double MODz = 90 + zz[4];

  ostringstream MODalignFileName; // output string stream

  MODalignFileName << "alignMOD_" << run << ".dat";

  ifstream iMODalignFile( MODalignFileName.str() );

  cout << endl;

  if( iMODalignFile.bad() || ! iMODalignFile.is_open() ) {
    cout << "no " << MODalignFileName.str() << ", will bootstrap" << endl;
  }
  else {

    cout << "read MODalignment from " << MODalignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string alignx( "alignx" );
    string aligny( "aligny" );
    string rot( "rot" );
    string tilt( "tilt" );
    string turn( "turn" );
    string dz( "dz" );

    while( ! iMODalignFile.eof() ) {

      string line;
      getline( iMODalignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) 
	tokenizer >> MODaligniteration;

      double val;
      tokenizer >> val;
      if(      tag == alignx )
	MODalignx = val;
      else if( tag == aligny )
	MODaligny = val;
      else if( tag == rot )
	MODrot = val;
      else if( tag == tilt )
	MODtilt = val;
      else if( tag == turn )
	MODturn = val;
      else if( tag == dz )
	MODz = val + zz[4];

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  iMODalignFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Mod gain:

  double m0[16][52][80];
  double m1[16][52][80];
  double m2[16][52][80];
  double m3[16][52][80];
  double m4[16][52][80];

  ifstream modgainFile( modgainFileName.c_str() );

  if(! modgainFile ) {
    cout << "modgain file " << modgainFileName << " not found" << endl;
    return 1;
  }
  else {
    cout << endl << "using MOD gain file " << modgainFileName << endl;
    int roc;
    int col;
    int row;
    double a0, a1, a2, a3, a4, a5;

    while( modgainFile >> roc ) {
      modgainFile >> col;
      modgainFile >> row;
      modgainFile >> a0;
      modgainFile >> a1;
      modgainFile >> a2;
      modgainFile >> a3;
      modgainFile >> a4;
      modgainFile >> a5;
      m0[roc][col][row] = a0;
      m1[roc][col][row] = a1;
      m2[roc][col][row] = a2;
      m3[roc][col][row] = a3;
      m4[roc][col][row] = a4;
    }

  } // modgainFile open

  double mke = 0.367; // [ke] to get mod q0 peak at 22 ke

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "scopes" << run << ".root";

  TFile * histoFile = new TFile( rootFileName.str(  ).c_str(  ), "RECREATE" );

  // book histos:

  double f = 5.6/pbeam;

  TH1I hdttlu( "dttlu", "TLU time between events;TLU time between events log_{10}(#Deltat [s]);events",
	       120, -4, 2 );
  TH1D hdttluw( "dttluw",
		"TLU time between events weighted;TLU time between events log_{10}(#Deltat [s]);#Deltat events [s]",
		120, -4, 2 );
  TH1I hdttlu1( "dttlu1", "TLU time between events;TLU time between events [s];events",
		120, 0, 0.06 );
  TH1I hdttlu2( "dttlu2", "TLU time between events;TLU time between events [s];events",
		500, 0, 5 );
  TH1I hdtdtb( "dtdtb", "DTB time between events;DTB time between events log_{10}(#Deltat [s]);events",
	       120, -4, 2 );

  TH1I hddt( "ddt", "#Deltadt TLU - DTB;TLU - DTB #Deltadt [ms];events", 200, -1, 1 );
  TProfile ddtvsdt( "ddtvsdt", "#Deltadt TLU - DTB;time to previous event [s];<#Deltadt TLU - DTB> [ms]",
		    500, 0, 5, -1e9, 1e9 );
  TProfile ddtvsev1( "ddtvsev1", "#Deltadt TLU - DTB;event;<#Deltadt TLU - DTB> [ms]",
		    100, 0, 50000, -1e9, 1e9 );
  TProfile ddtvsev2( "ddtvsev2", "#Deltadt TLU - DTB;event;<#Deltadt TLU - DTB> [ms]",
		    5000, 0, 5000*1000, -1e9, 1e9 );

  TH1I t1Histo( "t1", "event time;event time [s];events / 0.01 s", 100, 0, 1 );
  TH1I t2Histo( "t2", "event time;event time [s];events / s", 300, 0, 300 );
  TH1I t3Histo( "t3", "event time;event time [s];events / 10 s", 150, 0, 1500 );
  TH1I t4Histo( "t4", "event time;event time [s];events / 10 s", 600, 0, 6000 );
  TH1I t5Histo( "t5", "event time;event time [s];events / 100 s", 600, 0, 60000 );
  TH1I t6Histo( "t6", "event time;event time [h];events / 3 min", 1000, 0, 50 );

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
			max( 156, nx[ipl]/4 ), 0, nx[ipl] );
      hrow[ipl] = TH1I( Form( "row%i", ipl ),
			Form( "%i row;row;plane %i pixels", ipl, ipl ),
			max( 160, ny[ipl]/2 ), 0, ny[ipl] );
      hmap[ipl] = new TH2I( Form( "map%i", ipl ),
			    Form( "%i map;col;row;plane %i pixels", ipl, ipl ),
			    max( 156, nx[ipl]/4 ), 0, nx[ipl], max( 160, ny[ipl]/2 ), 0, ny[ipl] );
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

  // driplets:

  TH1I hdx35 = TH1I( "dx35", "3-5 dx;3-5 dx [mm];cluster pairs", 100, -f, f );
  TH1I hdy35 = TH1I( "dy35", "3-5 dy;3-5 dy [mm];cluster pairs", 100, -f, f );

  TH1I hdridx = TH1I( "dridx", "driplet dx;driplet dx [mm];driplets", 100, -0.1*f, 0.1*f );
  TH1I hdridy = TH1I( "dridy", "driplet dy;driplet dy [mm];driplets", 100, -0.1*f, 0.1*f );

  TH1I hdridxc = TH1I( "dridxc", "driplet dx;driplet dx [mm];driplets", 100, -0.1*f, 0.1*f );
  TH1I hdridyc = TH1I( "dridyc", "driplet dy;driplet dy [mm];driplets", 100, -0.1*f, 0.1*f );

  TProfile dridxvsy =
    TProfile( "dridxvsy",
	      "driplet dx vs y;driplet yB [mm];<driplets #Deltax> [mm]",
	      120, -6, 6, -0.05*f, 0.05*f );
  TProfile dridyvsx =
    TProfile( "dridyvsx",
	      "driplet dy vs x;driplet xB [mm];<driplets #Deltay> [mm]",
	      240, -12, 12, -0.05*f, 0.05*f );

  TProfile dridxvstx =
    TProfile( "dridxstx",
	      "driplet dx vs slope x;driplet slope x [rad];<driplets #Deltax> [mm]",
	      60, -0.003, 0.003, -0.05*f, 0.05*f );
  TProfile dridyvsty =
    TProfile( "dridysty",
	      "driplet dy vs slope y;driplet slope y [rad];<driplets #Deltay> [mm]",
	      60, -0.003, 0.003, -0.05*f, 0.05*f );

  TH1I drixHisto = TH1I( "drix", "driplets x;x [mm];driplets",
			  240, -12, 12 );
  TH1I driyHisto = TH1I( "driy", "driplets x;y [mm];driplets",
			  120, -6, 6 );
  TH2I * drixyHisto = new
    TH2I( "drixy", "driplets x-y;x [mm];y [mm];driplets",
	  240, -12, 12, 120, -6, 6 );
  TH1I dritxHisto = TH1I( "dritx", "driplet slope x;slope x [rad];driplets",
			    100, -0.005*f, 0.005*f );
  TH1I drityHisto = TH1I( "drity", "driplet slope y;slope y [rad];driplets",
			    100, -0.005*f, 0.005*f );

  TH1I ndriHisto = TH1I( "ndri", "driplets;driplets;events", 51, -0.5, 50.5 );

  TH1I drizixHisto = TH1I( "drizix",
			   "driplets x-z intersection;driplets intersect z [mm];driplet pairs",
			   140, -500, 200 );
  TH1I drizix2Histo = TH1I( "drizix2",
			    "driplets x-z intersection;driplets intersect z [mm];driplet pairs",
			    210, -10000, 500 );
  TH1I driziyHisto = TH1I( "driziy",
			   "driplets y-z intersection;driplets intersect z [mm];driplet pairs",
			   140, -500, 200 );

  TH1I dddmin1Histo = TH1I( "dddmin1",
			    "telescope driplets isolation;driplets min #Delta_{xy} [mm];driplet pairs",
			    100, 0, 1 );
  TH1I dddmin2Histo = TH1I( "dddmin2",
			    "telescope driplets isolation;driplets min #Delta_{xy} [mm];driplet pairs",
			    150, 0, 15 );

  // MOD vs driplets:

  TH1I modsxaHisto = TH1I( "modsxa",
			   "MOD + driplet x;MOD cluster + driplet #Sigmax [mm];MOD clusters",
			   1280, -32, 32 );
  TH1I moddxaHisto = TH1I( "moddxa",
			   "MOD - driplet x;MOD cluster - driplet #Deltax [mm];MOD clusters",
			   1280, -32, 32 );

  TH1I modsyaHisto = TH1I( "modsya",
			   "MOD + driplet y;MOD cluster + driplet #Sigmay [mm];MOD clusters",
			   320, -8, 8 );
  TH1I moddyaHisto = TH1I( "moddya",
			   "MOD - driplet y;MOD cluster - driplet #Deltay [mm];MOD clusters",
			   320, -8, 8 );

  TH1I moddxHisto = TH1I( "moddx",
			   "MOD - driplet x;MOD cluster - driplet #Deltax [mm];MOD clusters",
			   500, -2.5, 2.5 );
  TH1I moddxcHisto = TH1I( "moddxc",
			   "MOD - driplet x;MOD cluster - driplet #Deltax [mm];MOD clusters",
			   200, -0.5, 0.5 );
  TH1I moddxcqHisto = TH1I( "moddxcq",
			    "MOD - driplet x Landau peak;MOD cluster - driplet #Deltax [mm];Landau peak MOD clusters",
			    500, -0.5, 0.5 );
  TProfile moddxvsx =
    TProfile( "moddxvsx",
	      "MOD #Deltax vs x;x track [mm];<cluster - driplet #Deltax> [mm]",
	      216, -32.4, 32.4, -2.5, 2.5 );
  TProfile moddxvsy =
    TProfile( "moddxvsy",
	      "MOD #Deltax vs y;y track [mm];<cluster - driplet #Deltax> [mm]",
	      160, -8, 8, -2.5, 2.5 );
  TProfile moddxvstx =
    TProfile( "moddxvstx",
	      "MOD #Deltax vs #theta_{x};x track slope [rad];<cluster - driplet #Deltax> [mm]",
	      80, -0.002, 0.002, -2.5, 2.5 );

  TH1I moddyHisto = TH1I( "moddy",
			  "MOD - driplet y;MOD cluster - driplet #Deltay [mm];MOD clusters",
			  200, -0.5, 0.5 );
  TH1I moddycHisto = TH1I( "moddyc",
			   "MOD - driplet y;MOD cluster - driplet #Deltay [mm];MOD clusters",
			   200, -0.5, 0.5 );
  TH1I moddycqHisto = TH1I( "moddycq",
			    "MOD - driplet y Landau peak;MOD cluster - driplet #Deltay [mm];Landau peak MOD clusters",
			    500, -0.5, 0.5 );
  TProfile moddyvsx =
    TProfile( "moddyvsx",
	      "MOD #Deltay vs x;x track [mm];<cluster - driplet #Deltay> [mm]",
	      216, -32.4, 32.4, -0.5, 0.5 );
  TProfile moddyvsy =
    TProfile( "moddyvsy",
	      "MOD #Deltay vs y;y track [mm];<cluster - driplet #Deltay> [mm]",
	      160, -8, 8, -0.5, 0.5 );
  TProfile moddyvsty =
    TProfile( "moddyvsty",
	      "MOD #Deltay vs #theta_{y};y track slope [rad];<cluster - driplet #Deltay> [mm]",
	      80, -0.002, 0.002, -0.5, 0.5 );

  TH1I modnpxHisto =
    TH1I( "modnpx",
	  "MOD linked clusters;MOD cluster size [pixels];linked MOD cluster",
	  20, 0.5, 20.5 );

  TH1I modqHisto = TH1I( "modq",
			 "MOD linked clusters;MOD cluster charge [ke];linked MOD cluster",
			 80, 0, 80 );
  TH1I modq0Histo = TH1I( "modq0",
			 "MOD linked clusters;MOD normal cluster charge [ke];linked MOD cluster",
			 80, 0, 80 );

  TProfile2D * modnpxvsxmym = new
    TProfile2D( "modnpxvsxmym",
		"MOD cluster size vs xmod ymod;x track mod 300 [#mum];y track mod 200 [#mum];MOD <cluster size> [pixels]",
		120, 0, 300, 80, 0, 200, 0, 20 );

  TProfile2D * modqxvsxmym = new
    TProfile2D( "modqxvsxmym",
	      "MOD cluster charge vs xmod ymod;x track mod 300 [#mum];y track mod 200 [#mum];MOD <cluster charge> [ke]",
		120, 0, 300, 80, 0, 200, 0, 0.1 );

  TH1I modlkxBHisto = TH1I( "modlkxb",
			    "linked driplet at MOD x;driplet x at MOD [mm];linked driplets",
			    216, -32.4, 32.4 );
  TH1I modlkyBHisto = TH1I( "modlkyb",
			    "linked driplet at MOD y;driplet y at MOD [mm];linked driplets",
			    160, -8, 8 );
  TH1I modlkxHisto = TH1I( "modlkx",
			   "linked driplet at MOD x;driplet x at MOD [mm];linked driplets",
			   216, -32.4, 32.4 );
  TH1I modlkyHisto = TH1I( "modlky",
			   "linked driplet at MOD y;driplet y at MOD [mm];linked driplets",
			   160, -8, 8 );

  TH1I modlkcolHisto = TH1I( "modlkcol",
			     "MOD linked col;MOD linked col;linked MOD cluster",
			     216, 0, 432 );
  TH1I modlkrowHisto = TH1I( "modlkrow",
			     "MOD linked row;MOD linked row;linked MOD cluster",
			     182, 0, 182 );
  TProfile modlkvst1 =
    TProfile( "modlkvst1",
	      "driplet-MOD links vs time;time [s];driplets with MOD links",
	      300, 0, 300, -0.5, 1.5 );
  TProfile modlkvst2 =
    TProfile( "modlkvst2",
	      "driplet-MOD links vs time;time [s];driplets with MOD links",
	      150, 0, 1500, -0.5, 1.5 );
  TProfile modlkvst3 =
    TProfile( "modlkvst3",
	      "driplet-MOD links vs time;time [s];driplets with MOD links",
	      400, 0, 40000, -0.5, 1.5 );
  TProfile modlkvst6 =
    TProfile( "modlkvst6",
	      "driplet-MOD links vs time;time [h];driplets with MOD links",
	      1000, 0, 50, -0.5, 1.5 );

  TH1I ndrilkHisto = TH1I( "ndrilk", "driplet - MOD links;driplet - MOD links;events",
			    11, -0.5, 10.5 );

  // DUT:
  /*
  TH1I dutphHisto( "dutph",
		   "DUT pixel PH;pixel pulse height [ADC];DUT pixels",
		   200, -100, 900 );
  */
  TH1I dutpxq1stHisto( "dutpxq1st",
		       "DUT pixel charge 1st;1st pixel charge [ke];1st pixels",
		       100, 0, 25 );

  TH1I dutpxq2ndHisto( "dutpxq2nd",
		       "DUT pixel charge 2nd;2nd pixel charge [ke];2nd pixels",
		       100, 0, 25 );

  TH1I dutphHisto( "dutph", "DUT PH;ADC-PED [ADC];pixels", 500, -100, 900 );
  TH1I dutdphHisto( "dutdph", "DUT #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );
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
		   "normal cluster charge;normal cluster charge [ke];clusters",
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
		     120, -6, 6, -0.05, 0.05 );
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
		     240, -12, 12, -0.05, 0.05 );
  TProfile tridyvsty( "tridyvsty",
		      "triplet dy vs slope y;triplet slope y [rad];<triplet #Deltay> [mm]",
		      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile tridyvst3( "tridyvst3",
		      "triplet dy vs time;time [s];<triplet #Deltay> [mm]",
		      300, 0, 6000, -0.05, 0.05 );
  TProfile tridyvst6( "tridyvst6",
		      "triplet dy vs time;time [h];<triplet #Deltay> [mm]",
		      1000, 0, 50, -0.05, 0.05 );

  TH1I trix0Histo( "trix0", "triplets x at scint;x at scint [mm];triplets",
			  240, -12, 12 );
  TH1I triy0Histo( "triy0", "triplets y at scint;y at scint [mm];triplets",
			  120, -6, 6 );
  TH2I * trixy0Histo = new
    TH2I( "trixy0", "triplets x-y at scint;x at scint [mm];y at scint [mm];triplets",
	  240, -12, 12, 120, -6, 6 );

  TH1I tritxHisto( "tritx", "triplet slope x;slope x [rad];triplets",
			    100, -0.005, 0.005 );
  TH1I trityHisto( "trity", "triplet slope y;slope y [rad];triplets",
			    100, -0.005, 0.005 );

  TH1I ntriHisto( "ntri", "triplets;triplets;events", 51, -0.5, 50.5 );

  TH1I trixAHisto( "trixA", "triplets x at DUT;x at DUT [mm];triplets",
			  240, -12, 12 );
  TH1I triyAHisto( "triyA", "triplets y at DUT;y at DUT [mm];triplets",
			  120, -6, 6 );
  TH2I * trixyAHisto = new
    TH2I( "trixyA", "triplets x-y at DUT;x at DUT [mm];y at DUT [mm];triplets",
	  240, -12, 12, 120, -6, 6 );

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

  TH1I trixcHisto( "trixc", "triplets x in DUT;x in DUT [mm];triplets",
			  240, -12, 12 );
  TH1I triycHisto( "triyc", "triplets y in DUT;y in DUT [mm];triplets",
			  120, -6, 6 );
  TH2I * trixycHisto = new
    TH2I( "trixyc", "triplets x-y in DUT;x in DUT [mm];y in DUT [mm];triplets",
	  240, -12, 12, 120, -6, 6 );

  TProfile tridzcvsy( "tridzcvsy",
		   "z vs y;y track [mm];<z of intersect at DUT> [mm]",
		   60, -6, 6, -10, 102 );

  TH1I triz3Histo( "triz3",
		       "z3 should be zero;z3 [mm];triplets",
		       100, -0.01, 0.01 );

  // dripets - triplets:

  TH1I hsixdx = TH1I( "sixdx", "six dx;dx [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );
  TH1I hsixdy = TH1I( "sixdy", "six dy;dy [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );
  TH1I hsixdxc = TH1I( "sixdxc", "six dx;dx [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );
  TH1I hsixdxcsi = TH1I( "sixdxcsi", "six dx Si;#Deltax [mm];triplet-driplet pairs in Si", 200, -0.2*f, 0.2*f );
  TH1I hsixdxccu = TH1I( "sixdxccu", "six dx Cu;#Deltax [mm];triplet-driplet pairs in Cu", 200, -0.2*f, 0.2*f );
  TH1I hsixdxcsid = TH1I( "sixdxcsid", "six dx Si;#Deltax [mm];triplet-driplet pairs in Si", 200, -0.2*f, 0.2*f );

  TH1I hsixdyc = TH1I( "sixdyc", "six dy;dy [mm];triplet-driplet pairs", 200, -0.2*f, 0.2*f );
  TH1I hsixdycsi = TH1I( "sixdycsi", "six dy Si;#Deltay [mm];triplet-driplet pairs in Si", 200, -0.2*f, 0.2*f );
  TH1I hsixdyccu = TH1I( "sixdyccu", "six dy Cu;#Deltay [mm];triplet-driplet pairs in Cu", 200, -0.2*f, 0.2*f );

  TProfile sixdxvsx =
    TProfile( "sixdxvsx",
	      "six #Deltax vs x;xB [mm];<driplet - triplet #Deltax [mm]",
	      240, -12, 12, -0.1, 0.1 );
  TProfile sixmadxvsx =
    TProfile( "sixmadxvsx",
	      "six MAD x vs x;xB [mm];driplet - triplet MAD #Deltax [mm]",
	      240, -12, 12, 0, 0.1 );
  TProfile sixmadxvsy =
    TProfile( "sixmadxvsy",
	      "six MAD x vs y;yB [mm];driplet - triplet MAD #Deltax [mm]",
	      120, -6, 6, 0, 0.1 );
  TProfile sixmadxvstx =
    TProfile( "sixmadxvstx",
	      "six MAD x vs x;triplet #theta_{x} [rad];driplet - triplet MAD #Deltax [mm]",
	      80, -0.002, 0.002, 0, 0.1 );
  TProfile sixmadxvsdtx =
    TProfile( "sixmadxvsdtx",
	      "six MAD x vs x;driplet-triplet #Delta#theta_{x} [rad];driplet - triplet MAD #Deltax [mm]",
	      80, -0.002, 0.002, 0, 0.1 );
  TProfile sixdxvsy =
    TProfile( "sixdxvsy",
	      "six #Deltax vs y;yB [mm];<driplet - triplet #Deltax> [mm]",
	      100, -5, 5, -0.5, 0.5 );
  TProfile sixdxvstx =
    TProfile( "sixdxvstx",
	      "six #Deltax vs slope x;slope x [rad];<driplet - triplet #Deltax> [mm]",
	      100, -0.002, 0.002, -0.5, 0.5 );
  TProfile sixdxvsdtx =
    TProfile( "sixdxvsdtx",
	      "six #Deltax vs #Delta slope x;#Delta slope x [rad];<driplet - triplet #Deltax> [mm]",
	      100, -0.002, 0.002, -0.5, 0.5 );
  TProfile sixdxvst3 =
    TProfile( "sixdxvst3",
	      "sixplet dx vs time;time [s];<sixplet #Deltax> [mm]",
	      300, 0, 6000, -0.05, 0.05 );
  TProfile sixdxvst6 =
    TProfile( "sixdxvst6",
	      "sixplet dx vs time;time [h];<sixplet #Deltax> [mm]",
	      1000, 0, 50, -0.05, 0.05 );

  TProfile sixdyvsx =
    TProfile( "sixdyvsx",
	      "six #Deltay vs x;xB [mm];<driplet - triplet #Deltay> [mm]",
	      200, -10, 10, -0.5, 0.5 );
  TProfile sixdyvsy =
    TProfile( "sixdyvsy",
	      "six #Deltay vs y;yB [mm];<driplet - triplet #Deltay [mm]",
	      120, -6, 6, -0.1, 0.1 );
  TProfile sixdyvsty =
    TProfile( "sixdyvsty",
	      "six #Deltay vs slope y;slope y [rad];<driplet - triplet #Deltay> [mm]",
	      100, -0.002, 0.002, -0.5, 0.5 );
  TProfile sixdyvsdty =
    TProfile( "sixdyvsdty",
	      "six #Deltay vs #Delta slope y;#Delta slope y [rad];<driplet - triplet #Deltay> [mm]",
	      100, -0.002, 0.002, -0.5, 0.5 );
  TProfile sixdyvst3 =
    TProfile( "sixdyvst3",
	      "sixplet dy vs time;time [s];<sixplet #Deltay> [mm]",
	      300, 0, 6000, -0.05, 0.05 );
  TProfile sixdyvst6 =
    TProfile( "sixdyvst6",
	      "sixplet dy vs time;time [h];<sixplet #Deltay> [mm]",
	      1000, 0, 50, -0.05, 0.05 );
  TProfile sixmadyvsx =
    TProfile( "sixmadyvsx",
	      "six MAD y vs x;xB [mm];driplet - triplet MAD #Deltay [mm]",
	      240, -12, 12, 0, 0.1 );
  TProfile sixmadyvsy =
    TProfile( "sixmadyvsy",
	      "six MAD y vs y;yB [mm];driplet - triplet MAD #Deltay [mm]",
	      120, -6, 6, 0, 0.1 );
  TProfile sixmadyvsty =
    TProfile( "sixmadyvsty",
	      "six MAD y vs #theta_{y};triplet #theta_{y} [rad];driplet - triplet MAD #Deltay [mm]",
	      80, -0.002, 0.002, 0, 0.1 );
  TProfile sixmadyvsdty =
    TProfile( "sixmadyvsdty",
	      "six MAD y vs #Delta#theta_{y};driplet-triplet #Delta#theta_{y} [rad];driplet - triplet MAD #Deltay [mm]",
	      80, -0.002, 0.002, 0, 0.1 );

  TProfile2D * sixdxyvsxy = new
    TProfile2D( "sixdxyvsxy",
		"driplet - triplet #Delta_{xy} vs x-y;x_{mid} [mm];y_{mid} [mm];<sqrt(#Deltax^{2}+#Deltay^{2})> [rad]",
		120, -12, 12, 60, -6, 6, 0, 0.7 );

  TH1I hsixdtx =
    TH1I( "sixdtx",
	  "driplet slope x - triplet slope x;driplet slope x - triplet slope x;driplet-triplet pairs",
	  100, -0.005*f, 0.005*f );     
  TH1I hsixdty =
    TH1I( "sixdty",
	  "driplet slope y - triplet slope y;driplet slope y - triplet slope y;driplet-triplet pairs",
	  100, -0.005*f, 0.005*f );     
  TH1I hsixdtxsi =
    TH1I( "sixdtxsi",
	  "driplet triplet #Delta#theta_{x} Si;driplet - triplet #Delta#theta_{x} [rad];driplet-triplet pairs in Si",
	  100, -0.0025*f, 0.0025*f );
  TH1I hsixdtxcu =
    TH1I( "sixdtxcu",
	  "driplet triplet #Delta#theta_{x} Cu;driplet - triplet #Delta#theta_{x} [rad];driplet-triplet pairs in Cu",
	  100, -0.010*f, 0.010*f );

  TH1I hsixdtysi =
    TH1I( "sixdtysi",
	  "driplet triplet #Delta#theta_{y} Si;driplet - triplet #Delta#theta_{y} [rad];driplet-triplet pairs in Si",
	  100, -0.0025*f, 0.0025*f );
  TH1I hsixdtycu =
    TH1I( "sixdtycu",
	  "driplet triplet #Delta#theta_{y} Cu;driplet - triplet #Delta#theta_{y} [rad];driplet-triplet pairs in Cu",
	  100, -0.010*f, 0.010*f );

  TProfile sixdtvsxA =
    TProfile( "sixdtvsxA",
	      "driplet - triplet kink_{xy} vs x;x_{mid} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      240, -12, 12, 0, 0.1 );
  TProfile sixdtvsyA =
    TProfile( "sixdtvsyA",
	      "driplet - triplet kink_{xy} vs y;y_{mid} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      120, -6, 6, 0, 0.1 );
  TProfile2D * sixdtvsxyA = new
    TProfile2D( "sixdtvsxyA",
		"driplet - triplet kink_{xy} vs x-y;x_{mid} [mm];y_{mid} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		120, -12, 12, 60, -6, 6, 0, 0.1 );

  TProfile sixdtvsx =
    TProfile( "sixdtvsx",
	      "driplet - triplet kink_{xy} vs x;x_{DUT} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      240, -12, 12, 0, 0.1 );
  TProfile2D * sixdtvsxy = new
    TProfile2D( "sixdtvsxy",
		"driplet - triplet kink_{xy} vs x-y;x_{DUT} [mm];y_{DUT} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		120, -12, 12, 60, -6, 6, 0, 0.1 );

  TProfile sixdtvsxm =
    TProfile( "sixdtvsxm",
	      "driplet - triplet kink_{xy} vs xmod;track x mod 0.3 [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      60, 0, 0.3, 0, 0.1 );
  TProfile sixdtvsym =
    TProfile( "sixdtvsym",
	      "driplet - triplet kink_{xy} vs ymod;track y mod 0.2 [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      40, 0, 0.2, 0, 0.1 );
  TProfile2D * sixdtvsxmym = new
    TProfile2D( "sixdtvsxmym",
	      "driplet - triplet kink_{xy} vs xmod ymod;track x mod 300 [#mum];track y mod 200 [#mum];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      120, 0, 300, 80, 0, 200, 0, 0.1 );

  TH2I * sixxyHisto = new
    TH2I( "sixxy", "sixplets at z DUT;x [mm];y [mm];sixplets",
	  240, -12, 12, 120, -6, 6 );

  // DUT pixel vs triplets:

  TH1I cmssxaHisto( "cmssxa",
			   "DUT + Telescope x;cluster + triplet #Sigmax [mm];clusters",
			   440, -11, 11 );
  TH1I cmsdxaHisto( "cmsdxa",
			   "DUT - Telescope x;cluster - triplet #Deltax [mm];clusters",
			   440, -11, 11 );

  TH1I cmssyaHisto( "cmssya",
			   "DUT + Telescope y;cluster + triplet #Sigmay [mm];clusters",
			   440, -11, 11 ); // shallow needs wide range
  TH1I cmsdyaHisto( "cmsdya",
			   "DUT - Telescope y;cluster - triplet #Deltay [mm];clusters",
			   880, -22, 22 );

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
			   200, -0.25, 0.25 );

  TProfile cmsdxvsev =
    TProfile( "cmsdxvsev",
	      "DUT #Deltax vs time;time [events];<#Deltax> [mm]",
	      400, 0, 4000*1000, -0.2, 0.2 );

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
		     50, -3.75, 3.75, -0.2, 0.2 );
  TProfile cmsdxvsy( "cmsdxvsy",
		     "#Deltax vs y;y track [mm];<cluster - triplet #Deltax> [mm]",
		     76, -3.8, 3.8, -0.2, 0.2 );
  TProfile cmsdxvstx( "cmsdxvstx",
		      "#Deltax vs #theta_{x};x track slope [rad];<cluster - triplet #Deltax> [mm]",
		      80, -0.004, 0.004, -0.2, 0.2 );

  TH1I cmsdycHisto( "cmsdyc",
		    "#Deltay cut x;cluster - triplet #Deltay [mm];clusters",
		     200, -0.25, 0.25 );
  TH1I cmsdyciHisto( "cmsdyci",
		     "#Deltay cut x, isolated;cluster - triplet #Deltay [mm];isolated clusters",
		     200, -0.25, 0.25 );

  TH1I cmsdy8cHisto( "cmsdy8c",
		     "#Deltay cut x;cluster - sixplet #Deltay [mm];clusters",
		     200, -0.25, 0.25 );
  TH1I cmsdy8ciHisto( "cmsdy8ci",
		      "#Deltay cut x, isolated;cluster - sixplet #Deltay [mm];isolated clusters",
		      200, -0.25, 0.25 );

  TH1I cmsdyc3Histo( "cmsdyc3",
		     "#Deltay cut x, npx < 4;cluster - triplet #Deltay [mm];clusters, npx < 4",
		     200, -0.25, 0.25 );
  TH1I cmsdyc3iHisto( "cmsdyc3i",
		      "#Deltay cut x, npx < 4, isolated;cluster - triplet #Deltay [mm];isolated clusters, npx < 4",
		      200, -0.25, 0.25 );
  TH1I cmsdy8c3Histo( "cmsdy8c3",
		      "#Deltay cut x, npx < 4;cluster - sixplet #Deltay [mm];clusters, npx < 4",
		      200, -0.25, 0.25 );
  TH1I cmsdy8c3iHisto( "cmsdy8c3i",
		       "#Deltay cut x, npx < 4, isolated;cluster - sixplet #Deltay [mm];isolated clusters, npx < 4",
		       200, -0.25, 0.25 );

  TProfile cmsdyvsx( "cmsdyvsx",
		     "DUT #Deltay vs x;x track [mm];<cluster - triplet #Deltay> [mm]",
		     50, -3.75, 3.75, -0.2, 0.2 );
  TProfile cmsdyvsy( "cmsdyvsy",
		     "DUT #Deltay vs y;y track [mm];<cluster - triplet #Deltay> [mm]",
		     76, -3.8, 3.8, -0.2, 0.2 );
  TProfile cmsdyvsty( "cmsdyvsty",
		      "DUT #Deltay vs #theta_{y};y track slope [mrad];<cluster - triplet #Deltay> [mm]",
		      80, -2, 2, -0.2, 0.2 );

  TProfile cmsdyvsev =
    TProfile( "cmsdyvsev",
	      "DUT #Deltay vs time;time [events];<#Deltay> [mm]",
	      400, 0, 4000*1000, -0.2, 0.2 );
  TProfile cmsmadyvsev =
    TProfile( "cmsmadyvsev",
	      "DUT MAD(#Deltay) vs time;time [events];MAD(#Deltay) [mm]",
	      400, 0, 4000*1000, 0, 0.1 );

  TProfile cmsmadyvsx =
    TProfile( "cmsmadyvsx",
	      "DUT MAD(#Deltay) vs track x;track x [mm];MAD(#Deltay) [mm]",
	      80, -4, 4, 0, 0.1 );
  TProfile cmsmady8vsx =
    TProfile( "cmsmady8vsx",
	      "DUT MAD(#Deltay) vs six track x;track x [mm];MAD(#Deltay) six [mm]",
	      80, -4, 4, 0, 0.1 );

  TProfile cmsmadyvsy =
    TProfile( "cmsmadyvsy",
	      "DUT MAD(#Deltay) vs track y;track y [mm];MAD(#Deltay) [mm]",
	      80, -4, 4, 0, 0.1 );

  TProfile cmsmadyvsty =
    TProfile( "cmsmadyvsty",
	      "DUT MAD(#Deltay) vs track angle;y track angle [mrad];MAD(#Deltay) [mm]",
	      80, -2, 2, 0, 0.1 );

  TProfile cmsmadyvsq =
    TProfile( "cmsmadyvsq",
	      "DUT MAD(#Deltay) vs Q0;normal cluster charge [ke];MAD(#Deltay) [mm]",
	      100, 0, 100, 0, 0.1 );

  TProfile cmsmadyvsxm =
    TProfile( "cmsmadyvsxm",
	      "DUT y resolution vs xmod;x track mod 100 [#mum];MAD(#Deltay) [mm]",
	      100, 0, 100, 0, 0.2 );
  TProfile cmsmadyvsym =
    TProfile( "cmsmadyvsym",
	      "DUT y resolution vs ymod;y track mod 100 [#mum];MAD(#Deltay) [mm]",
	      100, 0, 100, 0, 0.2 );
  TProfile2D * cmsmadyvsxmym = new
    TProfile2D( "cmsmadyvsxmym",
		"DUT y resolution vs xmd, ymod;x track mod 100 [#mum];y track mod 100 [#mum];MAD(#Deltay) [mm]",
		50, 0, 100, 50, 0, 100, 0, 0.2 );

  TH1I cmsdycqHisto( "cmsdycq",
		    "#Deltay Landau peak;cluster - triplet #Deltay [mm];Landau peak clusters",
		     200, -0.25, 0.25 );
  TH1I cmsdycqiHisto( "cmsdycqi",
		      "#Deltay Landau peak, isolated;cluster - triplet #Deltay [mm];isolated Landau peak clusters",
		      200, -0.25, 0.25 );
  TH1I cmsdy8cqHisto( "cmsdy8cq",
		      "#Deltay Landau peak;cluster - sixplet #Deltay [mm];Landau peak clusters",
		      200, -0.25, 0.25 );
  TH1I cmsdy8cqiHisto( "cmsdy8cqi",
		       "#Deltay Landau peak, isolated;cluster - sixplet #Deltay [mm];isolated Landau peak clusters",
		       200, -0.25, 0.25 );
  TH1I cmsdycq2Histo( "cmsdycq2",
		      "#Deltay Landau peak;cluster - triplet #Deltay [mm];Landau peak clusters",
		      200, -0.25, 0.25 );

  TH1I trixAlkHisto( "trixAlk",
			   "linked triplet at DUT x;triplet x at DUT [mm];linked triplets",
			   240, -12, 12 );
  TH1I triyAlkHisto( "triyAlk",
			   "linked triplet at DUT y;triplet y at DUT [mm];linked triplets",
			   120, -6, 6 );
  TH2I * trixyAlkHisto = new
    TH2I( "trixyAlk", "linked triplets x-y at DUT;x at DUT [mm];y at DUT [mm];linked triplets",
	  240, -12, 12, 120, -6, 6 );

  TH1I trixclkHisto( "trixclk",
			   "linked triplet in DUT x;triplet x in DUT [mm];linked triplets",
			   240, -12, 12 );
  TH1I triyclkHisto( "triyclk",
			   "linked triplet in DUT y;triplet y in DUT [mm];linked triplets",
			   120, -6, 6 );
  TH2I * trixyclkHisto = new
    TH2I( "trixyclk", "linked triplets x-y in DUT;x in DUT [mm];y at DUT [mm];linked triplets",
	  240, -12, 12, 120, -6, 6 );

  TH1I cmsadcHisto( "cmsadc",
		    "DUT pixel ADC;pixel pulse height [ADC];linked pixels",
		    400, 0, 800 );
  TH1I cmscolHisto( "cmscol",
		    "DUT linked columns;DUT linked cluster column;linked pixels",
		    nx[iDUT], -0.5, nx[iDUT]-0.5 );
  TH1I cmsrowHisto( "cmsrow",
		    "DUT linked rows;DUT linked cluster row;linked pixels",
		    ny[iDUT], -0.5, ny[iDUT]-0.5 );
  TH1I cmsnpxHisto( "cmsnpx",
		     "linked DUT cluster size;cluster size [pixels];linked clusters",
		     80, 0.5, 80.5 );
  TH1I cmsncolHisto( "cmsncol",
		     "linked DUT cluster size;cluster size [columns];linked clusters",
		     52, 0.5, 52.5 );
  TH1I cmsnrowHisto( "cmsnrow",
		     "linked DUT cluster size;cluster size [rows];linked clusters",
		     80, 0.5, 80.5 );

  TProfile cmsncolvsxm( "cmsncolvsxm",
			"DUT cluster size vs xmod;x track mod 100 [#mum];<cluster size> [columns]",
			100, 0, 100, 0, 20 );
  TProfile cmsnrowvsxm( "cmsnrowvsxm",
			"DUT cluster size vs xmod;x track mod 100 [#mum];<cluster size> [rows]",
			100, 0, 100, 0, 80 );

  TProfile cmsncolvsym( "cmsncolvsym",
			"DUT cluster size vs ymod;y track mod 100 [#mum];<cluster size> [columns]",
			100, 0, 100, 0, 20 );
  TProfile cmsnrowvsym( "cmsnrowvsym",
			"DUT cluster size vs ymod;y track mod 100 [#mum];<cluster size> [rows]",
			100, 0, 100, 0, 80 );
  TProfile2D * cmsnpxvsxmym = new
    TProfile2D( "cmsnpxvsxmym",
		"DUT cluster size vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];<cluster size> [pixels]",
		50, 0, 100, 50, 0, 100, 0, 20 );
  TProfile2D * cmsnpxvsxmym2 = new
    TProfile2D( "cmsnpxvsxmym2",
		"DUT cluster size vs xmod ymod;x track mod 200 [#mum];y track mod 200 [#mum];<cluster size> [pixels]",
		100, 0, 200, 100, 0, 200, 0, 20 );

  TH1I cmsq0aHisto( "cmsq0a",
		    "normal cluster charge;normal cluster charge [ke];linked clusters",
		    160, 0, 80 );

  TH2I * cmsxmymq0 = new
    TH2I( "cmsxmymq0", "small cluster charge map;x track mod 100 [#mum];y track mod 100 [#mum];small charge tracks",
	  40, 0, 100, 40, 0, 100 );
  TH2I * cmsxyq0 = new
    TH2I( "cmsxyq0", "small cluster charge map;x track [mm];y track [mm];small charge tracks",
	  90, -4.5, 4.5, 90, -4.5, 4.5 );
  TH2I * cmsxmymq1 = new
    TH2I( "cmsxmymq1", "normal cluster charge map;x track mod 100 [#mum];y track mod 100 [#mum];normal charge tracks",
	  40, 0, 100, 40, 0, 100 );

  TH1I cmsph0Histo( "cmsph0",
		   "normal cluster PH;normal cluster pulse height [ADC];isolated linked clusters",
		   200, 0, 1000 );
  TH1I cmsq0Histo( "cmsq0",
		   "normal cluster charge;normal cluster charge [ke];isolated linked clusters",
		   160, 0, 80 );
  TH1I cmsq00Histo( "cmsq00",
		   "normal cluster charge;normal cluster charge [ke];isolated linked clusters",
		   160, 0, 40 );
  TH1I cmsq03Histo( "cmsq03",
		   "normal cluster charge, < 4 px;normal cluster charge [ke];linked clusters, < 4 px",
		    160, 0, 80 );
  TH1I cmsq0dHisto( "cmsq0d",
		    "cluster charge bias dot;cluster charge bias dot [ke];linked clusters",
		    160, 0, 80 );
  TH1I cmsq0nHisto( "cmsq0n",
		    "cluster charge no dot;cluster charge no dot [ke];linked clusters",
		    160, 0, 80 );

  double qxmax = 0.001; // with qwid = 1.1: 7.5 ke cutoff
  if( chip0 >= 300 ) qxmax = 0.007; // 3D with qwid = 1.6: 8 ke cutoff

  TProfile cmsqxvsx( "cmsqxvsx",
		     "DUT cluster charge vs x;x track [mm];<cluster charge> [ke]",
		     76, -3.8, 3.8, 0, qxmax ); // cutoff at 5 ke
  TProfile cmsqxvsy( "cmsqxvsy",
		     "DUT cluster charge vs y;y track [mm];<cluster charge> [ke]",
		     76, -3.8, 3.8, 0, qxmax );
  TProfile2D * cmsqxvsxy = new
		      TProfile2D( "cmsqxvsxy",
				  "DUT cluster charge vs xy;x track [mm];y track [mm];<cluster charge> [ke]",
				  76, -3.8, 3.8, 76, -3.8, 3.8, 0, qxmax );
  TProfile cmsqxvsxm( "cmsqxvsxm",
		      "DUT cluster charge vs xmod;x track mod 100 [#mum];<cluster charge> [ke]",
		      100, 0, 100, 0, qxmax );
  TProfile cmsqxvsym( "cmsqxvsym",
		      "DUT cluster charge vs ymod;y track mod 100 [#mum];<cluster charge> [ke]",
		      100, 0, 100, 0, qxmax );
  TProfile2D * cmsqxvsxmym = new
    TProfile2D( "cmsqxvsxmym",
	      "DUT cluster charge vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];<cluster charge> [ke]",
		50, 0, 100, 50, 0, 100, 0, qxmax );
  TProfile2D * cmsqxvsxmym2 = new
    TProfile2D( "cmsqxvsxmym2",
	      "DUT cluster charge vs xmod ymod;x track mod 200 [#mum];y track mod 200 [#mum];<cluster charge> [ke]",
		100, 0, 200, 100, 0, 200, 0, qxmax );

  TH1I cmspxqHisto( "cmspxq",
		    "DUT pixel charge linked;pixel charge [ke];linked pixels",
		    100, 0, 25 );

  TH1I cmsdminHisto( "cmsdmin",
		     "cluster - Telescope min dxy;cluster - triplet min #Deltaxy [mm];tracks",
		     200, 0, 1 );
  TH1I cmsdxminHisto( "cmsdxmin",
		     "cluster - Telescope min dx;cluster - triplet min #Deltax [mm];tracks",
		     200, -1, 1 );
  TH1I cmsdyminHisto( "cmsdymin",
		     "cluster - Telescope min dy;cluster - triplet min #Deltay [mm];tracks",
		     200, -1, 1 );
  TH1I cmspdxminHisto( "cmspdxmin",
		       "pixel - Telescope min dx;pixel - triplet min #Deltax [mm];tracks",
		       200, -1, 1 );
  TH1I cmspdyminHisto( "cmspdymin",
		       "pixel - Telescope min dy;pixel - triplet min #Deltay [mm];tracks",
		       200, -1, 1 );
  TH1I cmspdminHisto( "cmspdmin",
		      "pixel - Telescope min dxy;pixel - triplet min #Deltaxy [mm];tracks",
		      200, 0, 1 );

  TH2I * sixxylkHisto = new
    TH2I( "sixxylk",
	  "MOD-linked sixplets at z DUT;x_{mid} [mm];y_{mid} [mm];MOD-linked sixplets",
	  240, -12, 12, 120, -6, 6 );
  TH2I * sixxyeffHisto = new
    TH2I( "sixxyeff",
	  "MOD-linked sixplets with DUT;x_{mid} [mm];y_{mid} [mm];DUT-MOD-linked sixplets",
	  240, -12, 12, 120, -6, 6 );

  TProfile effvsdxy =
    TProfile( "effvsdxy",
	      "DUT efficiency vs dxy;xy search radius [mm];efficiency",
	      100, 0, 1, -1, 2 );

  TH1I effdminHisto = TH1I( "effdmin",
			 "min DUT - triplet distance;min DUT - triplet #Delta_{xy} [mm];MOD linked fiducial iso triplets",
			 200, 0, 20 );
  TH1I effdmin0Histo =
    TH1I( "effdmin0",
	  "min DUT - triplet distance;min DUT - triplet #Delta_{xy} [mm];inefficient triplets",
	  200, 0, 20 );
  TH1I effrxmin0Histo =
    TH1I( "effrxmin0",
	  "min DUT - triplet distance x;min DUT - triplet #Delta_{x}/#Delta_{xy};inefficient triplets",
	  200, -1, 1 );
  TH1I effrymin0Histo =
    TH1I( "effrymin0",
	  "min DUT - triplet distance y;min DUT - triplet #Delta_{y}/#Delta_{xy};inefficient triplets",
	  200, -1, 1 );
  TH1I effdxmin0Histo =
    TH1I( "effdxmin0",
	  "min DUT - triplet distance x;min DUT - triplet #Delta_{x} [mm];inefficient triplets",
	  120, -9, 9 );
  TH1I effdymin0Histo =
    TH1I( "effdymin0",
	  "min DUT - triplet distance y;min DUT - triplet #Delta_{y} [mm];inefficient triplets",
	  180, -9, 9 );

  TH1I effclq0Histo =
    TH1I( "effclq0",
	  "nearest cluster charge;DUT cluster charge [ke];nearest clusters",
	  160, 0, 80 );
  TH1I effclq1Histo =
    TH1I( "effclq1",
	  "nearest cluster charge dxy > 0.1;DUT cluster charge [ke];dxy > 0.1 nearest cluster",
	  160, 0, 80 );
  TH1I effclq2Histo =
    TH1I( "effclq2",
	  "nearest cluster charge dxy > 0.2;DUT cluster charge [ke];dxy > 0.2 nearest cluster",
	  160, 0, 80 );
  TH1I effclq3Histo =
    TH1I( "effclq3",
	  "nearest cluster charge dxy > 0.3;DUT cluster charge [ke];dxy > 0.3 nearest cluster",
	  160, 0, 80 );
  TH1I effclq4Histo =
    TH1I( "effclq4",
	  "nearest cluster charge dxy > 0.4;DUT cluster charge [ke];dxy > 0.4 nearest cluster",
	  160, 0, 80 );
  TH1I effclq5Histo =
    TH1I( "effclq5",
	  "nearest cluster charge dxy > 0.5;DUT cluster charge [ke];dxy > 0.5 nearest cluster",
	  160, 0, 80 );
  TH1I effclq6Histo =
    TH1I( "effclq6",
	  "nearest cluster charge dxy > 0.6;DUT cluster charge [ke];dxy > 0.6 nearest cluster",
	  160, 0, 80 );
  TH1I effclq7Histo =
    TH1I( "effclq7",
	  "nearest cluster charge dxy > 0.7;DUT cluster charge [ke];dxy > 0.7 mm nearest cluster",
	  160, 0, 80 );
  TH1I effclq8Histo =
    TH1I( "effclq8",
	  "nearest cluster charge dxy > 0.8;DUT cluster charge [ke];dxy > 0.8 mm nearest cluster",
	  160, 0, 80 );
  TH1I effclq9Histo =
    TH1I( "effclq9",
	  "nearest cluster charge dxy > 0.9;DUT cluster charge [ke];dxy > 0.9 mm nearest cluster",
	  160, 0, 80 );

  TH1I effclqrHisto =
    TH1I( "effclqr",
	  "nearest cluster charge, 1 driplet-MOD;DUT cluster charge [ke];mono nearest clusters",
	  160, 0, 80 );

  TProfile2D * effnpxvsxmym = new
    TProfile2D( "effnpxvsxmym",
		"DUT cluster size vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];<cluster size> [pixels]",
		50, 0, 100, 50, 0, 100, 0, 20 );
  TProfile2D * effq0vsxmym = new
    TProfile2D( "effq0vsxmym",
	      "DUT cluster charge vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];<cluster charge> [ke]",
		50, 0, 100, 50, 0, 100, 0, 50 );
  TH1I effq0dHisto( "effq0d",
		    "cluster charge bias dot;cluster charge bias dot [ke];tracks",
		    160, 0, 80 );
  TH1I effq0nHisto( "effq0n",
		    "cluster charge no dot;cluster charge outside bias dot [ke];tracks",
		    160, 0, 80 );

  TProfile2D * effvsxy = new
    TProfile2D( "effvsxy",
		"DUT efficiency vs x;x track at DUT [mm];y track at DUT [mm];efficiency",
		90, -4.5, 4.5, 90, -4.5, 4.5, -1, 2 ); // bin = pix
  TProfile effvsx =
    TProfile( "effvsx",
	      "DUT efficiency vs x;x track at DUT [mm];efficiency",
	      160, -4.0, 4.0, -1, 2 ); // bin = col
  TProfile effvsy =
    TProfile( "effvsy",
	      "DUT efficiency vs y;y track at DUT [mm];efficiency",
	      160, -4, 4, -1, 2 ); // bin = row

  TProfile effvsev1 =
    TProfile( "effvsev1",
	      "DUT efficiency vs events;events;efficiency",
	      300, 0, 30*1000, -1, 2 );
  TProfile effvsev2 =
    TProfile( "effvsev2",
	      "DUT efficiency vs events;events;efficiency",
	      160, 0, 160*1000, -1, 2 );

  TProfile effvst1 =
    TProfile( "effvst1",
	      "DUT efficiency vs time;time [s];<efficiency> / s",
	      300, 0, 300, -1, 2 );
  TProfile effvst2 =
    TProfile( "effvst2",
	      "DUT efficiency vs time;time [s];<efficiency> / 5 s",
	      200, 0, 1000, -1, 2 );
  TProfile effvst3 =
    TProfile( "effvst3",
	      "DUT efficiency vs time;time [s];<efficiency> / 10 s",
	      600, 0, 6000, -1, 2 );
  TProfile effvst4 =
    TProfile( "effvst4",
	      "DUT efficiency vs time;time [s];<efficiency> / 100 s",
	      600, 0, 60000, -1, 2 );
  TProfile effvst6 =
    TProfile( "effvst6",
	      "DUT efficiency vs time;time [h];efficiency",
	      1000, 0, 50, -1, 2 );

  TProfile effvst1t =
    TProfile( "effvst1t",
	      "DUT efficiency vs time;time [s];<efficiency> / s",
	      300, 0, 300, -1, 2 );
  TProfile effvst2t =
    TProfile( "effvst2t",
	      "DUT efficiency vs time;time [s];<efficiency> / 5 s",
	      200, 0, 1000, -1, 2 );
  TProfile effvst3t =
    TProfile( "effvst3t",
	      "DUT efficiency vs time;time [s];<efficiency> / 10 s",
	      600, 0, 6000, -1, 2 );
  TProfile effvst4t =
    TProfile( "effvst4t",
	      "DUT efficiency vs time;time [s];<efficiency> / 100 s",
	      600, 0, 60000, -1, 2 );
  TProfile effvst6t =
    TProfile( "effvst6t",
	      "DUT efficiency vs time;time [h];efficiency",
	      1000, 0, 50, -1, 2 );

  TProfile2D * effvsxt = new
    TProfile2D( "effvsxt",
	      "DUT efficiency vs time and x;time [s];x [mm];efficiency",
		100, 0, 1000, 50, -3.75, 3.75, -1, 2 );

  TProfile effvsntri =
    TProfile( "effvsntri",
	      "DUT efficiency vs triplets;triplets;efficiency",
	      20, 0.5, 20.5, -1, 2 );
  TProfile effvsndri =
    TProfile( "effvsndri",
	      "DUT efficiency vs driplets;driplets;efficiency",
	      20, 0.5, 20.5, -1, 2 );

  TProfile effvsxm =
    TProfile( "effvsxm",
	      "DUT efficiency vs xmod;x track mod 0.3 [mm];efficiency",
	      50, 0, 0.1, -1, 2 );
  TProfile effvsym =
    TProfile( "effvsym",
	      "DUT efficiency vs ymod;y track mod 0.2 [mm];efficiency",
	      50, 0, 0.1, -1, 2 );
  TProfile2D * effvsxmym = new
    TProfile2D( "effvsxmym",
		"DUT efficiency vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];efficiency",
		50, 0, 100, 50, 0, 100, -1, 2 );
  TProfile2D * effvsxmym2 = new
    TProfile2D( "effvsxmym2",
		"DUT efficiency vs xmod ymod;x track mod 200 [#mum];y track mod 200 [#mum];efficiency",
		510, 0, 200, 100, 0, 200, -1, 2 );

  TProfile effvstx =
    TProfile( "effvstx",
	      "DUT efficiency vs track slope x;x track slope [rad];efficiency",
	      100, -0.005, 0.005, -1, 2 );
  TProfile effvsty =
    TProfile( "effvsty",
	      "DUT efficiency vs track slope y;y track slope [rad];efficiency",
	      100, -0.005, 0.005, -1, 2 );
  TProfile effvstxy =
    TProfile( "effvstxy",
	      "DUT efficiency vs track slope;track slope [rad];efficiency",
	      80, 0, 0.004, -1, 2 );
  TProfile effvsdslp =
    TProfile( "effvsdslp",
	      "DUT efficiency vs kink;driplet - triplet kink angle [rad];efficiency",
	      100, 0, 0.01, -1, 2 );

  TProfile effvstmin =
    TProfile( "effvstmin",
	      "DUT efficiency vs triplet isolation;triplet isolation [mm];efficiency",
	      80, 0, 8, -1, 2 );
  TProfile effvsdmin =
    TProfile( "effvsdmin",
	      "DUT efficiency vs driplet isolation;driplet isolation [mm];efficiency",
	      80, 0, 8, -1, 2 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // data files:

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

  timespec ts;
  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s0 = ts.tv_sec; // seconds since 1.1.1970
  long f0 = ts.tv_nsec; // nanoseconds

  double zeit1 = 0; // eudaq event
  double zeit2 = 0; // eudaq planes
  double zeit3 = 0; // clus
  double zeit5 = 0; // DUT
  double zeit6 = 0; // track

  string E {"E"}; // empty  flag

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  list < vector <cluster> > clist[9];
  list < evInfo > infoList;
  int iev = 0;

  uint64_t tlutime0 = 0;
  uint64_t prevtlutime = 0;

  do {

    clock_gettime( CLOCK_REALTIME, &ts );
    time_t s1 = ts.tv_sec; // seconds since 1.1.1970
    long f1 = ts.tv_nsec; // nanoseconds

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

    if( iev <  0 )
      ldb = 1;

    if( lev < 10 )
      ldb = 1;

    if( ldb ) std::cout << "debug ev " << iev << endl << flush;

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

    evInfo evinf;
    evinf.tlutime = tlutime;
    evinf.dtbtime = tlutime; // init
    evinf.filled = E;
    infoList.push_back(evinf);

    double tludt = (tlutime - prevtlutime) / fTLU; // [s]
    if( iev > 1 && tludt > 1e-6 ) {
      hdttlu.Fill( log(tludt)/log10 );
      hdttluw.Fill( log(tludt)/log10, tludt );
    }
    hdttlu1.Fill( tludt );
    hdttlu2.Fill( tludt );

    prevtlutime = tlutime;

    if( iev < 10 )
      cout << "scopes processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "scopes processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "scopes processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 )
      cout << "scopes processing  " << run << "." << iev << "  taken " << evsec
	   << " dt " << tludt << endl;

    //if( tludt > 0.08 ) cout << "  ev " << iev << " dt " << tludt << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    clock_gettime( CLOCK_REALTIME, &ts );
    time_t s2 = ts.tv_sec; // seconds since 1.1.1970
    long f2 = ts.tv_nsec; // nanoseconds
    zeit1 += s2 - s1 + ( f2 - f1 ) * 1e-9; // read and clus

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      clock_gettime( CLOCK_REALTIME, &ts );
      time_t s1 = ts.tv_sec; // seconds since 1.1.1970
      long f1 = ts.tv_nsec; // nanoseconds

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldb ) std::cout << "PLANE " << plane.ID();

      // /home/pitzl/eudaq/main/include/eudaq/CMSPixelHelper.hh

      int ipl = plane.ID();

      if( run > 28000 && ipl > 0 && ipl < 7 ) // 2017, eudaq 1.6: Mimosa 1..6, DUT 7, REF 8, QAD 9
	ipl -= 1; // 0..5
      if( ipl > 8 ) ipl = 6; // QUAD = MOD

      if( ldb ) cout << " = ipl " << ipl << ", size " << pxl.size() << flush;

      vector <pixel> pb; // for clustering

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

	int xm = ix;
	int ym = iy;

	if( ipl == iMOD ) {

	  // leave space for big pixels:

	  int roc = ix / 52; // 0..7
	  int col = ix % 52; // 0..51
	  int row = iy;
	  int x = 1 + ix + 2*roc; // 1..52 per ROC with big pix
	  int y = iy;
	  if( iy > 79 ) y += 2;

	  ix = x;
	  iy = y;

	  // flip for upper ROCs into local addresses:

	  if( iy > 79 ) {
	    roc = 15 - roc; // 15..8
	    col = 51 - col; // 51..0
	    row = 159 - iy; // 79..0
	  }

	  if( adc > 0 &&
	      roc >= 0 && roc < 16 &&
	      col >= 0 && col < 52 &&
	      row >= 0 && row < 80 ) {

	    double Ared = adc - m4[roc][col][row]; // m4 is asymptotic maximum

	    if( Ared >= 0 )
	      Ared = -0.1; // avoid overflow

	    double a3 = m3[roc][col][row]; // positive
	    if( weib == 3 )
	      q = m1[roc][col][row] *
		( pow( -log( -Ared / a3 ), 1/m2[roc][col][row] ) - m0[roc][col][row] ) * mke;
	    // q = ( (-ln(-(A-m4)/m3))^1/m2 - m0 )*m1

	  } // valid

	} // MOD

	hcol[ipl].Fill( ix+0.5 );
	hrow[ipl].Fill( iy+0.5 );
	hmap[ipl]->Fill( ix+0.5, iy+0.5 );

	// fill pixel block for clustering:

	pixel px;
	px.col = ix; // col
	px.row = iy; // row
	px.adc = adc;
	px.q = q;
	px.ord = pb.size(); // readout order
	px.big = 0;
	pb.push_back(px);

	if( ipl == iMOD ) {

	  // double big pixels:
	  // 0+1
	  // 2..51
	  // 52+53

	  int col = xm % 52; // 0..51

	  if( col == 0 ) {
	    px.col = ix-1; // double
	    px.row = iy;
	    pb[pb.size()-1].adc *= 0.5;
	    pb[pb.size()-1].q *= 0.5;
	    px.adc = 0.5*adc;
	    px.q = 0.5*q;
	    px.big = 1;
	    pb.push_back(px);
	  }

	  if( col == 51 ) {
	    px.col = ix+1; // double
	    px.row = iy;
	    pb[pb.size()-1].adc *= 0.5;
	    pb[pb.size()-1].q *= 0.5;
	    px.adc = 0.5*adc;
	    px.q = 0.5*q;
	    px.big = 1;
	    pb.push_back(px);
	  }

	  if( ym == 79 ) {
	    px.col = ix; // double
	    px.row = 80;
	    pb[pb.size()-1].adc *= 0.5;
	    pb[pb.size()-1].q *= 0.5;
	    px.adc = 0.5*adc;
	    px.q = 0.5*q;
	    px.big = 1;
	    pb.push_back(px);
	  }

	  if( ym == 80 ) {
	    px.col = ix; // double
	    px.row = 81;
	    pb[pb.size()-1].adc *= 0.5;
	    pb[pb.size()-1].q *= 0.5;
	    px.adc = 0.5*adc;
	    px.q = 0.5*q;
	    px.big = 1;
	    pb.push_back(px);
	  }

	} // MOD

	if( pb.size() > 990 ) {
	  cout << "pixel buffer overflow in plane " << ipl
	       << ", event " << iev
	       << endl;
	  break;
	}

      } // pix

      hnpx[ipl].Fill( pb.size() );

      if( ldb ) std::cout << std::endl;

      clock_gettime( CLOCK_REALTIME, &ts );
      time_t s2 = ts.tv_sec; // seconds since 1.1.1970
      long f2 = ts.tv_nsec; // nanoseconds
      zeit2 += s2 - s1 + ( f2 - f1 ) * 1e-9; // read and clus

      // clustering:

      vector < cluster > vcl = getClus(pb);

      if( ldb ) cout << "clusters " << vcl.size() << endl;

      hncl[ipl].Fill( vcl.size() );

      for( vector<cluster>::iterator c = vcl.begin(); c != vcl.end(); ++c ) {

	hsiz[ipl].Fill( c->size );
	hncol[ipl].Fill( c->ncol );
	hnrow[ipl].Fill( c->nrow );

	c->vpix.clear(); // save space

      } // cl

      clist[ipl].push_back(vcl);

      clock_gettime( CLOCK_REALTIME, &ts );
      time_t s3 = ts.tv_sec; // seconds since 1.1.1970
      long f3 = ts.tv_nsec; // nanoseconds
      zeit3 += s3 - s2 + ( f3 - f2 ) * 1e-9; // read and clus

    } // eudaq planes

    ++iev;

  } while( reader->NextEvent() && iev < lev );

  delete reader;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read DUT stream:

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s4 = ts.tv_sec; // seconds since 1.1.1970
  long f4 = ts.tv_nsec; // nanoseconds
  
  string START {"START"};
  string hd;
  while( hd != START ) {
    getline( evFile, hd ); // read one line into string
    cout << "  " << hd << endl;
  }
  bool readnext = 1;

  string F {"F"}; // filled flag
  string A {"A"}; // added  flag

  int nresync = 0;

  uint64_t prevdtbtime = 0;

  bool ldbt = 0;

  iev = 0;

  for( list <evInfo>::iterator evinfo = infoList.begin(); evinfo != infoList.end() && iev < lev; ++evinfo ) {

    string DUTev;
    if( readnext )
      getline( evFile, DUTev );

    if( evFile.eof() ) {
      cout << evFileName << " EOF"
	   << " after event " << iev
	   << endl;
      break; // event loop
    }

    uint64_t tlutime = evinfo->tlutime;
    if( iev == 1  )
      tlutime0 = tlutime;
    double evsec = (tlutime - tlutime0) / fTLU;
    double tludt = (tlutime - prevtlutime) / fTLU; // [s]
    prevtlutime = tlutime;

    ++iev;

    bool ldb = 0;

    if( iev < 1 )
      ldb = 1;

    if( lev < 10 )
      ldb = 1;

    if( ldb ) std::cout << "debug ev " << iev << endl << flush;

    if( ldb ) cout << "  DUT ev " << DUTev << endl << flush;

    istringstream iss( DUTev ); // tokenize string

    int dut_ev;
    iss >> dut_ev; // DUT event

    if( dut_ev%1000 == 0 )
      cout << " " << dut_ev << flush;

    string filled;
    iss >> filled;

    int iblk; // event block number: 100, 200, 300...
    iss >> iblk;

    uint64_t dtbtime = prevdtbtime;

    if( filled == F )
      iss >> dtbtime; // from run 456 = 31093

    else if( filled == E )
      iss >> dtbtime; // from run 456 = 31093

    double dtbdt = ( dtbtime - prevdtbtime ) / fDTB;
    if( dtbdt > 1e-6 )
      hdtdtb.Fill( log(dtbdt)/log10 );

    hddt.Fill( (tludt - dtbdt)*1E3 ); // [ms]
    ddtvsdt.Fill( tludt, (tludt - dtbdt)*1E3 ); // [ms]
    ddtvsev1.Fill( iev, (tludt - dtbdt)*1E3 ); // [ms]
    ddtvsev2.Fill( iev, (tludt - dtbdt)*1E3 ); // [ms]

    if( iev > 65100 )
      ldbt = 0;

    if( ldbt && tludt > 0.5 ) cout << endl;
    if( ldbt )
      cout << "\t" << iev << " TLU " << tludt*1E3
	   << ", DTB " << dtbdt*1e3
	   << endl; // [ms]

    if( run == 31774 && iev > 0 )
      cout << "\t" << iev << " TLU " << tludt*1E3
	   << " - " << dtbdt*1e3 << " DTB"
	   << endl; // [ms]

    // large time gap = DTB readout of one data block

    while( iev > 88 && tludt > 0.5 && dtbdt < 0.2*tludt && run != 31166 ) {

      prevdtbtime = dtbtime;
      ++nresync;
      //if( ldbt )
      cout << "R4S resync " << nresync << endl;
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
	//if( ldbt )
	cout << "  ev " << iev << " TLUdt " << tludt*1E3
	     << ", DTBdt " << dtbdt*1e3
	     << endl; // [ms]
      }
      else { // added event
	if( ldbt )
	  cout << "\t DTB added" << endl;
	break;
      }

      break; // only once

    } // tludt

    evinfo->dtbtime = dtbtime;
    evinfo->filled = filled;

    vector <pixel> pb; // for clustering

    if( readnext && filled == F ) {

      string roi;
      getline( evFile, roi );
      istringstream roiss( roi ); // tokenize string

      int ipx = 0;
      vector <pixel> vpx;
      vpx.reserve(35);

      if( ldb ) cout << "  px";
      /* slow
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

	pixel px { col, row, ph, ph, ipx, 0 };
	vpx.push_back(px);
	++ipx;

      } // roi px
      */
      string BLANK{" "};
      size_t start = 0;
      size_t gap = 0;
      while( gap < roi.size()-1 ) { // data have trailing blank

	gap = roi.find( BLANK, start );
	string s1( roi.substr( start, gap - start ) );
	//cout << " " << s1 << "(" << gap << ")";
	//int col = stoi(s1);
	int col = atoi( s1.c_str() ); // 4% faster
	//int col = fast_atoi( s1.c_str() ); // another 6% faster
	start = gap + BLANK.size();

	gap = roi.find( BLANK, start );
	string s2( roi.substr( start, gap - start ) );
	//cout << " " << s2 << "(" << gap << ")";
	//int row = stoi(s2);
	int row = atoi( s2.c_str() );
	//int row = fast_atoi( s2.c_str() );
	start = gap + BLANK.size();

	gap = roi.find( BLANK, start );
	string s3( roi.substr( start, gap - start ) );
	//cout << " " << s1 << "(" << gap << ")";
	//double ph = stod(s3);
	double ph = atof(s3.c_str());
	start = gap + BLANK.size();

	hcol[iDUT].Fill( col+0.5 );
	hrow[iDUT].Fill( row+0.5 );

	pixel px { col, row, ph, ph, ipx, 0 };
	vpx.push_back(px);
	++ipx;

      } // roi

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

	dutphHisto.Fill( ph4 );
	dutdphHisto.Fill( dph ); // sig 2.7

	// r4scal.C

	double U = ( dph - p3[col4][row4] ) / p2[col4][row4];

	if( U >= 1 )
	  U = 0.9999999; // avoid overflow

	double vcal = p0[col4][row4] - p1[col4][row4] * log( (1-U)/U ); // inverse Fermi

	double q = ke*vcal;

	//double dphcut = 10; // dy8cq 4.80
	double dphcut = 12; // 31619 dy8cq 4.77
	//double dphcut = 15; // dy8cq 4.78

	if( chip0 > 300 ) dphcut = 20; // irradiated 3D is noisy
	if( chip0 == 122 ) dphcut = 24; // irrad, gain_2
	if( chip0 == 123 ) dphcut = 24; // irrad, gain_2
	if( chip0 == 128 ) dphcut = 24; // irrad, gain_2
	if( chip0 == 133 ) dphcut = 24; // irrad, gain_2

	if( dph > dphcut ) {

	//if( q > 0.8 ) { // 31166 cmsdycq 5.85
	//if( q > 0.9 ) { // 31166 cmsdycq 5.74
	//if( q > 1.0 ) { // 31166 cmsdycq 5.72
	//if( q > 1.1 ) { // 31166 cmsdycq 5.75
	//if( q > 1.2 ) { // 31166 cmsdycq 5.81
	//if( q > 1.5 ) { // 31166 cmsdycq 6.00
	//if( q > 2.0 ) { // 31166 cmsdycq 6.42
	//if( q > 2.5 ) { // 31166 cmsdycq 6.86  fittp0 +- 3 sig
	//if( q > 3.0 ) { // 31166 cmsdycq 7.37
	//if( q > 3.5 ) { // 31166 cmsdycq 8.01  eff 99.6
	//if( q > 4.0 ) { // 31166 cmsdycq 8.71  eff 99.4
	//if( q > 4.5 ) { // 31166 cmsdycq 9.61  eff 98.9
	//if( q > 5.0 ) { // 31166 cmsdycq 10.44  eff 98.0
	//if( q > 5.5 ) { // 31166 cmsdycq 11.3  eff 96.3
	//if( q > 6.0 ) { // 31166 cmsdycq 12.05  eff 93.5

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
	  px.adc = dph;
	  px.q = q;
	  //if( chip0 == 122 ) px.q = dph; // irrad, no gaincal
	  px.ord = pb.size(); // readout order, starts at zero
	  px.big = 0;
	  hmap[iDUT]->Fill( px.col+0.5, px.row+0.5 );
	  pb.push_back(px);

	  if( pb.size() > 990 ) {
	    cout << "R4S pixel buffer overflow in event " << iev
		 << endl;
	    break;
	  }

	} // dph cut

      } // roi px

      if( ldb ) cout << " npx " << pb.size() << endl << flush;

    } // filled

    readnext = 1;
    if( iev > 88 && dtbdt > 0.5 && tludt < 0.2*dtbdt && run != 31166 ) {
      readnext = 0;
      if( ldbt ) cout << "repeat DTB event " << DUTev << endl;
    }
    if( readnext )
      prevdtbtime = dtbtime;

    // DUT clustering:

    hnpx[iDUT].Fill( pb.size() );

    vector < cluster > vcl;
    vcl = getClus(pb);
    clist[iDUT].push_back(vcl);

    hncl[iDUT].Fill( vcl.size() );

    for( unsigned icl = 0; icl < vcl.size(); ++ icl ) {

      hsiz[iDUT].Fill( vcl[icl].size );
      //hclph.Fill( vcl[icl].sum );
      //hclmap->Fill( vcl[icl].col, vcl[icl].row );

    } // icl

    dutnpxvst2.Fill( evsec, pb.size() );
    dutnclvst2.Fill( evsec, vcl.size() );

    int DUTyld = 0;
    if( pb.size() ) DUTyld = 1; // no double counting: events with at least one px

    dutyldvst2.Fill( evsec, DUTyld );
    dutyldvst6.Fill( evsec/3600, DUTyld );

    if( ldb ) cout << "  DUT cl " << vcl.size() << endl << flush;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT:

    for( vector<cluster>::iterator c = vcl.begin(); c != vcl.end(); ++c ) {

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
      dutpxq2ndHisto.Fill( q2 ); // no noise ?

    } // DUT cl

    // tsunami corrected clusters:

    for( vector<cluster>::iterator c = vcl.begin(); c != vcl.end(); ++c ) {

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

	qcol[icol] += q; // project cluster onto cols
	qrow[irow] += q; // project cluster onto rows

      } // pix

      int ncol = colmax - colmin + 1;
      int nrow = rowmax - rowmin + 1;

      if( colmin > 0 && colmax < nx[iDUT]-1 && rowmin > 0 && rowmax < ny[iDUT]-1 ) {
	dutq0Histo.Fill( c->charge * norm );
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

  } // event loop

  cout << endl;
  cout << "R4S resyncs " << nresync << endl;

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s5 = ts.tv_sec; // seconds since 1.1.1970
  long f5 = ts.tv_nsec; // nanoseconds
  zeit5 += s5 - s4 + ( f5 - f4 ) * 1e-9; // read and clus

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT and MOD alignment iterations:

  int maxiter = DUTaligniteration + 1;
  if( maxiter < 4 ) maxiter = 4;
  unsigned niter = 0;

  for( int aligniteration = DUTaligniteration; aligniteration < maxiter; ++aligniteration ) {

    // normal vector on MOD surface:
    // N = ( 0, 0, -1 ) on MOD, towards -z
    // transform into tele system:
    // tilt alpha around x
    // turn omega around y

    double com = cos( MODturn*wt );
    double som = sin( MODturn*wt );
    double cam = cos( MODtilt*wt );
    double sam = sin( MODtilt*wt );
    double cfm = cos( MODrot );
    double sfm = sin( MODrot );

    double Nxm =-cam*som;
    double Nym = sam;
    double Nzm =-cam*com;

    double normM = cos( MODturn*wt ) * cos( MODtilt*wt ); // length of Nz
 
    double DUTalignx0 = DUTalignx; // at time 0
    double DUTaligny0 = DUTaligny;

    // normal vector on DUT surface:
    // N = ( 0, 0, -1 ) on DUT, towards -z
    // transform into tele system:
    // tilt alpha around x
    // turn omega around y

    double co = cos( DUTturn*wt );
    double so = sin( DUTturn*wt );
    double ca = cos( DUTtilt*wt );
    double sa = sin( DUTtilt*wt );
    double cf = cos( DUTrot );
    double sf = sin( DUTrot );

    double Nx =-ca*so;
    double Ny = sa;
    double Nz =-ca*co;

    norm = cos( DUTturn*wt ) * cos( DUTtilt*wt ); // length of Nz

    hdx35.Reset();
    hdy35.Reset();
    hdridx.Reset();
    hdridy.Reset();

    hdridxc.Reset();
    dridxvsy.Reset();
    dridxvstx.Reset();

    hdridyc.Reset();
    dridyvsx.Reset();
    dridyvsty.Reset();

    drixHisto.Reset();
    driyHisto.Reset();
    drixyHisto->Reset();
    dritxHisto.Reset();
    drityHisto.Reset();

    ndriHisto.Reset();

    drizixHisto.Reset();
    drizix2Histo.Reset();
    driziyHisto.Reset();
    dddmin1Histo.Reset();
    dddmin2Histo.Reset();

    modsxaHisto.Reset();
    moddxaHisto.Reset();
    modsyaHisto.Reset();
    moddyaHisto.Reset();

    moddxHisto.Reset();
    moddyHisto.Reset();

    moddycHisto.Reset();
    moddycqHisto.Reset();
    moddyvsx.Reset();
    moddyvsy.Reset();
    moddyvsty.Reset();

    moddxcHisto.Reset();
    moddxcqHisto.Reset();
    moddxvsx.Reset();
    moddxvsy.Reset();
    moddxvstx.Reset();

    modnpxHisto.Reset();
    modqHisto.Reset();
    modq0Histo.Reset();
    modnpxvsxmym->Reset();
    modqxvsxmym->Reset();

    modlkxBHisto.Reset();
    modlkyBHisto.Reset();
    modlkxHisto.Reset();
    modlkyHisto.Reset();
    modlkcolHisto.Reset();
    modlkrowHisto.Reset();

    modlkvst1.Reset();
    modlkvst2.Reset();
    modlkvst3.Reset();
    modlkvst6.Reset();
    ndrilkHisto.Reset();

    hdx02.Reset();
    hdy02.Reset();

    htridx.Reset();
    htridy.Reset();

    htridxc.Reset();
    tridxvsy.Reset();
    tridxvstx.Reset();
    tridxvst3.Reset();
    tridxvst6.Reset();

    htridxs1.Reset();
    htridxs2.Reset();
    htridxs3.Reset();
    htridxs4.Reset();
    htridxs5.Reset();

    htridxc1.Reset();
    htridxc2.Reset();
    htridxc3.Reset();
    htridxc4.Reset();
    htridxc5.Reset();

    htridyc.Reset();
    tridyvsx.Reset();
    tridyvsty.Reset();
    tridyvst3.Reset();
    tridyvst6.Reset();

    trix0Histo.Reset();
    triy0Histo.Reset();
    trixy0Histo->Reset();
    tritxHisto.Reset();
    trityHisto.Reset();

    ntriHisto.Reset();

    trixAHisto.Reset();
    triyAHisto.Reset();
    trixyAHisto->Reset();

    ttdxHisto.Reset();
    ttdx1Histo.Reset();

    ttdmin1Histo.Reset();
    ttdmin2Histo.Reset();

    trixcHisto.Reset();
    triycHisto.Reset();
    trixycHisto->Reset();

    tridzcvsy.Reset();
    triz3Histo.Reset();

    hsixdx.Reset();
    hsixdy.Reset();

    hsixdxc.Reset();
    hsixdxcsi.Reset();
    hsixdxccu.Reset();

    sixdxvsx.Reset();
    sixmadxvsx.Reset();
    sixdxvsy.Reset();
    sixdxvstx.Reset();
    sixdxvsdtx.Reset();
    sixdxvst3.Reset();
    sixdxvst6.Reset();
    sixmadxvsy.Reset();
    sixmadxvstx.Reset();
    sixmadxvsdtx.Reset();
    hsixdxcsid.Reset();

    hsixdyc.Reset();
    hsixdycsi.Reset();
    hsixdyccu.Reset();

    sixdyvsx.Reset();
    sixmadyvsx.Reset();
    sixdyvsy.Reset();
    sixdyvsty.Reset();
    sixdyvsdty.Reset();
    sixdyvst3.Reset();
    sixdyvst6.Reset();
    sixmadyvsy.Reset();
    sixmadyvsty.Reset();
    sixmadyvsdty.Reset();

    sixxyHisto->Reset();
    sixdxyvsxy->Reset();
    hsixdtx.Reset();
    hsixdtxsi.Reset();
    hsixdtysi.Reset();
    hsixdtxcu.Reset();
    hsixdtycu.Reset();
    hsixdty.Reset();
    sixdtvsxA.Reset();
    sixdtvsyA.Reset();
    sixdtvsxyA->Reset();
    sixdtvsx.Reset();
    sixdtvsxy->Reset();
    sixdtvsxm.Reset();
    sixdtvsym.Reset();
    sixdtvsxmym->Reset();

    cmsxvsx->Reset();
    cmsyvsy->Reset();

    cmssxaHisto.Reset();
    cmsdxaHisto.Reset();
    cmssyaHisto.Reset();
    cmsdyaHisto.Reset();

    cmsdxHisto.Reset();
    cmsdyHisto.Reset();
    cmsdxvsev.Reset();
    cmsdxvsev1->Reset();
    cmsdxvsev2->Reset();

    cmsdxcHisto.Reset();
    cmsdxvsx.Reset();
    cmsdxvsy.Reset();
    cmsdxvstx.Reset();

    cmsdyvsx.Reset();
    cmsmadyvsx.Reset();
    cmsmady8vsx.Reset();
    cmsdyvsy.Reset();
    cmsmadyvsy.Reset();

    cmsdycHisto.Reset();
    cmsdyciHisto.Reset();
    cmsdy8cHisto.Reset();
    cmsdy8ciHisto.Reset();

    cmsdyc3Histo.Reset();
    cmsdyc3iHisto.Reset();
    cmsdy8c3Histo.Reset();
    cmsdy8c3iHisto.Reset();
    cmsdyvsty.Reset();

    cmsdyvsev.Reset();
    cmsmadyvsev.Reset();
    cmsmadyvsty.Reset();
    cmsmadyvsq.Reset();

    cmsdycqHisto.Reset();
    cmsdycqiHisto.Reset();
    cmsdy8cqHisto.Reset();
    cmsdy8cqiHisto.Reset();

    cmsmadyvsxm.Reset();
    cmsmadyvsym.Reset();
    cmsmadyvsxmym->Reset();

    cmsdycq2Histo.Reset();

    cmsadcHisto.Reset();
    cmscolHisto.Reset();
    cmsrowHisto.Reset();

    cmsq0aHisto.Reset();
    cmsxmymq0->Reset();
    cmsxyq0->Reset();
    cmsxmymq1->Reset();

    trixAlkHisto.Reset();
    triyAlkHisto.Reset();
    trixyAlkHisto->Reset();
    trixclkHisto.Reset();
    triyclkHisto.Reset();
    trixyclkHisto->Reset();

    cmsnpxHisto.Reset();
    cmsncolHisto.Reset();
    cmsnrowHisto.Reset();

    cmsncolvsxm.Reset();
    cmsnrowvsxm.Reset();
    cmsncolvsym.Reset();
    cmsnrowvsym.Reset();
    cmsnpxvsxmym->Reset();
    cmsnpxvsxmym2->Reset();

    cmsph0Histo.Reset();
    cmsq0Histo.Reset();
    cmsq00Histo.Reset();
    cmsq03Histo.Reset();
    cmsq0dHisto.Reset();
    cmsq0nHisto.Reset();
    cmsqxvsx.Reset();
    cmsqxvsy.Reset();
    cmsqxvsxy->Reset();
    cmsqxvsxm.Reset();
    cmsqxvsym.Reset();
    cmsqxvsxmym->Reset();
    cmsqxvsxmym2->Reset();
    cmspxqHisto.Reset();

    effvsdmin.Reset();
    effvstmin.Reset();

    sixxylkHisto->Reset();
    sixxyeffHisto->Reset();
    effvsxy->Reset();
    effvsx.Reset();
    effvsy.Reset();

    effvsev1.Reset();
    effvsev2.Reset();
    effvst1.Reset();
    effvst2.Reset();
    effvst3.Reset();
    effvst4.Reset();
    effvst6.Reset();

    effvst1t.Reset();
    effvst2t.Reset();
    effvst3t.Reset();
    effvst4t.Reset();
    effvst6t.Reset();

    cmsdminHisto.Reset();
    cmsdxminHisto.Reset();
    cmsdyminHisto.Reset();
    cmspdxminHisto.Reset();
    cmspdyminHisto.Reset();
    cmspdminHisto.Reset();

    effvsdxy.Reset();
    effdminHisto.Reset();
    effdmin0Histo.Reset();
    effrxmin0Histo.Reset();
    effrymin0Histo.Reset();
    effdxmin0Histo.Reset();
    effdymin0Histo.Reset();

    effclq0Histo.Reset();
    effclqrHisto.Reset();
    effclq1Histo.Reset();
    effclq2Histo.Reset();
    effclq3Histo.Reset();
    effclq4Histo.Reset();
    effclq5Histo.Reset();
    effclq6Histo.Reset();
    effclq7Histo.Reset();
    effclq8Histo.Reset();
    effclq9Histo.Reset();

    effnpxvsxmym->Reset();
    effq0vsxmym->Reset();
    effq0dHisto.Reset();
    effq0nHisto.Reset();

    effvsxt->Reset();
    effvsntri.Reset();
    effvsndri.Reset();
    effvsxmym->Reset();
    effvsxmym2->Reset();
    effvsxm.Reset();
    effvsym.Reset();
    effvstx.Reset();
    effvsty.Reset();
    effvstxy.Reset();
    effvsdslp.Reset();

    prevtlutime = 0;
    prevdtbtime = 0;

    // loop over events, correlate planes:

    cout << endl << "tracking ev" << flush;

    list < vector <cluster> >::iterator evi[8];
    for( unsigned ipl = 0; ipl < 8; ++ipl )
      evi[ipl] = clist[ipl].begin();

    list <evInfo>::iterator evinfo = infoList.begin();

    iev = 0;

    for( ; evi[0] != clist[0].end() && evi[7] != clist[7].end() &&  evinfo != infoList.end();
	 ++evi[0], ++evi[1], ++evi[2], ++evi[3], ++evi[4], ++evi[5], ++evi[6], ++evi[7], ++evinfo ) {

      vector <cluster> cl[8]; // planes, one event
      for( unsigned ipl = 0; ipl < 8; ++ipl )
	cl[ipl] = *evi[ipl];

      uint64_t tlutime = evinfo->tlutime;
      if( iev == 0  )
	tlutime0 = tlutime;
      double evsec = (tlutime - tlutime0) / fTLU;
      double tludt = (tlutime - prevtlutime) / fTLU; // [s]
      prevtlutime = tlutime;

      uint64_t dtbtime = evinfo->dtbtime;
      double dtbdt = ( dtbtime - prevdtbtime ) / fDTB;
      prevdtbtime = dtbtime;
      /*
      cout << iev
	   << "  " << tlutime << "  " << tludt
	   << "  " << dtbtime << "  " << dtbdt
	   << "  " << tludt-dtbdt
	   << endl;
      */
      string filled = evinfo->filled;

      ++iev;
      if( iev%10000 == 0 )
	cout << " " << iev << flush;

      bool ldb = 0;

      if( iev <  0 )
	ldb = 1;

      if( lev < 10 )
	ldb = 1;

      if( ldb ) std::cout << "debug ev " << iev << endl << flush;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // make driplets 3+5-4:

      vector <triplet> driplets;

      double driCut = 0.1; // [mm]

      for( vector<cluster>::iterator cA = cl[3].begin(); cA != cl[3].end(); ++cA ) {

	double xA = cA->col*ptchx[3] - alignx[3];
	double yA = cA->row*ptchy[3] - aligny[3];
	double xmid = xA - midx[3];
	double ymid = yA - midy[3];
	xA = xmid - ymid*rotx[3];
	yA = ymid + xmid*roty[3];

	for( vector<cluster>::iterator cC = cl[5].begin(); cC != cl[5].end(); ++cC ) {

	  double xC = cC->col*ptchx[5] - alignx[5];
	  double yC = cC->row*ptchy[5] - aligny[5];
	  double xmid = xC - midx[5];
	  double ymid = yC - midy[5];
	  xC = xmid - ymid*rotx[5];
	  yC = ymid + xmid*roty[5];

	  double dx2 = xC - xA;
	  double dy2 = yC - yA;
	  double dz35 = zz[5] - zz[3]; // from 3 to 5 in z
	  hdx35.Fill( dx2 );
	  hdy35.Fill( dy2 );

	  if( fabs( dx2 ) > 0.005 * dz35 ) continue; // angle cut *f?
	  if( fabs( dy2 ) > 0.005 * dz35 ) continue; // angle cut

	  double avx = 0.5 * ( xA + xC ); // mid
	  double avy = 0.5 * ( yA + yC );
	  double avz = 0.5 * ( zz[3] + zz[5] ); // mid z
 
	  double slpx = ( xC - xA ) / dz35; // slope x
	  double slpy = ( yC - yA ) / dz35; // slope y

	  // middle plane B = 4:

	  for( vector<cluster>::iterator cB = cl[4].begin(); cB != cl[4].end(); ++cB ) {

	    double xB = cB->col*ptchx[4] - alignx[4];
	    double yB = cB->row*ptchy[4] - aligny[4];
	    double xmid = xB - midx[4];
	    double ymid = yB - midy[4];
	    xB = xmid - ymid*rotx[4];
	    yB = ymid + xmid*roty[4];

	    // interpolate track to B:

	    double dz = zz[4] - avz;
	    double xk = avx + slpx * dz; // driplet at k
	    double yk = avy + slpy * dz;

	    double dx3 = xB - xk;
	    double dy3 = yB - yk;
	    hdridx.Fill( dx3 );
	    hdridy.Fill( dy3 );

	    if( fabs( dy3 ) < 0.05 ) {
	      hdridxc.Fill( dx3 );
	      dridxvsy.Fill( yB, dx3 );
	      dridxvstx.Fill( slpx, dx3 );
	    }

	    if( fabs( dx3 ) < 0.05 ) {
	      hdridyc.Fill( dy3 );
	      dridyvsx.Fill( xB, dy3 );
	      dridyvsty.Fill( slpy, dy3 );
	    }

	    // telescope driplet cuts:

	    if( fabs(dx3) > driCut ) continue;
	    if( fabs(dy3) > driCut ) continue;

	    triplet dri;

	    // redefine triplet using planes 3 and 4 (A and B), avoiding MOD material:
	    avx = 0.5 * ( xB + xA ); // mid
	    avy = 0.5 * ( yB + yA );
	    avz = 0.5 * ( zz[4] + zz[3] ); // mid z
	    double dzAB = zz[4] - zz[3]; // from A to B in z
	    slpx = ( xB - xA ) / dzAB; // slope x
	    slpy = ( yB - yA ) / dzAB; // slope y

	    dri.xm = avx;
	    dri.ym = avy;
	    dri.zm = avz;
	    dri.sx = slpx;
	    dri.sy = slpy;
	    dri.lk = 0;
	    dri.ttdmin = 99.9; // isolation [mm]

	    driplets.push_back(dri);

	    drixHisto.Fill( avx );
	    driyHisto.Fill( avy );
	    drixyHisto->Fill( avx, avy );
	    dritxHisto.Fill( slpx );
	    drityHisto.Fill( slpy );

	  } // cl B

	} // cl C

      } // cl A

      ndriHisto.Fill( driplets.size() );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // driplets vs MOD:

      int nm = 0;
      int ndrilk = 0;

      double xcutMOD = 0.15;
      double ycutMOD = 0.15; // 502 in 25463 eff 99.87
      //double xcutMOD = 0.12;
      //double ycutMOD = 0.12; // 502 in 25463 eff 99.88
      //double xcutMOD = 0.10;
      //double ycutMOD = 0.10; // 502 in 25463 eff 99.88
      
      for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // jB = downstream

	double xmB = driplets[jB].xm;
	double ymB = driplets[jB].ym;
	double zmB = driplets[jB].zm;
	double sxB = driplets[jB].sx;
	double syB = driplets[jB].sy;

	double zB = MODz - zmB; // z MOD from mid of driplet
	double xB = xmB + sxB * zB; // driplet impact point on MOD
	double yB = ymB + syB * zB;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// dri vs dri: isolation at MOD

	double dddmin = 99.9;

	for( unsigned int jj = 0; jj < driplets.size(); ++jj ) {

	  if( jj == jB ) continue;

	  double xmj = driplets[jj].xm;
	  double ymj = driplets[jj].ym;
	  double sxj = driplets[jj].sx;
	  double syj = driplets[jj].sy;

	  double dz = MODz - driplets[jj].zm;
	  double xj = xmj + sxj * dz; // driplet impact point on DUT
	  double yj = ymj + syj * dz;

	  double dx = xB - xj;
	  double dy = yB - yj;
	  double dd = sqrt( dx*dx + dy*dy );
	  if( dd < dddmin )
	    dddmin = dd;

	  // intersection:

	  if( fabs( sxj - sxB ) > 0.002 ) {
	    double zi = (xmB-xmj)/(sxj - sxB);
	    drizixHisto.Fill( zi );
	    drizix2Histo.Fill( zi );
	  }
	  if( fabs( syj - syB ) > 0.002 ) {
	    double zi = (ymB-ymj)/(syj - syB);
	    driziyHisto.Fill( zi );
	  }

	} // jj

	dddmin1Histo.Fill( dddmin );
	dddmin2Histo.Fill( dddmin );
	driplets[jB].ttdmin = dddmin;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// intersect inclined track with tilted MOD plane:

	double zc = (Nzm*zB - Nym*ymB - Nxm*xmB) / (Nxm*sxB + Nym*syB + Nzm); // from zmB
	double yc = ymB + syB * zc;
	double xc = xmB + sxB * zc;

	double dzc = zc + zmB - MODz; // from MOD z0 [-8,8] mm

	// transform into MOD system: (passive).
	// large rotations don't commute: careful with order

	double x1 = com*xc - som*dzc; // turn o
	double y1 = yc;
	double z1 = som*xc + com*dzc;

	double x2 = x1;
	double y2 = cam*y1 + sam*z1; // tilt a

	double x3 = cfm*x2 + sfm*y2; // rot
	double y3 =-sfm*x2 + cfm*y2;

	double x4 =-x3 + MODalignx; // shift to mid
	double y4 = y3 + MODaligny; // invert y, shift to mid

	double xmod = fmod( 36.000 + x4, 0.3 ); // [0,0.3] mm, 2 pixel wide
	double ymod = fmod(  9.000 + y4, 0.2 ); // [0,0.2] mm

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// driplets vs MOD clusters:

	for( vector<cluster>::iterator c = cl[iMOD].begin(); c != cl[iMOD].end(); ++c ) {

	  double ccol = c->col;
	  double crow = c->row;
	  double modx = ( ccol + 0.5 - nx[iMOD]/2 ) * ptchx[iMOD]; // -33..33 mm
	  double mody = ( crow + 0.5 - ny[iMOD]/2 ) * ptchy[iMOD]; // -8..8 mm
	  double q = c->charge;
	  double q0 = q*normM;

	  bool lq = 1;
	  if( q0 < 18 ) lq = 0;
	  else if( q0 > 25 ) lq = 0;

	  double qx = exp( -q0 / qwid );

	  int npx = c->size;

	  // residuals for pre-alignment:

	  modsxaHisto.Fill( modx + x3 ); // peak
	  moddxaHisto.Fill( modx - x3 ); // 

	  modsyaHisto.Fill( mody + y3 ); // 
	  moddyaHisto.Fill( mody - y3 ); // peak

	  double moddx = modx - x4;
	  double moddy = mody - y4;

	  moddxHisto.Fill( moddx );
	  moddyHisto.Fill( moddy );

	  if( fabs( moddx ) < xcutMOD &&
	      c->big == 0 ) {

	    moddycHisto.Fill( moddy );
	    if( lq ) moddycqHisto.Fill( moddy );

	    //moddyvsx.Fill( x4, moddy ); // for rot
	    moddyvsx.Fill( -x3, moddy ); // for rot

	    //moddyvsy.Fill( y4, moddy ); // for tilt
	    moddyvsy.Fill( y2, moddy ); // for tilt

	    moddyvsty.Fill( syB, moddy );
	  }

	  if( fabs( moddy ) < ycutMOD &&
	      c->big == 0 ) {

	    moddxcHisto.Fill( moddx );
	    if( lq ) moddxcqHisto.Fill( moddx );

	    //moddxvsx.Fill( x4, moddx ); // for turn
	    moddxvsx.Fill( -x1, moddx ); // for turn

	    //moddxvsy.Fill( y4, moddx ); // for rot
	    moddxvsy.Fill( y3, moddx ); // for rot

	    moddxvstx.Fill( sxB, moddx );

	  }

	  if( fabs( moddx ) < xcutMOD &&
	      fabs( moddy ) < ycutMOD &&
	      c->big == 0 ) {
	    modnpxHisto.Fill( npx );
	    modqHisto.Fill( q );
	    modq0Histo.Fill( q0 );
	    modnpxvsxmym->Fill( xmod*1E3, ymod*1E3, npx );
	    modqxvsxmym->Fill( xmod*1E3, ymod*1E3, qx );
	  }

	  if( fabs( moddx ) < xcutMOD &&
	      fabs( moddy ) < ycutMOD ) {

	    modlkxBHisto.Fill( xB );
	    modlkyBHisto.Fill( yB );
	    modlkxHisto.Fill( x4 );
	    modlkyHisto.Fill( y4 );
	    modlkcolHisto.Fill( ccol );
	    modlkrowHisto.Fill( crow );

	    driplets[jB].lk = 1;
	    nm = 1; // we have a MOD-driplet match in this event
	    ++ndrilk;

	  } // MOD link x and y

	} // MOD

      } // driplets

      modlkvst1.Fill( evsec, nm ); // MOD yield vs time
      modlkvst2.Fill( evsec, nm );
      modlkvst3.Fill( evsec, nm );
      modlkvst6.Fill( evsec/3600, nm );
      ndrilkHisto.Fill( ndrilk );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // make triplets 2+0-1:

      vector <triplet> triplets;

      //double triCut = 0.1; // [mm]
      double triCut = 0.05; // [mm] 2.10.2017
      double zscint = -15; // [mm] scint

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

	    if( fabs( dy3 ) < triCut ) {

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

	    if( fabs( dx3 ) < triCut ) {
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

	    triplets.push_back(tri);

	    double dz0 = zscint - avz;
	    trix0Histo.Fill( avx + slpx*dz0 );
	    triy0Histo.Fill( avy + slpy*dz0 );
	    trixy0Histo->Fill( avx + slpx*dz0, avy + slpy*dz0 );
	    tritxHisto.Fill( slpx );
	    trityHisto.Fill( slpy );

	  } // cl B

	} // cl C

      } // cl A

      ntriHisto.Fill( triplets.size() );
      if( ldb ) cout << "  triplets " << triplets.size() << endl << flush;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT align vs time:

      if( run == 31153 ) {
	double p0 = -0.00467371;
	double p1 =  1.26074e-08;
	double p2 = -2.66354e-14;
	double p3 =  3.21548e-20;
	double p4 = -1.65260e-26;
	double p5 =  2.98921e-33;
	DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + p5 * iev ) * iev ) * iev ) * iev ) * iev;

	DUTaligny = DUTaligny0 - 0.0045 + 3.742e-9*iev; // trend correction
      }

      if( run == 31175 )
	DUTaligny = DUTaligny0 + 0.008 - 4.555e-9*iev; // trend correction

      if( run == 31315 ) { // cmsdyvsev->Fit("pol3")
	double p0 =   0.00136061;
	double p1 =  1.03449e-08;
	double p2 = -7.22758e-14;
	double p3 =  1.12509e-19;
	double p4 = -5.23639e-26;
	DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + p4 * iev ) * iev ) * iev ) * iev;

	p0 =  -0.00399797;
	p1 =  3.61043e-08;
	p2 = -1.47597e-13;
	p3 =  3.22517e-19;
	p4 = -3.39026e-25;
	double p5 =  1.34784e-31;
	DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + p5 * iev ) * iev ) * iev ) * iev ) * iev;

      }

      if( run == 31351 ) { // cmsdyvsev->Fit("pol3")
	DUTalignx = DUTalignx0 - 0.00218074 + ( 2.13713e-08 + ( -1.27865e-14 + 2.75008e-21*iev ) * iev ) * iev;
	DUTaligny = DUTaligny0 - 0.00195750 + ( 2.05995e-08 + ( -9.61125e-15 + 1.79382e-21*iev ) * iev ) * iev;
      }

      if( run == 31376 ) { // cmsdyvsev->Fit("pol3")
	DUTalignx = DUTalignx0 + 0.00141789 + (-3.52536e-08 + ( 5.36256e-14 -3.10804e-20*iev ) * iev ) * iev;
	DUTaligny = DUTaligny0 + 0.00126752 + (-4.88869e-08 + ( 2.63865e-13 + ( -6.82578e-19 + ( 8.13591e-25 -3.61649e-31*iev ) * iev ) * iev ) * iev ) * iev; // pol5
      }

      if( run == 31393 ) { // cmsdyvsev->Fit("pol9")
	double p0 = -0.000497698;
	double p1 =   7.4851e-09;
	double p2 = -1.18642e-13;
	double p3 =  3.84196e-19;
	double p4 = -5.30894e-25;
	double p5 =  4.00309e-31;
	double p6 = -1.78642e-37;
	double p7 =  4.73502e-44;
	double p8 = -6.91756e-51;
	double p9 =  4.30186e-58;
	DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;

	p0 =  -0.00294293;
	p1 =  1.41528e-08;
	p2 = -1.02278e-13;
	p3 =  3.71743e-19;
	p4 = -6.15689e-25;
	p5 =  5.48381e-31;
	p6 = -2.82389e-37;
	p7 =   8.4379e-44;
	p8 = -1.36005e-50;
	p9 =  9.15386e-58;
	DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;
      }

      if( run == 31473 ) { // cmsdyvsev->Fit("pol9")

	double p0 = -0.000563588;
	double p1 =  2.02659e-08;
	double p2 = -9.36375e-14;
	double p3 =  6.12014e-20;
	double p4 =  1.58543e-25;
	double p5 = -2.72422e-31;
	double p6 =  1.75228e-37;
	double p7 = -5.64298e-44;
	double p8 =  9.09636e-51;
	double p9 = -5.85518e-58;
	DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;

	p0 =  0.000341243;
	p1 =  -1.5951e-08;
	p2 =  9.43812e-14;
	p3 = -1.96965e-19;
	p4 =  1.70031e-25;
	p5 = -5.57192e-32;
	p6 = -5.40164e-39;
	p7 =  8.46203e-45;
	p8 = -2.06918e-51;
	p9 =  1.67139e-58;
	DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;
      }

      if( run == 31491 ) { // cmsdyvsev->Fit("pol9")
	double p0 = -0.000105696;
	double p1 =  1.59349e-08;
	double p2 = -4.04081e-13;
	double p3 =  3.01907e-18;
	double p4 = -9.66242e-24;
	double p5 =  1.53693e-29;
	double p6 = -1.10409e-35;
	double p7 =  5.69494e-43;
	double p8 =  3.44935e-48;
	double p9 = -1.31199e-54;
	DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;

	p0 = -0.000718211;
	p1 =  3.35502e-08;
	p2 = -6.88015e-13;
	p3 =   6.8559e-18;
	p4 = -3.28686e-23;
	p5 =  8.65805e-29;
	p6 = -1.32192e-34;
	p7 =  1.16459e-40;
	p8 =  -5.4897e-47;
	p9 =  1.07188e-53;
	DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;
      }

      if( run == 31610 ) { // cmsdyvsev->Fit("pol9")
      
	double p0 =   0.00879373;
	double p1 = -4.66823e-08;
	double p2 =  3.42111e-14;
	double p3 =  7.91648e-20;
	double p4 = -1.80482e-25;
	double p5 =  1.71464e-31;
	double p6 = -9.29736e-38;
	double p7 =  2.95246e-44;
	double p8 = -5.07344e-51;
	double p9 =  3.62521e-58;
	DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;

	p0 =  -0.00621059;
	p1 =  1.91826e-08;
	p2 = -1.50918e-14;
	p3 = -4.95739e-20;
	p4 =  1.26779e-25;
	p5 = -1.18179e-31;
	p6 =  5.50264e-38;
	p7 = -1.33154e-44;
	p8 =  1.52925e-51;
	p9 = -5.86036e-59;
	DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;
      }

      if( run == 31621 ) { // cmsdyvsev->Fit("pol9")
	double p0 =  -0.00082087;
	double p1 =  4.60501e-08;
	double p2 = -7.52843e-13;
	double p3 =  4.29444e-18;
	double p4 = -1.28897e-23;
	double p5 =  2.24598e-29;
	double p6 = -2.26304e-35;
	double p7 =  1.23123e-41;
	double p8 = -2.94649e-48;
	double p9 =  1.12485e-55;
	DUTalignx = DUTalignx0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;

	p0 =  2.53362e-05;
	p1 =   4.0548e-09;
	p2 = -5.29836e-13;
	p3 =  5.08858e-18;
	p4 = -2.19977e-23;
	p5 =  5.21882e-29;
	p6 = -7.17156e-35;
	p7 =  5.66452e-41;
	p8 = -2.37861e-47;
	p9 =  4.10141e-54;
	DUTaligny = DUTaligny0 + p0 + ( p1 + ( p2 + ( p3 + ( p4 + ( p5 + ( p6 + ( p7 + ( p8 + p9 * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;

      }

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // triplets at the DUT:

      double xcut = ptchx[iDUT];
      double ycut = ptchx[iDUT];
      if( rot90 ) {
	xcut = ptchy[iDUT];
	ycut = ptchx[iDUT];
      }
      if( fabs(DUTtilt) > 60 )
	ycut = 8;

      for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

	double xmA = triplets[iA].xm;
	double ymA = triplets[iA].ym;
	double zmA = triplets[iA].zm;
	double sxA = triplets[iA].sx;
	double syA = triplets[iA].sy;
	double txy = sqrt( sxA*sxA + syA*syA ); // angle

	double zA = DUTz - zmA; // z DUT from mid of triplet
	double xA = xmA + sxA * zA; // triplet at z(DUT)
	double yA = ymA + syA * zA;

	if( ldb ) cout << "  triplet " << iA << endl << flush;

	trixAHisto.Fill( xA ); // telescope coordinates
	triyAHisto.Fill( yA ); // telescope coordinates
	trixyAHisto->Fill( xA, yA ); // telescope coordinates

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
	//if( ttdmin > 0.33 ) liso = 1;
	if( ttdmin > 0.6 ) liso = 1; // harder cut = cleaner Q0

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// intersect inclined track with tilted DUT plane:

	double zc = (Nz*zA - Ny*ymA - Nx*xmA) / (Nx*sxA + Ny*syA + Nz); // from zmA
	double yc = ymA + syA * zc;
	double xc = xmA + sxA * zc;

	trixcHisto.Fill( xc ); // telescope coordinates
	triycHisto.Fill( yc ); // telescope coordinates
	trixycHisto->Fill( xc, yc ); // telescope coordinates

	double dzc = zc + zmA - DUTz; // from DUT z0 [-8,8] mm

	tridzcvsy.Fill( ymA, dzc );

	// transform into DUT system: (passive).
	// large rotations don't commute: careful with order

	double x1 = co*xc - so*dzc; // turn o
	double y1 = yc;
	double z1 = so*xc + co*dzc;

	double x2 = x1;
	double y2 = ca*y1 + sa*z1; // tilt a, sign cancels out
	double z2 =-sa*y1 + ca*z1; // should be zero (in DUT plane). is zero

	double x3 = cf*x2 + sf*y2; // rot
	double y3 =-sf*x2 + cf*y2;
	double z3 = z2; // should be zero (in DUT plane). is zero

	triz3Histo.Fill( z3 ); // is zero

	double x4 = upsignx*x3 + DUTalignx; // shift to mid
	double y4 = upsigny*y3 + DUTaligny; // invert y, shift to mid

	bool fiducial = 1;
	if( fabs( x4 ) > 3.9 ) fiducial = 0;
	if( fabs( y4 ) > 3.9 ) fiducial = 0;

	// reduce to 100x100 um region:

	double xmod = fmod( 9.000 + x4, 0.1 ); // [0,0.1] mm
	double ymod = fmod( 9.000 + y4, 0.1 ); // [0,0.1] mm
	if( chip0 == 142 || chip0 == 143 || chip0 == 144 )
	  ymod = fmod( 9.050 + y4, 0.1 ); // [0,0.1] mm
	double xmod2 = fmod( 9.000 + x4, 0.2 ); // [0,0.2] mm
	double ymod2 = fmod( 9.000 + y4, 0.2 ); // [0,0.2] mm

	double x8 = x4;
	double y8 = y4;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// match triplet and driplet for efficiency:

	bool lsixlk = 0;
	double sixcut = 0.1; // [mm] 502 in 25463 eff 99.87
	double dddmin = 99.9; // driplet isolation at MOD
	double sixdslp = 0.099; // [rad]

	for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // j = B = downstream

	  double xmB = driplets[jB].xm;
	  double ymB = driplets[jB].ym;
	  double zmB = driplets[jB].zm;
	  double sxB = driplets[jB].sx;
	  double syB = driplets[jB].sy;

	  // driplet at DUT:

	  double zB = zc + zmA - zmB; // z from mid of driplet to DUT intersect
	  double xB = xmB + sxB * zB; // driplet at DUT
	  double yB = ymB + syB * zB;

	  // driplet - triplet:

	  double dx = xB - xc; // at DUT intersect
	  double dy = yB - yc;
	  double dxy = sqrt( dx*dx + dy*dy );
	  double dtx = sxB - sxA;
	  double dty = syB - syA;
	  double dtxy = sqrt( dtx*dtx + dty*dty );

	  hsixdx.Fill( dx ); // for align fit
	  hsixdy.Fill( dy ); // for align fit

	  if( fabs(dy) < sixcut ) {

	    hsixdxc.Fill( dx );
	    if( x4 > xminCu && x4 < xmaxCu ) // no Cu
	      hsixdxcsi.Fill( dx );
	    else
	      hsixdxccu.Fill( dx );

	    sixdxvsx.Fill( x4, dx );
	    sixmadxvsx.Fill( x4, fabs(dx) );
	    if( x4 > xminCu && x4 < xmaxCu ) { // no Cu
	      sixdxvsy.Fill( yc, dx );
	      sixdxvstx.Fill( sxA, dx );
	      sixdxvsdtx.Fill( dtx, dx );
	      sixdxvst3.Fill( evsec, dx );
	      sixdxvst6.Fill( evsec/3600, dx );
	      sixmadxvsy.Fill( y4, fabs(dx) );
	      sixmadxvstx.Fill( sxA, fabs(dx) );
	      sixmadxvsdtx.Fill( dtx, fabs(dx) ); // U-shape
	      if( fabs( dtx ) < 0.0005 )
		hsixdxcsid.Fill( dx );
	    } // Si
	  } // dy

	  if( fabs(dx) < sixcut ) {

	    hsixdyc.Fill( dy );
	    if( x4 > xminCu && x4 < xmaxCu ) // no Cu
	      hsixdycsi.Fill( dy );
	    else
	      hsixdyccu.Fill( dy );

	    sixdyvsx.Fill( x4, dy );
	    sixmadyvsx.Fill( x4, fabs(dy) );
	    if( x4 > xminCu && x4 < xmaxCu ) { // no Cu
	      sixdyvsy.Fill( y4, dy );
	      sixdyvsty.Fill( syA, dy );
	      sixdyvsdty.Fill( dty, dy );
	      sixdyvst3.Fill( evsec, dy );
	      sixdyvst6.Fill( evsec/3600, dy );
	      sixmadyvsy.Fill( y4, fabs(dy) );
	      sixmadyvsty.Fill( syA, fabs(dy) );
	      sixmadyvsdty.Fill( dty, fabs(dy) ); // U-shape
	    }

	  } // dx

	  // match:

	  if( fabs(dx) < sixcut && fabs(dy) < sixcut ) {

	    sixxyHisto->Fill( xA, yA );

	    if( driplets[jB].lk ) { // driplet linked to MOD
	      lsixlk = 1;
	      dddmin = driplets[jB].ttdmin;
	      sixdslp = dtxy;
	    }

	    // compare slopes:

	    sixdxyvsxy->Fill( x4, y4, dxy );

	    hsixdtx.Fill( dtx );
	    if( x4 > xminCu && x4 < xmaxCu ) { // no Cu
	      hsixdtxsi.Fill( dtx );
	      hsixdtysi.Fill( dty );
	    }
	    else {
	      hsixdtxcu.Fill( dtx );
	      hsixdtycu.Fill( dty );
	    }
	    hsixdty.Fill( dty ); // width: 0.3 mrad
	    sixdtvsxA.Fill( xA, dtxy );
	    sixdtvsyA.Fill( yA, dtxy );
	    sixdtvsxyA->Fill( xA, yA, dtxy );
	    sixdtvsx.Fill( x4, dtxy );
	    sixdtvsxy->Fill( x4, y4, dtxy );
	    if( fiducial ) {
	      sixdtvsxm.Fill( xmod, dtxy );
	      sixdtvsym.Fill( ymod, dtxy );
	      sixdtvsxmym->Fill( xmod*1E3, ymod*1E3, dtxy );
	    }

	    // average driplet and triplet at DUT:

	    double xa = 0.5 * ( xB + xc );
	    double ya = 0.5 * ( yB + yc );

	    // transform into DUT system: (passive).

	    double dzc = zc + zmA - DUTz; // from DUT z0 [-8,8] mm

	    double x5 = co*xa - so*dzc; // turn o
	    double y5 = ya;
	    double z5 = so*xa + co*dzc;

	    double x6 = x5;
	    double y6 = ca*y5 + sa*z5; // tilt a

	    double x7 = cf*x6 + sf*y6; // rot
	    double y7 =-sf*x6 + cf*y6;

	    x8 = upsignx*x7 + DUTalignx; // shift to mid
	    y8 = upsigny*y7 + DUTaligny;

	    if( chip0 == 102 ) { // straight, no Cu behind DUT
	      x4 = x8; // sixplet interpol
	      y4 = y8;
	    }

	    if( chip0 == 123 ) { // straight, no Cu behind DUT
	      x4 = x8; // sixplet interpol
	      y4 = y8;
	    }

	    if( chip0 == 128 ) { // straight, no Cu behind DUT
	      x4 = x8; // sixplet interpol
	      y4 = y8;
	    }

	    if( chip0 == 133 ) { // straight, no Cu behind DUT
	      x4 = x8; // sixplet interpol
	      y4 = y8;
	    }

	    if( chip0 >= 300 ) { // straight, no Cu behind DUT
	      x4 = x8; // sixplet interpol
	      y4 = y8;
	    }

	    // update:

	    xmod = fmod( 9.000 + x4, 0.1 ); // [0,0.1] mm
	    ymod = fmod( 9.000 + y4, 0.1 ); // [0,0.1] mm
	    if( chip0 == 142 || chip0 == 143 || chip0 == 144 )
	      ymod = fmod( 9.050 + y4, 0.1 ); // [0,0.1] mm
	    xmod2 = fmod( 9.000 + x4, 0.2 ); // [0,0.2] mm
	    ymod2 = fmod( 9.000 + y4, 0.2 ); // [0,0.2] mm

	  } // six match

	} // driplets

	if( lsixlk && dddmin < 0.6 ) liso = 0; // require isolation at MOD, if present

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// DUT pixel clusters:

	int nm[99] = {0};
	double dmin = 19.9; // [mm]
	double dxmin = 9;
	double dymin = 9;
	double pdxmin = 9;
	double pdymin = 9;
	double pdmin = 19;
	double clQ0 = 0;
	int clsz0 = 0;

	for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

	  double ccol = c->col;
	  double crow = c->row;

	  // cluster isolation:

	  bool isoc = 1;
	  for( vector<cluster>::iterator c2 = cl[iDUT].begin(); c2 != cl[iDUT].end(); ++c2 ) {
	    if( c2 == c ) continue;
	    if( fabs( c2->col - ccol ) < 8 &&
		fabs( c2->row - crow ) < 8 )
	      isoc = 0;
	  }

	  if( chip0 > 300 ) // noisy 3D
	    isoc = 1;
	  if( chip0 == 122 ) // noisy irrad
	    isoc = 1;
	  if( chip0 == 123 ) // noisy irrad
	    isoc = 1;
	  if( chip0 == 128 ) // noisy irrad
	    isoc = 1;
	  if( chip0 == 133 ) // noisy irrad
	    isoc = 1;

	  double Q0 = c->charge * norm; // cluster charge normalized to vertical incidence

	  double Qx = exp( -Q0 / qwid ); // Moyal weight

	  int colmin = 999;
	  int colmax = -1;
	  int rowmin = 999;
	  int rowmax = -1;

	  double qcol[nx[iDUT]];
	  for( int icol = 0; icol < nx[iDUT]; ++icol ) qcol[icol] = 0;

	  double qrow[ny[iDUT]];
	  for( int irow = 0; irow < ny[iDUT]; ++irow ) qrow[irow] = 0;

	  double sumph = 0;

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

	    qcol[icol] += q; // project cluster onto cols
	    qrow[irow] += q; // project cluster onto rows

	    sumph += px->adc;

	  } // pix

	  int ncol = colmax - colmin + 1;
	  int nrow = rowmax - rowmin + 1;

	  // eta-algo in rows:

	  double q1 = 0; // highest charge
	  double q2 = 0; // 2nd highest
	  int i1 = 0;
	  int i2 = 0;
	  double sumq = 0;
	  double sumrow = 0;
	  double sumrow2 = 0;
	  double sumrow3 = 0;

	  for( int irow = rowmin; irow <= rowmax; ++irow ) {

	    double q = qrow[irow];
	    sumq += q;
	    sumrow += irow*q;
	    double drow = irow - crow; // for central moments
	    sumrow2 += drow*drow*q;
	    sumrow3 += drow*drow*drow*q;

	    if( q > q1 ) {
	      q2 = q1;
	      q1 = q;
	      i2 = i1;
	      i1 = irow;
	    }
	    else if( q > q2 ) {
	      q2 = q;
	      i2 = irow;
	    }

	  } // rows

	  double q12 = q1 + q2;
	  double eta = 0;
	  if( q12 > 1 ) eta = ( q1 - q2 ) / q12;
	  if( i2 > i1 ) eta = -eta;

	  // column cluster:

	  sumq = 0;
	  double sumcol = 0;
	  double sumcol2 = 0;
	  double sumcol3 = 0;
	  double p1 = 0; // highest charge
	  double p2 = 0; // 2nd highest
	  int j1 = 0;
	  int j2 = 0;

	  for( int icol = colmin; icol <= colmax; ++icol ) {

	    double q = qcol[icol];
	    if( q > p1 ) {
	      p2 = p1;
	      p1 = q;
	      j2 = j1;
	      j1 = icol;
	    }
	    else if( q > p2 ) {
	      p2 = q;
	      j2 = icol;
	    }

	    //Tue 21.7.2015 if( q > 17 ) q = 17; // truncate Landau tail [ke]
	    sumq += q;
	    sumcol += icol*q;
	    double dcol = icol - ccol; // distance from COG
	    sumcol2 += dcol*dcol*q; // 2nd central moment
	    sumcol3 += dcol*dcol*dcol*q; // 3rd central moment

	  } // cols

	  double p12 = p1 + p2;
	  double uta = 0;
	  if( p12 > 1 ) uta = ( p1 - p2 ) / p12;
	  if( j2 > j1 ) uta = -uta;

	  if( rot90 )
	    eta = uta;

	  bool lq = 1; // Landau peak
	  if( ( Q0 < 8 || Q0 > 16 ) ) lq = 0; // r102 50x50
	  if( chip0 >= 300 ) { // 3D
	    lq = 1;
	    if( ( Q0 < 11 || Q0 > 22 ) ) lq = 0; // 3D
	  }

	  // DUT - triplet:

	  double cmsx = ( ccol + 1.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	  double cmsy = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm

	  if( rot90 ) {
	    cmsx = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	    cmsy = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	  }

	  if( liso &&  isoc ) {
	    cmsxvsx->Fill( x4, cmsx );
	    cmsyvsy->Fill( y4, cmsy );
	  }

	  // residuals for pre-alignment:

	  if( liso &&  isoc ) {

	    cmssxaHisto.Fill( cmsx + x3 ); // rot, tilt and turn but no shift
	    cmsdxaHisto.Fill( cmsx - x3 ); // peak

	    cmssyaHisto.Fill( cmsy + y3 );
	    cmsdyaHisto.Fill( cmsy - y3 );

	  }

	  // residuals:

	  double cmsdx = cmsx - x4; // triplet extrapol
	  double cmsdy = cmsy - y4;

	  //double cmsdx8 = cmsx - x8; // sixplet interpol
	  double cmsdy8 = cmsy - y8;

	  double dxy = sqrt( cmsdx*cmsdx + cmsdy*cmsdy );

	  if( ldb ) cout << "    DUT dxy " << dxy << endl << flush;

	  // for eff: nearest pixel

	  for( unsigned ipx = 0; ipx < c->vpix.size(); ++ipx ) {

	    double px = ( c->vpix[ipx].col + 1.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	    double py = ( c->vpix[ipx].row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	    if( rot90 ) {
	      px = ( c->vpix[ipx].row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	      py = ( c->vpix[ipx].col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	    }
	    double pdx = px - x4; // triplet extrapol
	    double pdy = py - y4;
	    double pdxy = sqrt( pdx*pdx + pdy*pdy );
	    if( fabs(pdx) < fabs(pdxmin) ) pdxmin = pdx;
	    if( fabs(pdy) < fabs(pdymin) ) pdymin = pdy;
	    if( fabs(pdxy) < fabs(pdmin) ) pdmin = pdxy;
	    for( int iw = 1; iw < 99; ++iw )
	      if( pdxy < iw*0.010 ) // 10 um bins
		nm[iw] = 1; // eff

	  } // pix

	  if( liso &&  isoc ) {

	    cmsdxHisto.Fill( cmsdx );
	    cmsdyHisto.Fill( cmsdy );
	    cmsdxvsev.Fill( iev, cmsdx ); // align stability
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

	  } // iso

	  // for dy:

	  if( fabs(cmsdx) < xcut ) {

	    if( liso && isoc && fabs(y4) < 3.9 ) { // fiducial y
	      cmsdyvsx.Fill( x4, cmsdy );
	      cmsmadyvsx.Fill( x4, fabs(cmsdy) );
	      cmsmady8vsx.Fill( x4, fabs(cmsdy8) );
	    }

	    if( liso && isoc && fabs(x4) < 3.9 ) { // fiducial x
	      cmsdyvsy.Fill( y4, cmsdy );
	      cmsmadyvsy.Fill( y4, fabs(cmsdy) );
	    }

	    if( fabs(x4) < 3.9 && fabs(y4) < 3.9 ) { // fiducial x and y

	      cmsdycHisto.Fill( cmsdy );
	      if( lsixlk && liso && isoc )
		cmsdyciHisto.Fill( cmsdy );

	      if( x4 < 1.4 ) { // Cu cutout rot90
		cmsdy8cHisto.Fill( cmsdy8 );
		if( lsixlk && liso &&  isoc )
		  cmsdy8ciHisto.Fill( cmsdy8 );
	      }

	      if( ncol < 3 && nrow < 3 ) {
		cmsdyc3Histo.Fill( cmsdy ); // 31166: side peaks at +-0.125 mm
		if( lsixlk && liso &&  isoc )
		  cmsdyc3iHisto.Fill( cmsdy ); // 31166: side peaks eliminated
		if( x4 < 1.4 ) { // Cu cutout rot90
		  cmsdy8c3Histo.Fill( cmsdy8 );
		  if( lsixlk && liso &&  isoc )
		    cmsdy8c3iHisto.Fill( cmsdy8 );
		} // Cu
	      } // ncol

	      if( liso &&  isoc ) {
		cmsdyvsty.Fill( syA*1E3, cmsdy );
		cmsdyvsev.Fill( iev, cmsdy ); // trend?

		cmsmadyvsev.Fill( iev, fabs(cmsdy) );
		cmsmadyvsty.Fill( syA*1E3, fabs(cmsdy) );
		cmsmadyvsq.Fill( Q0, fabs(cmsdy) ); // resolution vs charge
	      }

	      if( lq ) {

		cmsdycqHisto.Fill( cmsdy ); // 3D 31215: 5.9
		if( lsixlk && liso &&  isoc )
		  cmsdycqiHisto.Fill( cmsdy );

		if( x4 < 1.4 ) { // Cu cutout rot90
		  cmsdy8cqHisto.Fill( cmsdy8 );
		  if( lsixlk && liso &&  isoc )
		    cmsdy8cqiHisto.Fill( cmsdy8 );
		}

		if( liso &&  isoc ) {
		  cmsmadyvsxm.Fill( xmod*1E3, fabs(cmsdy) ); // within pixel
		  cmsmadyvsym.Fill( ymod*1E3, fabs(cmsdy) ); // within pixel
		  cmsmadyvsxmym->Fill( xmod*1E3, ymod*1E3, fabs(cmsdy) ); // within pixel
		}

	      } // Q0

	      if( Q0 > 9 && Q0 < 14 )
		cmsdycq2Histo.Fill( cmsdy );

	    } // fiducial

	  } // dx

	  // xy cuts:

	  if( fabs(cmsdx) < xcut &&
	      fabs(cmsdy) < ycut ) {

	    for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {
	      cmsadcHisto.Fill( px->adc );
	      cmscolHisto.Fill( px->col );
	      cmsrowHisto.Fill( px->row );
	    }

	  } // dxy

	  if( fabs(x4) < 3.9 && fabs(y4) < 3.9 && // fiducial
	      fabs(cmsdx) < xcut &&
	      fabs(cmsdy) < ycut ) {

	    cmsq0aHisto.Fill( Q0 ); // Landau
	    if( Q0 < 7.5 ) {
	      cmsxmymq0->Fill( xmod*1E3, ymod*1E3 ); // cluster size map
	      cmsxyq0->Fill( x4, y4 ); // cluster size map
	    }
	    else
	      cmsxmymq1->Fill( xmod*1E3, ymod*1E3 ); // cluster size map

	  } // fiducial linked

	  if( liso && isoc &&
	      fabs(cmsdx) < xcut &&
	      fabs(cmsdy) < ycut ) {

	    trixAlkHisto.Fill( xA ); // telescope coordinates at DUT
	    triyAlkHisto.Fill( yA );
	    trixyAlkHisto->Fill( xA, yA );
	    trixclkHisto.Fill( xc ); // telescope coordinates crossing DUT
	    triyclkHisto.Fill( yc );
	    trixyclkHisto->Fill( xc, yc );

	    if( fabs(x4) < 3.9 && fabs(y4) < 3.9 ) { // isolated fiducial

	      cmsnpxHisto.Fill( c->size );
	      cmsncolHisto.Fill( ncol );
	      cmsnrowHisto.Fill( nrow );

	      cmsncolvsxm.Fill( xmod*1E3, ncol );
	      cmsnrowvsxm.Fill( xmod*1E3, nrow );

	      cmsncolvsym.Fill( ymod*1E3, ncol ); // within pixel
	      cmsnrowvsym.Fill( ymod*1E3, nrow ); // within pixel

	      cmsnpxvsxmym->Fill( xmod*1E3, ymod*1E3, c->size ); // cluster size map
	      cmsnpxvsxmym2->Fill( xmod2*1E3, ymod2*1E3, c->size ); // cluster size map

	      cmsph0Histo.Fill( sumph*norm ); // Landau
	      cmsq0Histo.Fill( Q0 ); // Landau
	      cmsq00Histo.Fill( Q0 ); // Landau

	      if( ncol < 3 && nrow < 3 )
		cmsq03Histo.Fill( Q0 ); // Landau
 
	      if( sqrt( pow( xmod-0.050, 2 ) + pow( ymod-0.050, 2 ) ) < 0.010 ) // 50x50 bias dot
		cmsq0dHisto.Fill( Q0 );

	      if( sqrt( pow( xmod-0.050, 2 ) + pow( ymod-0.050, 2 ) ) > 0.020 ) // no dot
		cmsq0nHisto.Fill( Q0 );

	      // cluster charge profiles, exponential weighting Qx:

	      cmsqxvsx.Fill( x4, Qx );
	      cmsqxvsy.Fill( y4, Qx );
	      cmsqxvsxy->Fill( x4, y4, Qx );
	      cmsqxvsxm.Fill( xmod*1E3, Qx ); // Q within pixel, depends on ph cut
	      cmsqxvsym.Fill( ymod*1E3, Qx ); // Q within pixel
	      cmsqxvsxmym->Fill( xmod*1E3, ymod*1E3, Qx ); // cluster charge profile
	      cmsqxvsxmym2->Fill( xmod2*1E3, ymod2*1E3, Qx ); // cluster charge profile

	      for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px )
		cmspxqHisto.Fill( px->q );

	    } // fiducial

	  } // cut xy

	  if( dxy < dmin ) {
	    dmin = dxy;
	    clQ0 = Q0;
	    clsz0 = c->size;
	    dxmin = cmsdx;
	    dymin = cmsdy;
	  }

	} // loop DUT clusters

	if( ldb ) cout << "    eff " << nm[49] << endl << flush;

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// DUT efficiency vs isolated MOD-linked fiducial tracks:

	if( lsixlk
	    //&& cl[iDUT].size() < 2 // empty or single cluster, same eff
	    ) {

	  double fidx0 =-3.8;
	  double fidx9 = 3.8;
	  double fidy0 =-3.8;
	  double fidy9 = 3.8;

	  bool lddt = 1;
	  if( fabs(tludt - dtbdt) > 0.02E-3 ) lddt = 0; // NED TO REMEMBER

	  if( x4 > fidx0 && x4 < fidx9 &&
	      y4 > fidy0 && y4 < fidy9 &&
	      lddt ) { // fiducial

	    effvsdmin.Fill( dddmin, nm[49] ); // at MOD, small effect

	    if( dddmin > 0.4 )
	      effvstmin.Fill( ttdmin, nm[49] ); // at DUT, flat

	  } // fid

	  //if( dddmin > 0.4 && ttdmin > 1.4 ) { // stronger iso [mm] 99.94
	  //if( dddmin > 0.4 && ttdmin > 0.6 ) { // iso [mm] 99.93
	  if( iev > 100 &&
	      ( run != 31147 || iev < 6400 || iev > 8900 ) &&
	      filled != A &&
	      dddmin > 0.6 ) { // iso [mm] 99.94

	    if( lddt ) {

	      sixxylkHisto->Fill( xA, yA );
	      if( nm[49] ) sixxyeffHisto->Fill( xA, yA );

	      effvsxy->Fill( x4, y4, nm[49] ); // map

	      if( y4 > fidy0 && y4 < fidy9 )
		effvsx.Fill( x4, nm[49] );

	      if( x4 > fidx0 && x4 < fidx9 )
		effvsy.Fill( y4, nm[49] );

	    } // ddt

	    bool fiducial = 0;

	    if( x4 > fidx0 && x4 < fidx9 &&
		y4 > fidy0 && y4 < fidy9 ) { // fiducial
	      fiducial = 1;
	    }

	    if( chip0 == 142 && y4 < -1.5 && y4 > -1.6 )
	      fiducial = 0;

	    if( fiducial ) {

	      effvsev1.Fill( iev, nm[49] );
	      effvsev2.Fill( iev, nm[49] );
	      effvst1.Fill( evsec, nm[49] );
	      effvst2.Fill( evsec, nm[49] );
	      effvst3.Fill( evsec, nm[49] );
	      effvst4.Fill( evsec, nm[49] );
	      effvst6.Fill( evsec/3600, nm[49] );

	      if( lddt ) {

		effvst1t.Fill( evsec, nm[49] );
		effvst2t.Fill( evsec, nm[49] );
		effvst3t.Fill( evsec, nm[49] );
		effvst4t.Fill( evsec, nm[49] );
		effvst6t.Fill( evsec/3600, nm[49] );

		cmsdminHisto.Fill( dmin );
		cmsdxminHisto.Fill( dxmin );
		cmsdyminHisto.Fill( dymin );
		cmspdxminHisto.Fill( pdxmin );
		cmspdyminHisto.Fill( pdymin );
		cmspdminHisto.Fill( pdmin );

		for( int iw = 1; iw < 99; ++iw )
		  effvsdxy.Fill( iw*0.010-0.001, nm[iw] );

		effdminHisto.Fill( dmin );
		if( nm[49] == 0 ) {
		  effdmin0Histo.Fill( dmin );
		  effrxmin0Histo.Fill( dxmin/dmin );
		  effrymin0Histo.Fill( dymin/dmin );
		  effdxmin0Histo.Fill( dxmin );
		  effdymin0Histo.Fill( dymin );
		}

		effclq0Histo.Fill( clQ0 );
		if( ndrilk == 1 ) effclqrHisto.Fill( clQ0 );
		if( dmin > 0.1 ) effclq1Histo.Fill( clQ0 );
		if( dmin > 0.2 ) effclq2Histo.Fill( clQ0 ); // Landau tail
		if( dmin > 0.3 ) effclq3Histo.Fill( clQ0 );
		if( dmin > 0.4 ) effclq4Histo.Fill( clQ0 );
		if( dmin > 0.5 ) effclq5Histo.Fill( clQ0 );
		if( dmin > 0.6 ) effclq6Histo.Fill( clQ0 );
		if( dmin > 0.7 ) effclq7Histo.Fill( clQ0 );
		if( dmin > 0.8 ) effclq8Histo.Fill( clQ0 );
		if( dmin > 0.9 ) effclq9Histo.Fill( clQ0 );

		effnpxvsxmym->Fill( xmod*1E3, ymod*1E3, clsz0 ); // cluster size map
		effq0vsxmym->Fill( xmod*1E3, ymod*1E3, clQ0 ); // cluster charge profile
		if( sqrt( pow( xmod-0.050, 2 ) + pow( ymod-0.050, 2 ) ) < 0.010 ) // 50x50 bias dot
		  effq0dHisto.Fill( clQ0 );
		if( sqrt( pow( xmod-0.050, 2 ) + pow( ymod-0.050, 2 ) ) > 0.020 ) // no dot
		  effq0nHisto.Fill( clQ0 );

		effvsxt->Fill( evsec, x4, nm[49] );
		effvsntri.Fill( triplets.size(), nm[49] ); // flat
		effvsndri.Fill( driplets.size(), nm[49] ); // flat
		effvsxmym->Fill( xmod*1E3, ymod*1E3, nm[49] );
		effvsxmym2->Fill( xmod2*1E3, ymod2*1E3, nm[49] );
		effvsxm.Fill( xmod, nm[49] ); // bias dot
		effvsym.Fill( ymod, nm[49] ); // bias dot
		effvstx.Fill( sxA, nm[49] );
		effvsty.Fill( syA, nm[49] );
		effvstxy.Fill( txy, nm[49] ); // no effect
		effvsdslp.Fill( sixdslp, nm[49] ); // no effect

	      } // ddt

	    } // fiducial

	  } // iso

	} // six

      } // loop triplets iA

      if( ldb ) cout << "done ev " << iev << endl << flush;

    } // events

    cout << endl;

    if( aligniteration+1 == maxiter ) {
        histoFile->Write();
	histoFile->Close();
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // MOD alignment:

    cout << endl
	 << "MOD alignment fits for iteration " << MODaligniteration << endl;

    if( moddxaHisto.GetEntries() > 9999 ) {

      double newMODalignx = MODalignx;
      double newMODaligny = MODaligny;

      if( moddxaHisto.GetMaximum() > modsxaHisto.GetMaximum() ) {
	cout << endl << moddxaHisto.GetTitle()
	     << " bin " << moddxaHisto.GetBinWidth(1)
	     << endl;
	TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
	double xpk = moddxaHisto.GetBinCenter( moddxaHisto.GetMaximumBin() );
	fgp0->SetParameter( 0, moddxaHisto.GetMaximum() ); // amplitude
	fgp0->SetParameter( 1, xpk );
	fgp0->SetParameter( 2, moddxaHisto.GetBinWidth(1) ); // sigma
	fgp0->SetParameter( 3, moddxaHisto.GetBinContent( moddxaHisto.FindBin(xpk-1) ) ); // BG
	moddxaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0->GetParameter(0)
	     << endl << "mid " << fgp0->GetParameter(1)
	     << endl << "sig " << fgp0->GetParameter(2)
	     << endl << " BG " << fgp0->GetParameter(3)
	     << endl;
	newMODalignx = fgp0->GetParameter(1);
	delete fgp0;
      }
      else {
	cout << endl << modsxaHisto.GetTitle()
	     << " bin " << modsxaHisto.GetBinWidth(1)
	     << endl;
	TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
	double xpk = modsxaHisto.GetBinCenter( modsxaHisto.GetMaximumBin() );
	fgp0->SetParameter( 0, modsxaHisto.GetMaximum() ); // amplitude
	fgp0->SetParameter( 1, xpk );
	fgp0->SetParameter( 2, modsxaHisto.GetBinWidth(1) ); // sigma
	fgp0->SetParameter( 3, modsxaHisto.GetBinContent( modsxaHisto.FindBin(xpk-1) ) ); // BG
	modsxaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1  );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0->GetParameter(0)
	     << endl << "mid " << fgp0->GetParameter(1)
	     << endl << "sig " << fgp0->GetParameter(2)
	     << endl << " BG " << fgp0->GetParameter(3)
	     << endl;
	newMODalignx = fgp0->GetParameter(1);
	delete fgp0;
      }

      if( moddyaHisto.GetMaximum() > modsyaHisto.GetMaximum() ) {
	cout << endl << moddyaHisto.GetTitle()
	     << " bin " << moddyaHisto.GetBinWidth(1)
	     << endl;
	TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
	double xpk = moddyaHisto.GetBinCenter( moddyaHisto.GetMaximumBin() );
	fgp0->SetParameter( 0, moddyaHisto.GetMaximum() ); // amplitude
	fgp0->SetParameter( 1, xpk );
	fgp0->SetParameter( 2, moddyaHisto.GetBinWidth(1) ); // sigma
	fgp0->SetParameter( 3, moddyaHisto.GetBinContent( moddyaHisto.FindBin(xpk-1) ) ); // BG
	moddyaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0->GetParameter(0)
	     << endl << "mid " << fgp0->GetParameter(1)
	     << endl << "sig " << fgp0->GetParameter(2)
	     << endl << " BG " << fgp0->GetParameter(3)
	     << endl;
	newMODaligny = fgp0->GetParameter(1);
	delete fgp0;
      }
      else {
	cout << endl << modsyaHisto.GetTitle()
	     << " bin " << modsyaHisto.GetBinWidth(1)
	     << endl;
	TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
	double xpk = modsyaHisto.GetBinCenter( modsyaHisto.GetMaximumBin() );
	fgp0->SetParameter( 0, modsyaHisto.GetMaximum() ); // amplitude
	fgp0->SetParameter( 1, xpk );
	fgp0->SetParameter( 2, modsyaHisto.GetBinWidth(1) ); // sigma
	fgp0->SetParameter( 3, modsyaHisto.GetBinContent( modsyaHisto.FindBin(xpk-1) ) ); // BG
	modsyaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0->GetParameter(0)
	     << endl << "mid " << fgp0->GetParameter(1)
	     << endl << "sig " << fgp0->GetParameter(2)
	     << endl << " BG " << fgp0->GetParameter(3)
	     << endl;
	newMODaligny = fgp0->GetParameter(1);
	delete fgp0;
      }

      // finer alignment:

      if( MODaligniteration > 0 && fabs( newMODalignx - MODalignx ) < 0.1 ) {

	cout << endl << moddxcHisto.GetTitle()
	     << " bin " << moddxcHisto.GetBinWidth(1)
	     << endl;
	TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
	fgp0->SetParameter( 0, moddxcHisto.GetMaximum() ); // amplitude
	fgp0->SetParameter( 1, moddxcHisto.GetBinCenter( moddxcHisto.GetMaximumBin() ) );
	fgp0->SetParameter( 2, 8*moddxcHisto.GetBinWidth(1) ); // sigma
	fgp0->SetParameter( 3, moddxcHisto.GetBinContent(1) ); // BG
	moddxcHisto.Fit( "fgp0", "q" );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0->GetParameter(0)
	     << endl << "mid " << fgp0->GetParameter(1)
	     << endl << "sig " << fgp0->GetParameter(2)
	     << endl << " BG " << fgp0->GetParameter(3)
	     << endl;
	newMODalignx = MODalignx + fgp0->GetParameter(1);
	delete fgp0;

	// dxvsx -> turn:

	if( fabs(som) > 0.01 ) {
	  moddxvsx.Fit( "pol1", "q", "", -midx[iMOD]+0.2, midx[iMOD]-0.2 );
	  TF1 * fdxvsx = moddxvsx.GetFunction( "pol1" );
	  cout << endl << moddxvsx.GetTitle()
	       << ": slope " << fdxvsx->GetParameter(1)
	       << ", extra turn " << fdxvsx->GetParameter(1)/wt/som
	       << " deg"
	       << endl;
	  MODturn += fdxvsx->GetParameter(1)/wt/som; // [deg] min 0.6 deg
	  delete fdxvsx;
	}

      } // iter

      if( MODaligniteration > 0 && fabs( newMODaligny - MODaligny ) < 0.1 ) {

	cout << endl << moddycHisto.GetTitle()
	     << " bin " << moddycHisto.GetBinWidth(1)
	     << endl;
	TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
	fgp0->SetParameter( 0, moddycHisto.GetMaximum() ); // amplitude
	fgp0->SetParameter( 1, moddycHisto.GetBinCenter( moddycHisto.GetMaximumBin() ) );
	fgp0->SetParameter( 2, 5*moddycHisto.GetBinWidth(1) ); // sigma
	fgp0->SetParameter( 3, moddycHisto.GetBinContent(1) ); // BG
	moddycHisto.Fit( "fgp0", "q" );
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0->GetParameter(0)
	     << endl << "mid " << fgp0->GetParameter(1)
	     << endl << "sig " << fgp0->GetParameter(2)
	     << endl << " BG " << fgp0->GetParameter(3)
	     << endl;
	newMODaligny = MODaligny + fgp0->GetParameter(1);
	delete fgp0;

	// dyvsx -> rot

	moddyvsx.Fit( "pol1", "q", "", -midx[iMOD]+0.2, midx[iMOD]-0.2 );
	TF1 * fdyvsx = moddyvsx.GetFunction( "pol1" );
	cout << endl << moddyvsx.GetTitle()
	     << ": extra rot " << fdyvsx->GetParameter(1) << endl;
	MODrot += fdyvsx->GetParameter(1);
	delete fdyvsx;

	// dyvsy -> tilt:

	if( fabs( sam ) > 0.01 ) {
	  moddyvsy.Fit( "pol1", "q", "", -midy[iMOD]+0.2, midy[iMOD]-0.2 );
	  TF1 * fdyvsy = moddyvsy.GetFunction( "pol1" );
	  cout << endl << moddyvsy.GetTitle()
	       << ": slope " << fdyvsy->GetParameter(1)
	       << ", extra tilt " << fdyvsy->GetParameter(1)/wt/sam
	       << " deg"
	       << endl;
	  MODtilt += fdyvsy->GetParameter(1)/wt/sam; // [deg] min 0.6 deg
	  delete fdyvsy;
	}

	// dyvsty -> dz:

	moddyvsty.Fit( "pol1", "q", "", -0.002, 0.002 );
	TF1 * fdyvsty = moddyvsty.GetFunction( "pol1" );
	cout << endl << moddyvsty.GetTitle()
	     << ": z shift " << fdyvsty->GetParameter(1)
	     << " mm"
	     << endl;
	MODz += fdyvsty->GetParameter(1);
	delete fdyvsty;

      }

      MODalignx = newMODalignx;
      MODaligny = newMODaligny;
      ++MODaligniteration;

    } // MOD
    else
      cout << "  not enough   data" << endl;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT alignment:

    cout << endl
	 << "DUT alignment fits for iteration " << DUTaligniteration << endl;

    double newDUTalignx = DUTalignx0; // time 0
    double newDUTaligny = DUTaligny0;

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
	cmsdxaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 ); // fit range around peak
	cout << "Fit Gauss + BG:"
	     << endl << "  A " << fgp0->GetParameter(0)
	     << endl << "mid " << fgp0->GetParameter(1)
	     << endl << "sig " << fgp0->GetParameter(2)
	     << endl << " BG " << fgp0->GetParameter(3)
	     << endl;
	newDUTalignx = fgp0->GetParameter(1);
	delete fgp0;
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
	delete fgp0;
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
	delete fgp0;
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
	delete fgp0;
      }

    } // cmsdxa
    else
      cout << "  not enough   data" << endl;

    // finer alignment:

    if( DUTaligniteration > 0 && fabs( newDUTalignx - DUTalignx0 ) < 0.1 &&
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
      newDUTalignx = DUTalignx0 + fgp0->GetParameter(1);
      delete fgp0;

    }

    if( DUTaligniteration > 0 && fabs( newDUTaligny - DUTaligny0 ) < 0.1 &&
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
      newDUTaligny = DUTaligny0 + fgp0->GetParameter(1);
      delete fgp0;

      // dyvsx -> -rot

      if( cmsdyvsx.GetEntries() > 999 ) {
	cmsdyvsx.Fit( "pol1", "q", "", -midx[iDUT]+0.2, midx[iDUT]-0.2 );
	TF1 * fdyvsx = cmsdyvsx.GetFunction( "pol1" );
	cout << endl << cmsdyvsx.GetTitle();
	if( rot90 )
	  cout << ": extra rot " << upsigny*fdyvsx->GetParameter(1);
	else
	  cout << ": extra rot " << -upsigny*fdyvsx->GetParameter(1);
	cout << endl;
	if( rot90 )
	  DUTrot += upsigny*fdyvsx->GetParameter(1);
	else
	  DUTrot -= upsigny*fdyvsx->GetParameter(1);
	delete fdyvsx;
      }

      // dyvsy -> tilt:

      if( fabs(DUTtilt0) > 0.4 && cmsdyvsy.GetEntries() > 999 ) {
	cmsdyvsy.Fit( "pol1", "q", "", -midy[iDUT]+0.2, midy[iDUT]-0.2 );
	TF1 * fdyvsy = cmsdyvsy.GetFunction( "pol1" );
	cout << endl << cmsdyvsy.GetTitle()
	     << ": slope " << fdyvsy->GetParameter(1)
	     << ", extra tilt " << fdyvsy->GetParameter(1)/wt/sa
	     << " deg"
	     << endl;
	DUTtilt += fdyvsy->GetParameter(1)/wt/sa;
	delete fdyvsy;
      }

      // dyvsty -> dz:

      if( cmsdyvsty.GetEntries() > 999 ) {
	cmsdyvsty.Fit( "pol1", "q", "", -2, 2 );
	TF1 * fdyvsty = cmsdyvsty.GetFunction( "pol1" );
	cout << endl << cmsdyvsty.GetTitle()
	     << ": z shift " << upsigny*fdyvsty->GetParameter(1)*1E3
	     << " mm"
	     << endl;
	DUTz += upsigny*fdyvsty->GetParameter(1)*1E3;
	delete fdyvsty;
      }

    } // iter

    DUTalignx = newDUTalignx; // time 0
    DUTaligny = newDUTaligny;
    ++DUTaligniteration;

    //cout << "Enter any letter to continue;" << endl; string any; cin >> any;

    ++niter;

  } // aligniterations

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  // write new MOD alignment:

  ofstream MODalignFile( MODalignFileName.str() );

  MODalignFile << "# MOD alignment for run " << run << endl;
  MODalignFile << "iteration " << MODaligniteration << endl;
  MODalignFile << "alignx " << MODalignx << endl;
  MODalignFile << "aligny " << MODaligny << endl;
  MODalignFile << "rot " << MODrot << endl;
  MODalignFile << "tilt " << MODtilt << endl;
  MODalignFile << "turn " << MODturn << endl;
  MODalignFile << "dz " << MODz - zz[4] << endl;

  MODalignFile.close();

  cout << endl << "wrote MOD alignment iteration " << MODaligniteration
       << " to " << MODalignFileName.str() << endl
       << "  alignx " << MODalignx << endl
       << "  aligny " << MODaligny << endl
       << "  rot    " << MODrot << endl
       << "  tilt   " << MODtilt << endl
       << "  turn   " << MODturn << endl
       << "  dz     " << MODz - zz[4] << endl
    ;

  // write new DUT alignment:

  cout << "update DUT alignment file? (y/n)" << endl;
  string ans;
  cin >> ans;
  string YES{"y"};
  if( ans == YES ) {

    ofstream DUTalignFile( DUTalignFileName.str() );

    DUTalignFile << "# DUT alignment for run " << run << endl;
    DUTalignFile << "iteration " << DUTaligniteration << endl;
    DUTalignFile << "alignx " << DUTalignx << endl;
    DUTalignFile << "aligny " << DUTaligny << endl;
    DUTalignFile << "rot " << DUTrot << endl;
    DUTalignFile << "tilt " << DUTtilt << endl;
    DUTalignFile << "turn " << DUTturn << endl;
    DUTalignFile << "dz " << DUTz - zz[2] << endl;

    DUTalignFile.close();

    cout << endl << "wrote DUT alignment iteration " << DUTaligniteration
	 << " to " << DUTalignFileName.str() << endl
	 << "  alignx " << DUTalignx << endl
	 << "  aligny " << DUTaligny << endl
	 << "  rot    " << DUTrot << endl
	 << "  tilt   " << DUTtilt << endl
	 << "  turn   " << DUTturn << endl
	 << "  dz     " << DUTz - zz[2] << endl
      ;
  }

  cout << endl
       << "DUT efficiency " << 100*effvst6t.GetMean(2) << "%"
       << " from " << effvst6t.GetEntries() << " events"
       << endl;

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s6 = ts.tv_sec; // seconds since 1.1.1970
  long f6 = ts.tv_nsec; // nanoseconds
  zeit6 += s6 - s5 + ( f6 - f5 ) * 1e-9; // track

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s9 = ts.tv_sec; // seconds since 1.1.1970
  long f9 = ts.tv_nsec; // nanoseconds

  cout << endl
       << "done after " << infoList.size() << " events"
       << " in " << s9 - s0 + ( f9 - f0 ) * 1e-9 << " s"
       << endl << "eudaq event  " << zeit1 << " s"
       << endl << "eudaq planes " << zeit2 << " s"
       << endl << "tele clus    " << zeit3 << " s"
       << endl << "R4S DUT      " << zeit5 << " s"
       << endl << "tracking     " << zeit6 << " s"
       << " (with " << niter << " alignment iterations)"
       << endl << "timed sum    " << zeit1+zeit2+zeit3+zeit5+zeit6 << " s"
       << endl;

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  /*
  delete histoFile;
  for( int ipl = 0; ipl < 9; ++ipl ) delete hmap[ipl];
  delete drixyHisto;
  delete modnpxvsxmym;
  delete modqxvsxmym;
  delete trixy0Histo;
  delete trixyAHisto;
  delete trixycHisto;
  delete sixdxyvsxy;
  delete sixdtvsxyA;
  delete sixdtvsxy;
  delete sixdtvsxmym;
  delete sixxyHisto;
  delete cmsxvsx;
  delete cmsyvsy;
  delete cmsdxvsev1;
  delete cmsdxvsev2;
  delete cmsmadyvsxmym;
  delete trixyAlkHisto;
  delete trixyclkHisto;
  delete cmsnpxvsxmym;
  delete cmsnpxvsxmym2;
  delete cmsxmymq0;
  delete cmsxyq0;
  delete cmsxmymq1;
  delete cmsqxvsxy;
  delete cmsqxvsxmym;
  delete cmsqxvsxmym2;
  delete sixxylkHisto;
  delete sixxyeffHisto;
  delete effnpxvsxmym;
  delete effq0vsxmym;
  delete effvsxy;
  delete effvsxt;
  delete effvsxmym;
  delete effvsxmym2;
  */
  return 0;
}

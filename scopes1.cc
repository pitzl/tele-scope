
// Daniel Pitzl, DESY, Sep 2017
// telescope analysis with eudaq and ROC4Sens and MOD
// one event at a time - saves memory

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
#include <TH2.h>
#include <TH3.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <set> // for hotset
#include <map> // colpx
#include <cmath>
#include <sys/time.h> // gettimeofday, timeval

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  double adc;
  double ph;
  double q;
  int ord;
  bool big;
  double noise;
  bool operator < (const pixel & pxObj ) const
  {
    return row < pxObj.row;
  }
};

struct cluster {
  vector <pixel> vpix; // Armin Burgmeier: list
  int size;
  int ncol, nrow;
  double col, row;
  double charge;
  bool big;
  bool fake = true;
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

struct cpixel // for hot pixel counting
{
  int col;
  int row;
  int cnt;
  bool operator < (const cpixel & pxObj ) const
  {
    return cnt > pxObj.cnt;
  }
};

//------------------------------------------------------------------------------
vector < cluster > getClus( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with local coordinates
  // next-neighbour topological clustering (allows fCluCut-1 empty pixels)

  vector < cluster > v;
  if( pb.size() == 0 ) return v;

  vector <bool> gone( pb.size() ); // initialized to zero

  unsigned seed = 0;

  while( seed < pb.size() ) {

    // start a new cluster

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
      //cout << "GetClus: cluster with non-positive charge" << endl;
    }

    c.charge = sumQ;
    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;

    v.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  // nothing left, return clusters

  return v;
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
  bool syncmod = 0; // re-sync required ?
  bool writeMODalign = 1;
  bool writeDUTalign = 1;
  double forcedphcut = -1; // negative -> not used
  
  for( int i = 1; i < argc; ++i ) {
    
    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event
      
      if( !strcmp( argv[i], "-s" ) )
        syncmod = 1;
      
      if( !strcmp( argv[i], "-m" ) )
        writeMODalign = 0;
      
      if( !strcmp( argv[i], "-u" ) )
        writeDUTalign = 0;
      
      if( !strcmp( argv[i], "-c" ) )
        forcedphcut = atoi( argv[++i] );
        
  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  string geoFileName( "geo.dat" );
  double DUTtilt0 = 19.3;
  double DUTturn0 = 19.3;
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
    string TURN( "turn" );
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

      if( tag == TURN ) {
	tokenizer >> DUTturn0;
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
	<< "  nominal DUT turn " << DUTturn0 << " deg" << endl
	<< "  DUT chip " << chip0 << endl
	<< "  DUT gain file " << gainFileName << endl
	<< "  MOD gain file " << modgainFileName << endl
	<< "  Weibull version " << weib << endl
	;
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
  if( run >= 32032 )  // Apr 2018
    fDTB = 39.997110E6; // 40 MHz DTB clock from ddt and ddtvsdt: fDTB*(1-slope*1E-3)
  if( run >= 33511 )  // Sep 2018
    fDTB = 39.996810E6; // 40 MHz DTB clock from ddt and ddtvsdt: fDTB*(1-slope*1E-3)

  double xbeam = 2.0; // 133i
  double ybeam = 0.0;
  double rbeam = 4.0;
  /* rlq scopes33576.root or rlq scopes33610.root
    .x qx2.C+("cmsqxvsxy");
    qmap->Draw("colz");
    TArc c( -1.5, 0.9, 4.4 ); // circle
    c.SetFillStyle(0);
    c.SetLineWidth(3);
    c.SetX1(-1.5);
    c.SetY1( 0.9);
    c.SetR1(4.4);
    c.SetR2(4.4);
    c.Draw("same");
    TArc m( -1.5, 0.9, 0.1 ); // mid
    m.SetFillStyle(1000);
    m.Draw("same");
  */
  if( run >= 33568 ) { // 174i, 179i at 125 mA
    xbeam =-1.5;
    ybeam = 0.9;
    rbeam = 4.4;
  }
  /* rlq scopes33701.root or rlq scopes33712.root
    .x qx2.C+("cmsqxvsxy");
    qmap->Draw("colz");
    TArc c( -2.2, 0.5, 3.2 ); // circle
    c.SetFillStyle(0);
    c.SetLineWidth(3);
    c.SetX1(-2.2);
    c.SetR1(3.2);
    c.SetR2(3.2);
    c.Draw("same");
    TArc m( -2.2, 0.5, 0.1 ); // mid
    m.SetFillStyle(1000);
    m.Draw("same");
  */
  if( run >= 33626 ) { // 179i 50 mA
    xbeam =-2.2;
    ybeam = 0.5;
    rbeam = 3.2;
  }

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

  int aligniteration = 0;
  double alignx[9];
  double aligny[9];
  double rotx[9];
  double roty[9];

  ostringstream alignFileName; // output string stream

  alignFileName << "align/align_" << run << ".dat";

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

  hotFileName << "hot/hot_" << run << ".dat";

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
  // dut alignment over time
  
  double p0x = 0;
  double p1x = 0;
  double p2x = 0;
  double p3x = 0;
  double p4x = 0;
  double p5x = 0;
  double p6x = 0;
  double p7x = 0;
  double p8x = 0;
  double p9x = 0;
  
  double p0y = 0;
  double p1y = 0;
  double p2y = 0;
  double p3y = 0;
  double p4y = 0;
  double p5y = 0;
  double p6y = 0;
  double p7y = 0;
  double p8y = 0;
  double p9y = 0;
  
  // open infile
  string filename( "alignDUT/inputAlignVsT.dat" );
  ifstream alignvstFile( filename );
  if( alignvstFile.bad() || ! alignvstFile.is_open() ) {
    cout << "Error opening " << filename << endl;
    return 1;
  }
  else{
    cout << "Start to read data from " << filename << endl;
    
    string hash( "#" );
    string RUN( "run" );
    
    while( ! alignvstFile.eof() ) {
      
      string line;
      getline( alignvstFile, line );
      
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
          
          string word;
          double numb;
          do{
            
            getline( alignvstFile, line );
            if( line.empty() ) continue;
            
            stringstream tokenizer2( line );
            tokenizer2 >> word;
            tokenizer2 >> numb;
            
            if( word == "p0x" ) p0x = numb;
            if( word == "p1x" ) p1x = numb;
            if( word == "p2x" ) p2x = numb;
            if( word == "p3x" ) p3x = numb;
            if( word == "p4x" ) p4x = numb;
            if( word == "p5x" ) p5x = numb;
            if( word == "p6x" ) p6x = numb;
            if( word == "p7x" ) p7x = numb;
            if( word == "p8x" ) p8x = numb;
            if( word == "p9x" ) p9x = numb;
            if( word == "p0y" ) p0y = numb;
            if( word == "p1y" ) p1y = numb;
            if( word == "p2y" ) p2y = numb;
            if( word == "p3y" ) p3y = numb;
            if( word == "p4y" ) p4y = numb;
            if( word == "p5y" ) p5y = numb;
            if( word == "p6y" ) p6y = numb;
            if( word == "p7y" ) p7y = numb;
            if( word == "p8y" ) p8y = numb;
            if( word == "p9y" ) p9y = numb;
            
          } while( !( word == "run" ) );
          
        }
      } // wanted run
      
    } // while not EOF
  } // file open
  alignvstFile.close();
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT:

  const double log10 = log(10);
  const double wt = atan(1.0) / 45.0; // pi/180 deg

  double qwid = 1.1; // [ke] for Moyal in 150 um
  if( chip0 >= 119 && chip0 < 138 ) qwid = 0.8; // irrad gain_1
  if( chip0 >= 119 && chip0 <= 138 && run >= 33841 ) qwid = 1.6; // irrad_1 icy
  if( chip0 >= 300 ) qwid = 1.6; // 230 um 3D

  double qxmax = 0.001; // with qwid = 1.1: 7.5 ke cutoff
  if( chip0 >= 119 && chip0 <= 138 ) qxmax = 0.10; // irrad_1 = 1.8 ke cutoff
  if( chip0 >= 119 && chip0 <= 138 && run >= 33841 ) qxmax = 0.03; // irrad_1 icy 5.5 ke cutoff
  if( chip0 == 174 || chip0 == 179 || chip0 == 193 || chip0 == 109 || chip0 == 155 || chip0 == 115 )
    qxmax = 0.10; // irrad = 2.5 ke cutoff
  if( chip0 >= 300 ) qxmax = 0.007; // 3D with qwid = 1.6: 8 ke cutoff

  double pwid = 25; // 193i
    
  if ( run == 34411 ) pwid = 25.148;
  if ( run == 34418 ) pwid = 23.824;
  if ( run == 34419 ) pwid = 22.438;
  if ( run == 34420 ) pwid = 20.674;
  if ( run == 34421 ) pwid = 18.750;
  if ( run == 34422 ) pwid = 15.930;
  if ( run == 34423 ) pwid = 13.176;
  if ( run == 34424 ) pwid = 11.629;
  if ( run == 34458 ) pwid = 36.050;
  if ( run == 34459 ) pwid = 34.603;
  if ( run == 34460 ) pwid = 33.047;
  if ( run == 34461 ) pwid = 30.845;
  if ( run == 34462 ) pwid = 26.436;
  if ( run == 34463 ) pwid = 23.812;
  if ( run == 34464 ) pwid = 22.063;
  if ( run == 34465 ) pwid = 18.907;
  if ( run == 34466 ) pwid = 16.726;
  if ( run == 34467 ) pwid = 14.607;
  if ( run == 34468 ) pwid = 11.229;
  if ( run == 34469 ) pwid =  9.881;
  if ( run == 34470 ) pwid =  8.974;
  if ( run == 34475 ) pwid = 24.921;
  if ( run == 34477 ) pwid = 24.209;
  if ( run == 34478 ) pwid = 22.523;
  if ( run == 34479 ) pwid = 19.785;
  if ( run == 34480 ) pwid = 18.500;
  if ( run == 34483 ) pwid = 16.139;
  if ( run == 34484 ) pwid = 15.335;
  if ( run == 34485 ) pwid = 12.847;
  if ( run == 34486 ) pwid = 11.504;
  if ( run == 34487 ) pwid =  8.458;
  if ( run == 34490 ) pwid = 25.736;
  if ( run == 34491 ) pwid = 22.763;
  if ( run == 34492 ) pwid = 18.027;
  if ( run == 34493 ) pwid = 16.570;
  if ( run == 34494 ) pwid = 14.620;
  if ( run == 34495 ) pwid = 12.490;
  if ( run == 34496 ) pwid = 11.805;
  if ( run == 34497 ) pwid =  8.744;
  if ( run == 34498 ) pwid = 25.838;
  if ( run == 34500 ) pwid = 26.925;
  if ( run == 34501 ) pwid = 25.327;
  if ( run == 34503 ) pwid = 20.991;
  if ( run == 34504 ) pwid = 18.687;
  if ( run == 34505 ) pwid = 15.548;
  if ( run == 34506 ) pwid = 23.988;
  if ( run == 34507 ) pwid = 14.538;
  if ( run == 34508 ) pwid = 12.685;
  if ( run == 34509 ) pwid = 10.677;
  if ( run == 34510 ) pwid =  8.473;
  if ( run == 34511 ) pwid = 26.412;
  if ( run == 35215 ) pwid = 17.664;
  if ( run == 35216 ) pwid = 16.259;
  if ( run == 35217 ) pwid = 14.863;
  if ( run == 35218 ) pwid = 13.446;
  if ( run == 35219 ) pwid = 13.767;
  if ( run == 35220 ) pwid = 12.428;
  if ( run == 35221 ) pwid = 11.930;
  if ( run == 35222 ) pwid = 12.123;
  if ( run == 35223 ) pwid =  9.913;
  if ( run == 35224 ) pwid =  8.747;
  if ( run == 35609 ) pwid = 14.986;
  if ( run == 35610 ) pwid = 14.025;
  if ( run == 35611 ) pwid = 12.257;
  if ( run == 35615 ) pwid = 11.516;
  if ( run == 35617 ) pwid =  9.679;
  if ( run == 35623 ) pwid = 16.610;
  if ( run == 35624 ) pwid =  7.515;
  if ( run == 35625 ) pwid =  6.850;
  if ( run == 35626 ) pwid =  6.568;
  if ( run == 36210 ) pwid = 15.810;
  if ( run == 36211 ) pwid = 14.201;
  if ( run == 36212 ) pwid = 12.503;
  if ( run == 36213 ) pwid = 11.723;
  if ( run == 36214 ) pwid = 10.197;
  if ( run == 36215 ) pwid =  8.771;
  if ( run == 36216 ) pwid =  7.797;
  if ( run == 36217 ) pwid =  7.336;
  if ( run == 36218 ) pwid = 17.522;
  if ( run == 36230 ) pwid = 23.508;
  if ( run == 36231 ) pwid = 21.439;
  if ( run == 36232 ) pwid = 20.524;
  if ( run == 36233 ) pwid = 19.779;
  if ( run == 36234 ) pwid = 14.063;
  if ( run == 36235 ) pwid =  8.540;
  if ( run == 36236 ) pwid = 10.958;
  if ( run == 36237 ) pwid = 17.888;
  if ( run == 36238 ) pwid = 23.439;
  if ( run == 36239 ) pwid = 26.593;
  if ( run == 36240 ) pwid = 27.195;
  if ( run == 36241 ) pwid = 27.662;
  if ( run == 36242 ) pwid = 28.335;
  if ( run == 36252 ) pwid = 18.775;
  if ( run == 36253 ) pwid = 16.545;
  if ( run == 36254 ) pwid = 12.765;
  if ( run == 36255 ) pwid =  9.550;
  if ( run == 36257 ) pwid = 20.541;
  if ( run == 36258 ) pwid = 22.697;
  if ( run == 36259 ) pwid = 23.770;
  if ( run == 36260 ) pwid = 25.752;
  if ( run == 36263 ) pwid = 20.982;
  if ( run == 36264 ) pwid = 20.500;
  if ( run == 36265 ) pwid = 20.180;
  if ( run == 36266 ) pwid = 19.645;
  if ( run == 36267 ) pwid = 15.355;
  if ( run == 36268 ) pwid = 11.097;
  if ( run == 36269 ) pwid = 22.749;
  if ( run == 36270 ) pwid = 24.285;
  if ( run == 36271 ) pwid = 25.129;
  if ( run == 36272 ) pwid = 25.586;
  if ( run == 36276 ) pwid = 23.516;
  if ( run == 36277 ) pwid = 21.975;
  if ( run == 36281 ) pwid = 21.059;
  if ( run == 36282 ) pwid = 20.268;
  if ( run == 36283 ) pwid = 18.391;
  if ( run == 36284 ) pwid = 17.702;
  if ( run == 36285 ) pwid = 16.523;
  if ( run == 36286 ) pwid = 16.402;
  if ( run == 36287 ) pwid = 13.208;
  if ( run == 36288 ) pwid =  9.670;
  if ( run == 36289 ) pwid =  7.156;
  if ( run == 36309 ) pwid = 19.595;
  if ( run == 36312 ) pwid = 18.690;
  if ( run == 36313 ) pwid = 17.160;
  if ( run == 36314 ) pwid = 17.267;
  if ( run == 36315 ) pwid = 15.363;
  if ( run == 36316 ) pwid = 14.491;
  if ( run == 36317 ) pwid = 13.554;
  if ( run == 36318 ) pwid = 11.057;
  if ( run == 36319 ) pwid =  9.076;
  if ( run == 36320 ) pwid =  7.645;
  
  double pxmax = exp(-10./pwid);
  // does this matter? linked and in time es noise free?
  // what about non-irradiated?
  

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
  if( chip0 == 151 ) rot90 = 1;
  if( chip0 == 155 ) rot90 = 1;
  if( chip0 == 156 ) rot90 = 1;
  if( chip0 == 157 ) rot90 = 1;

  if( chip0 == 163 ) rot90 = 1;
  if( chip0 == 164 ) rot90 = 1;
  if( chip0 == 165 ) rot90 = 1;

  if( chip0 == 174 ) rot90 = 1;
  if( chip0 == 179 ) rot90 = 1;
  if( chip0 == 193 ) rot90 = 1;
  if( chip0 == 195 ) rot90 = 1;
  
  if( chip0 == 227 ) rot90 = 1;
  
  if( chip0 == 196 ) rot90 = 1;
  if( chip0 == 194 ) rot90 = 1;
  if( chip0 == 206 ) rot90 = 1;
  if( chip0 == 197 ) rot90 = 1;
  if( chip0 == 207 ) rot90 = 1;
  
  bool fifty = 0;
  if( chip0 == 102 ) fifty = 1;
  if( chip0 == 106 ) fifty = 1;
  if( chip0 == 111 ) fifty = 1;
  if( chip0 == 117 ) fifty = 1;
  if( chip0 == 118 ) fifty = 1;
  if( chip0 == 119 ) fifty = 1; // irr
  if( chip0 == 122 ) fifty = 1; // irr
  if( chip0 == 123 ) fifty = 1; // irr
  if( chip0 == 126 ) fifty = 1; // irr
  if( chip0 == 127 ) fifty = 1; // irr
  if( chip0 == 132 ) fifty = 1; // irr
  if( chip0 == 133 ) fifty = 1; // irr
  if( chip0 == 134 ) fifty = 1; // irr
  if( chip0 == 135 ) fifty = 1; // irr
  if( chip0 == 139 ) fifty = 1;
  if( chip0 == 140 ) fifty = 1;
  if( chip0 == 142 ) fifty = 1;
  if( chip0 == 143 ) fifty = 1;
  if( chip0 == 144 ) fifty = 1;
  if( chip0 == 147 ) fifty = 1;
  if( chip0 == 149 ) fifty = 1;
  if( chip0 == 151 ) fifty = 1;
  if( chip0 == 152 ) fifty = 1; // FDD
  if( chip0 == 155 ) fifty = 1; // FDD
  if( chip0 == 158 ) fifty = 1; // FDD
  if( chip0 == 159 ) fifty = 1;
  if( chip0 == 160 ) fifty = 1;
  if( chip0 == 161 ) fifty = 1; // poly
  if( chip0 == 166 ) fifty = 1;
  if( chip0 == 172 ) fifty = 1;
  if( chip0 == 173 ) fifty = 1;
  if( chip0 == 175 ) fifty = 1;
  if( chip0 == 176 ) fifty = 1;
  if( chip0 == 177 ) fifty = 1;
  if( chip0 == 181 ) fifty = 1;
  if( chip0 == 182 ) fifty = 1;
  if( chip0 == 183 ) fifty = 1;
  if( chip0 == 184 ) fifty = 1;
  if( chip0 == 185 ) fifty = 1;
  if( chip0 == 186 ) fifty = 1;
  if( chip0 == 187 ) fifty = 1;
  if( chip0 == 191 ) fifty = 1;
  if( chip0 == 225 ) fifty = 1;
  if( chip0 == 226 ) fifty = 1;
  if( chip0 == 229 ) fifty = 1;
  if( chip0 == 230 ) fifty = 1;
  if( chip0 >= 300 ) fifty = 1; // 3D
  
  if( chip0 == 198 ) fifty = 1;
  
  double upsignx =  1; // w.r.t. telescope
  double upsigny =  1;

  if( rot90 ) {
    upsignx = -1;
    upsigny =  1;
  }
  if( chip0 == 108 && run >= 32194 ) upsigny = -1; // Apr 2018
  if( chip0 == 109 && run >= 32195 ) upsigny = -1; // Apr 2018
  if( chip0 == 111 && run >= 32214 ) upsigny = -1; // Apr 2018
  if( chip0 == 114 && run >= 32196 ) upsigny = -1; // Apr 2018
  if( chip0 == 115 && run >= 32199 ) upsigny = -1; // Apr 2018
  if( chip0 == 116 && run >= 32200 ) upsigny = -1; // Apr 2018
  if( chip0 == 118 && run >= 32170 ) upsigny = -1; // Apr 2018
  if( chip0 == 120 ) upsigny = -1;
  if( chip0 == 124 ) upsigny = -1;
  if( chip0 == 128 ) upsigny = -1;
  if( chip0 == 130 ) upsigny = -1;
  if( chip0 == 137 ) upsigny = -1;
  if( chip0 == 139 && run >= 32285 ) upsigny = -1; // Apr 2018
  if( chip0 == 142 && run >= 32287 ) upsigny = -1; // Apr 2018
  if( chip0 == 144 && run >= 32289 ) upsigny = -1; // Apr 2018
  if( chip0 == 145 && run >= 32290 ) upsigny = -1; // Apr 2018
  if( chip0 == 146 ) upsigny = -1; // Apr 2018
  if( chip0 == 147 ) upsigny = -1; // Apr 2018
  if( chip0 == 148 && run >= 32292 ) upsigny = -1; // Apr 2018
  if( chip0 == 149 && run >= 32293 ) upsigny = -1; // Apr 2018
  if( chip0 == 150 && run >= 32294 ) upsigny = -1; // Apr 2018
  if( chip0 == 151 && run >= 32295 ) upsigny = -1; // Apr 2018
  if( chip0 == 152 && run >= 32296 ) upsigny = -1; // Apr 2018
  if( chip0 == 153 ) upsigny = -1; // Apr 2018
  if( chip0 == 155 && run >= 32217 ) upsigny = -1; // Apr 2018
  if( chip0 == 156 ) upsigny = -1; // Apr 2018
  if( chip0 == 157 ) upsigny = -1; // Apr 2018
  if( chip0 == 158 && run >= 32307 ) upsigny = -1; // Apr 2018
  if( chip0 == 159 && run >= 32308 ) upsigny = -1; // Apr 2018
  if( chip0 == 160 && run >= 32309 ) upsigny = -1; // Apr 2018
  if( chip0 == 164 ) upsigny = -1; // Apr 2018
  if( chip0 == 165 ) upsigny = -1; // Apr 2018
  if( chip0 == 166 ) upsigny = -1; // Apr 2018
  if( chip0 == 172 ) upsigny = -1; // Apr 2018
  if( chip0 == 173 ) upsigny = -1; // Apr 2018
  if( chip0 == 174 && run < 33000 ) upsigny = -1; // Apr 2018
  if( chip0 == 175 ) upsigny = -1; // Apr 2018
  if( chip0 == 176 ) upsigny = -1; // Apr 2018
  if( chip0 == 179 ) upsigny = -1; // Apr 2018
  if( chip0 == 191 ) upsigny = -1; // Apr 2018
  if( chip0 == 193 ) upsigny = -1; // Apr 2018

  if( ( chip0 == 174 || chip0 == 179 || chip0 == 193 || chip0 == 109 || chip0 == 115 || chip0 == 120 ) &&
      run > 33500 ) { 
    upsignx = 1;
    upsigny = 1; // Sep 2018
  }
  if( chip0 == 133 && run >= 33500 ) { // 50x50 on straight
    upsignx =-1;
    upsigny = 1; // Sep 2018
  }
  if( chip0 == 155 && run >= 33737 ) { // 50x50 on rot90
    upsignx = 1;
    upsigny = 1; // Sep 2018
  }
  if( ( chip0 == 166 || chip0 == 173 || chip0 == 175 || chip0 == 192 ) &&
      run >= 33740 ) { // 50x50 on straight, cable from below
    upsignx =-1;
    upsigny = 1;
  }
  if( run >= 34246 ) { // Oct 2018
    upsignx = 1;
    upsigny = 1;
  }
  if( !rot90 && run >= 34246 ) { // Oct 2018
    upsignx =-1;
    upsigny = 1;
  }
  if( run >= 35215 ) { // Feb 2019
    upsignx = -1;
    upsigny = -1;
  }
  if( run >= 35609 ) { // Mar 2019
    upsignx =  1;
    upsigny =  1;
  }
  if( run >= 36209) { // Apr 2019
    upsignx = -1;
    upsigny =  1;
  }
  if( run >= 36230) { // Mar 2019
    upsignx =  1;
    upsigny =  1;
  }
  if( run >= 36500) { // May 2019
    upsignx = -1;
    upsigny =  1;
  }
  if( run >= 36783) { // July 2019 straight
    upsignx = -1;
    upsigny =  1;
  }
  if( run >= 36924) { // July 2019 rot90
    upsignx =  1;
    upsigny =  1;
  }
  if( run >= 36979) { // July 2019 straight again
    upsignx = -1;
    upsigny =  1;
  }
  if( run >= 37015) { // July 2019 rot90 again
    upsignx =  1;
    upsigny =  1;
  }
  
  cout << "upsignx " << upsignx << endl;
  cout << "upsigny " << upsigny << endl;

  int iDUT = 7;

  // error handling for pixel geometry
  if( fifty ){
    if( nx[iDUT] == 156 && ny[iDUT] == 160 ){
      cout << "fifty flag matches DUT pixel sizes in geometry file." << endl << flush;      
    }
    else{
      cout << endl << "ERROR -- " 
      << "fifty flag mismatches DUT pixel sizes in geometry file. Is this 100x25?" << endl << endl << flush;
      return 0;
    }
  }
  else{
    if( nx[iDUT] == 78 && ny[iDUT] == 320 ){
      cout << "fifty flag matches DUT pixel sizes in geometry file." << endl << flush;      
    }
    else{
      cout << endl << "ERROR -- "
      << "fifty flag mismatches DUT pixel sizes in geometry file. Is this 50x50?" << endl <<endl << flush;
      return 0;
    }
  }
  
  int DUTaligniteration = 0;
  double DUTalignx = 0.0;
  double DUTaligny = 0.0;
  double DUTrot = 0.0;
  double DUTtilt = DUTtilt0; // [deg]
  double DUTturn = DUTturn0; // [deg]
  double DUTz = 0.5*( zz[2] + zz[3] );

  ostringstream DUTalignFileName; // output string stream

  DUTalignFileName << "alignDUT/alignDUT_" << run << ".dat";

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

  if( DUTaligniteration <= 1 ) {
    DUTtilt = DUTtilt0; // from runs.dat
    DUTturn = DUTturn0; // from runs.dat
  }

  double DUTalignx0 = DUTalignx; // at time 0
  double DUTaligny0 = DUTaligny;

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
  //if( run >= 31148 ) ke = 0.039; // chip 111  1004.dat
  if( run >= 31148 ) ke = 0.042; // chip 111 -v0
  if( run >= 31163 ) ke = 0.035; // chip 117  1006.dat
  //if( run >= 31173 ) ke = 0.037; // chip 117
  if( run >= 31173 ) ke = 0.0400; // chip 117 -v0
  if( run >= 31210 ) ke = 0.068; // chip 332  3D 230 um at 17.2
  if( run >= 31237 ) ke = 0.050; // chip 352  3D 230 um at 17.2
  if( run >= 31295 ) ke = 0.040; // chip 143  at 11 ke
  if( run >= 31304 ) ke = 0.0385; // chip 144  at 11 ke
  if( run >= 31309 ) ke = 0.0274; // chip 142  at 11 ke
  //if( run >= 31351 ) ke = 0.025; // chip 142  at 11 ke
  if( run >= 31351 ) ke = 0.0282; // chip 142  at 11 ke -V0
  //if( run >= 31355 ) ke = 0.039; // chip 144  at 11 ke
  if( run >= 31355 ) ke = 0.0434; // chip 144  at 11 ke -V0
  if( run >= 31390 ) ke = 0.031; // chip 102  at 11 ke
  if( run >= 31425 ) ke = 0.036; // chip 144  at 11 ke
  //if( run >= 31473 ) ke = 0.0296; // chip 144  at 11 ke
  if( run >= 31473 ) ke = 0.0319; // chip 144  at 11 ke -V0
  if( run >= 31488 ) ke = 0.035; // chip 160  at 11 ke
  if( run >= 31512 ) ke = 0.068; // chip 332  3D 230 um
  if( run >= 31562 ) ke = 0.050; // chip 352  3D 230 um
  if( run >= 31592 ) ke = 0.031; // chip 152  deep diff 180 um at 11 ke
  if( run >= 31610 ) ke = 0.035; // chip 160  at 11 ke
  //if( run >= 31614 ) ke = 0.030; // chip 152  deep diff 180 um at 11 ke
  //if( run >= 31614 ) ke = 0.0332; // chip 152  deep diff 180 um at 11 ke -v0
  if( run >= 31614 ) ke = 0.0329; // chip 152  deep diff 180 um at 11 ke -v0
  if( run >= 31632 ) ke = 0.035; // irrad default

  if( run >= 32170 ) ke = 0.0313; // 118
  if( run >= 32186 ) ke = 0.0388; // 193
  if( run >= 32189 ) ke = 0.0357; // 179
  if( run >= 32191 ) ke = 0.0393; // 174
  if( run >= 32193 ) ke = 0.0379; // 165
  if( run >= 32194 ) ke = 0.0369; // 108
  if( run >= 32195 ) ke = 0.0371; // 109
  if( run >= 32196 ) ke = 0.0386; // 114
  if( run >= 32197 ) ke = 0.0369; // 108
  if( run >= 32199 ) ke = 0.0377; // 115
  if( run >= 32200 ) ke = 0.0372; // 116
  if( run >= 32201 ) ke = 0.0375; // 146

  if( run >= 32214 ) ke = 0.0344; // 111
  if( run >= 32215 ) ke = 0.0347; // 147
  if( run >= 32219 ) ke = 0.0350; // 155 DD 13.2
  if( run >= 32236 ) ke = 0.0347; // 166 at 11 ke
  if( run >= 32246 ) ke = 0.0333; // 173 at 11 ke
  if( run >= 32247 ) ke = 0.0343; // 175 at 11 ke
  if( run >= 32248 ) ke = 0.0336; // 176 at 11 ke
  if( run >= 32249 ) ke = 0.0334; // 191 at 11 ke

  if( run >= 32268 ) ke = 0.0346; // 118 at 11 ke gain_1
  if( run >= 32269 ) ke = 0.0365; // 111 at 11 ke gain_1
  if( run >= 32270 ) ke = 0.0380; // 147 at 11 ke gain_1
  if( run >= 32271 ) ke = 0.0358; // 155 DD 13.2 gain_1
  if( run >= 32275 ) ke = 0.0338; // 155 DD 13.2 gain_2
  if( run >= 32277 ) ke = 0.0357; // 109 gain_1
  if( run >= 32278 ) ke = 0.0376; // 114 gain_1
  if( run >= 32279 ) ke = 0.0364; // 115 gain_1
  if( run >= 32280 ) ke = 0.0358; // 116 gain_1
  if( run >= 32281 ) ke = 0.0366; // 165 gain_1
  if( run >= 32282 ) ke = 0.0378; // 174 gain_1
  if( run >= 32283 ) ke = 0.0342; // 179 gain_1
  if( run >= 32284 ) ke = 0.0380; // 193 gain_1
  if( run >= 32285 ) ke = 0.0369; // 139 gain_1
  if( run >= 32287 ) ke = 0.0272; // 142 gain_1  noisy
  if( run >= 32288 ) ke = 0.0352; // 108 gain_1
  if( run >= 32289 ) ke = 0.0371; // 144 gain_1
  if( run >= 32290 ) ke = 0.0369; // 145 gain_1
  if( run >= 32291 ) ke = 0.0360; // 146 gain_1
  if( run >= 32292 ) ke = 0.0394; // 148 gain_1
  if( run >= 32293 ) ke = 0.0380; // 149 gain_1
  if( run >= 32294 ) ke = 0.0370; // 150 gain_1
  if( run >= 32295 ) ke = 0.0384; // 151 gain_1
  if( run >= 32298 ) ke = 0.0360; // 152 gain_1 FDD 13.2
  if( run >= 32299 ) ke = 0.0375; // 166 gain_1
  if( run >= 32302 ) ke = 0.0371; // 172 gain_1
  if( run >= 32303 ) ke = 0.0361; // 173 gain_1
  if( run >= 32304 ) ke = 0.0373; // 175 gain_1
  //if( run >= 32305 ) ke = 0.0366; // 176 gain_1
  if( run >= 32305 ) ke = 0.0390; // 176 gain_1 v0
  if( run >= 32306 ) ke = 0.0359; // 192 gain_1
  if( run >= 32307 ) ke = 0.0348; // 158 gain_1 FDD 13.2
  if( run >= 32308 ) ke = 0.0373; // 159 gain_1
  if( run >= 32309 ) ke = 0.0360; // 153 gain_1 FDD 13.2
  if( run >= 32311 ) ke = 0.0384; // 156 gain_1 FDD 13.2
  if( run >= 32312 ) ke = 0.0344; // 157 gain_1 FDD 13.2
  if( run >= 32313 ) ke = 0.0366; // 164 gain_1
  if( run >= 32315 ) ke = 0.0375; // 160 gain_1
  if( run >= 32687 ) ke = 0.0395; // 192 at 11.3 ke
  // Sep 2018: irrad_2, using fresh gains
  if( run >= 33500 ) ke = 0.0367; // default
  if( run >= 33568 ) ke = 0.0378; // 174 gain_1 fresh
  if( run >= 33609 ) ke = 0.0342; // 179 gain_1 fresh
  if( run >= 33664 ) ke = 0.0367; // irrad_1 default for 133i annealed
  if( run >= 33686 ) ke = 0.0378; // 174 gain_1 fresh for irrad_2 annealed
  if( run >= 33692 ) ke = 0.0357; // 109 gain_1
  if( run >= 33695 ) ke = 0.0380; // 193 gain_1
  if( run >= 33708 ) ke = 0.0342; // 179 gain_1 fresh
  if( run >= 33737 ) ke = 0.0358; // 155 DD 13.2 gain_1
  if( run >= 33740 ) ke = 0.0375; // 166 gain_1
  if( run >= 33752 ) ke = 0.0361; // 173 gain_1
  if( run >= 33754 ) ke = 0.0373; // 175 gain_1
  if( run >= 33761 ) ke = 0.0377; // 115
  if( run >= 33794 ) ke = 0.0367; // default
  if( chip0 == 102 && run >= 34345 ) ke = 0.0356;

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

  DUThotFileName.seekp(0);
  DUThotFileName << "hotTrash/hotDUT_" << run << ".datx"; // for filling later
  
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

  MODalignFileName << "alignMOD/alignMOD_" << run << ".dat";

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

  // normal vector on MOD surface:
  // N = ( 0, 0, -1 ) on MOD, towards -z
  // transform into tele system:
  // tilt alpha around x
  // turn omega around y

  const double com = cos( MODturn*wt );
  const double som = sin( MODturn*wt );
  const double cam = cos( MODtilt*wt );
  const double sam = sin( MODtilt*wt );
  const double cfm = cos( MODrot );
  const double sfm = sin( MODrot );

  const double Nxm =-cam*som;
  const double Nym = sam;
  const double Nzm =-cam*com;

  const double normm = cos( MODturn*wt ) * cos( MODtilt*wt ); // length of Nz

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

  if( forcedphcut > 0 )
    rootFileName << "scopes1" << run << "_th" << 10*forcedphcut << ".root";
  else
    rootFileName << "scopes1" << run << ".root";

  TFile* histoFile = new TFile( rootFileName.str(  ).c_str(  ), "RECREATE" );
  //TFile * histoFile = new TFile( rootFileName.str(  ).c_str(  ), "UPDATE" );
  
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
	      110, -5.5, 5.5, -0.05*f, 0.05*f );
  TProfile dridyvsx =
    TProfile( "dridyvsx",
	      "driplet dy vs x;driplet xB [mm];<driplets #Deltay> [mm]",
	      110, -11, 11, -0.05*f, 0.05*f );

  TProfile dridxvstx =
    TProfile( "dridxstx",
	      "driplet dx vs slope x;driplet slope x [rad];<driplets #Deltax> [mm]",
	      60, -0.003, 0.003, -0.05*f, 0.05*f );
  TProfile dridyvsty =
    TProfile( "dridysty",
	      "driplet dy vs slope y;driplet slope y [rad];<driplets #Deltay> [mm]",
	      60, -0.003, 0.003, -0.05*f, 0.05*f );
  TProfile dridxvst3( "dridxvst3",
		      "driplet dx vs time;time [s];<driplet #Deltax> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile dridxvst5( "dridxvst5",
		      "driplet dx vs time;time [s];<driplet #Deltax> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );
  TProfile dridyvst3( "dridyvst3",
		      "driplet dy vs time;time [s];<driplet #Deltay> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile dridyvst5( "dridyvst5",
		      "driplet dy vs time;time [s];<driplet #Deltay> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

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
  
  TProfile modlkvsev =
  TProfile( "modlkvsev",
            "driplet-MOD links vs event;event [k];driplets with MOD links",
            500, 0, 500, -0.5, 1.5 );

  TH1I ndrilkHisto = TH1I( "ndrilk", "driplet - MOD links;driplet - MOD links;events",
			    11, -0.5, 10.5 );

  // DUT:
  /*
  TH1I dutphHisto( "dutph",
		   "DUT pixel PH;pixel pulse height [ADC];DUT pixels",
		   200, -100, 900 );
  */

  TH1I dutaHisto( "duta", "DUT ADC for Vcal 300 mV;pulse height [ADC];pixels", 300, 0, 600 );

  for( int col = 0; col < 155; ++col )
    for( int row = 0; row < 160; ++row ) {
      double t = ( 300 - p0[col][row] ) / p1[col][row]; // from Vcal
      double a = p3[col][row] + p2[col][row] / ( 1 + exp(-t) ); // [ADC]
      dutaHisto.Fill( a );
    }

  // #################################################################################  
  // finn for tsunami and noise analysis and eta telescope calculation
  
  TH1I dutadcAll( "dutadcAll", "DUT PH;ADC-PED [ADC];pixels", 500, -100, 900 );
  TH1I dutadcPix( "dutadcPix", "DUT PH;ADC-PED [ADC];pixels", 500, -100, 900 );
  
  // 1 event display
  TProfile2D * evtdspl = new TProfile2D( "evtdspl", "Display last Event;col;row; #DeltaPH signi.",
                                         //PH <ADC>",
                                         78, 0, 78, 320, 0, 320, -2222, 2222 );
  //155, 0, 155, 160, 0, 160, -2222, 2222 );
  // 1 event telescope hit position
  TProfile2D * hitposi = new TProfile2D( "hitposi", "Telescope hit position;col;row;PH <ADC>", 
                                         78, 0, 78, 320, 0, 320, -2222, 2222 );
  //155, 0, 155, 160, 0, 160, -2222, 2222 );
  
  // 100 event dispalys
  //TProfile2D* evtdspls[100];
  //for( int i = 0; i < 100; i++ ) 
  //  evtdspls[i] = new TProfile2D( Form( "evtdspl%i", i ), Form( "Display Event 401 + %i;col;row;PH <ADC>", i ), 
  //                                      155, 0, 155, 160, 0, 160, -2222, 2222 );
 
  //plots for correlation studies
  // 09 10 11 12 13
  // 14  1  2  3 15
  // 16  4  x  5 17
  // 18  6  7  8 19
  // 20 21 22 23 24
  
  TH1I posxdph ( "posxdph", "Hit Pixel #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );
  TH1I pos2dph ( "pos2dph", "Pos 2 Pixel #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );
  TH1I pos7dph ( "pos7dph", "Pos 7 Pixel #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );
  
  // encoding neighbours
  int corPosx[24] = {-1,+0,+1,-1,+1,-1,+0,+1,
                     -2,-1,+0,+1,+2,-2,+2,-2,+2,-2,+2,-2,-1,+0,+1,+2};
  int corPosy[24] = {+1,+1,+1,+0,+0,-1,-1,-1,
                     +2,+2,+2,+2,+2,+1,+1,+0,+0,-1,-1,-2,-2,-2,-2,-2};
  
  //TGraph* gphCor [24]; // correlations of dph
  //TGraph* gdphCor[24]; // correlations of dph
  //TGraph* gq0Cor [24]; // correlations of q
  //int gphCorPtnr [24];
  //int gdphCorPtnr[24];
  //int gq0CorPtnr [24];
  //TProfile pphCor [24];
  //TProfile pq0Cor [24];
  TProfile pdphCor[24];
  for( int i = 0; i < 24; i++ ){
    // cout << "at i " << i+1 << " look at " << corPosx[i] << " -- " << corPosy[i] << endl << flush;
    //gphCor[i] = new TGraph();
    //gphCor[i]->SetName( Form( "gphCor%i", i+1 ) );
    //gphCor[i]->SetTitle( Form( "PH Correlations %i", i+1 ) );
    //gphCor[i]->GetXaxis()->SetTitle( "PH Center" );
    //gphCor[i]->GetYaxis()->SetTitle( Form( "PH Pos %i", i+1 ) );
    //gphCorPtnr[i]= 0;
    //
    //gdphCor[i] = new TGraph();
    //gdphCor[i]->SetName( Form( "gdphCor%i", i+1 ) );
    //gdphCor[i]->SetTitle( Form( "DPH Correlations %i", i+1 ) );
    //gdphCor[i]->GetXaxis()->SetTitle( "DPH Center" );
    //gdphCor[i]->GetYaxis()->SetTitle( Form( "DPH Pos %i", i+1 ) );
    //gdphCorPtnr[i]= 0;
    //
    //gq0Cor[i] = new TGraph();
    //gq0Cor[i]->SetName( Form( "gq0Cor%i", i+1 ) );
    //gq0Cor[i]->SetTitle( Form( "q0 Correlations %i", i+1 ) );
    //gq0Cor[i]->GetXaxis()->SetTitle( "Q Center" );
    //gq0Cor[i]->GetYaxis()->SetTitle( Form( "Q Pos %i", i+1 ) );
    //gq0CorPtnr[i]= 0;
    //
    //pphCor[i] = TProfile( Form( "pphCor%i", i+1 ),
    //                      Form( "PH Correlations %i;PH Center; PH Pos %i", i+1, i+1 ),
    //                      180, -200, 1600, -200, 1600 );
    //pq0Cor[i] = TProfile( Form( "pq0Cor%i", i+1 ),
    //                      Form( "Q0 Correlations %i;Q Center; Q Pos %i", i+1, i+1 ),
    //                      180, -5, 35, -5, 35 );
    pdphCor[i] = TProfile( Form( "pdphCor%i", i+1 ),
                           Form( "DPH Correlations %i;DPH Center; DPH Pos %i", i+1, i+1 ),
                           180, -200, 1600, -200, 1600 );
  }
  
  TH1I dutcoldist( "dutcoldist", "Distance to next hit column in ROI; ncol; pixels", 1100, -10, 100 );  
  TH1I dutdphCntrl( "dutdphCntrl", "Common mode control plot ;ADC-PED [ADC];pixels", 500, -100, 900 );  
  
  TProfile2D* cntrlvsxmym = new TProfile2D( "cntrlvsxmym","accepted vs xmod, ymod;x track mod 50 [#mum];y track mod     50 [#mum];central pixel <charge> [ke]",
                          50, 0, 100, 50, 0, 100, 0, 20 );
  
  // plots for noise studies
  double avgdph[155][160] = {0}; // storage for dph average calculation
  double rmsdph[155][160] = {0}; // storage for dph rms calcultation
  double cntdph[155][160] = {0}; // storage for number of dph stored for average and rms
  TProfile2D* pavgdph = new TProfile2D( "pavgdph", "Average pulse height if not hit;col;row;<ADC>",
                      //155, 0, 155, 160, 0, 160, -2222, 2222 );
                      77, 0, 154, 80, 0, 160, -2222, 2222 );
  TProfile2D* prmsdph = new TProfile2D( "prmsdph", "RMS pulse height if not hit;col;row;<ADC>", 
                      155, 0, 155, 160, 0, 160, -2222, 2222 );
                      //77, 0, 154, 80, 0, 160, -2222, 2222 );
  TProfile2D* pcntdph = new TProfile2D( "pcntdph", "N entries for avgdph and rmsdph;col;row;<ADC>", 
                      //155, 0, 155, 160, 0, 160, -2222, 2222 );
                      77, 0, 154, 80, 0, 160, -2222, 2222 );
  
  TH1I hrmsdph( "hrmsdph", "Dph RMS histogram ;ADC-PED [ADC];pixels", 500, -10, 90 );
  TH1I hrmsdphEff1( "hrmsdphEff1", "Dph RMS histogram for efficient pixels;ADC-PED [ADC];pixels", 100, -10, 90 );
  TH1I hrmsdphEff0( "hrmsdphEff0", "Dph RMS histogram for inefficient pixels;ADC-PED [ADC];pixels", 100, -10, 90 );
  
  TProfile rmsvseff( "rmsvseff",
                     "Pixel DPH RMS vs Efficiency; DPH RMS; efficiency",
                     100 , 0, 200, -1, 2 );
  
  // eta plots
  TProfile2D* cntrletacol1 = new TProfile2D( "cntrletacol1","accepted eta col 1 vs xmod, ymod;x track mod 50 [#mum];y track mod 50 [#mum];central pixel <charge> [ke]",
                                50, 0, 100, 50, 0, 100, 0, 20 );
  TProfile2D* cntrletacol2 = new TProfile2D( "cntrletacol2","accepted eta col 2vs xmod, ymod;x track mod 50 [#mum];y track mod 50 [#mum];central pixel <charge> [ke]",
                                50, 0, 100, 50, 0, 100, 0, 20 );
  
  // etacol pixel selection grid
  int etacolPosx[6] = {-1,+0,+1,-1,+0,+1};
  int etacolPosy[6] = {+0,+0,+0,+1,+1,+1};
  // 0 1 2
  // 3 4 5
  
  TH1I etaCol2( "etaCol2",
                    "Two Pixel Eta in Columns; #eta;pixel",
                    200, -0.5, 1.5 );
  TH1I etaCol5( "etaCol5",
                    "Five.5 Pixel Eta in Columns; #eta;pixel",
                    200, -0.5, 1.5 );
  TH1I etaCol6( "etaCol6",
                    "Six Pixel Eta in Columns; #eta;pixel",
                    200, -0.5, 1.5 );
  TH1I dphEtaCol( "dphEtaCol", "Eta Col #DeltaPH;#DeltaPH [ADC];pixels", 1000, -100, 900 );
  
  TProfile2D* cntrletarow1 = new TProfile2D( "cntrletarow1","accepted eta row 1 vs xmod, ymod;x track mod 50 [#mum];y track mod 50 [#mum];central pixel <charge> [ke]",
                                             50, 0, 100, 50, 0, 100, 0, 20 );
  TProfile2D* cntrletarow2 = new TProfile2D( "cntrletarow2","accepted eta row 2vs xmod, ymod;x track mod 50 [#mum];y track mod 50 [#mum];central pixel <charge> [ke]",
                                             50, 0, 100, 50, 0, 100, 0, 20 );
  
  // etarow pixel selection grid
  int etarowPosx[6] = {+0,+0,+0,+1,+1,+1};
  int etarowPosy[6] = {-1,+0,+1,-1,+0,+1};
  // 0 3
  // 1 4
  // 2 5
  
  TH1I etaRow2( "etaRow2",
                    "Two Pixel Eta in Rows; #eta;pixel",
                    200, -0.5, 1.5 );
  TH1I etaRow5( "etaRow5",
                    "Five.5 Pixel Eta in Rows; #eta;pixel",
                    200, -0.5, 1.5 );
  TH1I etaRow6( "etaRow6",
                    "Six Pixel Eta in Rows; #eta;pixel",
                    200, -0.5, 1.5 );
  TH1I dphEtaRow( "dphEtaRow", "Eta Row #DeltaPH;#DeltaPH [ADC];pixels", 1000, -100, 900 );
  
  TH1I dutq0allHisto( "dutq0all",
                      "normal pixel charge;normal pixel charge [ke];pixel",
                      500, -5.0, 45.0 ); // 10.75/222 = 0.048423 ==>> -100 -> -5 and 900 -> 45
  
  TH1I dutq0cutHisto( "dutq0cut",
                      "normal pixel charge dph < 10;normal pixel charge [ke];pixel",
                      500, -5.0, 45.0 );
  
  TH1I dutbigdph( "dutbigdph", "Central Single Pixel Hit Seed #DeltaPH Distribution; #DeltaPH signi.; pixels",
                  //#DeltaPH [ADC];pixels",
                  //2000, -1000, 1000 );
                  2000, -100, 100 ); // for signi
  
  TH1I pingpong( "pingpong", "Control Online like selection #DeltaPH Distribution;#DeltaPH [ADC];pixels",
                 2000, -1000, 1000 );
  
  TH1I fandbdifph( "fandbdifph", "Control Online like selection #DeltaPH Distribution;#DeltaPH [ADC];pixels",
                   2000, -1000, 1000 );
  
  TH1I fandbsigph( "fandbsigph", "Control Online like selection #DeltaPH Signi. Distribution;Significance ;pixels",
                   2000, -200, 200 );
  
  TH1I fandbnoiph( "fandbnoiph", "Noise Only Online #DeltaPH Distribution;#DeltaPH [ADC];pixels",
                   2000, -1000, 1000 );
  
  TH1I fandbnoisi( "fandbnoisi", "Noise Only Online #DeltaPH Signifi. Distribution;Significance ;pixels",
                   2000, -200, 200 );
  
  // #################################################################################  
    
  TH1I dutpxq1stHisto( "dutpxq1st",
		       "DUT pixel charge 1st;1st pixel charge [ke];1st pixels",
		       100, 0, 25 );

  TH1I dutpxq2ndHisto( "dutpxq2nd",
		       "DUT pixel charge 2nd;2nd pixel charge [ke];2nd pixels",
		       100, 0, 25 );

  TProfile dutphvsprev( "dutphvsprev", "Tsunami;previous PH [ADC];<PH> [ADC]", 90, -100, 800, -999, 1999 );
  TProfile dutdphvsprev( "dutdphvsprev", "Tsunami;previous #DeltaPH [ADC];<#DeltaPH> [ADC]", 90, -100, 300, -999, 1999 );
  TH1I dutphHisto( "dutph", "DUT PH;ADC-PED [ADC];pixels", 500, -100, 900 );
  TH1I dutdphHisto( "dutdph", "DUT #DeltaPH;#DeltaPH [ADC];pixels", 1000, -100, 900 );
  TH1I dutdpiHisto( "dutdpi", "DUT -#DeltaPH;-#DeltaPH [ADC];pixels", 1000, -100, 900 );
  TProfile dutqvsdph( "dutqvsdph", "gain;#DeltaPH [ADC];q [ke]", 500, -100, 900, -11, 99 );
  TH1I dutpxqHisto( "dutpxq", "DUT pixel charge;pixel charge [ke];all ROI pixels", 100, -2, 8 );

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
                   480, 0, 120 );

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
		      "triplet dx vs time;time [s];<triplet #Deltax> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile tridxvst5( "tridxvst5",
		      "triplet dx vs time;time [s];<triplet #Deltax> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

  TProfile tridyvsx( "tridyvsx",
		     "triplet dy vs x;triplet xB [mm];<triplet #Deltay> [mm]",
		     110, -11, 11, -0.05, 0.05 );
  TProfile tridyvsty( "tridyvsty",
		      "triplet dy vs slope y;triplet slope y [rad];<triplet #Deltay> [mm]",
		      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile tridyvst3( "tridyvst3",
		      "triplet dy vs time;time [s];<triplet #Deltay> [mm] / 10s",
		      300, 0, 3000, -0.05, 0.05 );
  TProfile tridyvst5( "tridyvst5",
		      "triplet dy vs time;time [s];<triplet #Deltay> [mm] / min",
		      1100, 0, 66000, -0.05, 0.05 );

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

  TProfile tridzcvsx( "tridzcvsx",
		   "z vs x;x track [mm];<z of intersect at DUT> [mm]",
		   60, -6, 6, -10, 102 );
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
	      220, -11, 11, -0.1, 0.1 );
  TProfile sixmadxvsx =
    TProfile( "sixmadxvsx",
	      "six MAD x vs x;xB [mm];driplet - triplet MAD #Deltax [mm]",
	      220, -11, 11, 0, 0.1 );
  TProfile sixmadxvsy =
    TProfile( "sixmadxvsy",
	      "six MAD x vs y;yB [mm];driplet - triplet MAD #Deltax [mm]",
	      110, -5.5, 5.5, 0, 0.1 );
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
	      "sixplet dx vs time;time [s];<sixplet #Deltax> [mm] / 10s",
	      300, 0, 3000, -0.05, 0.05 );
  TProfile sixdxvst5 =
    TProfile( "sixdxvst5",
	      "sixplet dx vs time;time [s];<sixplet #Deltax> [mm] / min",
	      1100, 0, 66000, -0.05, 0.05 );

  TProfile sixdyvsx =
    TProfile( "sixdyvsx",
	      "six #Deltay vs x;xB [mm];<driplet - triplet #Deltay> [mm]",
	      200, -10, 10, -0.5, 0.5 );
  TProfile sixdyvsy =
    TProfile( "sixdyvsy",
	      "six #Deltay vs y;yB [mm];<driplet - triplet #Deltay [mm]",
	      110, -5.5, 5.5, -0.1, 0.1 );
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
	      "sixplet dy vs time;time [s];<sixplet #Deltay> [mm] . 10s",
	      300, 0, 3000, -0.05, 0.05 );
  TProfile sixdyvst5 =
    TProfile( "sixdyvst5",
	      "sixplet dy vs time;time [s];<sixplet #Deltay> [mm] / min",
	      1100, 0, 66000, -0.05, 0.05 );
  TProfile sixmadyvsx =
    TProfile( "sixmadyvsx",
	      "six MAD y vs x;xB [mm];driplet - triplet MAD #Deltay [mm]",
	      220, -11, 11, 0, 0.1 );
  TProfile sixmadyvsy =
    TProfile( "sixmadyvsy",
	      "six MAD y vs y;yB [mm];driplet - triplet MAD #Deltay [mm]",
	      110, -5.5, 5.5, 0, 0.1 );
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
		110, -11, 11, 55, -5.5, 5.5, 0, 0.7 );

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

  TProfile sixdtvsxav =
    TProfile( "sixdtvsxav",
	      "driplet - triplet kink_{xy} vs x;x_{avg} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      240, -12, 12, 0, 0.1 );
  TProfile sixdtvsyav =
    TProfile( "sixdtvsyav",
	      "driplet - triplet kink_{xy} vs y;y_{avg} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      120, -6, 6, 0, 0.1 );
  TProfile2D * sixdtvsxyav = new
    TProfile2D( "sixdtvsxyav",
		"driplet - triplet kink_{xy} vs x-y;x_{avg} [mm];y_{avg} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		120, -12, 12, 60, -6, 6, 0, 0.1 );

  TProfile sixdtvsx =
    TProfile( "sixdtvsx",
	      "driplet - triplet kink_{xy} vs x;x_{DUT} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      110, -11, 11, 0, 0.1 );
  TProfile2D * sixdtvsxy = new
    TProfile2D( "sixdtvsxy",
		"driplet - triplet kink_{xy} vs x-y;x_{DUT} [mm];y_{DUT} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		110, -11, 11, 55, -5.5, 5.5, 0, 0.1 );

  TProfile sixdtvsxm =
    TProfile( "sixdtvsxm",
	      "driplet - triplet kink_{xy} vs xmod;track x mod 0.1 [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      50, 0, 0.1, 0, 0.1 );
  TProfile sixdtvsym =
    TProfile( "sixdtvsym",
	      "driplet - triplet kink_{xy} vs ymod;track y mod 0.1 [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      50, 0, 0.1, 0, 0.1 );
  TProfile2D * sixdtvsxmym = new
    TProfile2D( "sixdtvsxmym",
	      "driplet - triplet kink_{xy} vs xmod ymod;track x mod 100 [#mum];track y mod 100 [#mum];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		50, 0, 100, 50, 0, 100, 0, 0.1 );

  TH2I * sixxyHisto = new
    TH2I( "sixxy", "sixplets at z DUT;x [mm];y [mm];sixplets",
	  240, -12, 12, 120, -6, 6 );

  // ROI:

  TH1I dutgHisto( "dutg", "DUT gain;1/gain [mv/ADC];pixels", 200, 0, 2 );

  TProfile roiqvsdph( "roiqvsdph", "gain;#DeltaPH [ADC];q [ke]", 500, -100, 900, -11, 99 );

  TH1I roiq1Histo( "roiq1",
		   "1st ring pixel charge;1st ring pixel charge [ke];1st ring pixels", 180, -5, 40 );
  TH1I roiq2Histo( "roiq2",
		   "2nd ring pixel charge;2nd ring pixel charge [ke];2nd ring pixels", 180, -5, 40 );

  TProfile roiq0vsxm( "roiq0vsxm",
		      "central pixel vs xmod;x track mod 50 [#mum];central pixel <charge> [ke]",
		      50, 0, 50, -5, 20 );
  TProfile roiq1vsxm( "roiq1vsxm",
		      "left pixel vs xmod;x track mod 50 [#mum];left pixel <charge> [ke]",
		      50, 0, 50, -5, 20 );
  TProfile roiq2vsxm( "roiq2vsxm",
		      "right pixel vs xmod;x track mod 50 [#mum];right pixel <charge> [ke]",
		      50, 0, 50, -5, 20 );

  TProfile roiq1vsq0( "roiq1vsq0", "tsunami left;q_{mid} [ke];<q_{left}> [ke]", 25, 0, 25, -99, 99 );
  TProfile roiq2vsq0( "roiq2vsq0", "tsunami right;q_{mid} [ke];<q_{right}> [ke]", 25, 0, 25, -99, 99 );
  TProfile roiq3vsq0( "roiq3vsq0", "tsunami up;q_{mid} [ke];<q_{up}> [ke]", 25, 0, 25, -99, 99 );
  TProfile roiq4vsq0( "roiq4vsq0", "tsunami down;q_{mid} [ke];<q_{down}> [ke]", 25, 0, 25, -99, 99 );

  TH1I roiphHisto( "roiph", "ROI charge;ROI charge [ADC];fiducial ROI", 100, -200, 800 );
  TH1I roip0Histo( "roip0",
		   "central pixel charge;central pixel charge [ADC];fiducial ROI", 100, -200, 800 );
  TH1I roippHisto( "roipp",
		   "central + 2nd pixel charge;central + 2nd pixel charge [ADC];fiducial ROI",
		   100, -200, 800 );
  TH1I roipcHisto( "roipc",
		   "ROI cluster charge;ROI cluster charge [ADC];fiducial ROI", 100, -200, 800 );
  TH1I roipAHisto( "roipA",
 		   "1st ring charge;1st ring charge [ADC];fiducial ROI", 100, -200, 800 );
  TH1I roipBHisto( "roipB",
		   "2nd ring charge;2nd ring charge [ADC];fiducial ROI", 100, -200, 800 );

  TH1I roiqHisto( "roiq", "ROI charge;ROI charge [ke];fiducial ROI", 180, -5, 40 );
  TH1I roiq0Histo( "roiq0",
		   "central pixel charge;central pixel charge [ke];fiducial ROI", 180, -5, 40 );
  TH1I roiq012Histo( "roiq012",
		   "central 3 pixel charge;central 3 pixel charge [ke];fiducial ROI", 180, -5, 40 );
  TH1I roiqqHisto( "roiqq",
		   "central + 2nd pixel charge;central + 2nd pixel charge [ke];fiducial ROI",
		   180, -5, 40 );
  TH1I roiqcHisto( "roiqc", "ROI cluster charge;ROI cluster charge [ke];fiducial ROI", 180, -5, 40 );
  TH1I roiqAHisto( "roiqA",
		   "1st ring charge;1st ring charge [ke];fiducial ROI", 180, -5, 40 );
  TH1I roiqBHisto( "roiqB",
		   "2nd ring charge;2nd ring charge [ke];fiducial ROI", 180, -5, 40 );

  TProfile2D * roiq0vsxmym = new
    TProfile2D( "roiq0vsxmym",
		"central pixel vs xmod, ymod;x track mod 50 [#mum];y track mod 50 [#mum];central pixel <charge> [ke]",
		50, 0, 50, 50, 0, 50, -5, 20 );
  TProfile2D * roiqqvsxmym = new
    TProfile2D( "roiqqvsxmym",
		"central + 2nd pixel vs xmod, ymod;x track mod 50 [#mum];y track mod 50 [#mum];central +2nd pixel <charge> [ke]",
		50, 0, 50, 50, 0, 50, -5, 20 );
  TProfile2D * roiqcvsxmym = new
    TProfile2D( "roiqcvsxmym",
		"ROI cluster vs xmod, ymod;x track mod 50 [#mum];y track mod 50 [#mum];ROI cluster <charge> [ke]",
		50, 0, 50, 50, 0, 50, -5, 20 );
  TProfile2D * roiqAvsxmym = new
    TProfile2D( "roiqAvsxmym",
		"1st ring vs xmod, ymod;x track mod 50 [#mum];y track mod 50 [#mum];1st ring <charge> [ke]",
		50, 0, 50, 50, 0, 50, -5, 20 );
  TProfile2D * roiqBvsxmym = new
    TProfile2D( "roiqBvsxmym",
		"2nd ring vs xmod, ymod;x track mod 50 [#mum];y track mod 50 [#mum];2nd ring <charge> [ke]",
		50, 0, 50, 50, 0, 50, -5, 20 );

  TH1I roipxq1Histo( "roipxq1", "ROI 1st pixel charge;ROI 1st  pixel charge [ke];fiducial ROI", 180, -5, 40 );
  TH1I roipxq2Histo( "roipxq2", "ROI 2nd pixel charge;ROI 2nd  pixel charge [ke];fiducial ROI", 180, -5, 40 );
  TH1I roipxq3Histo( "roipxq3", "ROI 3rd pixel charge;ROI 3rd  pixel charge [ke];fiducial ROI", 180, -5, 40 );

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
			   440, -11, 11 );

  TH2I * cmsxvsx = new TH2I( "cmsxvsx",
			     "DUT vs Telescope x;track x [mm];DUT x [mm];track-cluster combinations",
			     160, -4, 4, 160, -4, 4 );
  TH2I * cmsyvsy = new TH2I( "cmsyvsy",
			     "DUT vs Telescope y;track y [mm];DUT y [mm];track-cluster combinations",
			     160, -4, 4, 160, -4, 4 );

  TH1I cmsxHisto("cmsx", "DUT cluster x-position; x [mm]; clusters",
                 168, -4.2, 4.2);
  TH1I cmsyHisto("cmsy", "DUT cluster y-position; y [mm]; clusters",
                 168, -4.2, 4.2);
  
  TH1I cmsdxHisto( "cmsdx",
		   "DUT - Telescope x;cluster - triplet #Deltax [mm];clusters",
		   200, -0.5, 0.5 );
  TH1I cmsdyHisto( "cmsdy",
		   "DUT - Telescope y;cluster - triplet #Deltay [mm];clusters",
		   200, -0.25, 0.25 );

  TProfile cmsdxvst5 =
    TProfile( "cmsdxvst5",
	      "DUT #Deltax vs time;time [s];<#Deltax> [mm] / min",
	      1100, 0, 66000, -0.2, 0.2 );
  TProfile cmsdxvsev =
    TProfile( "cmsdxvsev",
	      "DUT #Deltax vs time;time [events];<#Deltax> [mm]",
	      1400, 0, 14000*1000, -0.2, 0.2 );

  TH2I * cmsdxvsev1 = new TH2I( "cmsdxvsev1",
				"DUT - Telescope x;event;#Deltax [mm]",
				100, 0, 50000, 100, -5, 5 );
  TH2I * cmsdxvsev2 = new TH2I( "cmsdxvsev2",
				"DUT - Telescope x;event;#Deltax [mm]",
				1000, 0, 1000*1000, 50, -5, 5 );

  TH1I cmsdxcHisto( "cmsdxc",
		    "DUT - Telescope x, cut dy;cluster - triplet #Deltax [mm];clusters",
		    400, -0.2, 0.2 );
  TH1I cmsdxcinfHisto( "cmsdxcinf",
                       "DUT - Telescope x, cut dy, isolated, no fiducial;cluster - triplet #Deltax [mm];isolated clusters",
                       400, -0.2, 0.2 );
  TH1I cmsdxciHisto( "cmsdxci",
		     "DUT - Telescope x, cut dy, isolated;cluster - triplet #Deltax [mm];isolated clusters",
		     400, -0.2, 0.2 );
  TH1I cmsdxcqHisto( "cmsdxcq",
		     "DUT - Telescope x, cut dy, Landau peak;cluster - triplet #Deltax [mm];Landau peak clusters",
		     400, -0.2, 0.2 );

  TProfile cmsdxvsx( "cmsdxvsx",
		     "#Deltax vs x;x track [mm];<cluster - triplet #Deltax> [mm]",
		     50, -3.75, 3.75, -0.2, 0.2 );

  TH2I * cmsdxvsxmq = new
    TH2I( "cmsdxvsxmq",
	  "#Deltax vs xmod;x track mod 50 [#mum];Landau peak cluster - triplet #Deltax [#mum]",
	  50, 0, 50, 100, -50, 50 );

  TProfile cmsdxvsy( "cmsdxvsy",
		     "#Deltax vs y;y track [mm];<cluster - triplet #Deltax> [mm]",
                     76, -3.8, 3.8, -0.2, 0.2 );
  TProfile cmsdxvsyc( "cmsdxvsyc",
                      "#Deltax vs y;y track [mm];<cluster - triplet #Deltax> [mm]",
                      76, -3.8, 3.8, -0.04, 0.04 );
  TProfile cmsdxvstx( "cmsdxvstx",
		      "#Deltax vs #theta_{x};x track slope [mrad];<cluster - triplet #Deltax> [mm]",
		      80, -2, 2, -0.2, 0.2 );
  TProfile cmsmadxvsx =
    TProfile( "cmsmadxvsx",
	      "DUT MAD(#Deltax) vs track x;track x [mm];MAD(#Deltax) [mm]",
              80, -4, 4, 0, 0.1 );
    TProfile cmsmadxvsy =
    TProfile( "cmsmadxvsy",
              "DUT MAD(#Deltax) vs track y;track y [mm];MAD(#Deltax) [mm]",
              80, -4, 4, 0, 0.1 );
  TProfile cmsmadxvsxm =
    TProfile( "cmsmadxvsxm",
	      "DUT x resolution vs xmod;x track mod 50 [#mum];MAD(#Deltax) [#mum]",
	      50, 0, 50, 0, 100 );
  TProfile cmsmadxvsq =
    TProfile( "cmsmadxvsq",
              "DUT MAD(#Deltax) vs Q0;normal cluster charge [ke];MAD(#Deltax) [mm]",
              100, 0, 100, 0, 0.1 );
    
  TH2I * cmsdxvsxm2 = new
    TH2I( "cmsdxvsxm2",
	  "#Deltax vs xmod;x track mod 50 [#mum];2-row cluster - triplet #Deltax [#mum]",
	  50, 0, 50, 100, -50, 50 );

  TH2I * cmsduvsxm2 = new
    TH2I( "cmsduvsxm2",
	  "#Deltau vs xmod;x track mod 50 [#mum];2-row cluster - triplet #Deltau [#mum]",
	  50, 0, 50, 100, -50, 50 );

  TH1I cmsdxc2Histo( "cmsdxc2",
		     "DUT - Telescope x, cut dy, 2-row;cluster - triplet #Deltax [mm];2-row clusters",
		     400, -0.2, 0.2 );
  TH1I cmsduc2Histo( "cmsduc2",
		     "DUT - Telescope x, cut dy, 2-row;cluster - triplet #Deltau [mm];2-row clusters",
		     400, -0.2, 0.2 );

  TProfile cmsdyvsx( "cmsdyvsx",
		     "DUT #Deltay vs x;x track [mm];<cluster - triplet #Deltay> [mm]",
		     50, -3.75, 3.75, -0.2, 0.2 );
  TProfile cmsdyvsxc( "cmsdyvsxc",
                      "DUT #Deltay vs x;x track [mm];<cluster - triplet #Deltay> [mm]",
                      50, -3.75, 3.75, -0.1, 0.1 );
  
  TProfile cmsmadyvsx =
    TProfile( "cmsmadyvsx",
	      "DUT MAD(#Deltay) vs track x;track x [mm];MAD(#Deltay) [mm]",
	      80, -4, 4, 0, 0.1 );
  TProfile cmsmady8vsx =
    TProfile( "cmsmady8vsx",
	      "DUT MAD(#Deltay) vs six track x;track x [mm];MAD(#Deltay) six [mm]",
	      80, -4, 4, 0, 0.1 );

  TProfile cmsdyvsy( "cmsdyvsy",
		     "DUT #Deltay vs y;y track [mm];<cluster - triplet #Deltay> [mm]",
		     76, -3.8, 3.8, -0.2, 0.2 );

  TProfile cmsmadyvsy =
    TProfile( "cmsmadyvsy",
	      "DUT MAD(#Deltay) vs track y;track y [mm];MAD(#Deltay) [mm]",
	      80, -4, 4, 0, 0.1 );

  TH1I cmsdycHisto( "cmsdyc",
		    "#Deltay cut x;cluster - triplet #Deltay [mm];clusters",
		     200, -0.25, 0.25 );
  TH1I cmsdycinfHisto( "cmsdycinf",
                       "#Deltay cut x, isolated, no fiducial;cluster - triplet #Deltay [mm];clusters",
                       400, -0.2, 0.2 );
  TH1I cmsdyciHisto( "cmsdyci",
		     "#Deltay cut x, isolated;cluster - triplet #Deltay [mm];isolated clusters",
		     200, -0.25, 0.25 );
  TH1I cmsdyci6Histo( "cmsdyci6",
		     "#Deltay cut x, isolated six;cluster - triplet #Deltay [mm];isolated clusters",
		     400, -0.2, 0.2 );

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

  TProfile cmsmadyvsq =
    TProfile( "cmsmadyvsq",
	      "DUT MAD(#Deltay) vs Q0;normal cluster charge [ke];MAD(#Deltay) [mm]",
	      100, 0, 100, 0, 0.1 );

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

  TProfile cmsdyvsym =
    TProfile( "cmsdyvsym",
	      "DUT y #Deltay vs ymod;y track mod 100 [#mum];<cluster - triplet #Deltay> [mm]",
	      100, 0, 100, -0.1, 0.1 );
  TProfile cmsdyvsty( "cmsdyvsty",
		      "DUT #Deltay vs #theta_{y};y track slope [mrad];<cluster - triplet #Deltay> [mm]",
		      80, -2, 2, -0.2, 0.2 );

  TProfile cmsdyvst5 =
    TProfile( "cmsdyvst5",
	      "DUT #Deltay vs time;time [s];<#Deltay> [mm] / min",
	      1100, 0, 66000, -0.2, 0.2 );
  TProfile cmsdyvsev =
    TProfile( "cmsdyvsev",
	      "DUT #Deltay vs time;time [events];<#Deltay> [mm]",
	      1400, 0, 14000*1000, -0.2, 0.2 );
  TProfile cmsmadyvsev =
    TProfile( "cmsmadyvsev",
	      "DUT MAD(#Deltay) vs time;time [events];MAD(#Deltay) [mm]",
	      920, 0, 9200*1000, 0, 0.1 );

  TProfile cmsmadyvsty =
    TProfile( "cmsmadyvsty",
	      "DUT MAD(#Deltay) vs track angle;y track angle [mrad];MAD(#Deltay) [mm]",
	      80, -2, 2, 0, 0.1 );

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

  TH2I * cmsdyvsym2 = new
    TH2I( "cmsdyvsym2",
	  "#Deltay vs ymod;y track mod 50 [#mum];2-col cluster - triplet #Deltay [#mum]",
	  50, 0, 50, 100, -50, 50 );

  TH2I * cmsdvvsym2 = new
    TH2I( "cmsdvvsym2",
	  "#Deltav vs ymod;y track mod 50 [#mum];2-col cluster - triplet #Deltav [#mum]",
	  50, 0, 50, 100, -50, 50 );

  TH1I cmsdyc2Histo( "cmsdyc2",
		     "DUT - Telescope y, cut dx, 2-col;cluster - triplet #Deltay [mm];2-col clusters",
		     400, -0.2, 0.2 );
  TH1I cmsdvc2Histo( "cmsdvc2",
		     "DUT - Telescope y, cut dx, 2-col;cluster - triplet #Deltav [mm];2-col clusters",
		     400, -0.2, 0.2 );

  TH1I cmsdycq2Histo( "cmsdycq2",
		      "#Deltay Landau peak;cluster - triplet #Deltay [mm];Landau peak clusters",
		      400, -0.2, 0.2 );

  TH1I cmsadcHisto( "cmsadc",
		    "DUT pixel ADC;pixel pulse height [ADC];linked pixels",
		    400, 0, 800 );
  TH1I cmscolHisto( "cmscol",
		    "DUT linked columns;DUT linked cluster column;linked pixels",
		    nx[iDUT], -0.5, nx[iDUT]-0.5 );
  TH1I cmsrowHisto( "cmsrow",
		    "DUT linked rows;DUT linked cluster row;linked pixels",
		    ny[iDUT], -0.5, ny[iDUT]-0.5 );

  TH1I cmsq0aHisto( "cmsq0a",
		    "normal cluster charge;normal cluster charge [ke];linked clusters",
		    160, 0, 80 );

  TH2I * cmsxmymq0 = new
    TH2I( "cmsxmymq0", "small cluster charge map;x track mod 100 [#mum];y track mod 100 [#mum];small charge tracks",
	  40, 0, 100, 40, 0, 100 );
  TH2I * cmsxyq0 = new
    TH2I( "cmsxyq0", "small cluster charge map;x track [mm];y track [mm];small charge tracks",
          //90, -4.5, 4.5, 90, -4.5, 4.5 );
          360, -4.5, 4.5, 360, -4.5, 4.5 ); // #rb4all
  TH2I * cmsxmymq1 = new
    TH2I( "cmsxmymq1", "normal cluster charge map;x track mod 100 [#mum];y track mod 100 [#mum];normal charge tracks",
	  40, 0, 100, 40, 0, 100 );
  
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
			   "linked triplet at DUT x;triplet x at DUT [mm];linked triplets",
			   240, -12, 12 );
  TH1I triyclkHisto( "triyclk",
			   "linked triplet at DUT y;triplet y at DUT [mm];linked triplets",
			   120, -6, 6 );
  TH2I * trixyclkHisto = new
    TH2I( "trixyclk", "linked triplets x-y at DUT;x at DUT [mm];y at DUT [mm];linked triplets",
	  240, -12, 12, 120, -6, 6 );

  TH1I cmsnpxHisto( "cmsnpx",
		     "linked DUT cluster size;cluster size [pixels];linked clusters",
		     80, 0.5, 80.5 );
  TH1I cmsncolHisto( "cmsncol",
		     "linked DUT cluster size;cluster size [columns];linked clusters",
		     52, 0.5, 52.5 );
  TH1I cmsnrowHisto( "cmsnrow",
		     "linked DUT cluster size;cluster size [rows];linked clusters",
		     80, 0.5, 80.5 );

  TH1I cmsetaHisto( "cmseta",
		     "linked DUT 2-row cluster eta;eta;linked 2-row clusters",
		     100, -1, 1 );
  TH1I cmsetaqHisto( "cmsetaq",
		     "linked Landau peak DUT 2-row cluster eta;eta;linked Landau peak 2-row clusters",
		     100, -1, 1 );

  TH1I cmsetapHisto( "cmsetap",
		     "linked DUT 2-row cluster PH eta;PH eta;linked 2-row clusters",
		     100, -1, 1 );
  TH1I cmsetapqHisto( "cmsetapq",
		     "linked Landau peak DUT 2-row cluster PH eta;PH eta;linked Landau peak 2-row clusters",
		     100, -1, 1 );
  TH2I * cmsetavsxm2 = new
    TH2I( "cmsetavsxm2",
	  "eta vs xmod;x track mod 50 [#mum];2-row cluster eta",
	  50, 0, 50, 100, -1, 1 );
  TProfile cmsetavsym( "cmsetavsym",
			"DUT eta vs ymod;y track mod 100 [#mum];2-row cluster <eta>",
			100, 0, 100, -1, 1 );

  TH1I cmsutaHisto( "cmsuta",
		     "linked DUT 2-col cluster uta;uta;linked 2-col clusters",
		     100, -1, 1 );
  TH1I cmsutaqHisto( "cmsutaq",
		     "linked Landau peak DUT 2-col cluster uta;uta;linked Landau peak 2-col clusters",
		     100, -1, 1 );

  TH1I cmsutapHisto( "cmsutap",
		     "linked DUT 2-col cluster PH uta;PH uta;linked 2-col clusters",
		     100, -1, 1 );
  TH1I cmsutapqHisto( "cmsutapq",
		     "linked Landau peak DUT 2-col cluster PH uta;PH uta;linked Landau peak 2-col clusters",
		     100, -1, 1 );

  TProfile cmsutavsxm( "cmsutavsxm",
			"DUT uta vs xmod;x track mod 100 [#mum];2-col cluster <uta>",
			100, 0, 100, -1, 1 );
  TProfile cmsutavsym( "cmsutavsym",
			"DUT uta vs ymod;y track mod 100 [#mum];2-col cluster <uta>",
			100, 0, 100, -1, 1 );

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
		"DUT cluster size vs xmod ymod;x track mod 200 [#mum];y track mod 100 [#mum];<cluster size> [pixels]",
		100, 0, 200, 50, 0, 100, 0, 20 );

  TProfile cmsncolvsxmq( "cmsncolvsxmq",
			 "DUT cluster size vs xmod;x track mod 100 [#mum];Landau peak <cluster size> [columns]",
			 100, 0, 100, 0, 20 );
  TProfile cmsnrowvsxmq( "cmsnrowvsxmq",
			 "DUT cluster size vs xmod;x track mod 100 [#mum];Landau peak <cluster size> [rows]",
			 100, 0, 100, 0, 80 );

  TProfile cmsncolvsymq( "cmsncolvsymq",
			 "DUT cluster size vs ymod;y track mod 100 [#mum];Landau peak <cluster size> [columns]",
			 100, 0, 100, 0, 20 );
  TProfile cmsnrowvsymq( "cmsnrowvsymq",
			 "DUT cluster size vs ymod;y track mod 100 [#mum];Landau peak <cluster size> [rows]",
			 100, 0, 100, 0, 80 );
  TProfile2D * cmsnpxvsxmymq = new
    TProfile2D( "cmsnpxvsxmymq",
		"DUT cluster size vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];Landau peak <cluster size> [pixels]",
		50, 0, 100, 50, 0, 100, 0, 20 );
  TProfile2D * cmsnpxvsxmym2q = new
    TProfile2D( "cmsnpxvsxmym2q",
		"DUT cluster size vs xmod ymod;x track mod 200 [#mum];y track mod 100 [#mum];Landau peak <cluster size> [pixels]",
		100, 0, 200, 50, 0, 100, 0, 20 );

  TH1I cmsad0Histo( "cmsad0",
		    "normal cluster ADC;normal cluster pulse height [ADC];isolated linked clusters",
		    320, 0, 1600 );
  TH1I cmsph0Histo( "cmsph0",
		    "normal cluster PH;normal cluster pulse height [ADC];isolated linked clusters",
		    320, 0, 1600 );
  TH1I cmsv0Histo( "cmsv0",
		   "normal cluster Vcal;normal cluster pulse height [mV];isolated linked clusters",
		   320, 0, 2400 );
  TH1I cmsq0Histo( "cmsq0",
		   "normal cluster charge;normal cluster charge [ke];isolated linked clusters",
		   480, 0, 120 );
  TH1I dutclq( "dutclq",
                   "dut cluster charge;cluster charge [ke];all clusters",
                    480, 0, 120 );
  TH1I cmsq00Histo( "cmsq00",
		   "normal cluster charge;normal cluster charge [ke];isolated linked clusters",
		   160, 0, 40 );
  TH1I cmsq0dtHisto( "cmsq0dt",
		     "normal cluster charge;normal cluster charge [ke];isolated linked in-time clusters",
		     480, 0, 120 );
  TH1I cmsq0sixHisto( "cmsq0six",
		      "normal cluster charge;normal cluster charge [ke];isolated six-linked clusters",
		      480, 0, 120 );
  TH1I cmsq0bHisto( "cmsq0b",
		   "normal cluster charge beam spot;normal cluster charge [ke];beam spot linked clusters",
		   480, 0, 120 );
  TH1I cmsq0hHisto( "cmsq0h",
		   "normal cluster charge beam halo;normal cluster charge [ke];beam halo linked clusters",
		   480, 0, 120 );
  TH1I cmsq03Histo( "cmsq03",
		   "normal cluster charge, < 4 px;normal cluster charge [ke];linked clusters, < 4 px",
		    480, 0, 120 );
  TH1I cmsq0dHisto( "cmsq0d",
		    "cluster charge bias dot;cluster charge bias dot [ke];linked clusters",
		    480, 0, 120 );
  TH1I cmsq0nHisto( "cmsq0n",
		    "cluster charge no dot;cluster charge no dot [ke];linked clusters",
		    480, 0, 120 );
  TProfile cmsqxvst5 =
    TProfile( "cmsqxvst5",
	      "cluster charge vs time;time [s];<cluster charge> [ke]",
	      1100, 0, 66000, 0, qxmax );

  TH2I * cmslqxyHisto = new
    TH2I( "cmslqxy",
	  "DUT low Q map;track x [mm];track y [mm];low Q clusters",
	  160, -4, 4, 160, -4, 4 );
  TH2I * cmslqxmymHisto = new
    TH2I( "cmslqxmym",
	  "DUT low Q map;track x mod 100 [#mum];track y mod 100 [#mum];low Q clusters",
	  50, 0, 100, 50, 0, 100 );

  TH1I cmslqdxHisto( "cmslqdx",
		   "DUT - Telescope x;cluster - triplet #Deltax [mm];low Q clusters",
		   100, -0.1, 0.1 );
  TH1I cmslqdyHisto( "cmslqdy",
		   "DUT - Telescope y;cluster - triplet #Deltay [mm];low Q clusters",
		   100, -0.1, 0.1 );
  TH1I cmslqt4Histo( "cmslqt4", "event time;event time [s];low Q events / 10 s",
		     600, 0, 6000 );
  TH1I cmslqt5Histo( "cmslqt5", "event time;event time [s];low Q events / 100 s",
		     600, 0, 60000 );
  TH1I cmslqnpxHisto( "cmslqnpx", "DUT pixels;DUT event size [pixels];low Q events",
		     101, -0.5, 100.5 );
  TH1I cmslqnclHisto( "cmslqncl", "DUT clusters;DUT event size [clusters];low Q events",
		     21, -0.5, 20.5 );
  TH1I cmslqntriHisto( "cmslqntri", "triplet tracks;triplet tracks;low Q events",
		     21, -0.5, 20.5 );
  TH1I cmslqndriHisto( "cmslqndri", "driplet tracks;driplet tracks;low Q events",
		     21, -0.5, 20.5 );

  TProfile cmsqxvsx( "cmsqxvsx",
		     "DUT cluster charge vs x;x track [mm];<cluster charge> [ke]",
		     76, -3.8, 3.8, 0, qxmax ); // cutoff at 5 ke
  TProfile cmsqxvsy( "cmsqxvsy",
		     "DUT cluster charge vs y;y track [mm];<cluster charge> [ke]",
		     76, -3.8, 3.8, 0, qxmax );
  TProfile cmsqxvsr( "cmsqxvsr",
		     "DUT cluster charge vs beam;beam distance [mm];<cluster charge> [ke]",
		     80, 0, 8, 0, qxmax );
  TProfile2D * cmsqxvsxy = new
		      TProfile2D( "cmsqxvsxy",
				  "DUT cluster charge vs xy;x track [mm];y track [mm];<cluster charge> [ke]",
                                  //76, -3.8, 3.8, 76, -3.8, 3.8, 0, qxmax );
                                  304, -4.5, 4.5, 304, -4.5, 4.5 ); // #rb4all

  TProfile cmspxvsx( "cmspxvsx",
		     "DUT cluster charge vs x;x track [mm];<cluster charge> [ADC]",
		     76, -3.8, 3.8, 0, pxmax ); // cutoff at 5 ke
  TProfile cmspxvsy( "cmspxvsy",
		     "DUT cluster charge vs y;y track [mm];<cluster charge> [ADC]",
		     76, -3.8, 3.8, 0, pxmax );
  TProfile cmspxvsr( "cmspxvsr",
		     "DUT cluster charge vs beam;beam distance [mm];<cluster charge> [ADC]",
		     80, 0, 8, 0, pxmax );
  TProfile2D * cmspxvsxy = new
		      TProfile2D( "cmspxvsxy",
				  "DUT cluster charge vs xy;x track [mm];y track [mm];<cluster charge> [ADC]",
                                  //76, -3.8, 3.8, 76, -3.8, 3.8, 0, qxmax );
                                  304, -4.5, 4.5, 304, -4.5, 4.5 ); // #rb4all

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
	      "DUT cluster charge vs xmod ymod;x track mod 200 [#mum];y track mod 100 [#mum];<cluster charge> [ke]",
		100, 0, 200, 50, 0, 100, 0, qxmax );

  TH1I cmspxphHisto( "cmspxph",
		    "DUT pixel PH linked;pixel PH [ADC];linked pixels",
		    500, 0, 1000 );
  TH1I cmspxqHisto( "cmspxq",
		    "DUT pixel charge linked;pixel charge [ke];linked pixels",
		    250, 0, 25 );

  TH1I cmssdpHisto( "cmssdp",
		    "DUT seed pixel PH linked;seed pixel PH [ADC];linked clusters",
		    500, 0, 1000 );
  TH1I cmssdqHisto( "cmssdq",
		    "DUT seed pixel charge linked;seed pixel charge [ke];linked clusters",
		    100, 0, 20 );

  TProfile effvsdmin =
    TProfile( "effvsdmin",
	      "DUT efficiency vs driplet isolation;driplet isolation [mm];efficiency",
	      80, 0, 8, -1, 2 );
  TProfile effvstmin =
    TProfile( "effvstmin",
	      "DUT efficiency vs triplet isolation;triplet isolation [mm];efficiency",
	      80, 0, 8, -1, 2 );

  TH2I * sixxylkHisto = new
    TH2I( "sixxylk",
	  "MOD-linked sixplets at z DUT;x [mm];y [mm];MOD-linked sixplets",
	  240, -12, 12, 120, -6, 6 );
  TH2I * sixxyeffHisto = new
    TH2I( "sixxyeff",
	  "MOD-linked sixplets with DUT;x [mm];y [mm];DUT-MOD-linked sixplets",
	  240, -12, 12, 120, -6, 6 );

  TProfile2D * effvsxy = new
    TProfile2D( "effvsxy",
		"DUT efficiency vs xy;x track at DUT [mm];y track at DUT [mm];efficiency",
                //90, -4.5, 4.5, 90, -4.5, 4.5, -1, 2 ); // bin = pix
                //180, -4.5, 4.5, 180, -4.5, 4.5, -1, 2 ); // bin = pix
                360, -4.5, 4.5, 360, -4.5, 4.5 ); // #rb4all
  TProfile2D * effvsxyf = new
  TProfile2D( "effvsxyf",
                "DUT efficiency vs xy;x track at DUT [mm] (fine) ;y track at DUT [mm];efficiency",
              //80, 0, 160, 80, 0, 160, -1, 2 ); // combine 4
              //160, 0, 160, 160, 0, 160, -1, 2 );
              360, -4.5, 4.5, 360, -4.5, 4.5 ); // #rb4all
  
  TProfile2D * effvsxyfid = new
  TProfile2D( "effvsxyfid",
              "DUT fid efficiency vs xy;x track at DUT [mm];y track at DUT [mm];efficiency",
              //90, -4.5, 4.5, 90, -4.5, 4.5, -1, 2 ); // bin = pix
              360, -4.5, 4.5, 360, -4.5, 4.5 ); // #rb4all
  
  // plots for studies differential in fluence
  // vs x and y
  TProfile2D * cmsnoivsxy = new
  TProfile2D( "cmsnoivsxy",
              "DUT Noise vs xy;x track at DUT [mm];y track at DUT [mm]; Pixel PH RMS [ADC]",
              //90, -4.5, 4.5, 90, -4.5, 4.5, 0, 20 ); // bin = pix
              360, -4.5, 4.5, 360, -4.5, 4.5 ); // #rb4all
  
  TProfile2D * cmsclpxvsxy = new
  TProfile2D( "cmsclpxvsxy",
              "Estimator for Cluster PH vs xy;x track at DUT [mm];y track at DUT [mm]; Est. Cluster PH",
              //90, -4.5, 4.5, 90, -4.5, 4.5, 0, pxmax ); // bin = pix
              360, -4.5, 4.5, 360, -4.5, 4.5 ); // #rb4all
  
  TH2D * cmsclphvsxy = new
  TH2D( "cmsclphvsxy",
        "DUT Cluster PH vs xy;x track at DUT [mm];y track at DUT [mm]; Cluster PH [ADC]",
        //90, -4.5, 4.5, 90, -4.5, 4.5 ); // bin = pix
        360, -4.5, 4.5, 360, -4.5, 4.5 ); // #rb4all
  
  TH2D * cmsstonvsxy = new
  TH2D( "cmsstonvsxy",
        "DUT S/N vs xy;x track at DUT [mm];y track at DUT [mm]; Cluster PH / Pixel PH RMS",
        //90, -4.5, 4.5, 90, -4.5, 4.5 ); // bin = pix
        360, -4.5, 4.5, 360, -4.5, 4.5 ); // #rb4all
  
  // 3D histograms
  TH3D* cmsclph3D = new 
  TH3D( "cmsclph3D",
        "Cluster PH histogram vs xy; x track at DUT [mm];y track at DUT [mm]; Cluster PH",
        //90, -4.5, 4.5, 90, -4.5, 4.5, 320, 0, 1600 ); // bin = pix
        360, -4.5, 4.5, 360, -4.5, 4.5, 320, 0, 1600 ); // #rb4all  
  
  TProfile2D * cmsmadxvsxy = new
  TProfile2D( "cmsmadxvsxy",
              "DUT x resolution vs xy;x track at DUT [mm];y track at DUT [mm]; MAD(#Deltax) [#mum]",
              //90, -4.5, 4.5, 90, -4.5, 4.5, 0, 100 ); // bin = pix
              360, -4.5, 4.5, 360, -4.5, 4.5 ); // #rb4all
  
  TProfile2D * cmsmadyvsxy = new
  TProfile2D( "cmsmadyvsxy",
              "DUT y resolution vs xy;x track at DUT [mm];y track at DUT [mm]; MAD(#Deltay) [#mum]",
              //90, -4.5, 4.5, 90, -4.5, 4.5, 0, 100 ); // bin = pix
              360, -4.5, 4.5, 360, -4.5, 4.5 ); // #rb4all
  
  // control histograms
  TH1I * cmsnoih = new
  TH1I( "cmsnoih",
        "DUT Noise Distribution; Pixel PH RMS [ADC]; entries", 
        100, 0, 20 );
  
  TH1I* cmsclpxh = new
  TH1I( "cmsclpxh",
        "Estimator for Cluster PH Distribution; Est. Cluster PH; entries",
        100, 0, pxmax+10 );
  
  TH1I* cmsclphh = new
  TH1I( "cmsclphh",
        "DUT Cluster PH Distribution; Cluster PH [ADC]; entries",
        1700, -100, 1600 );
  
  TH1I * cmsstonh = new
  TH1I( "cmsstonh",
        "DUT S/N Distribution; Cluster PH / Pixel PH RMS",
        250, 0, 100 );
  
  TH1I* cmsmadxh = new
  TH1I( "cmsmadxh",
        "DUT x resolution Distribution; MAD(#Deltax) [#mum]; entries",
        100, -50, 50 );
  
  TH1I* cmsmadyh = new
  TH1I( "cmsmadyh",
        "DUT y resolution Distribution; MAD(#Deltay) [#mum]; entries",
        100, -50, 50  );
  
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
    
    TProfile effvsdr =
    TProfile( "effvsdr",
              "DUT efficiency vs beam spot;beam spot distance [mm];efficiency",
              120, 0, 12, -1, 2 );
    
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

  TH1I effq0Histo =
    TH1I( "effq0",
	  "nearest cluster charge;DUT cluster charge [ke];cluster nearest to track",
	  480, 0, 120 );
  TH1I effqrHisto =
    TH1I( "effqr",
	  "nearest cluster charge, 1 driplet-MOD;DUT cluster charge [ke];mono nearest clusters",
	  480, 0, 120 );
  TH1I effq1Histo =
    TH1I( "effq1",
	  "nearest cluster charge dxy > 0.1;DUT cluster charge [ke];dxy > 0.1 nearest cluster",
	  480, 0, 120 );
  TH1I effq2Histo =
    TH1I( "effq2",
	  "nearest cluster charge dxy > 0.2;DUT cluster charge [ke];dxy > 0.2 nearest cluster",
	  480, 0, 120 );
  TH1I effq3Histo =
    TH1I( "effq3",
	  "nearest cluster charge dxy > 0.3;DUT cluster charge [ke];dxy > 0.3 nearest cluster",
	  480, 0, 120 );
  TH1I effq4Histo =
    TH1I( "effq4",
	  "nearest cluster charge dxy > 0.4;DUT cluster charge [ke];dxy > 0.4 nearest cluster",
	  480, 0, 120 );
  TH1I effq5Histo =
    TH1I( "effq5",
	  "nearest cluster charge dxy > 0.5;DUT cluster charge [ke];dxy > 0.5 nearest cluster",
	  480, 0, 120 );
  TH1I effq6Histo =
    TH1I( "effq6",
	  "nearest cluster charge dxy > 0.6;DUT cluster charge [ke];dxy > 0.6 nearest cluster",
	  480, 0, 120 );
  TH1I effq7Histo =
    TH1I( "effq7",
	  "nearest cluster charge dxy > 0.7;DUT cluster charge [ke];dxy > 0.7 mm nearest cluster",
	  480, 0, 120 );
  TH1I effq8Histo =
    TH1I( "effq8",
	  "nearest cluster charge dxy > 0.8;DUT cluster charge [ke];dxy > 0.8 mm nearest cluster",
	  480, 0, 120 );
  TH1I effq9Histo =
    TH1I( "effq9",
	  "nearest cluster charge dxy > 0.9;DUT cluster charge [ke];dxy > 0.9 mm nearest cluster",
	  480, 0, 120 );

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
		    480, 0, 120 );
  TH1I effq0nHisto( "effq0n",
		    "cluster charge no dot;cluster charge outside bias dot [ke];tracks",
		    480, 0, 120 );

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

  TProfile2D * effvsxmym = new
    TProfile2D( "effvsxmym",
		"DUT efficiency vs xmod ymod;x track mod 100 [#mum];y track mod 100 [#mum];efficiency",
		50, 0, 100, 50, 0, 100, -1, 2 );
  TProfile2D * effvsxmym2 = new
    TProfile2D( "effvsxmym2",
		"DUT efficiency vs xmod ymod;x track mod 200 [#mum];y track mod 100 [#mum];efficiency",
                510, 0, 200, 50, 0, 100, -1, 2 ); 
  TProfile2D * effvsxmym5 = new
    TProfile2D( "effvsxmym5",
                "DUT efficiency vs xmod ymod;x track mod 50 [#mum];y track mod 50 [#mum];efficiency",
                50, 0, 50, 50, 0, 50, -1, 2 );
    
  TProfile effvsxm =
    TProfile( "effvsxm",
	      "DUT efficiency vs xmod;x track mod 0.3 [mm];efficiency",
	      50, 0, 0.1, -1, 2 );
  TProfile effvsym =
    TProfile( "effvsym",
	      "DUT efficiency vs ymod;y track mod 0.2 [mm];efficiency",
	      50, 0, 0.1, -1, 2 );

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
  // fake rate plots
  TH1I allncl = TH1I( "allncl",
		      "number of all clusters; # cluster; events",
		      51, -0.5, 50.5 );
  TH1I fksncl = TH1I( "fksncl",
		      "number of fake clusters; # cluster; events",
		      51, -0.5, 50.5 );
  TH1I truncl = TH1I( "truncl",
		      "number of true clusters; # cluster; events",
		      51, -0.5, 50.5 );
  TH1I allnpx = TH1I( "allnpx",
		      "cluster size for all clusters; # pixel; events",
		      51, -0.5, 50.5 );
  TH1I fksnpx = TH1I( "fksnpx",
		      "cluster size for fake clusters; # pixel; events",
		      51, -0.5, 50.5 );
      
  TProfile2D * fksvsxy = new
    TProfile2D( "fksvsxy",
		"DUT fake rate; cluster x [mm];cluster y DUT [mm]; fake rate",
		//90, -4.5, 4.5, 90, -4.5, 4.5, -1, 2 ); // bin = pix
                360, -4.5, 4.5, 360, -4.5, 4.5 ); // #rb4all
    
  TH1I fksperevnt = TH1I( "fksperevnt",
		      "fake rate distribution; fake rate; events",
		      100, 0, 1.01 );
  
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

  string evFileName = Form( "data/roi%06i.txt", run );
  cout << "try to open  " << evFileName;
  ifstream evFile( evFileName.c_str() );
  if( !evFile ) {
    cout << " : failed " << endl;
    return 2;
  }
  cout << " : succeed " << endl;

  // try to read back online file for pedestal and pixel noise
  cout << "try to open online rootfile " << Form( "data/roi%06i.root", run );
  TFile* frootroi = new TFile( Form( "data/roi%06i.root", run ) );
  if( frootroi->IsOpen() ){
    cout << "Succeed!" << endl;
  }
  else{
    cout << "Online file not found, ABORT!" << endl;
    return 2;
  }
  //TProfile2D* pedxy = (TProfile2D*) frootroi->Get( "pedxy" );
  TProfile2D* pxrms = (TProfile2D*) frootroi->Get( "phCRMS" );
  if( pxrms ){
    cout << "Online file intact" << endl;
  }
  else{
    cout << "Online file corrupt -- simlink other file, ABORT!" << endl;
    return 2;
  }
  histoFile->cd(); // reactivate rootfile for output
    
  string START {"START"};
  string hd;
  while( hd != START ) {
    getline( evFile, hd ); // read one line into string
    cout << "  " << hd << endl;
  }
  bool readnext = 1;

  string F {"F"}; // filled flag
  string E {"E"}; // empty  flag
  string A {"A"}; // added  flag

  int iev = 0;
  int nresync = 0;

  vector < cluster > cl0[9]; // remember from previous event
  vector < cluster > cl1[9]; // remember from previous event

  uint64_t tlutime0 = 0;
  uint64_t prevtlutime = 0;
  uint64_t prevdtbtime = 0;

  bool ldbt = 0;

  timeval tv;
  gettimeofday( &tv, NULL );
  long s0 = tv.tv_sec; // seconds since 1.1.1970
  long u0 = tv.tv_usec; // microseconds

  double zeit1 = 0; // tele
  double zeit2 = 0; // DUT
  double zeit3 = 0; // tracking

  // calculate fake rate
  double totalncl = 0;
  double totalfks = 0;
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  do {

    gettimeofday( &tv, NULL );
    long s1 = tv.tv_sec; // seconds since 1.1.1970
    long u1 = tv.tv_usec; // microseconds

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

    double tludt = (tlutime - prevtlutime) / fTLU; // [s]
    if( iev > 1 && tludt > 1e-6 ) {
      hdttlu.Fill( log(tludt)/log10 );
      hdttluw.Fill( log(tludt)/log10, tludt );
    }
    hdttlu1.Fill( tludt );
    hdttlu2.Fill( tludt );

    prevtlutime = tlutime;

    if( iev < 10 )
      cout << "scopes1 processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "scopes1 processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "scopes1 processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 )
      cout << "scopes1 processing  " << run << "." << iev << "  taken " << evsec << " dt " << tludt << endl;

    //if( tludt > 0.08 ) cout << "  ev " << iev << " dt " << tludt << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    vector < cluster > cl[9];

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldb ) std::cout << "PLANE " << plane.ID();

      // /home/pitzl/eudaq/main/include/eudaq/CMSPixelHelper.hh

      int ipl = plane.ID();

      if( ( ( run > 28000 && run <= 32516 ) // 2017, eudaq 1.7: Mimosa 1..6, DUT 7, REF 8, QAD 9
	    || run >= 32612 ) // Jun 2018
	    && ipl > 0 && ipl < 7 )
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
	int adc = plane.GetPixel(ipix); // ADC 0..255

	// skip hot pixels:

	int ipx = ix*ny[ipl] + iy;
	if( hotset[ipl].count(ipx) ) {
	  if( ldb ) cout << " hot" << flush;
	  continue;
	}

	double q = adc;

	int xm = ix; // remember original
	int ym = iy;

	if( ipl == iMOD ) {

	  int roc = ix / 52; // 0..7
	  int col = ix % 52; // 0..51
	  int row = iy;

	  // leave space for big pixels:

	  ix = 1 + ix + 2*roc; // 1..52 per ROC with big pix
	  if( iy > 79 ) iy += 2; // 0..79, 82..161

	  // flip for upper ROCs into local addresses:

	  if( ym > 79 ) {
	    roc = 15 - roc; // 15..8
	    col = 51 - col; // 51..0
	    row = 159 - ym; // 79..0
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
	px.ph = adc;
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
	    pb[pb.size()-1].ph *= 0.5;
	    pb[pb.size()-1].q *= 0.5;
	    px.ph = 0.5*adc;
	    px.q = 0.5*q;
	    px.big = 1;
	    pb.push_back(px);
	  }

	  if( col == 51 ) {
	    px.col = ix+1; // double
	    px.row = iy;
	    pb[pb.size()-1].ph *= 0.5;
	    pb[pb.size()-1].q *= 0.5;
	    px.ph = 0.5*adc;
	    px.q = 0.5*q;
	    px.big = 1;
	    pb.push_back(px);
	  }

	  if( ym == 79 ) {
	    px.col = ix; // double
	    px.row = 80;
	    pb[pb.size()-1].ph *= 0.5;
	    pb[pb.size()-1].q *= 0.5;
	    px.ph = 0.5*adc;
	    px.q = 0.5*q;
	    px.big = 1;
	    pb.push_back(px);
	  }

	  if( ym == 80 ) {
	    px.col = ix; // double
	    px.row = 81;
	    pb[pb.size()-1].ph *= 0.5;
	    pb[pb.size()-1].q *= 0.5;
	    px.ph = 0.5*adc;
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

    gettimeofday( &tv, NULL );
    long s2 = tv.tv_sec; // seconds since 1.1.1970
    long u2 = tv.tv_usec; // microseconds
    zeit1 += s2 - s1 + ( u2 - u1 ) * 1e-6; // read and cluster

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // read R4S DUT stream:

    string DUTev;
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

    if( run == -31610 && iev > 1972100 )
      //if( run == 31610 && iev > 1202900 )
      cout << "\t" << iev << " TLU " << tludt*1E3
	   << " - " << dtbdt*1e3 << " DTB"
	   << endl; // [ms]

    // large time gap = DTB readout of one data block

    while( iev > 88 && tludt > 0.5 && dtbdt < 0.2*tludt && run != 31166 ) {

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
	//if( ldbt )
	cout << "\t" << iev << " TLU " << tludt*1E3
	     << ", DTB " << dtbdt*1e3
	     << endl; // [ms]
      }
      else { // added event
	if( ldbt )
	  cout << "\t DTB added" << endl;
	break;
      }

      break; // only once

    } // tludt

    vector <pixel> pb; // for clustering

    map < int, set <pixel> > colpxmap; // per col, sorted along row

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

	if( roiss.eof() ) 
	  continue;
	if( ldb )
	  cout << " " << col << " " << row << " " << ph << endl;
	hcol[iDUT].Fill( col+0.5 );
	hrow[iDUT].Fill( row+0.5 );

	//pixel px { col, row, ph, ph, ph, ipx, 0 };
        pixel px { col, row, ph, ph, ph, ipx, 0, 999 }; // include noise
        vpx.push_back(px);
        dutadcAll.Fill( ph );
	++ipx;

      } // roi px

      unsigned sizecut = 999;
      if( vpx.size() > sizecut ) {
        cout << endl << "R4S vpx " << vpx.size() << " cleared in event " << iev;
        vpx.clear();
        filled = A;
      }
      
      // column-wise common mode correction:

      set <pixel> compx[156]; // per column, sorted along row

      for( unsigned ipx = 0; ipx < vpx.size(); ++ipx ) {
	int col = vpx[ipx].col;
	//pixel px { col, vpx[ipx].row, vpx[ipx].adc, vpx[ipx].ph, vpx[ipx].q, vpx[ipx].ord, 0 };
        pixel px { col, vpx[ipx].row, vpx[ipx].adc, vpx[ipx].ph, vpx[ipx].q, vpx[ipx].ord, 0, 999 };
        compx[col].insert(px); // sorted: by row
      }

      // find seed pixel of "one central roi events"
      int colcent = -17;
      int rowcent = -17;
      if( vpx.size() == 49 ){ // one central roi
        colcent = vpx[0].col + 3; // we start readout botom left, this should be central (seed)
        rowcent = vpx[0].row + 3;
        //cout << iev << " found central roi event with center in " << colcent << " " << rowcent << endl;
      }
            
      for( unsigned col = 0; col < 155; ++col ) {

	if( compx[col].size() < 2 ) continue;

	set <pixel> colpx; // per col, sorted along row

	auto px1 = compx[col].begin();
	auto px7 = compx[col].end(); --px7; // last row

	// cut off noisy pixel for simulation studies
	if( run == 32277 ){ // need to change fiducial for eff
          if( col      < 14 ) continue;
          if( px1->row <  2 ) continue;
        }
        	
	int row1 = px1->row;
	//int row7 = px7->row;
        int rowprv = row1;
	
        double ph1 = px1->ph;
	double ph7 = px7->ph;
        double dphprev = 0;
        double phprev = ph1;
        double phprv = ph1;
        
	auto px4 = px1;
        ++px4;
        
        // store seed pixel candidates
        double posseed = 0;
        double negseed = 0;
        
	for( ; px4 != px7; ++px4 ) { // between 1 and 7, exclusively

	  int col4 = px4->col;
	  int row4 = px4->row;
	  double ph4 = px4->ph;

          // reconstruct differential as online
          const double bdiffph = ph4 - phprv;
          ++px4;
          const double fdiffph = px4->ph - ph4;
          //int col5 = px4->col;
          //int row5 = px4->row;
          --px4;
          //cout << "here " << px4->col << " " << col5 << " " << px4->row << " " << row5 << " " << endl;
          
          // correction for x-talk -- switched OFF for nw
          double cx = 0.; // x-talk correction factor
          if( run == 32171 ) cx =  0.146; // c118 50x50  p-stop default -- gain = 2
          if( run == 32197 ) cx =  0.076; // c108 100x25 p-stop default -- gain = 2
          if( run == 32214 ) cx =  0.082; // c111 50x50  p-stop open    -- gain = 2
          if( run == 32246 ) cx =  0.094; // c173 50x50  p-stop open    -- gain = 2
          if( run == 32216 ) cx =  0.097; // c147 50x50  p-stop default -- gain = 2
          if( run == 32268 ) cx =  0.040; // c118 50x50  p-stop default -- gain = 1
          if( run == 32269 ) cx = -0.016; // c111 50x50  p-stop open    -- gain = 1
          if( run == 32288 ) cx = -0.020; // c108 100x25 p-stop default -- gain = 1
          
          // Wed May 29 12:43:22 CEST 2019 -- retry
          if( run >= 36252 && run <= 36260 ) cx = -0.09; // c179 for proton irrad 6.6e15 from run 36260
          if( run >= 36263 && run <= 36273 ) cx = -0.16; // c174 for proton irrad 6.6e15 from run 36272
          if( run >= 36276 && run <= 36301 ) cx = -0.18; // c193 for proton irrad 6.6e15 from run 36276
          if( run >= 36309 && run <= 36321 ) cx = -0.10; // c115 for proton irrad 6.6e15 from run 36309
          
          // Fri Aug  9 16:11:49 CEST 2019
          if( run >= 36786 && run <= 36797 ) cx =  0.36; // c119 for proton irrad 3.3e15 from run 
          if( run >= 36803 && run <= 36829 ) cx = -0.03; // c133 for proton irrad 3.3e15 from run 
          if( run >= 36831 && run <= 36849 ) cx = -0.16; // c126 for proton irrad 3.3e15 from run 
           
          //cx = 0.; // SWITCH
          if ( row4 - rowprv == 1 )
            ph4 = ph4 - cx * phprv;
          phprv = ph4;
          
          // cross hairs -- this is manipulation to get missing columns
          //if( col4 ==  20 || row4 ==  20 ||
          //    col4 ==  21 || row4 ==  21 ||
          //    col4 == 140 || row4 == 140 ||
          //    col4 == 141 || row4 == 141
          //  )
          //ph4 = 0;s
          
          double dph;
          // -> (take whats closer -- works fine but mean is better)
          //if( row4 - row1 < row7 - row4 )
          //  dph = ph4 - ph1;
          //else
          //  dph = ph4 - ph7;
          
          // -> (linear interpolation -- bad idea)
          //dph = ph4 - ((ph7-ph1)/(row7-row1))*(row4-row1)+ph1;
          
          // -> (mean -- good idea )
          dph = ph4 - ( ph1 + ph7 ) / 2;

          pixel px; // sensor px
          
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
	  
	  // fill differential seed for central single seed events
	  if( col4 == colcent && row4 == rowcent ){
            //posseed = bdiffph;
            //negseed = fdiffph;
            // signi instead adc
            posseed = bdiffph / pxrms->GetBinContent( col4+1 , row4+1 );
            negseed = fdiffph / pxrms->GetBinContent( col4+1 , row4+2 );
            //cout << iev <<" posseed " << posseed << " negseed " << negseed << endl;
          }
          //evtdspl->Fill( col4, row4, dph );
          
	  bool cool = 1;
	  int hpx = col4 * 160 + row4;
	  if( hotset[iDUT].count(hpx) ) {
	    if( ldb ) cout << " hot" << flush;
	    cool = 0;
	  }
	  
	  dutadcPix.Fill( px4->adc );
	  
	  if( cool ) {
            dutphHisto.Fill( ph4 );
            dutphvsprev.Fill( phprev, ph4 );
            dutdphvsprev.Fill( dphprev, dph );
            dutdphHisto.Fill( dph );
            dutdpiHisto.Fill( -dph );
          }
          
          rowprv = px4->row;
          phprev = ph4;
          dphprev = dph;
          

	  // r4scal.C

	  double U = ( dph - p3[col4][row4] ) / p2[col4][row4];

	  if( U >= 1 )
	    U = 0.9999999; // avoid overflow

	  double vcal = p0[col4][row4] - p1[col4][row4] * log( (1-U)/U ); // inverse Fermi

	  // subtract Vcal offset:

	  double U0 = -p3[col4][row4] / p2[col4][row4]; // dph = 0
	  double v0 = p0[col4][row4] - p1[col4][row4] * log( (1-U0)/U0 ); // inverse Fermi

	  double q = ke * ( vcal - v0 );

	  // linear gain from Vcal:
	  /*
	  double v = 300; // fresh
	  if( chip0 >= 119 && chip0 <= 138 )
	    v = 200; // irrad
	  double t = ( v - p0[col4][row4] ) / p1[col4][row4]; // from Vcal
	  double a = p3[col4][row4] + p2[col4][row4] / ( 1 + exp(-t) ); // [ADC]
	  double g = v / a; // mv/ADC
	  q = ke * dph * g; // overwrite! cmsdy8cq 5.09 instead 4.94
	  */

	  if( col4 > 19 && cool ) {
	    dutqvsdph.Fill( dph, q );
	  }
	  if( cool )
	    dutpxqHisto.Fill( q );

	  px.adc = ph4;
	  px.ph = dph;
	  px.q = q;
	  px.ord = pb.size(); // readout order, starts at zero
	  px.big = 0;
          // matching dph assuming same noise for all pixels
          px.noise = pxrms->GetBinContent( col4+1 , row4+1 ) * sqrt(3) / 2.; 
          //px.noise = bdiffph / pxrms->GetBinContent( col4+1 , row4+1 ); // significance
          
          
          // fill
          colpx.insert(px); // sorted along col

          dutq0allHisto.Fill( q * norm );
          
	  double dphcut = 12; // gain_1 2017  144 31376 cmsdy8cq 4.98
	  //dphcut = 8; // 144 31376 cmsdy8cq 4.94
	  //dphcut = 16; // 144 31376 cmsdy8cq 5.16

	  if( chip0 == 144 ) dphcut = 10; // for pixelav
	  if( chip0 == 160 ) dphcut = 10; // for RD53A

	  if( run >= 31635 ) dphcut = 24; // gain_2 2018
	  if( run == 32263 ) dphcut = 12; // gain_1 test
	  if( run == 32264 ) dphcut = 12; // gain_1 test
	  if( run >= 32267 && run <= 32273 ) dphcut = 12; // gain_1
	  if( run >= 32277 ) dphcut = 12; // gain_1 2018
	  if( run == 32301 ) dphcut = 24; // gain_2 20 MHz
	  if( run >= 33511 ) dphcut = 16; // irrad Sep 2018

	  if( chip0 == 118 && run <= 32266 ) dphcut = 20; // gain_2
	  if( chip0 == 109 && run >= 32277 ) dphcut = 20; // gain_1, noisy
	  if( chip0 == 139 && run >= 32285 ) dphcut = 16; // gain_1, noisy

	  if( chip0 == 120 ) dphcut = 33; // irrad, gain_2
	  if( chip0 == 120 && run >= 32582 ) dphcut = 24; // irrad, gain_1
	  if( chip0 == 122 ) dphcut = 33; // irrad, gain_2
	  if( chip0 == 123 ) dphcut = 33; // irrad, gain_2
	  if( chip0 == 123 && run >= 32564 ) dphcut = 28; // irrad, gain_1
	  if( chip0 == 126 ) dphcut = 24; // irrad, gain_2, noise
	  if( chip0 == 128 ) dphcut = 33; // irrad, gain_2
	  if( chip0 == 130 ) dphcut = 33; // irrad, gain_2
	  //if( chip0 == 133 ) dphcut = 55; // irrad, gain_2, online 4 sigma
	  if( chip0 == 133 ) dphcut = 33; // irrad, gain_2
	  //if( chip0 == 133 ) dphcut = 22; // irrad, gain_1 March
	  if( chip0 == 133 && run >= 32498 ) dphcut = 28; // irrad, gain_1
	  if( chip0 == 134 ) dphcut = 40; // irrad, gain_2 rms 13
	  if( chip0 == 137 ) dphcut = 33; // irrad, gain_2
	  if( chip0 == 156 ) dphcut = 16; // gain_1 
	  if( chip0 > 300 ) dphcut = 20; // irradiated 3D is noisy
	  //if( chip0 == 197 ) dphcut = 50;
	  
	  // dph cut in sigma
	  double sgmct = 4;
          
          if( forcedphcut > 0 )   // force dphcut (also pass it as variable into scopes)
            sgmct = forcedphcut;
          
          // control plots
          fandbdifph.Fill( bdiffph );
          fandbdifph.Fill( fdiffph );
          fandbsigph.Fill( bdiffph / pxrms->GetBinContent( col4+1 , row4+1 ) );
          if( row4 < 158 )
            fandbsigph.Fill( fdiffph / pxrms->GetBinContent( col4+1 , row4+2 ) );
          
          if(  abs( col4 - colcent ) > 1 
            && abs( row4 - rowcent ) > 1 
            && colcent != -17 
            && rowcent != -17
          ){
            fandbnoiph.Fill( bdiffph );
            fandbnoiph.Fill( fdiffph );
            fandbnoisi.Fill( bdiffph / pxrms->GetBinContent( col4+1 , row4+1 ) );
            if( row4 < 158 )
              fandbnoisi.Fill( fdiffph / pxrms->GetBinContent( col4+1 , row4+2 ) );
          }
          
          
          if( bdiffph > + ( sgmct * pxrms->GetBinContent( col4+1 , row4+1 ) ) )
            pingpong.Fill( bdiffph );
          if( ( fdiffph < - ( sgmct * pxrms->GetBinContent( col4+1 , row4+2 ) ) && row4 < 158 ) )
            pingpong.Fill( fdiffph );
          
          // make a simple ph cut instead  
          // if( dph > 9 ) {  
          //if( dph > 40 ) {  
          
          // make differential forward backwards selection as online
          if(        bdiffph > + ( sgmct * pxrms->GetBinContent( col4+1 , row4+1 ) ) 
            || ( ( fdiffph < - ( sgmct * pxrms->GetBinContent( col4+1 , row4+2 ) ) && row4 < 158 ) )
          ) {
            
            dutq0cutHisto.Fill( q * norm );
            hmap[iDUT]->Fill( col4+0.5, row4+0.5 ); // chip coordinates
            // skip hot pixels:
	    if( cool ) {
              pb.push_back(px);
              //hmap[iDUT]->Fill( col4+0.5, row4+0.5 ); // chip coordinates
            }
            
            double bigger = bdiffph / pxrms->GetBinContent( col4+1 , row4+1 );
            if( ( - fdiffph / pxrms->GetBinContent( col4+1 , row4+2 ) ) > bigger ) 
              bigger = - fdiffph / pxrms->GetBinContent( col4+1 , row4+2 );
            
            //evtdspl->Fill( px.col , px.row, dph );
            evtdspl->Fill( px.col , px.row, bigger );
            
            if( pb.size() > 990 ) {
              cout << "R4S pixel buffer overflow in event " << iev
              << endl;
              break;
            }
            
	  } // dph cut

	} // px4

	colpxmap[col] = colpx;
        
        if( pb.size() > 990 ) {
          cout << "R4S pixel buffer overflow in event " << iev
          << endl;
          break;
        }
        
      } // cols

      if( ldb ) cout << " npx " << pb.size() << endl << flush;

    } // filled

    readnext = 1;
    if( iev > 88 && dtbdt > 0.5 && tludt < 0.2*dtbdt && run != 31166 ) {
      readnext = 0;
      if( ldbt ) cout << "repeat DTB event " << DUTev << endl;
    }
    if( readnext )
      prevdtbtime = dtbtime;

    bool lddt = 1;
    if( fabs(tludt - dtbdt) > 0.02E-3 ) lddt = 0;

    // DUT clustering:

    hnpx[iDUT].Fill( pb.size() );

    cl[iDUT] = getClus(pb);

    hncl[iDUT].Fill( cl[iDUT].size() );

    for( unsigned icl = 0; icl < cl[iDUT].size(); ++ icl ) {

      hsiz[iDUT].Fill( cl[iDUT][icl].size );
      dutclq.Fill(     cl[iDUT][icl].charge );
      
      //hclph.Fill( cl[iDUT][icl].sum );
      //hclmap->Fill( cl[iDUT][icl].col, cl[iDUT][icl].row );

    } // icl

    dutnpxvst2.Fill( evsec, pb.size() );
    dutnclvst2.Fill( evsec, cl[iDUT].size() );

    int DUTyld = 0;
    if( pb.size() ) DUTyld = 1; // no double counting: events with at least one px

    dutyldvst2.Fill( evsec, DUTyld );
    dutyldvst6.Fill( evsec/3600, DUTyld );

    if( ldb ) cout << "  DUT cl " << cl[iDUT].size() << endl << flush;

    gettimeofday( &tv, NULL );
    long s3 = tv.tv_sec; // seconds since 1.1.1970
    long u3 = tv.tv_usec; // microseconds
    zeit2 += s3 - s2 + ( u3 - u2 ) * 1e-6; // read and cluster DUT

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT align vs time: parameters are read from file now
    DUTalignx = DUTalignx0 + p0x + ( p1x + ( p2x + ( p3x + ( p4x + ( p5x + ( p6x + ( p7x + ( p8x + p9x
    * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;
    DUTaligny = DUTaligny0 + p0y + ( p1y + ( p2y + ( p3y + ( p4y + ( p5y + ( p6y + ( p7y + ( p8y + p9y
    * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev ) * iev;
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // syncronisation

    if( ! syncmod )
      for( int ipl = 0; ipl < 8; ++ ipl )
	cl0[ipl] = cl[ipl];
    else
      cl0[iMOD] = cl[iMOD]; // shift all but MOD

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make driplets 3+5-4:

    vector <triplet> driplets;

    double driCut = 0.1; // [mm]

    for( vector<cluster>::iterator cA = cl0[3].begin(); cA != cl0[3].end(); ++cA ) {

      double xA = cA->col*ptchx[3] - alignx[3];
      double yA = cA->row*ptchy[3] - aligny[3];
      double xmid = xA - midx[3];
      double ymid = yA - midy[3];
      xA = xmid - ymid*rotx[3];
      yA = ymid + xmid*roty[3];

      for( vector<cluster>::iterator cC = cl0[5].begin(); cC != cl0[5].end(); ++cC ) {

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

	for( vector<cluster>::iterator cB = cl0[4].begin(); cB != cl0[4].end(); ++cB ) {

	  double xB = cB->col*ptchx[4] - alignx[4];
	  double yB = cB->row*ptchy[4] - aligny[4];
	  double xmid = xB - midx[4];
	  double ymid = yB - midy[4];
	  xB = xmid - ymid*rotx[4];
	  yB = ymid + xmid*roty[4];

	  // interpolate track to B:

	  double dz = zz[4] - avz;
	  double xm = avx + slpx * dz; // driplet at m
	  double ym = avy + slpy * dz;

	  double dxm = xB - xm;
	  double dym = yB - ym;
	  hdridx.Fill( dxm );
	  hdridy.Fill( dym );

	  if( fabs( dym ) < 0.05 ) {
	    hdridxc.Fill( dxm );
	    dridxvsy.Fill( yB, dxm );
	    dridxvstx.Fill( slpx, dxm );
	    dridxvst3.Fill( evsec, dxm );
	    dridxvst5.Fill( evsec, dxm );
	  }

	  if( fabs( dxm ) < 0.05 ) {
	    hdridyc.Fill( dym );
	    dridyvsx.Fill( xB, dym );
	    dridyvsty.Fill( slpy, dym );
	    dridyvst3.Fill( evsec, dym );
	    dridyvst5.Fill( evsec, dym );
	  }

	  // telescope driplet cuts:

	  if( fabs(dxm) > driCut ) continue;
	  if( fabs(dym) > driCut ) continue;

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

      if( run == -31175 && evsec > 35000 ) // iev 1'610'100
	cout << evsec
	     << "  " << iev
	     << "  " << x4
	     << "  " << y4
	     << endl;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // driplets vs MOD clusters:

      for( vector<cluster>::iterator c = cl0[iMOD].begin(); c != cl0[iMOD].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;
	double modx = ( ccol + 0.5 - nx[iMOD]/2 ) * ptchx[iMOD]; // -33..33 mm
	double mody = ( crow + 0.5 - ny[iMOD]/2 ) * ptchy[iMOD]; // -8..8 mm
	double q = c->charge;
	double q0 = q*normm;

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

	if( run == -31175 && evsec > 35000 ) {// iev 1'612'000
	  cout << "\t\t\t\t"
	       << "  " << modx
	       << "  " << mody;
	  if( fabs( moddx ) < xcutMOD &&
	      fabs( moddy ) < ycutMOD )
	    cout << " lk";
	  cout << endl;
	}

      } // MOD

    } // driplets

    modlkvst1.Fill( evsec, nm ); // MOD yield vs time
    modlkvst2.Fill( evsec, nm );
    modlkvst3.Fill( evsec, nm );
    modlkvst6.Fill( evsec/3600, nm );
    modlkvsev.Fill( iev/1000. , nm );
    ndrilkHisto.Fill( ndrilk );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT:

    for( vector<cluster>::iterator c = cl0[iDUT].begin(); c != cl0[iDUT].end(); ++c ) {

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

    for( vector<cluster>::iterator c = cl0[iDUT].begin(); c != cl0[iDUT].end(); ++c ) {

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

	dutadcHisto.Fill( px->ph );
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
    
    if( iev != lev ) 
      evtdspl->Reset();
    
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make triplets 2+0-1:

    vector <triplet> triplets;

    //double triCut = 0.1; // [mm]
    double triCut = 0.05; // [mm] 2.10.2017
    double zscint = -15; // [mm] scint

    for( vector<cluster>::iterator cA = cl0[0].begin(); cA != cl0[0].end(); ++cA ) {

      double xA = cA->col*ptchx[0] - alignx[0];
      double yA = cA->row*ptchy[0] - aligny[0];
      double xmid = xA - midx[0];
      double ymid = yA - midy[0];
      xA = xmid - ymid*rotx[0];
      yA = ymid + xmid*roty[0];

      for( vector<cluster>::iterator cC = cl0[2].begin(); cC != cl0[2].end(); ++cC ) {

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

	for( vector<cluster>::iterator cB = cl0[1].begin(); cB != cl0[1].end(); ++cB ) {

	  double xB = cB->col*ptchx[1] - alignx[1];
	  double yB = cB->row*ptchy[1] - aligny[1];
	  double xmid = xB - midx[1];
	  double ymid = yB - midy[1];
	  xB = xmid - ymid*rotx[1];
	  yB = ymid + xmid*roty[1];

	  // interpolate track to B:

	  double dz = zz[1] - avz;
	  double xm = avx + slpx * dz; // triplet at mid
	  double ym = avy + slpy * dz;

	  double dxm = xB - xm;
	  double dym = yB - ym;
	  htridx.Fill( dxm );
	  htridy.Fill( dym );

	  if( fabs( dym ) < triCut ) {

	    htridxc.Fill( dxm );
	    tridxvsy.Fill( yB, dxm );
	    tridxvstx.Fill( slpx, dxm );
	    tridxvst3.Fill( evsec, dxm );
	    tridxvst5.Fill( evsec, dxm );

	    if(      cB->size == 1 )
	      htridxs1.Fill( dxm ); // 4.2 um
	    else if( cB->size == 2 )
	      htridxs2.Fill( dxm ); // 4.0 um
	    else if( cB->size == 3 )
	      htridxs3.Fill( dxm ); // 3.8 um
	    else if( cB->size == 4 )
	      htridxs4.Fill( dxm ); // 4.3 um
	    else
	      htridxs5.Fill( dxm ); // 3.6 um

	    if(      cB->ncol == 1 )
	      htridxc1.Fill( dxm ); // 4.0 um
	    else if( cB->ncol == 2 )
	      htridxc2.Fill( dxm ); // 4.1 um
	    else if( cB->ncol == 3 )
	      htridxc3.Fill( dxm ); // 3.6 um
	    else if( cB->ncol == 4 )
	      htridxc4.Fill( dxm ); // 3.5 um
	    else
	      htridxc5.Fill( dxm ); // 4.1 um

	  } // dy

	  if( fabs( dxm ) < triCut ) {
	    htridyc.Fill( dym );
	    tridyvsx.Fill( xB, dym );
	    tridyvsty.Fill( slpy, dym );
	    tridyvst3.Fill( evsec, dym );
	    tridyvst5.Fill( evsec, dym );
	  }

	  // telescope triplet cuts:

	  if( fabs(dxm) > triCut ) continue;
	  if( fabs(dym) > triCut ) continue;

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
    // triplets at the DUT:

    double xcut = ptchx[iDUT];
    double ycut = ptchy[iDUT];
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
      double xA = xmA + sxA * zA; // triplet impact point on DUT
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

      tridzcvsx.Fill( xmA, dzc );
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
      double xmod5 = fmod( 9.000 + x4, 0.05 );
      double ymod5 = fmod( 9.000 + y4, 0.05 );

      // DUT row and column from telescope
      double xcol = (x4 + 3.85 ) / 0.05; 
      double yrow = (y4 + 4.00 ) / 0.05;
      if( !fifty ){
        xcol = (x4 + 3.85 ) / 0.025; 
        yrow = (y4 + 4.00 ) / 0.100;
      }
      
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
	    sixdxvst5.Fill( evsec, dx );
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
	    sixdyvst5.Fill( evsec, dy );
	    sixmadyvsy.Fill( y4, fabs(dy) );
	    sixmadyvsty.Fill( syA, fabs(dy) );
	    sixmadyvsdty.Fill( dty, fabs(dy) ); // U-shape
	  }

	} // dx

	// match:

	if( fabs(dx) < sixcut && fabs(dy) < sixcut ) {

	  // average driplet and triplet at DUT:

	  double xa = 0.5 * ( xB + xc );
	  double ya = 0.5 * ( yB + yc );

	  sixxyHisto->Fill( xa, ya );

	  if( driplets[jB].lk ) { // driplet linked to MOD
	    lsixlk = 1;
	    dddmin = driplets[jB].ttdmin;
	    sixdslp = dtxy;
	  }

	  // compare slopes:

	  sixdxyvsxy->Fill( xa, ya, dxy );

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
	  sixdtvsxav.Fill( xa, dtxy );
	  sixdtvsyav.Fill( ya, dtxy );
	  sixdtvsxyav->Fill( xa, ya, dtxy );
	  sixdtvsx.Fill( x4, dtxy );
	  sixdtvsxy->Fill( x4, y4, dtxy );
	  if( fiducial ) {
	    sixdtvsxm.Fill( xmod, dtxy );
	    sixdtvsym.Fill( ymod, dtxy );
	    sixdtvsxmym->Fill( xmod*1E3, ymod*1E3, dtxy ); // gStyle->SetPalette(55) 105 107
	  }

	  // transform into DUT system: (passive)

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

	  if( chip0 == 126 ) { // straight, no Cu behind DUT
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

	  if( chip0 == 134 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	  if( chip0 == 152 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	  if( chip0 == 158 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	  if( chip0 == 159 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	  if( chip0 == 160 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	  if( chip0 == 166 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	  if( chip0 == 172 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	  if( chip0 == 173 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	  if( chip0 == 175 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	  if( chip0 == 176 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	  if( chip0 == 191 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	  if( chip0 == 192 ) { // straight, no Cu behind DUT
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
	  xmod5 = fmod( 9.000 + x4, 0.05 );
	  ymod5 = fmod( 9.000 + y4, 0.05 );

	} // six match

      } // driplets

      if( lsixlk && dddmin < 0.6 ) liso = 0; // require isolation at MOD, if present

      double dxbeam = x4 - xbeam;
      double dybeam = y4 - ybeam;
      double drbeam = sqrt( dxbeam*dxbeam + dybeam*dybeam );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT ROI:

      double fidx = 3.7;
      double fidy = 3.8;
      if( rot90 ) {
	fidx = 3.8;
	fidy = 3.7;
      }

      if( iev > 400 &&
	  filled != A &&
	  lddt &&
	  lsixlk &&
	  liso &&
	  fabs(x4) < fidx &&
	  fabs(y4) < fidy ) { // fiducial x and y

	double sump0 = 0;
	double sumpA = 0;
	double sumpB = 0;
	double p2nd = 0;

	double sumq0 = 0;
	double sumqA = 0;
	double sumqB = 0;
	double q2nd = 0;

	double dmin = 99.9;
	double ph0 = 0;
	double ph1 = 0;
	double ph2 = 0;
	double ph3 = 0;
	double ph4 = 0;
	double ph5 = 0;
	double ph6 = 0;
	double ph7 = 0;
	double ph8 = 0;
	double q0 = 0;
	double q1 = 0;
	double q2 = 0;
	double q3 = 0;
	double q4 = 0;
	double q5 = 0;
	double q6 = 0;
	double q7 = 0;
	double q8 = 0;

	set <double> qset;

	double ptx = ptchx[iDUT];
	double pty = ptchy[iDUT];
	if( rot90 ) {
	  ptx = ptchy[iDUT];
	  pty = ptchx[iDUT];
	}
	
	vector<pixel> hitList; // storing hit pixels (flowing some hit criteria)
	vector<vector<pixel>> coletaList;
        vector<vector<pixel>> rowetaList;
        bool etadb = false;

	for( int col = 0; col < 155; ++col ) { // x

	  if( colpxmap.count(col) < 1 ) continue; // empty

	  for( auto px = colpxmap[col].begin(); px != colpxmap[col].end(); ++px ) {

            // 1 event display
            //evtdspl->Fill( px->col, px->row, px->ph );
            // 100 event dispalys
            //if( iev < 501 ) 
            //  evtdspls[iev-401]->Fill( px->col, px->row, px->adc );
            
            //-------------------------- #tsunamifinn - estimate tsunami correction --------------------------
            //   coments
            // find hit pixels
            // later draw the neighbor pixels pulse hight as function of hit pixels pulse hight
            // this is also helpfull to veto hit pixels for noise analysis
            
            //pixel hitPix { px->col, px->row, px->adc, px->ph, px->q, px->ord, px->big };
            //pixel eptPix { -17, -17, 0, 0, 0, 0, 0}; // empty pixel for initialisation
            pixel hitPix { px->col, px->row, px->adc, px->ph, px->q, px->ord, px->big, 999 };
            pixel eptPix { -17, -17, 0, 0, 0, 0, 0, 999}; // empty pixel for initialisation
            
            // using pixel pulse height to find hit pixels
            // loose and tight comparable in all quantities
            //if( px->adc > 100 ){ // tight threshold
            //if( px->adc >  50 ){ // loose threshold
            //if( px->ph  > 100 ){ // tight threshold
            //if( px->ph  >  50 ){ // loose threshold
            //if( px->q   > 3.0 ){ // tight threshold
            //if( px->q   > 1.5 ){ // loose threshold
            //  hitList.push_back( hitPix );
            //  cntrlvsxmym->Fill( xmod * 1e3, ymod * 1e3, 7 );  
            //}
            
            // using telescope pointing  to find hit pixels
            double xpix = ( px->col + 1.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
            double ypix = ( px->row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT];
            if( rot90 ) {
              xpix = ( px->row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
              ypix = ( px->col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
            }
            double dx = xpix - x4;
            double dy = ypix - y4;
            
            if ( fifty ) { // for 50x50
              if( fabs( dx ) < 0.010 && fabs( dy ) < 0.010 ){ // nptm (new pointing tele medium)
                
                hitList.push_back( hitPix );
                cntrlvsxmym->Fill( xmod * 1e3, ymod * 1e3, 7 );
                posxdph.Fill( hitPix.ph );
                
              }
            }
            else { // for 25x100
              if( fabs( dx ) < 0.0025 && fabs( dy ) < 0.040 ){ // nptm (pointing tele tight)
                
                hitList.push_back( hitPix );
                cntrlvsxmym->Fill( xmod * 1e3, ymod * 1e3, 7 );
                posxdph.Fill( hitPix.ph );
                
              }
            }
            //-------------------------- #tsunamifinn - estimate tsunami correction -------------------------- end
            
            // fill pixel noise plot to telescope xy-map
            cmsnoivsxy->Fill( x4, y4, px->noise );
            cmsnoih->Fill( px->noise );
                        
            //-------------------------- #etafinn - calculate "telescope eta" -----------------------------
            //   coments
            // usually we fix telescope coordinates and rotate the sensor.
            // since i want to calculate eta for columns and rows, the senor needs to stay fixed
            // and i have to rotate the telescope coordinates.
            
            // dutdx in column direction
            // dutdy in row    direction
            double coldx = dx;
            double rowdx = dy;
            if( rot90 ) {
              coldx = dy;
              rowdx = dx;
            }
            
            // fill columns to coletaList ------------------------------------------------------
            // hit one half 
            if( fabs( coldx - 0.0 )           < ptchx[iDUT]/2 && // allows /pm half pitch
                fabs( rowdx - ptchy[iDUT]/4 ) < ptchy[iDUT]/4    // allows lower half off the pixel 
            ){
              
              cntrletacol1->Fill( xmod * 1e3, ymod * 1e3, 7 );
              
              vector<pixel> coleta;
              coleta.push_back( eptPix );
              coleta.push_back( eptPix );
              coleta.push_back( eptPix );
              coleta.push_back( eptPix );
              coleta.push_back( hitPix ); // hit in the lower part, want to store pixels below
              coleta.push_back( eptPix );
              
              coletaList.push_back( coleta );
              
              if( etadb )
                cout << "found lower hitPix at col : row " << hitPix.col << " : " << hitPix.row << endl;
              
            }
            
            // hit the other half
            if( fabs( coldx - 0.0 )           < ptchx[iDUT]/2 && // allows /pm half pitch
                fabs( rowdx + ptchy[iDUT]/4 ) < ptchy[iDUT]/4    // allows higher half off the pixel
            ){
              
              cntrletacol2->Fill( xmod * 1e3, ymod * 1e3, 7 );
              
              vector<pixel> coleta;
              coleta.push_back( eptPix );
              coleta.push_back( hitPix ); // hit in the upper part, want to store pixels above
              coleta.push_back( eptPix );
              coleta.push_back( eptPix );
              coleta.push_back( eptPix );
              coleta.push_back( eptPix );
              
              coletaList.push_back( coleta );
              
              if( etadb )
                cout << "found upper hitPix at col : row " << hitPix.col << " : " << hitPix.row << endl;
              
            }
            // fill columns to coletaList ------------------------------------------------------ end
            
            //fill rows to rowetaList ------------------------------------------------------
            // hit one half
            if( fabs( rowdx - 0.0 )           < ptchy[iDUT]/2 && // allows /pm half pitch
                fabs( coldx - ptchx[iDUT]/4 ) < ptchx[iDUT]/4    // allows left half off the pixel
            ){
              
              cntrletarow1->Fill( xmod * 1e3, ymod * 1e3, 7 );
              
              vector<pixel> roweta;
              roweta.push_back( eptPix );
              roweta.push_back( eptPix );
              roweta.push_back( eptPix );
              roweta.push_back( eptPix );
              roweta.push_back( hitPix ); // hit on the left side, want to store pixels on the left
              roweta.push_back( eptPix );
              
              rowetaList.push_back( roweta );
              
              if( etadb )
                cout << "found left hitPix at evt col : row " << iev << " "
                << hitPix.col << " : " << hitPix.row << endl;
              
            }
            
            // hit the other half
            if( fabs( rowdx - 0.0 )           < ptchy[iDUT]/2 && // allows /pm half pitch
                fabs( coldx + ptchx[iDUT]/4 ) < ptchx[iDUT]/4    // allows right half off the pixel
            ){
              
              cntrletarow2->Fill( xmod * 1e3, ymod * 1e3, 7 );
              
              vector<pixel> roweta;
              roweta.push_back( eptPix );
              roweta.push_back( hitPix ); // hit on the right side, want to store pixels on the right
              roweta.push_back( eptPix );
              roweta.push_back( eptPix );
              roweta.push_back( eptPix );
              roweta.push_back( eptPix );
              
              rowetaList.push_back( roweta );
              
              if( etadb )
                cout << "found right hitPix at evt col : row " << iev << " " 
                << hitPix.col << " : " << hitPix.row << endl;
              
            }
            //fill rows to rowetaList ------------------------------------------------------ end
            //-------------------------- #etafinn - calculate "telescope eta" ----------------------------- end
            
          } //pixel
        } // columns
        
        // for debugging
        if( etadb ) {
          // test: printout col and roweta
          cout << "---------------------------------------" << endl;
          cout << endl;
          for( unsigned int ihl = 0; ihl < coletaList.size(); ihl++){ // coletaList
            
            vector<pixel> coleta = coletaList.at( ihl );
            for( int i = 0; i < 6; i++ )
              cout << "coletaList i:" << i 
              << " col : row " << coleta.at(i).col << " : " << coleta.at(i).row << endl;
          }
          cout << endl;
          for( unsigned int ihl = 0; ihl < rowetaList.size(); ihl++){ // rowetaList
            
            vector<pixel> roweta = rowetaList.at( ihl );
            for( int i = 0; i < 6; i++ )
              cout << "rowetaList i:" << i 
              << " col : row " << roweta.at(i).col << " : " << roweta.at(i).row << endl;
            
          }
          cout << endl;
          cout << "---------------------------------------" << endl;
        }
        
        for( int col = 0; col < 155; ++col ) { // columns
          
          if( colpxmap.count(col) < 1 ) continue; // empty
          
          for( auto px = colpxmap[col].begin(); px != colpxmap[col].end(); ++px ) {
            
            //pixel thisPx { px->col, px->row, px->adc, px->ph, px->q, px->big, px->ord };
            pixel thisPx { px->col, px->row, px->adc, px->ph, px->q, px->big, px->ord, 999 };
            
            // distance to next hit column for common mode control plot -- initialize
            double dHitCol = 500;
            
            // fill plots for tsunami analysis
            for( unsigned int ihl = 0; ihl < hitList.size(); ihl++){ // hitList
              pixel hitPix = hitList.at( ihl );
              
              // compare dph distribution for different positions 
              // relative to hit pixel (-1 for c++ arrays)
              if( px->col == ( hitPix.col + corPosx[2-1] ) && px->row == ( hitPix.row + corPosy[2-1] ) )
                pos2dph.Fill( px->ph );
              if( px->col == ( hitPix.col + corPosx[7-1] ) && px->row == ( hitPix.row + corPosy[7-1] ) )
                pos7dph.Fill( px->ph );
              
              //-------------------------- #tsunamifinn - estimate tsunami correction --------------------------
              // fill plots for tsunami analysis
              
              // 09 10 11 12 13
              // 14  1  2  3 15
              // 16  4  0  5 17
              // 18  6  7  8 19
              // 20 21 22 23 24
              
              // position condition
              for( int i = 0; i < 24; i++ ){
                if( px->col == ( hitPix.col + corPosx[i] ) && px->row == ( hitPix.row + corPosy[i] ) ){
                  
                  // fill graph and increment point counter
                  //gphCor[i]->SetPoint( gphCorPtnr[i], hitPix.adc, px->adc );
                  //gdphCor[i]->SetPoint( gdphCorPtnr[i], hitPix.ph, px->ph );
                  //gq0Cor[i]->SetPoint( gq0CorPtnr[i], hitPix.q*norm, px->q*norm );
                  
                  //pphCor[i]. Fill( hitPix.adc, px->adc );
                  //pq0Cor[i]. Fill( hitPix.q*norm, px->q*norm );
                  pdphCor[i].Fill( hitPix.ph, px->ph );
                  
                  //++gphCorPtnr[i];
                  //++gdphCorPtnr[i];
                  //++gq0CorPtnr[i];
                  
                } // position conition
              } // pixels for scatter plot
              //-------------------------- #tsunamifinn - estimate tsunami correction -------------------------- end
              
              
              // distance to next hit column for common mode control plot -- calculate
              if( dHitCol > fabs( px->col - hitPix.col ) )
                dHitCol = fabs( px->col - hitPix.col );
              
            } // hitList
            
            //-------------------------- #etafinn - calculate "telescope eta" -----------------------------
            // get the other needed pixels for eta telescope calculation
            
            // #coleta #finn #eta
            // 0 1 2
            // 3 4 5
            // etacol pixel selection grid
            // etacolPosx[6] = {-1,+0,+1,-1,+0,+1};
            // etacolPosy[6] = {+0,+0,+0,+1,+1,+1};
            
            for( unsigned int ihl = 0; ihl < coletaList.size(); ihl++){ // coletaList
              
              vector<pixel> coleta = coletaList.at( ihl );
              
              if( coleta.at(4).row != -17 && coleta.at(4).col != -17 ){ // if down
                for( int i = 0; i < 6; i++ ){
                  if( i == 4 ) continue;
                  if( px->col == ( coleta.at(4).col + etacolPosx[i] ) && 
                    px->row == ( coleta.at(4).row + etacolPosy[i] - 1 ) ){ 
                    
                    // fill the other 5 pixels in coletaList.at( ihl )
                    coleta.at(i) = thisPx;
                  
                  if( etadb )  
                    cout << "doloop " << iev << " -- " << i << " col : row " 
                    << px->col << " : " << px->row << endl;
                  
                    } // position condition
                } // pixels for col eta
              } // if down
              
              if( coleta.at(1).row != -17 && coleta.at(1).col != -17 ){ // if up
                for( int i = 0; i < 6; i++ ){
                  if( i == 1 ) continue;
                  if( px->col == ( coleta.at(1).col + etacolPosx[i] ) && 
                    px->row == ( coleta.at(1).row + etacolPosy[i] - 0 ) ){
                    
                    // fill the other 5 pixels in coletaList.at( ihl )
                    coleta.at(i) = thisPx;
                  
                  if( etadb )
                    cout << "uploop " << iev << " -- " << i << " col : row " 
                    << px->col << " : " << px->row << endl;
                  
                    } // position condition
                } // pixels for col eta
              } // if up
              
              // fill back the selected pixels
              coletaList.at( ihl ) = coleta;
              
            } // coletaList
            
            // #roweta #finn #eta
            // 0 3
            // 1 4
            // 2 5
            // etarow pixel selection grid
            // int etarowPosx[6] = {+0,+0,+0,+1,+1,+1};
            // int etarowPosy[6] = {-1,+0,+1,-1,+0,+1};
            
            for( unsigned int ihl = 0; ihl < rowetaList.size(); ihl++){ // rowetaList
              
              vector <pixel> roweta = rowetaList.at( ihl );
              
              if( roweta.at(4).row != -17 && roweta.at(4).col != -17 ){
                for( int i = 0; i < 6; i++ ){
                  if( i == 4 ) continue;
                  if( px->col == ( roweta.at(4).col + etarowPosx[i] - 1 ) && 
                    px->row == ( roweta.at(4).row + etarowPosy[i] ) ){ 
                    
                    // fill the other 5 pixels in rowetaList.at( ihl )
                    roweta.at(i) = thisPx;
                  
                  if( etadb )  
                    cout << "left_loop " << iev << " -- " << i << " col : row " 
                    << px->col << " : " << px->row << endl;
                  
                    } // position condition
                } // pixels for col eta
              } // if left
              
              if( roweta.at(1).row != -17 && roweta.at(1).col != -17 ){
                for( int i = 0; i < 6; i++ ){
                  if( i == 1 ) continue;
                  if( px->col == ( roweta.at(1).col + etarowPosx[i] - 0 ) && 
                    px->row == ( roweta.at(1).row + etarowPosy[i] ) ){
                    
                    // fill the other 5 pixels in rowetaList.at( ihl )
                    roweta.at(i) = thisPx;
                  
                  if( etadb )
                    cout << "rightloop " << iev << " -- " << i << " col : row " 
                    << px->col << " : " << px->row << endl;
                  
                    } // position condition
                } // pixels for row eta
              } // if right
              
              // fill back the selected pixels
              rowetaList.at( ihl ) = roweta;
              
            } // rowetaList
            
            dutcoldist.Fill( dHitCol );
            //-------------------------- #etafinn - calculate "telescope eta" ----------------------------- end
            
            // fill plots for offline noise analysis
            if( dHitCol > 0 ){ // one is to hard for loose selection
              
              // mask for noisy pixel (run 32171)
              //if( 
              //  px->col  >  10 && px->col  < 145 && // fiducial col
              //  px->row  >  10 && px->row  < 150    // fiducial row
              //){  
              
              dutdphCntrl.Fill( px->ph );
              
              avgdph[px->col][px->row] += px->ph;
              rmsdph[px->col][px->row] += px->ph * px->ph;
              cntdph[px->col][px->row] += 1;
              
              //}
            }            
          
	    // match pixel and track:

	    double xpix = ( px->col + 1.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // mm
	    double ypix = ( px->row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT];
	    if( rot90 ) {
	      xpix = ( px->row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	      ypix = ( px->col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	    }
	    double dx = xpix - x4;
	    double dy = ypix - y4;
	    double dd = sqrt( dx*dx + dy*dy );

	    int ring = 9;

	    if(      fabs( dx ) < 0.5 * ptx && // central pixel
		     fabs( dy ) < 0.5 * pty )
	      ring = 0;

	    else if( fabs( dx ) < 1.5 * ptx &&
		     fabs( dy ) < 1.5 * pty )
	      ring = 1;

	    else if( fabs( dx ) < 2.5 * ptx &&
		     fabs( dy ) < 2.5 * pty )
	      ring = 2;

	    // linear gain from Vcal:
	    /*
	    double v = 300; // fresh
	    if( chip0 >= 119 && chip0 <= 138 )
	      v = 200; // irrad
	    double t = ( v - p0[px->col][px->row] ) / p1[px->col][px->row]; // from Vcal
	    double a = p3[px->col][px->row] + p2[px->col][px->row] / ( 1 + exp(-t) ); // [ADC]
	    double g = v / a; // mv/ADC
	    dutgHisto.Fill( g );
	    */

	    double ph = px->ph;
	    //double q = ke * ph * g;
	    double q = px->q;

	    if( ring < 3 ) {
	      qset.insert( -q); // negative, sort increasing
	      roiqvsdph.Fill( ph, q );
	    }

	    if(      ring == 0 ) { // central pixel
	      sump0 += ph;
	      sumq0 += q;
	    }
	    else if( ring == 1 ) {
	      roiq1Histo.Fill( q );
	      sumpA += ph;
	      sumqA += q;
	      if( dd < dmin ) {
		dmin = dd;
		p2nd = ph;
		q2nd = q;
	      }
	    }
	    else if( ring == 2 ) {
	      roiq2Histo.Fill( q );
	      sumpB += ph;
	      sumqB += q;
	    }

	    // 5 3 6
	    // 1 0 2
	    // 7 4 8

	    if(      fabs( dy ) < 0.5 * pty ) { // central row
	      if(      fabs( dx ) < 0.5 * ptx ) {
		ph0 = ph;
		q0 = q;
	      }
	      else if( dx > -1.5 * ptx && dx < -0.5 * ptx ) {
		ph1 = ph;
		q1 = q;
	      }
	      else if( dx >  0.5 * ptx && dx <  1.5 * ptx ) {
		ph2 = ph;
		q2 = q;
	      }
	    }
	    else if( dy > 0.5 * pty && dy < 1.5 * pty) { // top row
	      if(      fabs( dx ) < 0.5 * ptx ) { // mid {
		ph3 = ph;
		q3 = q;
	      }
	      else if( dx > -1.5 * ptx && dx < -0.5 * ptx ) {
		ph5 = ph;
		q5 = q;
	      }
	      else if( dx >  0.5 * ptx && dx <  1.5 * ptx ) {
		ph6 = ph;
		q6 = q;
	      }
	    }
	    else if( dy > -1.5 * pty && dy < -0.5 * pty) { // bot row
	      if(      fabs( dx ) < 0.5 * ptx ) {
		ph4 = ph;
		q4 = q;
	      }
	      else if( dx > -1.5 * ptx && dx < -0.5 * ptx ) {
		ph7 = ph;
		q7 = q;
	      }
	      else if( dx >  0.5 * ptx && dx <  1.5 * ptx ) {
		ph8 = ph;
		q8 = q;
	      }
	    }

	  } // colpx

	} // cols
	
	//-------------------------- #etafinn - calculate "telescope eta" -----------------------------
	// calculate eta and fill the plots
	
	// for debugging
	if( etadb ) {
          // test: printout col and roweta
          cout << "---------------------------------------" << endl;
          cout << endl;
          for( unsigned int ihl = 0; ihl < coletaList.size(); ihl++){ // coletaList
            
            vector<pixel> coleta = coletaList.at( ihl );
            for( int i = 0; i < 6; i++ )
              cout << "coletaList i:" << i 
              << " col : row " << coleta.at(i).col << " : " << coleta.at(i).row << endl;
          }
          cout << endl;
          for( unsigned int ihl = 0; ihl < rowetaList.size(); ihl++){ // rowetaList
            
            vector<pixel> roweta = rowetaList.at( ihl );
            for( int i = 0; i < 6; i++ )
              cout << "rowetaList i:" << i 
              << " col : row " << roweta.at(i).col << " : " << roweta.at(i).row << endl;
            
          }
          cout << endl;
          cout << "---------------------------------------" << endl;
          cout << endl;
          cout << endl;
        }
        
        // calculate and fill
        
        for( unsigned int ihl = 0; ihl < coletaList.size(); ihl++){ // coletaList
          
          if( etadb )
            cout << "this is coleta hit nr " << ihl << endl;
          
          vector<pixel> coleta = coletaList.at( ihl );
          
          bool allFill = true; // not always all needed pixels in roi
          for(int i = 0; i < 6; i++){ // pixels in coleta
            
            if( etadb )
              cout << "i | row : col | ph ---- " << i << " | " 
              << coleta.at(i).row << " : " 
              << coleta.at(i).col << " | " 
              << coleta.at(i).ph  << endl; 
            
            if( coleta[i].row == -17 ) allFill = false; // 
            
          }
          
          // central only
          double loval = coleta.at(1).ph;
          double upval = coleta.at(4).ph;
          etaCol2.Fill( ( upval ) / ( loval + upval ) );
          
          // add neighbours
          loval += coleta.at(0).ph;
          loval += coleta.at(2).ph;
          upval += coleta.at(3).ph;
          upval += coleta.at(5).ph;
          etaCol5.Fill( ( upval ) / ( loval + upval ) );
          
          if( allFill )
            etaCol6.Fill( ( upval ) / ( loval + upval ) );
          
          if( etadb )
            cout << "loval, upval, eta -- " << loval << ", "
            << upval << ", "
            << ( upval ) / ( loval + upval ) << endl;
          
          dphEtaCol.Fill( loval + upval ); // charge control plot
          
        } // coletaList
        
        for( unsigned int ihl = 0; ihl < rowetaList.size(); ihl++){ // rowetaList
          
          if( etadb )
            cout << "this is roweta hit nr " << ihl << endl;
          
          vector<pixel> roweta = rowetaList.at( ihl );
          
          bool allFill = true; // not always all needed pixels in roi
          for(int i = 0; i < 6; i++){ // pixels in roweta
            
            if( etadb )
              cout << "i | row : col | ph ---- " << i << " | " 
              << roweta.at(i).row << " : " 
              << roweta.at(i).col << " | " 
              << roweta.at(i).ph  << endl; 
            
            if( roweta[i].row == -17 ) allFill = false; // 
            
          }
          
          // central only
          double loval = roweta.at(1).ph;
          double upval = roweta.at(4).ph;
          etaRow2.Fill( ( upval ) / ( loval + upval ) );
          
          // add neighbours
          loval += roweta.at(0).ph;
          loval += roweta.at(2).ph;
          upval += roweta.at(3).ph;
          upval += roweta.at(5).ph;
          etaRow5.Fill( ( upval ) / ( loval + upval ) );
          
          if( allFill )
            etaRow6.Fill( ( upval ) / ( loval + upval ) );
          
          if( etadb )
            cout << "loval, upval, eta -- " << loval << ", "
            << upval << ", "
            << ( upval ) / ( loval + upval ) << endl;
          
          dphEtaRow.Fill( loval + upval ); // charge control plot
          
        } // rowetaList

	// for 25x10 on rot90:

	if( fifty && ymod > 0.010 && ymod < 0.040 ) { // central row
          roiq0vsxm.Fill( xmod5*1E3, q0 );
          roiq1vsxm.Fill( xmod5*1E3, q1 );
          roiq2vsxm.Fill( xmod5*1E3, q2 );
        }
        else if( !fifty && ymod > 0.010 && ymod < 0.090 ) { // central row
          roiq0vsxm.Fill( xmod5*1E3, q0 );
          roiq1vsxm.Fill( xmod5*1E3, q1 );
          roiq2vsxm.Fill( xmod5*1E3, q2 );
        }
	// 5 3 6
	// 1 0 2
	// 7 4 8

	double sumpc = ph0; // 2x2 mask
	double sumqc = q0;

	if( xmod5 < 0.025 )  { // 50x50
	  sumpc += ph1;
	  sumqc += q1;
	  roiq2vsq0.Fill( q0+q3+q4, q2+q6+q8 ); // Tsunami right
	  if( ymod5 > 0.025 ) {
	    sumpc += ph5;
	    sumqc += q5;
	  }
	  if( ymod5 < 0.025 ) {
	    sumpc += ph7;
	    sumqc += q7;
	  }
	}
	if( xmod5 > 0.025 )  {
	  sumpc += ph2;
	  sumqc += q2;
	  roiq1vsq0.Fill( q0+q3+q4, q1+q5+q7 ); // Tsunami left
	  if( ymod5 > 0.025 ) {
	    sumpc += ph6;
	    sumqc += q6;
	  }
	  if( ymod5 < 0.025 ) {
	    sumpc += ph8;
	    sumqc += q8;
	  }
	}
	if( ymod5 < 0.025 ) {
	  sumpc += ph4;
	  sumqc += q4;
	  roiq3vsq0.Fill( q0+q1+q2, q3+q5+q6 ); // Tsunami up
	}
	if( ymod5 > 0.025 ) {
	  sumpc += ph3;
	  sumqc += q3;
	  roiq4vsq0.Fill( q0+q1+q2, q4+q7+q8 ); // Tsunami down
	}

	roiphHisto.Fill( ( sump0+sumpA+sumpB) * norm );
	roip0Histo.Fill( sump0 * norm );
	roippHisto.Fill( (sump0+p2nd) * norm );
	roipcHisto.Fill( sumpc * norm );
	roipAHisto.Fill( sumpA * norm );
	roipBHisto.Fill( sumpB * norm );

	roiqHisto.Fill( (sumq0+sumqA+sumqB) * norm );
	roiq0Histo.Fill( sumq0 * norm );
	roiq012Histo.Fill( (sumq0+q1+q2) * norm );
	roiqqHisto.Fill( (sumq0+q2nd) * norm );
	roiqcHisto.Fill( sumqc * norm );
	roiqAHisto.Fill( sumqA * norm );
	roiqBHisto.Fill( sumqB * norm );

	roiq0vsxmym->Fill( xmod5*1E3, ymod5*1E3, sumq0 * norm ); // edge deficit
	roiqqvsxmym->Fill( xmod5*1E3, ymod5*1E3, (sumq0+q2nd) * norm ); // corner deficit
	roiqcvsxmym->Fill( xmod5*1E3, ymod5*1E3, sumqc * norm ); // corner deficit
	roiqAvsxmym->Fill( xmod5*1E3, ymod5*1E3, sumqA * norm );
	roiqBvsxmym->Fill( xmod5*1E3, ymod5*1E3, sumqB * norm );

	if( qset.size() > 2 ) {
	  auto qit = qset.begin();
	  roipxq1Histo.Fill( -(*qit++) ); // highest
	  roipxq2Histo.Fill( -(*qit++) ); // 2nd
	  roipxq3Histo.Fill( -(*qit) ); // 3rd
	}

      } // fid

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT pixel clusters:

      int nm[99] = {0};
      double minnm = -1;
      double dmin = 19.9; // [mm]
      double dxmin = 9;
      double dymin = 9;
      double pdxmin = 9;
      double pdymin = 9;
      double pdmin = 19;
      double clQ0 = 0;
      //double clQ0 = -1; // hide ineficiency
      double clpx = 0;
      double clph = 0;
      int clsz0 = 0;

      for( vector<cluster>::iterator c = cl0[iDUT].begin(); c != cl0[iDUT].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;

	// cluster isolation:

	bool isoc = 1;
        //for( vector<cluster>::iterator c2 = cl0[iDUT].begin(); c2 != cl0[iDUT].end(); ++c2 ) {
        //  if( c2 == c ) continue;
        //  if( fabs( c2->col - ccol ) < 8 &&
        //	fabs( c2->row - crow ) < 8 &&
        //	c2->charge > c->charge ) // mask close-by small clusters
        //    isoc = 0;
        //}
        
        /*
         ( *chip0 == 120 ) // noisy irrad
         isoc = 1;
         if( chip0 == 122 ) // noisy irrad
           isoc = 1;
         if( chip0 == 123 ) // noisy irrad
           isoc = 1;
         if( chip0 == 126 ) // noisy irrad
           isoc = 1;
         if( chip0 == 128 ) // noisy irrad
           isoc = 1;
         if( chip0 == 130 ) // noisy irrad
           isoc = 1;
         if( chip0 == 133 ) // noisy irrad
           isoc = 1;
         if( chip0 == 134 ) // noisy irrad
           isoc = 1;
         if( chip0 == 137 ) // noisy irrad
           isoc = 1;
         if( chip0 == 174 ) // noisy irrad
           isoc = 1;
         if( chip0 == 179 ) // noisy irrad
           isoc = 1;
         if( chip0 == 109 ) // noisy irrad
           isoc = 1;
         if( chip0 == 115 ) // noisy irrad
           isoc = 1;
         if( chip0 == 193 ) // noisy irrad
           isoc = 1;
         if( chip0 == 191 ) // noisy irrad
           isoc = 1;
         if( chip0 == 119 ) // noisy irrad
           isoc = 1;
         if( chip0 > 300 ) // noisy 3D
           isoc = 1;
         */

	double Q0 = c->charge * norm; // cluster charge normalized to vertical incidence

	double Qx = exp( -Q0 / qwid ); // Moyal weight

	int colmin = 999;
	int colmax = -1;
	int rowmin = 999;
	int rowmax = -1;

	double pcol[nx[iDUT]];
	for( int icol = 0; icol < nx[iDUT]; ++icol ) pcol[icol] = 0;

	double prow[ny[iDUT]];
	for( int irow = 0; irow < ny[iDUT]; ++irow ) prow[irow] = 0;

	double qcol[nx[iDUT]];
	for( int icol = 0; icol < nx[iDUT]; ++icol ) qcol[icol] = 0;

	double qrow[ny[iDUT]];
	for( int irow = 0; irow < ny[iDUT]; ++irow ) qrow[irow] = 0;

	double sumad = 0;
	double sumph = 0;
	double qseed = 0;
	double pseed = 0;

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

	  pcol[icol] += px->ph; // project cluster onto cols
	  prow[irow] += px->ph; // project cluster onto rows
	  qcol[icol] += q; // project cluster onto cols
	  qrow[irow] += q; // project cluster onto rows

	  sumad += px->adc;
	  sumph += px->ph;

	  if( q > qseed ) {
	    qseed = q;
	    pseed = px->ph;
	  }

	} // pix

	int ncol = colmax - colmin + 1;
	int nrow = rowmax - rowmin + 1;

        double PH0 = sumph*norm;
        //double Px = exp( -sumph / pwid );
        double Px = exp( - PH0  / pwid ); // normalized

	// eta-algo in rows:

	double p1 = 0; // highest charge
	double p2 = 0; // 2nd highest
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
	    p2 = p1;
	    p1 = prow[irow];
	    i2 = i1;
	    i1 = irow;
	  }
	  else if( q > q2 ) {
	    q2 = q;
	    i2 = irow;
	    p2 = prow[irow];
	  }

	} // rows

	double q12 = q1 + q2;
	double eta = 0;
	if( q12 > 1 ) eta = ( q2 - q1 ) / q12; // always negative
	if( i1 > i2 ) eta = -eta; // right-left

	double p12 = p1 + p2;
	double etap = 0;
	if( p12 > 1 ) etap = ( p2 - p1 ) / p12; // always negative
	if( i1 > i2 ) etap = -etap; // right-left

	double ceta = i1; // eta from diffusion: limited range
	//if( q12 > 1 ) ceta = 0.5*(i1+i2) + 0.30*(i2-i1)*(q2-q1)/q12; // 322291 cmsduc2 5.58
	//if( q12 > 1 ) ceta = 0.5*(i1+i2) + 0.25*(i2-i1)*(q2-q1)/q12; // 322291 cmsduc2 5.39
	if( q12 > 1 ) ceta = 0.5*(i1+i2) + 0.20*(i2-i1)*(q2-q1)/q12; // 322291 cmsduc2 5.31
	//if( q12 > 1 ) ceta = 0.5*(i1+i2) + 0.15*(i2-i1)*(q2-q1)/q12; // 322291 cmsduc2 5.35
	//if( q12 > 1 ) ceta = 0.5*(i1+i2) + 0.10*(i2-i1)*(q2-q1)/q12; // 322291 cmsduc2 5.53
	//if( q12 > 1 ) ceta = 0.5*(i1+i2) + 0.05*(i2-i1)*(q2-q1)/q12; // 322291 cmsduc2 5.81

	// column cluster:

	sumq = 0;
	double sumcol = 0;
	double sumcol2 = 0;
	double sumcol3 = 0;
	p1 = 0; // highest charge
	p2 = 0; // 2nd highest
	q1 = 0; // highest charge
	q2 = 0; // 2nd highest
	int j1 = 0;
	int j2 = 0;

	for( int jcol = colmin; jcol <= colmax; ++jcol ) {

	  double q = qcol[jcol];
	  if( q > p1 ) {
	    q2 = q1;
	    q1 = q;
	    p2 = p1;
	    p1 = pcol[jcol];
	    j2 = j1;
	    j1 = jcol;
	  }
	  else if( q > p2 ) {
	    q2 = q;
	    j2 = jcol;
	    p2 = pcol[jcol];
	  }

	  sumq += q;
	  sumcol += jcol*q;
	  double dcol = jcol - ccol; // distance from COG
	  sumcol2 += dcol*dcol*q; // 2nd central moment
	  sumcol3 += dcol*dcol*dcol*q; // 3rd central moment

	} // cols

	q12 = q1 + q2;
	double uta = 0;
	if( q12 > 1 ) uta = ( q2 - q1 ) / q12;
	if( j1 > j2 ) uta = -uta;

	p12 = p1 + p2;
	double utap = 0;
	if( p12 > 1 ) utap = ( p1 - p2 ) / p12; // always negative
	if( j1 > j2 ) utap = -utap; // right-left

	double cuta = j1;
	//if( q12 > 1 ) cuta = 0.5*(j1+j2) + 0.10*(j2-j1)*(q2-q1)/q12; // cmsdvc2 13.6
	if( q12 > 1 ) cuta = 0.5*(j1+j2) + 0.20*(j2-j1)*(q2-q1)/q12; // cmsdvc2 13.7

	bool lq = 1; // Landau peak
	if( ( Q0 < 8 || Q0 > 16 ) ) lq = 0; // r102 50x50
	//if( ( Q0 < 8 || Q0 > 14 ) ) lq = 0; // r144

	if( run >= 31632 && chip0 >= 119 && chip0 <= 138 ) { // Feb-Mar 2018 irrad gain_2
	  lq = 1;
	  if( ( Q0 < 3 || Q0 > 12 ) ) lq = 0;
	}

	if( run >= 33500 && // Sep 2018
	    ( chip0 == 174 || chip0 == 179 || chip0 == 193 || chip0 == 109 || chip0 == 155 || chip0 == 115 ) ) {
	  lq = 1;
	  if( ( Q0 < 3 || Q0 > 10 ) ) lq = 0;
	}

	if( chip0 >= 300 ) { // 3D
	  lq = 1;
	  if( ( Q0 < 11 || Q0 > 22 ) ) lq = 0; // 3D
	}

	// DUT - triplet:

	double cmsx = ( ccol + 1.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	double cmsy = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	double cmsu = ( cuta + 1.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // x from uta
	double cmsv = ( ceta + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // y from eta

	if( rot90 ) {
	  cmsx = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	  cmsy = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	  cmsu = ( ceta + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // x from eta
	  cmsv = ( cuta + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // y from uta
	}
	
	cmsxHisto.Fill( cmsx ); // check positioning
        cmsyHisto.Fill( cmsy );

	// residuals for pre-alignment:

	if( liso && isoc ){//&& lddt ) {
          
	  cmsxvsx->Fill( x4, cmsx );
	  cmsyvsy->Fill( y4, cmsy );

	  cmssxaHisto.Fill( cmsx + x3 ); // rot, tilt and turn but no shift
	  cmsdxaHisto.Fill( cmsx - x3 ); // peak

	  cmssyaHisto.Fill( cmsy + y3 );
	  cmsdyaHisto.Fill( cmsy - y3 );

	}

	// residuals:

	double cmsdx = cmsx - x4; // triplet extrapol
	double cmsdy = cmsy - y4;
	double cmsdu = cmsu - x4; // triplet extrapol
	double cmsdv = cmsv - y4;

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
	      
	  if( pdxy < 0.2 ){ // same as iw = 20
            c->fake = false;
            //cout << "event " << iev << endl
            //     << "  track x, y " <<   x4 << ", " <<   y4 << endl
            //     << "  clust x, y " << cmsx << ", " << cmsy << endl
            //     << "  =>> NO FAKE" << endl
            //     << "  col, row " << c->vpix[ipx].col << " " << c->vpix[ipx].row << endl;
          }
          //cout << "  " << c->fake << endl;
          
          //if( c->vpix[ipx].significance < 3.5 || !c->fake ){
          //  cout << "event " << iev << endl
          //       << "  col, row " << c->vpix[ipx].col << " " << c->vpix[ipx].row << endl
          //       << "  signi " << c->vpix[ipx].significance << endl
          //       << "  pdxy  " << pdxy << endl
          //       << "  fake? " << c->fake << endl;
          //}
              
	} // pix
	
	//for( int iw = 1; iw < 99; ++iw ){ // minimum efficient radius
	//  if( nm[iw] == 0 ){
	//    minnm = iw*0.010;
	//  }
	//} // this is inside cluster loop, just want to know once per triplet
	// anyhow i don't realy need this any more, i know my fake events now
	
	//hitposi->Fill( xcol, yrow, 1 ); // somehow shifted
	
	
	cmsdxHisto.Fill( cmsdx );
	cmsdyHisto.Fill( cmsdy );

        if( fabs(cmsdy) < ycut ) { // cut y draw x
          
          cmsdxcHisto.Fill( cmsdx ); // align: same sign
          
          if( liso && isoc && lddt ) { // triplet isolation, cluster isolation, timing
            
            cmsdxcinfHisto.Fill( cmsdx );
            
            //if( fabs(x4) < 3.8 ) { // fiducial x
            //
            //  if( lq ) { // landau
            //    
            //  } // landau
            //
            //} // fiducial x
            //
            //if( fabs(y4) < 3.9 ) { // fiducial y
            //  
            //  if( lq ) { // landau
            //    
            //  } // landau
            //
            //} // fiducial y
            
            if( fabs(x4) < 3.8 && fabs(y4) < 3.9 ) { // fiducial xy
              
              cmsdxciHisto.Fill( cmsdx );         // align: same sign
              cmsdxvsev.Fill( iev, cmsdx );       // align stability
              cmsdxvsev1->Fill( iev, cmsdx );     // sync stability
              cmsdxvsev2->Fill( iev, cmsdx );     // sync stability
              cmsmadxvsq.Fill( Q0, fabs(cmsdx) ); // resolution vs charge
              
              if( lq ) { // landau
                
                cmsdxcqHisto.Fill( cmsdx );
                cmsdxvsx.Fill( x4, cmsdx ); // align: same sign
                cmsdxvsxmq->Fill( xmod5*1E3, cmsdx*1E3 );
                cmsdxvsy.Fill( y4, cmsdx ); // align: opposite sign, depends on turn
                cmsdxvsyc.Fill( y4, cmsdx ); // tighter dx
                cmsdxvstx.Fill( sxA*1E3, cmsdx );
                cmsmadxvsx.Fill( x4, fabs(cmsdx) );
                cmsmadxvsy.Fill( y4, fabs(cmsdx) );
                cmsmadxvsxm.Fill( xmod*1E3, fabs(cmsdx)*1E3 ); // within pixel
                
              } // landau
              
              if( rot90 && nrow == 2 ) {  // cluster size
                cmsdxvsxm2->Fill( xmod5*1E3, cmsdx*1E3 );
                cmsduvsxm2->Fill( xmod5*1E3, cmsdu*1E3 );
                cmsdxc2Histo.Fill( cmsdx );
                cmsduc2Histo.Fill( cmsdu );
              } // cluster size
              
              else if( !rot90 && ncol == 2 ) { // cluster size
                
                cmsdxvsxm2->Fill( xmod5*1E3, cmsdx*1E3 );
                cmsduvsxm2->Fill( xmod5*1E3, cmsdu*1E3 );
                cmsdxc2Histo.Fill( cmsdx );
                cmsduc2Histo.Fill( cmsdu );
                
              } // cluster size
              
            } // fiducial xy
            
          } // triplet isolation, cluster isolation, timing
          
        } // cut y draw x

        if( fabs(cmsdx) < xcut ) { // cut x draw y
          
          cmsdycHisto.Fill( cmsdy );
          
          if( liso && isoc && lddt ) { // triplet isolation, cluster isolation, timing
            
            cmsdycinfHisto.Fill( cmsdy );
            
            if( fabs(x4) < 3.8 ) { // fiducial x
              
              if( lq ) { // landau
                
                cmsdyvsy.Fill( y4, cmsdy );
                cmsmadyvsy.Fill( y4, fabs(cmsdy) );
                
              } // landau
              
            } // fiducial x
            if( fabs(y4) < 3.9 ) { // fiducial y
              
              if( lq ) { // landau
                
                cmsdyvsx.Fill( x4, cmsdy );
                cmsdyvsxc.Fill( x4, cmsdy ); // tighter dy range
                cmsmadyvsx.Fill( x4, fabs( cmsdy ) );
                cmsmady8vsx.Fill( x4, fabs( cmsdy8 ) );
                
              } // landau
              
            } // fiducial y
            if( fabs(x4) < 3.8 && fabs(y4) < 3.9 ) { // fiducial xy
              
              cmsdyciHisto.Fill( cmsdy );
              cmsmadyvsq.Fill( Q0, fabs(cmsdy) ); // resolution vs charge
              
              if( lsixlk ) { // modlink
                cmsdyci6Histo.Fill( cmsdy );
              } // modlink
              
              if( x4 < 1.4 ) { // Cu cutout rot90
                cmsdy8cHisto.Fill( cmsdy8 );
                if( lsixlk && liso &&  isoc )
                  cmsdy8ciHisto.Fill( cmsdy8 );
              }
              
              if( ncol < 3 && nrow < 3 ) { // ncol
                
                cmsdyc3Histo.Fill( cmsdy ); // 31166: side peaks at +-0.125 mm
                
                if( lsixlk ) { // modlink
                  cmsdyc3iHisto.Fill( cmsdy ); // 31166: side peaks eliminated
                } // modlink
                
                if( x4 < 1.4 ) { // Cu cutout rot90
                  cmsdy8c3Histo.Fill( cmsdy8 );
                  if( lsixlk && liso &&  isoc )
                    cmsdy8c3iHisto.Fill( cmsdy8 );
                } // Cu
                
              } // ncol
              
              if( lq ) { // landau
                
                cmsdycqHisto.Fill( cmsdy ); // 3D 31215: 5.9
                cmsdyvsym.Fill( ymod*1E3, cmsdy );
                cmsdyvsty.Fill( syA*1E3, cmsdy );
                cmsdyvsev.Fill( iev, cmsdy ); // trend?
                cmsmadyvsev.Fill( iev, fabs(cmsdy) );
                cmsmadyvsty.Fill( syA*1E3, fabs(cmsdy) );
                cmsmadyvsxm.Fill( xmod*1E3, fabs(cmsdy)*1E3 ); // within pixel
                cmsmadyvsym.Fill( ymod*1E3, fabs(cmsdy)*1E3 ); // within pixel
                cmsmadyvsxmym->Fill( xmod*1E3, ymod*1E3, fabs(cmsdy)*1E3 ); // within pixel
                
                if( lsixlk ) { // modlink
                  cmsdycqiHisto.Fill( cmsdy );
                }
                
                if( x4 < 1.4 ) { // Cu cutout rot90
                  cmsdy8cqHisto.Fill( cmsdy8 );
                  if( lsixlk && liso &&  isoc )
                    cmsdy8cqiHisto.Fill( cmsdy8 );
                }
                
                if( rot90 && ncol == 2 ) { // cluster size
                  cmsdyvsym2->Fill( ymod5*1E3, cmsdy*1E3 );
                  cmsdvvsym2->Fill( ymod5*1E3, cmsdv*1E3 );
                  cmsdyc2Histo.Fill( cmsdy );
                  cmsdvc2Histo.Fill( cmsdv );
                } // cluster size
                
                if( !rot90 && nrow == 2 ) { // cluster size
                  
                  cmsdyvsym2->Fill( ymod5*1E3, cmsdy*1E3 );
                  cmsdvvsym2->Fill( ymod5*1E3, cmsdv*1E3 );
                  cmsdyc2Histo.Fill( cmsdy );
                  cmsdvc2Histo.Fill( cmsdv );
                  
                } // cluster size
                
              } // landau
              
              if( Q0 > 9 && Q0 < 14 ) { // landau2
                cmsdycq2Histo.Fill( cmsdy );
              }
              
            } // fiducial xy
            
          } // triplet isolation, cluster isolation, timing
          
        } // cut x draw y
        
        // xy cuts:
        
        if( fabs(cmsdx) < xcut &&
            fabs(cmsdy) < ycut && lddt ) {
          
          for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {
            cmsadcHisto.Fill( px->ph );
            cmscolHisto.Fill( px->col );
            cmsrowHisto.Fill( px->row );
          }
          
          } // dxy

	if( fabs(x4) < 3.8 && fabs(y4) < 3.9 && // fiducial
	    fabs(cmsdx) < xcut &&
	    fabs(cmsdy) < ycut ) {

	  cmsq0aHisto.Fill( Q0 ); // Landau
	  cmsqxvst5.Fill( evsec, Qx );

	  if( Q0 < 7.5 ) {
	    cmsxmymq0->Fill( xmod*1E3, ymod*1E3 ); // cluster size map
	    cmsxyq0->Fill( x4, y4 ); // cluster size map
	  }
	  else
	    cmsxmymq1->Fill( xmod*1E3, ymod*1E3 ); // cluster size map

	} // fiducial linked

	if( liso && isoc &&  // cut43D
	    fabs(cmsdx) < xcut &&
	    fabs(cmsdy) < ycut ) {

	  trixAlkHisto.Fill( xA ); // telescope coordinates at DUT
	  triyAlkHisto.Fill( yA );
	  trixyAlkHisto->Fill( xA, yA );
	  trixclkHisto.Fill( xc ); // telescope coordinates
	  triyclkHisto.Fill( yc );
	  trixyclkHisto->Fill( xc, yc ); // telescope coordinates

	  if( fabs(x4) < 3.8 && fabs(y4) < 3.9 ) { // isolated fiducial

	    cmsnpxHisto.Fill( c->size );
	    cmsncolHisto.Fill( ncol );
	    cmsnrowHisto.Fill( nrow );

	    if( nrow == 2 ) {

	      cmsetaHisto.Fill( eta );
	      cmsetapHisto.Fill( etap );
	      if( lq ) {
		cmsetaqHisto.Fill( eta );
		cmsetapqHisto.Fill( etap );
		cmsetavsxm2->Fill( xmod5*1E3, eta ); // within pixel, rot90
		cmsetavsym.Fill( ymod*1E3, eta ); // within pixel
	      }

	    } // 2-row

	    if( ncol == 2 ) {

	      cmsutaHisto.Fill( uta );
	      cmsutapHisto.Fill( utap );
	      if( lq ) {
		cmsutaqHisto.Fill( uta );
		cmsutapqHisto.Fill( utap );
		cmsutavsxm.Fill( xmod*1E3, uta ); // within pixel
		cmsutavsym.Fill( ymod*1E3, uta ); // within pixel
	      }

	    } // 2-col

	    cmsncolvsxm.Fill( xmod*1E3, ncol );
	    cmsnrowvsxm.Fill( xmod*1E3, nrow );

	    cmsncolvsym.Fill( ymod*1E3, ncol ); // within pixel
	    cmsnrowvsym.Fill( ymod*1E3, nrow ); // within pixel

	    cmsnpxvsxmym->Fill( xmod*1E3, ymod*1E3, c->size ); // cluster size map
	    cmsnpxvsxmym2->Fill( xmod2*1E3, ymod*1E3, c->size ); // cluster size map

	    if( lq ) {

	      cmsncolvsxmq.Fill( xmod*1E3, ncol );
	      cmsnrowvsxmq.Fill( xmod*1E3, nrow );

	      cmsncolvsymq.Fill( ymod*1E3, ncol ); // within pixel
	      cmsnrowvsymq.Fill( ymod*1E3, nrow ); // within pixel

	      cmsnpxvsxmymq->Fill( xmod*1E3, ymod*1E3, c->size ); // cluster size map
	      cmsnpxvsxmym2q->Fill( xmod2*1E3, ymod*1E3, c->size ); // cluster size map

	    } // Landau peak

	    cmsad0Histo.Fill( sumad*norm ); // adc with common mode
	    cmsph0Histo.Fill( sumph*norm ); // dph common mode subtracted
	    cmsv0Histo.Fill( Q0/ke ); // Landau [mV]
	    cmsq0Histo.Fill( Q0 ); // Landau
	    cmsq00Histo.Fill( Q0 ); // Landau
            // 3D for differential
            cmsclph3D->Fill( x4, y4, PH0 );
            
            if( lddt ) cmsq0dtHisto.Fill( Q0 ); // low Q
	    if( lsixlk ) cmsq0sixHisto.Fill( Q0 ); // clean

	    if( drbeam < rbeam )
	      cmsq0bHisto.Fill( Q0 ); // beam spot
	    else
	      cmsq0hHisto.Fill( Q0 ); // beam halo

	    if( ncol < 3 && nrow < 3 )
	      cmsq03Histo.Fill( Q0 ); // Landau
 
	    if( sqrt( pow( xmod-0.050, 2 ) + pow( ymod-0.050, 2 ) ) < 0.010 ) // 50x50 bias dot
	      cmsq0dHisto.Fill( Q0 );

	    if( sqrt( pow( xmod-0.050, 2 ) + pow( ymod-0.050, 2 ) ) > 0.020 ) // no dot
	      cmsq0nHisto.Fill( Q0 );

	    if( Q0 < 7 ) { // SCM144 in 31425, 31473: low Q peak: out-of-time crap
	      cmslqxyHisto->Fill( x4, y4 ); // everywhere
	      cmslqxmymHisto->Fill( xmod*1E3, ymod*1E3 ); // everywhere
	      cmslqdxHisto.Fill( cmsdx ); // symmetric, centered
	      cmslqdyHisto.Fill( cmsdy ); // symmetric, centered
	      cmslqt4Histo.Fill( evsec ); // anytime
	      cmslqt5Histo.Fill( evsec ); // 
	      cmslqnpxHisto.Fill( pb.size() );
	      cmslqnclHisto.Fill( cl[iDUT].size() ); // 1.55
	      cmslqntriHisto.Fill( triplets.size() ); // 2.8
	      cmslqndriHisto.Fill( driplets.size() );
	    }

	    // cluster charge profiles, exponential weighting Qx:

	    cmsqxvsx.Fill( x4, Qx );
	    cmsqxvsy.Fill( y4, Qx );
	    cmsqxvsr.Fill( drbeam, Qx );
	    cmsqxvsxy->Fill( x4, y4, Qx );

	    cmspxvsx.Fill( x4, Px );
	    cmspxvsy.Fill( y4, Px );
	    cmspxvsr.Fill( drbeam, Px );
	    cmspxvsxy->Fill( x4, y4, Px );

	    cmsqxvsxm.Fill( xmod*1E3, Qx ); // Q within pixel, depends on ph cut
	    cmsqxvsym.Fill( ymod*1E3, Qx ); // Q within pixel
	    cmsqxvsxmym->Fill( xmod*1E3, ymod*1E3, Qx ); // cluster charge profile
	    cmsqxvsxmym2->Fill( xmod2*1E3, ymod*1E3, Qx ); // cluster charge profile

            for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {
	      cmspxphHisto.Fill( px->ph );
              cmspxqHisto.Fill( px->q );
	    }

	    cmssdpHisto.Fill( pseed );
	    cmssdqHisto.Fill( qseed );

	  } // fiducial

	} // cut xy

	if( dxy < dmin ) {
	  dmin = dxy;
	  clQ0 = Q0;
	  clsz0 = c->size;
	  dxmin = cmsdx;
          dymin = cmsdy;
          clpx = Px;
          clph = PH0;
	}

      } // loop DUT clusters

      if( ldb ) cout << "    eff " << nm[20] << endl << flush;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT efficiency vs isolated MOD-linked fiducial tracks:

      if( lsixlk
	  //&& cl0[iDUT].size() < 2 // empty or single cluster, same eff
	  ) {

	double fidx0 =-3.7;
	double fidx9 = 3.7;
	double fidy0 =-3.8;
	double fidy9 = 3.8;

	if( rot90 ) {
	  fidx0 =-3.8;
	  fidx9 = 3.8;
	  fidy0 =-3.7;
	  fidy9 = 3.7;
	}
	
	// cut off corner for 137
	if( chip0 == 137 ){
          fidy9 = 2.5 - 0.33 *x4;
          if( fidy9 > 3.1 ) 
            fidy9 = 3.1;
        }
        
        // beam spot cut --------------> 
        // circle parameters
        double xcnt = 0;
        double ycnt = 0;
        double rads = 1000;
        
        // for each sensor
        //if( run >= 36252 && run <= 36260 ) { xcnt = -2.56; ycnt = 1.71; rads = 1.5; }
        //if( run >= 36263 && run <= 36273 ) { xcnt = -2.38; ycnt = 2.03; rads = 1.5; }
        //if( run >= 36276 && run <= 36301 ) { xcnt = -2.55; ycnt = 0.74; rads = 1.5; }
        //if( run >= 36309 && run <= 36321 ) { xcnt = -1.94; ycnt = 1.88; rads = 1.5; }
        
        // radius cut
        double fidx0r = - ( sqrt( rads * rads - ( y4 - ycnt ) * ( y4 - ycnt ) ) - xcnt );
        double fidx9r =   ( sqrt( rads * rads - ( y4 - ycnt ) * ( y4 - ycnt ) ) + xcnt );
        double fidy0r = - ( sqrt( rads * rads - ( x4 - xcnt ) * ( x4 - xcnt ) ) - ycnt );
        double fidy9r =   ( sqrt( rads * rads - ( x4 - xcnt ) * ( x4 - xcnt ) ) + ycnt );
        
        // if sqrt argument negative -> outside of circle. yields nan, should be 0 instead
        if( 0 > ( rads * rads - ( y4 - ycnt ) * ( y4 - ycnt ) ) ) fidx0r = 0;
        if( 0 > ( rads * rads - ( y4 - ycnt ) * ( y4 - ycnt ) ) ) fidx9r = 0;
        if( 0 > ( rads * rads - ( x4 - xcnt ) * ( x4 - xcnt ) ) ) fidy0r = 0;
        if( 0 > ( rads * rads - ( x4 - xcnt ) * ( x4 - xcnt ) ) ) fidy9r = 0;
        
        // always use harder cut
        if( fidx0r > fidx0 ) fidx0 = fidx0r;
        if( fidx9r < fidx9 ) fidx9 = fidx9r;
        if( fidy0r > fidy0 ) fidy0 = fidy0r;
        if( fidy9r < fidy9 ) fidy9 = fidy9r;
        // <-------------- beam spot cut
        
        if( run == 32277 ){ // cut off noisy pixel for simulation studies -> adapt fiducial
          fidx0 =-3.6;
          fidy0 =-3.0;
        }

	if( x4 > fidx0 && x4 < fidx9 &&
	    y4 > fidy0 && y4 < fidy9 &&
	    lddt ) { // fiducial

	  effvsdmin.Fill( dddmin, nm[20] ); // at MOD, small effect

	  if( dddmin > 0.4 )
	    effvstmin.Fill( ttdmin, nm[20] ); // at DUT, flat

	} // fid

	//if( dddmin > 0.4 && ttdmin > 1.4 ) { // stronger iso [mm] 99.94
	//if( dddmin > 0.4 && ttdmin > 0.6 ) { // iso [mm] 99.93
	if( iev > 400 &&
	    ( run != 31147 || iev < 6400 || iev > 8900 ) &&
	    filled != A &&
	    dddmin > 0.6 ) { // iso [mm] 99.94

	  if( lddt ) {

	    sixxylkHisto->Fill( xA, yA );
            if( nm[20] ) sixxyeffHisto->Fill( xA, yA );
            
            effvsxy->Fill( x4, y4, nm[20] ); // map
            effvsxyf->Fill( xcol, yrow, 1-nm[20] ); // map; now really in pixel size
            
            // plots for studies differential in fluence
            
            cmsclpxvsxy->Fill( x4, y4, clpx ); 
            cmsmadxvsxy->Fill( x4, y4, fabs( dxmin )*1E3 );
            cmsmadyvsxy->Fill( x4, y4, fabs( dymin )*1E3 );
            cmsclpxh->Fill( clpx );
            cmsmadxh->Fill( dxmin*1E3 );
            cmsmadyh->Fill( dymin*1E3 );
            // 3D for differential
            //cmsclph3D->Fill( x4, y4, clph );
            
	    if( y4 > fidy0 && y4 < fidy9 )
	      effvsx.Fill( x4, nm[20] );

	    if( x4 > fidx0 && x4 < fidx9 )
	      effvsy.Fill( y4, nm[20] );

	  } // ddt

	  bool fiducial = 0;

	  if( x4 > fidx0 && x4 < fidx9 &&
	      y4 > fidy0 && y4 < fidy9 ) { // fiducial
	    fiducial = 1;
	  }

	  if( chip0 == 142 && y4 < -1.5 && y4 > -1.6 )
	    fiducial = 0;

          if( fiducial ) {
            
            effvsxyfid->Fill( x4, y4, nm[20] ); // map
            
            effvsev1.Fill( iev, nm[20] );
            effvsev2.Fill( iev, nm[20] );
            effvst1.Fill( evsec, nm[20] );
            effvst2.Fill( evsec, nm[20] );
            effvst3.Fill( evsec, nm[20] );
            effvst4.Fill( evsec, nm[20] );
            effvst6.Fill( evsec/3600, nm[20] );
            
            if( lddt ) {
              
              effvst1t.Fill( evsec, nm[20] );
              effvst2t.Fill( evsec, nm[20] );
              effvst3t.Fill( evsec, nm[20] );
              effvst4t.Fill( evsec, nm[20] );
              effvst6t.Fill( evsec/3600, nm[20] );
              
              effvsdr.Fill( drbeam, nm[20] );

	      cmsdminHisto.Fill( dmin );
	      cmsdxminHisto.Fill( dxmin );
	      cmsdyminHisto.Fill( dymin );
	      cmspdxminHisto.Fill( pdxmin );
	      cmspdyminHisto.Fill( pdymin );
	      cmspdminHisto.Fill( pdmin );

	      for( int iw = 1; iw < 99; ++iw )
		effvsdxy.Fill( iw*0.010-0.001, nm[iw] );

	      effdminHisto.Fill( dmin );
	      if( nm[20] == 0 ) {
		effdmin0Histo.Fill( dmin );
		effrxmin0Histo.Fill( dxmin/dmin );
		effrymin0Histo.Fill( dymin/dmin );
		effdxmin0Histo.Fill( dxmin );
		effdymin0Histo.Fill( dymin );
	      }

	      effq0Histo.Fill( clQ0 );
	      if( ndrilk == 1 ) effqrHisto.Fill( clQ0 );
	      if( dmin > 0.1 ) effq1Histo.Fill( clQ0 );
	      if( dmin > 0.2 ) effq2Histo.Fill( clQ0 ); // Landau tail
	      if( dmin > 0.3 ) effq3Histo.Fill( clQ0 );
	      if( dmin > 0.4 ) effq4Histo.Fill( clQ0 );
	      if( dmin > 0.5 ) effq5Histo.Fill( clQ0 );
	      if( dmin > 0.6 ) effq6Histo.Fill( clQ0 );
	      if( dmin > 0.7 ) effq7Histo.Fill( clQ0 );
	      if( dmin > 0.8 ) effq8Histo.Fill( clQ0 );
	      if( dmin > 0.9 ) effq9Histo.Fill( clQ0 );

	      effnpxvsxmym->Fill( xmod*1E3, ymod*1E3, clsz0 ); // cluster size map
	      effq0vsxmym->Fill( xmod*1E3, ymod*1E3, clQ0 ); // cluster charge profile
	      if( sqrt( pow( xmod-0.050, 2 ) + pow( ymod-0.050, 2 ) ) < 0.010 ) // 50x50 bias dot
		effq0dHisto.Fill( clQ0 );
	      if( sqrt( pow( xmod-0.050, 2 ) + pow( ymod-0.050, 2 ) ) > 0.020 ) // no dot
		effq0nHisto.Fill( clQ0 );

              effvsxt->Fill( evsec, x4, nm[20] );
              effvsntri.Fill( triplets.size(), nm[20] ); // flat
              effvsndri.Fill( driplets.size(), nm[20] ); // flat
              
              effvsxmym->Fill( xmod*1E3, ymod*1E3, nm[20] );
              if( rot90 )
                effvsxmym2->Fill( xmod*1E3, ymod2*1E3, nm[20] );
              else
                effvsxmym2->Fill( xmod2*1E3, ymod*1E3, nm[20] );
              effvsxmym5->Fill( xmod5*1E3, ymod5*1E3, nm[20] );
              effvsxm.Fill( xmod, nm[20] ); // bias dot
              effvsym.Fill( ymod, nm[20] ); // bias dot
              
              effvstx.Fill( sxA, nm[20] );
              effvsty.Fill( syA, nm[20] );
              effvstxy.Fill( txy, nm[20] ); // no effect
              effvsdslp.Fill( sixdslp, nm[20] ); // no effect

	    } // ddt

	  } // fiducial

	} // iso

      } // six

    } // loop triplets iA

    gettimeofday( &tv, NULL );
    long s4 = tv.tv_sec; // seconds since 1.1.1970
    long u4 = tv.tv_usec; // microseconds
    zeit3 += s4 - s3 + ( u4 - u3 ) * 1e-6; // tracking

    // fill fake rate plots
    double nfks   = 0;
    double ncl    = 0;
    double isFake = 0;
    
    for( vector<cluster>::iterator c = cl0[iDUT].begin(); c != cl0[iDUT].end(); ++c ){
      
      totalncl++;
      
      double px = ( c->col + 1.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
      double py = ( c->row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
      if( rot90 ) {
        px = ( c->row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
        py = ( c->col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
      }
      
      isFake = 0;
      if( c->fake ){
        isFake = 1;
        nfks++;
        totalfks++;
        fksnpx.Fill( c->vpix.size() );
      }
      
      fksvsxy->Fill( px, py, isFake ); // map
      allnpx.Fill( c->vpix.size() );
      
    }
    
    ncl = (double)cl[iDUT].size();
    if( ncl > 0 ){
      allncl.Fill( ncl );
      fksncl.Fill( nfks );
      truncl.Fill( ncl - nfks );
      fksperevnt.Fill( nfks / ncl );
    }
    
    //cout << "event " << iev << endl
    //     << "nfks " << nfks << endl
    //     << "clsz " << cl[iDUT].size() << endl
    //     << "fake rate " << (double)nfks / ncl << endl;
    
    if( ldb ) cout << "done ev " << iev << endl << flush;

    if( iev != lev ){ 
      //evtdspl->Reset();
      hitposi->Reset();
    }

    if( syncmod ) { // shift all but MOD

      for( int ipl = 0; ipl < 6; ++ipl ) {
	cl0[ipl] = cl1[ipl]; // remember
	cl1[ipl] = cl[ipl]; // remember
      }

      cl0[iDUT] = cl1[iDUT]; // remember
      cl1[iDUT] = cl[iDUT]; // remember

    }

    ++iev;
  } while( reader->NextEvent() && iev <= lev ); // <= for same as scopes

  delete reader;

  gettimeofday( &tv, NULL );
  long s9 = tv.tv_sec; // seconds since 1.1.1970
  long u9 = tv.tv_usec; // microseconds

  cout << endl
       << "done after " << iev << " events"
       << " in " << s9 - s0 + ( u9 - u0 ) * 1e-6 << " s"
       << " (read and cluster tele " << zeit1 << " s"
       << ", DUT " << zeit2 << " s"
       << ", tracking " << zeit3 << " s)"
       << endl
       << "resyncs " << nresync
       << endl;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // hot pixels:

  cout << endl << "DUT hot pixel list for run " << run << endl;

  int hotcut = iev / 1024;

  ofstream oDUThotFile( DUThotFileName.str() );

  oDUThotFile << "# DUT hot pixel list for run " << run
              << " with " << iev << " events"
              << ", hot cut " << hotcut
              << endl;

  multiset <cpixel> pxset;

  cout << hmap[iDUT]->GetTitle() << endl;
  cout << hmap[iDUT]->GetEntries() << endl;
  cout << hmap[iDUT]->GetNbinsX() << endl;
  cout << hmap[iDUT]->GetNbinsY() << endl << flush;

  for( int ii = 1; ii <= hmap[iDUT]->GetNbinsX(); ++ii )

    for( int jj = 1; jj <= hmap[iDUT]->GetNbinsY(); ++jj ) {

      //cout << ii << " " << jj << endl;
        
      int n = hmap[iDUT]->GetBinContent(ii,jj);

      if( n > hotcut ) {
        cpixel px{ ii-1, jj-1, n };
        pxset.insert(px); // sorted by descending cnt
      }
  }

  cout << "hot " << pxset.size() << endl;

  for( auto px = pxset.begin(); px != pxset.end(); ++px ) {
    //cout << setw(3) << px->col << ", "
    //     << setw(3) << px->row << ":  "
    //     << px->cnt << endl;
    oDUThotFile << "pix "
              << setw(3) << px->col
              << setw(5) << px->row
              << "  " << px->cnt
              << endl;
  }
  cout << "DUT hot pixel list written to " << DUThotFileName.str() << endl;

  oDUThotFile.close();
    
  //for( int i = 0; i < 24; i++){ tgraphs not used now
  //  gphCor[i]->Write();
  //  gdphCor[i]->Write();
  //  gq0Cor[i]->Write();
  //}
  
  // finish noise efficiency analysis
  for(int col = 0; col < 155; col++){
    for(int row = 0; row < 160; row++){
      pavgdph->Fill( col+0.5, row+0.5, avgdph[col][row] / cntdph[col][row] );
      pcntdph->Fill( col+0.5, row+0.5, cntdph[col][row] );
      
      if( cntdph[col][row] > 5 ){ // want at least 5
        prmsdph->Fill( col+0.5, row+0.5, sqrt( rmsdph[col][row] / cntdph[col][row] ) );
        hrmsdph.Fill( sqrt( rmsdph[col][row] / cntdph[col][row] ) );
      }
    }
  }
  for(int col = 5; col < 155; col++){
    for(int row = 5; row < 155; row++){
      if( effvsxyf->GetBinContent( col, row ) < 0.01 ) { // effvsxyf is inverted efficiency
        hrmsdphEff1.Fill( prmsdph->GetBinContent( col, row ) );
      }
      else {
        hrmsdphEff0.Fill( prmsdph->GetBinContent( col, row ) );
      }
      rmsvseff.Fill( prmsdph->GetBinContent( col, row ) , 1 - effvsxyf->GetBinContent( col, row ) ); //re-invert
    }
  }  
  
  // convert px to ph
  for( int ii = 1; ii <= cmsclpxvsxy->GetNbinsX(); ++ii )
    for( int jj = 1; jj <= cmsclpxvsxy->GetNbinsY(); ++jj ) {
      double qx = cmsclpxvsxy->GetBinContent(ii,jj);
      double ns = cmsnoivsxy->GetBinContent(ii,jj);
      if( qx > 0 ) {
        double q = -pwid*log(qx);
        cmsclphvsxy->SetBinContent( ii, jj, q );
        cmsclphh->Fill( q );
        if( ns > 0 ){
          cmsstonvsxy->SetBinContent( ii, jj, q/ns );
          cmsstonh->Fill( q/ns );
        }
      }
      
    } // jj
    
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  histoFile->Write();
  histoFile->Close();
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOD alignment:

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
	//delete fdxvsx;
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
      //delete fdyvsx;

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
	//delete fdyvsy;
      }

      // dyvsty -> dz:

      moddyvsty.Fit( "pol1", "q", "", -0.002, 0.002 );
      TF1 * fdyvsty = moddyvsty.GetFunction( "pol1" );
      cout << endl << moddyvsty.GetTitle()
	   << ": z shift " << fdyvsty->GetParameter(1)
	   << " mm"
	   << endl;
      MODz += fdyvsty->GetParameter(1);
      //delete fdyvsty;

    }

    // write new MOD alignment:
    if( writeMODalign ) {
      ofstream MODalignFile( MODalignFileName.str() );
    
      MODalignFile << "# MOD alignment for run " << run << endl;
      ++MODaligniteration;
      MODalignFile << "iteration " << MODaligniteration << endl;
      MODalignFile << "alignx " << newMODalignx << endl;
      MODalignFile << "aligny " << newMODaligny << endl;
      MODalignFile << "rot " << MODrot << endl;
      MODalignFile << "tilt " << MODtilt << endl;
      MODalignFile << "turn " << MODturn << endl;
      MODalignFile << "dz " << MODz - zz[4] << endl;
    
      MODalignFile.close();
    
      cout << endl << "wrote MOD alignment iteration " << MODaligniteration
	   << " to " << MODalignFileName.str() << endl
	   << "  alignx " << newMODalignx << endl
	   << "  aligny " << newMODaligny << endl
	   << "  rot    " << MODrot << endl
	   << "  tilt   " << MODtilt << endl
	   << "  turn   " << MODturn << endl
	   << "  dz     " << MODz - zz[4] << endl
        ;
    
    } // MOD
  }
  
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT alignment:

  cout << endl
       << "DUT alignment fits for iteration " << DUTaligniteration << endl;

  double newDUTalignx = DUTalignx0;
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
      cmsdxcHisto.GetEntries() > 999 ) {

    cout << endl << cmsdxcHisto.GetTitle()
	 << " bin " << cmsdxcHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -0.5, 0.5 );
    fgp0->SetParameter( 0, cmsdxcHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, cmsdxcHisto.GetBinCenter( cmsdxcHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 8*cmsdxcHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, cmsdxcHisto.GetBinContent(1) ); // BG
    cmsdxcHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newDUTalignx = DUTalignx0 + fgp0->GetParameter(1);
    delete fgp0;

    // dxvsy -> -rot for 25x100 on rot90:

    if( cmsdxvsy.GetEntries() > 999 && !fifty && rot90 ) {
      cmsdxvsy.Fit( "pol1", "q", "", -midy[iDUT]+0.2, midy[iDUT]-0.2 );
      TF1 * fdxvsy = cmsdxvsy.GetFunction( "pol1" );
      cout << endl << cmsdxvsy.GetTitle();
      cout << ": extra rot " << upsignx*upsigny*fdxvsy->GetParameter(1); // converges
      cout << endl;
      DUTrot += upsignx*upsigny*fdxvsy->GetParameter(1);
      delete fdxvsy;
    }

    // dxvsx -> turn:

    if( fabs(DUTturn) > 0.4 && cmsdxvsx.GetEntries() > 999 ) {

      cmsdxvsx.Fit( "pol1", "q", "", -midx[iDUT]+0.2, midx[iDUT]-0.2 );
      TF1 * fdxvsx = cmsdxvsx.GetFunction( "pol1" );
      cout << endl << cmsdxvsx.GetTitle()
	   << ": slope " << fdxvsx->GetParameter(1)
	   << ", extra turn " << fdxvsx->GetParameter(1)/wt/so
	   << " deg"
	   << endl;
      DUTturn += fdxvsx->GetParameter(1)/wt/so;
      //delete fdxvsx;

    }

    // dxvstx -> dz:

    if( cmsdxvstx.GetEntries() > 999 && !fifty && rot90 ) {
      cmsdxvstx.Fit( "pol1", "q", "", -2, 2 );
      TF1 * fdxvstx = cmsdxvstx.GetFunction( "pol1" );
      cout << endl << cmsdxvstx.GetTitle()
	   << ": z shift " << upsignx*fdxvstx->GetParameter(1)*1E3 // converges
	   << " mm"
	   << endl;
      DUTz += upsignx*fdxvstx->GetParameter(1)*1E3;
      delete fdxvstx;
    }

  } // good x

  // fine y:

  if( DUTaligniteration > 0 && fabs( newDUTaligny - DUTaligny0 ) < 0.1 &&
      cmsdycHisto.GetEntries() > 999 ) {

    cout << endl << cmsdycHisto.GetTitle()
	 << " bin " << cmsdycHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -0.5, 0.5 );
    fgp0->SetParameter( 0, cmsdycHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, cmsdycHisto.GetBinCenter( cmsdycHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 5*cmsdycHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, cmsdycHisto.GetBinContent(1) ); // BG
    cmsdycHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newDUTaligny = DUTaligny0 + fgp0->GetParameter(1);
    delete fgp0;

    // dyvsx -> -rot

    if( cmsdyvsx.GetEntries() > 999 && ( fifty || ! rot90 ) ) {
      cmsdyvsx.Fit( "pol1", "q", "", -midx[iDUT]+0.2, midx[iDUT]-0.2 );
      TF1 * fdyvsx = cmsdyvsx.GetFunction( "pol1" );
      cout << endl << cmsdyvsx.GetTitle();
      cout << ": extra rot " << upsigny*fdyvsx->GetParameter(1);
      cout << endl;
      DUTrot -= upsigny*fdyvsx->GetParameter(1);
      delete fdyvsx;
    }
    
    // dyvsy -> tilt:

    if( fabs(DUTtilt) > 0.4 && cmsdyvsy.GetEntries() > 999 ) {
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

    if( cmsdyvsty.GetEntries() > 999 && ( fifty || ! rot90 ) ) {
      cmsdyvsty.Fit( "pol1", "q", "", -2, 2 );
      TF1 * fdyvsty = cmsdyvsty.GetFunction( "pol1" );
      cout << endl << cmsdyvsty.GetTitle()
	   << ": z shift " << upsigny*fdyvsty->GetParameter(1)*1E3
	   << " mm"
	   << endl;
      DUTz += upsigny*fdyvsty->GetParameter(1)*1E3;
      //delete fdyvsty;
    }

  } // good y

  cout << endl
       << "run " << run << endl
       << "  DUT efficiency " << 100*effvst6.GetMean(2)  << " %" << endl
       << "  DUT fake rate  " << 100*(totalfks/totalncl) << " %" << endl
       << "  from " << effvst6.GetEntries() << " events"
       << endl;

  cout << endl
       << "DUT alignment iteration " << DUTaligniteration << endl
       << "  alignx " << newDUTalignx << endl
       << "  aligny " << newDUTaligny << endl
       << "  rot    " << DUTrot << endl
       << "  tilt   " << DUTtilt << endl
       << "  turn   " << DUTturn << endl
       << "  dz     " << DUTz - zz[2] << endl
    ;

  // write new DUT alignment:

  //cout << "update DUT alignment file? (y/n)" << endl;
  //string ans{"y"};
  //string YES{"y"};
  //cin >> ans;
  //if( ans == YES ) {
  if( writeDUTalign ) {
    cout << "writing DUT alignment" << endl;
      
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

    cout << " to " << DUTalignFileName.str() << endl;
  }  
  else
    cout << "NOT writing DUT alignment" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  return 0;
}

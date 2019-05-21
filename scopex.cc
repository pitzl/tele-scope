
// Daniel Pitzl, DESY, Sep 2017, May 2019
// telescope analysis with eudaq and ROC4Sens and MOD

// event display

// make scopex
// needs runs.dat
// needs align_24500.dat from tele

// scopex 24500
// scopex -l 99999 24500
// scopex -t 2.2 25200
// scopex -s 24093
// scopex -l 99999 28027
// scopex -s 28037
// scopex 31210 # 3D
// scopex -l 312100 31425

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TApplication.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <set> // for hotset
#include <map> // colpx
#include <cmath>
#include <sys/time.h> // gettimeofday, timeval
#include <sys/ioctl.h>

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  double adc; // online ph
  double ph;  // offline dph
  double q;   // calibrated [ke]
  int ord;    // readout order
  bool big;   // module big pixel flag
  bool operator < (const pixel & pxObj ) const // row-wise ordering
  {
    return row < pxObj.row;
  }
};

struct cluster {
  vector <pixel> vpix;
  int size;
  int ncol, nrow;
  double col, row;
  double charge;
  bool big;
  bool iso;
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
    c.iso = 1;
    int minx = 9999;
    int maxx = 0;
    int miny = 9999;
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

} // getClus

//------------------------------------------------------------------------------
bool kbhit()
{
  usleep(1000); // [us]
  int byteswaiting;
  ioctl( 0, FIONREAD, &byteswaiting );
  return byteswaiting > 0;
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

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

    if( !strcmp( argv[i], "-m" ) )
      syncmod = 1;

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  int r4srun{0};
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
    string R4SRUN( "r4srun" );
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

      if( tag == R4SRUN ) {
	tokenizer >> r4srun;
	continue;
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

    } // while getline runs.dat

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

  if( run >= 31148 )  // 4.10.2017
    fDTB = 39.99712E6; // ddtvsdt->Fit("pol1")  fDTB*(1-slope*1E-3)

  if( run >= 31592 )  // 16.12.2017 different DTB
    fDTB = 39.996943E6; // ddtvsdt->Fit("pol1")  fDTB*(1-slope*1E-3)

  if( run >= 31632 )  // Feb 2018
    fDTB = 39.996984E6; // ddtvsdt->Fit("pol1")  fDTB*(1-slope*1E-3)

  if( run >= 32032 )  // Apr 2018
    fDTB = 39.997110E6; // 40 MHz DTB clock from ddt and ddtvsdt: fDTB*(1-slope*1E-3)

  if( run >= 33511 )  // Sep 2018
    fDTB = 39.996810E6; // 40 MHz DTB clock from ddt and ddtvsdt: fDTB*(1-slope*1E-3)

  if( run >= 34248 )  // Oct 2018
    fDTB = 39.996899E6; // ddtvsdt->Fit("pol1")  fDTB*(1-slope*1E-3)

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

  } // telescope alignFile

  ialignFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Mimosa hot pixels:

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

  } // telescope hotFile

  ihotFile.close();

  for( int ipl = 0; ipl < 6; ++ipl )
    cout << ipl << ": hot " << hotset[ipl].size() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT:

  const double wt = atan(1.0) / 45.0; // pi/180 deg
    
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

  if( chip0 == 194 ) rot90 = 1; // 4E15 n
  if( chip0 == 195 ) rot90 = 1; // 0.6E15 n
  if( chip0 == 196 ) rot90 = 1; // 8E15 n
  if( chip0 == 197 ) rot90 = 1; // 16E15 n
  if( chip0 == 200 ) rot90 = 1; // 16E15 n
  if( chip0 == 206 ) rot90 = 1; // p-spray max imp 8E15 n
  if( chip0 == 207 ) rot90 = 1; // 16E15 n
  if( chip0 == 210 ) rot90 = 1; // FDD 4E15 n
  if( chip0 == 211 ) rot90 = 1; // FDD 0.6E15 n
  if( chip0 == 216 ) rot90 = 1; // FDD 8E15 n

  if( chip0 == 227 ) rot90 = 1;

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
  if( chip0 == 152 ) fifty = 1;
  if( chip0 == 155 ) fifty = 1;
  if( chip0 == 158 ) fifty = 1;
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

  if( chip0 == 198 ) fifty = 1; // 8E15 n
  if( chip0 == 199 ) fifty = 1; // 0.6E15 n
  if( chip0 == 200 ) fifty = 1; // 16E15 n
  if( chip0 == 202 ) fifty = 1; // 4E15 n
  if( chip0 == 203 ) fifty = 1; // 0.6E15 n
  if( chip0 == 204 ) fifty = 1; // 8E15 n
  if( chip0 == 212 ) fifty = 1; // FDD 4E15 n
  if( chip0 == 214 ) fifty = 1; // FDD 0.6E15 n
  if( chip0 == 215 ) fifty = 1; // FDD 8E15 n

  if( chip0 == 225 ) fifty = 1;
  if( chip0 == 226 ) fifty = 1;
  if( chip0 == 229 ) fifty = 1;
  if( chip0 == 230 ) fifty = 1;

  if( chip0 >= 300 ) fifty = 1; // 3D

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

  if( ( chip0 == 174 || chip0 == 179 || chip0 == 193 || chip0 == 109 ||
	chip0 == 115 || chip0 == 120 || chip0 == 108 ) &&
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
  if( run >= 35200 ) { // Feb 2019, cable from above
    upsignx = -1;
    upsigny = -1;
  }
  if( run >= 35600 ) { // Mar 2019: cable from below
    upsignx = 1;
    upsigny = 1;
  }
  if( !rot90 && run >= 36200 ) { // Apr 2019: cable from below
    upsignx =-1; // why?
    upsigny = 1;
  }

  cout << endl;
  cout << "upsignx " << upsignx << endl;
  cout << "upsigny " << upsigny << endl;

  int iDUT = 7;

  int DUTaligniteration = 0;
  double DUTalignx = 0.0;
  double DUTaligny = 0.0;
  double DUTrot = 0.0;
  double DUTtilt = DUTtilt0; // [deg]
  double DUTturn = DUTturn0; // [deg]
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

  if( DUTaligniteration <= 1 ) {
    DUTtilt = DUTtilt0; // from runs.dat
    DUTturn = DUTturn0; // from runs.dat
  }

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

  double ke = 0.037; // Landau peak at 11 ke  chip 102  1002.dat
  //if( run >= 31148 ) ke = 0.039; // chip 111  1004.dat
  if( run >= 31148 ) ke = 0.0426; // chip 111  1004.dat -V0
  if( run >= 31163 ) ke = 0.035; // chip 117  1006.dat
  if( run >= 31173 ) ke = 0.037; // chip 117
  if( run >= 31210 ) ke = 0.068; // chip 332  3D 230 um at 17.2
  if( run >= 31237 ) ke = 0.050; // chip 352  3D 230 um at 17.2
  if( run >= 31295 ) ke = 0.040; // chip 143  at 11 ke
  if( run >= 31304 ) ke = 0.0385; // chip 144  at 11 ke
  if( run >= 31309 ) ke = 0.0274; // chip 142  at 11 ke
  //if( run >= 31351 ) ke = 0.025; // chip 142  at 11 ke
  if( run >= 31351 ) ke = 0.0282; // chip 142  at 11 ke -V0
  //if( run >= 31355 ) ke = 0.039; // chip 144  at 11 ke
  if( run >= 31355 ) ke = 0.0438; // chip 144  at 11 ke -V0
  if( run >= 31390 ) ke = 0.031; // chip 102  at 11 ke
  if( run >= 31425 ) ke = 0.036; // chip 144  at 11 ke
  //if( run >= 31473 ) ke = 0.0296; // chip 144  at 11 ke
  if( run >= 31473 ) ke = 0.0319; // chip 144  at 11 ke -V0
  if( run >= 31488 ) ke = 0.035; // chip 160  at 11 ke
  if( run >= 31512 ) ke = 0.068; // chip 332  3D 230 um
  if( run >= 31562 ) ke = 0.050; // chip 352  3D 230 um
  if( run >= 31592 ) ke = 0.031; // chip 152  deep diff 180 um at 11 ke
  if( run >= 31610 ) ke = 0.035; // chip 160  at 11 ke
  if( run >= 31614 ) ke = 0.030; // chip 152  deep diff 180 um at 11 ke
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
  if( run >= 32463 ) ke = 0.068; // chip 332  3D 230 um at 17.2
  if( run >= 32498 ) ke = 0.0367; // default
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
  // ROOT:

  gStyle->SetTextFont(62); // 62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetTitleFont( 62, "XYZ" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.5, "y" );
  gStyle->SetTitleOffset( 1.7, "z" );

  gStyle->SetLabelFont( 62, "XYZ" );

  gStyle->SetLabelOffset( 0.022, "xyz" );

  gStyle->SetTickLength( -0.02, "xyz" ); // tick marks outside

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  // statistics box:

  gStyle->SetOptStat(10);
  gStyle->SetStatFont(42); // 42 = Helvetica normal
  gStyle->SetStatBorderSize(1); // no 'shadow'
  gStyle->SetStatX(0.82);
  gStyle->SetStatY(0.92);

  gStyle->SetPalette(55); // sunset
  //gStyle->SetPalette(56); // white to red
  //gStyle->SetPalette(73); // blue
  //gStyle->SetPalette(90); // green to magenta
  //gStyle->SetPalette(109); // sea blue to magenta
  //gStyle->SetNumberContours(32); // -20..300

  TApplication app( "app", 0, 0 );
  TCanvas c1( "c1", "pixel event display", 900, 800 ); // square
  c1.SetTopMargin( 0.12 );
  c1.SetRightMargin( 0.18 );
  gPad->SetGrid();

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

  string START {"START"};
  string hd;
  while( hd != START ) {
    getline( evFile, hd ); // read one line into string
    cout << "  " << hd << endl;
  }
  bool readnext = 1;
  string DUTev; // May 2019: must be declared outside event loop, else empty reads

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

  bool more = 1; // event displays
  string Q{"q"};

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

    if( lev < 212 )
      ldb = 1;

    uint64_t tlutime = evt.GetTimestamp(); // 384 MHz = 2.6 ns
    if( iev < 2  )
      tlutime0 = tlutime;
    double evsec = (tlutime - tlutime0) / fTLU;

    double tludt = (tlutime - prevtlutime) / fTLU; // [s]

    prevtlutime = tlutime;

    if( iev < 10 )
      cout << "scopex processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "scopex processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "scopex processing  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%1000 == 0 )
      cout << "scopex processing  " << run << "." << iev << "  taken " << evsec << " dt " << tludt << endl;

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

	// fill pixel block for clustering:

	pixel px;
	px.col = ix; // col
	px.row = iy; // row
	px.adc = adc;
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

      if( ldb ) std::cout << std::endl;

      // clustering:

      cl[ipl] = getClus(pb);

      if( ldb ) cout << "clusters " << cl[ipl].size() << endl;

    } // eudaq planes

    gettimeofday( &tv, NULL );
    long s2 = tv.tv_sec; // seconds since 1.1.1970
    long u2 = tv.tv_usec; // microseconds
    zeit1 += s2 - s1 + ( u2 - u1 ) * 1e-6; // read and cluster

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // read R4S DUT stream:

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
      cout << endl
	   << iev
	   << " DUT " << DUTev
	   << " TLU " << tludt*1E3
	   << " - DTB " << dtbdt*1e3 << " ms"
	   << " = " << (tludt-dtbdt)*1e6 << " us"
	   << ", R4S resync " << nresync
	   << endl;
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
	cout << iev
	     << " DUT " << DUTev
	     << " TLU " << tludt*1E3
	     << " - DTB " << dtbdt*1e3 << " ms"
	     << " = " << (tludt-dtbdt)*1e6 << " us"
	     << endl;
      }
      else { // added event
	if( ldbt )
	  cout << "\t DTB added" << endl;
	break;
      }

      //break; // only once

    } // tludt

    vector <pixel> pb; // for clustering

    map < int, set <pixel> > colpxmap; // per col, sorted along row

    if( readnext && filled == F ) {

      string roi;
      getline( evFile, roi );
      istringstream roiss( roi ); // tokenize string

      int ipx = 0;
      vector <pixel> vpx;
      vpx.reserve(49);

      if( ldb ) cout << "  px";

      while( ! roiss.eof() ) {

	int col;
	int row;
	double ph;
	roiss >> col; // ROC
	roiss >> row; // ROC
	roiss >> ph;
	if( ldb ) cout << " " << col << " " << row << " " << ph;

	pixel px { col, row, ph, ph, ph, ipx, 0 };
	vpx.push_back(px);
	++ipx;

      } // roi px

      // column-wise common mode correction:

      set <pixel> compx[156]; // per column, sorted along row

      for( unsigned ipx = 0; ipx < vpx.size(); ++ipx ) {
	int col = vpx[ipx].col;
	pixel px { col, vpx[ipx].row,
		   vpx[ipx].adc, vpx[ipx].ph, vpx[ipx].q,
		   vpx[ipx].ord, 0 };
	compx[col].insert(px); // sorted: by row
      }

      for( unsigned col = 0; col < 155; ++col ) {

	if( compx[col].size() < 2 ) continue;

	set <pixel> colpx; // per col, sorted along row

	auto px1 = compx[col].begin();
	auto px7 = compx[col].end(); --px7; // last row

	double ph1 = px1->ph;
	double ph7 = px7->ph;

	auto px4 = px1;
	double phprev = px4->ph;
	++px4;

	for( ; px4 != px7; ++px4 ) { // between 1 and 7, exclusively

	  int col4 = px4->col;
	  int row4 = px4->row;
	  double ph4 = px4->ph;

	  if( chip0 == 118 && run <= 32262 ) // gain2
	    ph4 -= 0.146*phprev; // Tsunami

	  if( chip0 == -198 )
	    ph4 -= 0.13*phprev; // Tsunami, worse resolution

	  if( chip0 == 179 )
	    ph4 += 0.03*phprev; // anti-Tsunami, see roipvsdxdy

	  double dph = ph4 - 0.5*(ph1+ph7); // 20.11.2018 less noise

	  phprev = ph4;

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
	  px.adc = ph4;
	  px.ph = dph;
	  px.q = q;
	  px.ord = pb.size(); // readout order, starts at zero
	  px.big = 0;

	  colpx.insert(px); // sorted along col

	  double dphcut = 12; // gain_1 2017

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
	  if( chip0 == 137 && run >= 36223 ) dphcut = 20; // irrad, gain_1 rms 7
	  if( chip0 == 156 ) dphcut = 16; // gain_1 
	  if( chip0 == 198 ) dphcut = 20; // irrad, gain_1
	  if( chip0 == 201 ) dphcut = 20; // irrad, RMS 6 ADC
	  if( chip0 > 300 ) dphcut = 20; // irradiated 3D is noisy

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

	    pb.push_back(px);

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

    // DUT clustering:

    cl[iDUT] = getClus(pb);

    if( ldb ) cout << "  DUT cl " << cl[iDUT].size() << endl << flush;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT:

     for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

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

      // cluster isolation:

      vector<cluster>::iterator d = c;
      ++d;
      for( ; d != cl[iDUT].end(); ++d ) {

	bool done = 0;

	for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	  for( vector<pixel>::iterator qx = d->vpix.begin(); qx != d->vpix.end(); ++qx ) {

	    if( abs( px->col - qx->col ) < 4 &&
		abs( px->row - qx->row ) < 4 ) {

	      if( c->charge < d->charge )
		c->iso = 0; // flag smaller cluster
	      else
		d->iso = 0;

	      done = 1; // once is enough
	      break; // qx

	    }

	  } // qx

	  if( done ) break;

	} // px

      } // d

    } // DUT clusters

    gettimeofday( &tv, NULL );
    long s3 = tv.tv_sec; // seconds since 1.1.1970
    long u3 = tv.tv_usec; // microseconds
    zeit2 += s3 - s2 + ( u3 - u2 ) * 1e-6; // read and cluster DUT

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if( run == 31175 && iev >= 1612000 )
      syncmod = 1;

    if( run == 31351 && evsec >= 17100 )
      syncmod = 1; // does not help

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

	} // cl B

      } // cl C

    } // cl A

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // driplets vs MOD:

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

      } // jj

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

	int npx = c->size;

	double moddx = modx - x4;
	double moddy = mody - y4;

	if( fabs( moddx ) < xcutMOD &&
	    fabs( moddy ) < ycutMOD ) {

	  driplets[jB].lk = 1;

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

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make triplets 2+0-1:

    vector <triplet> triplets;

    //double triCut = 0.1; // [mm]
    double triCut = 0.05; // [mm] 2.10.2017

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

	} // cl B

      } // cl C

    } // cl A

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

      } // jj

      triplets[iA].ttdmin = ttdmin;

      if( ttdmin < 0.2 ) continue; // track isolation at DUT

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

      double x2 = x1;
      double y2 = ca*y1 + sa*z1; // tilt a, sign cancels out
      double z2 =-sa*y1 + ca*z1; // should be zero (in DUT plane). is zero

      double x3 = cf*x2 + sf*y2; // rot
      double y3 =-sf*x2 + cf*y2;
      double z3 = z2; // should be zero (in DUT plane). is zero

      double x4 = upsignx*x3 + DUTalignx; // shift to mid
      double y4 = upsigny*y3 + DUTaligny; // invert y, shift to mid

      double x8 = x4;
      double y8 = y4;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // match triplet and driplet:

      bool lsixlk = 0;
      double sixcut = 0.1; // [mm] 502 in 25463 eff 99.87
      double dddmin = 99.9; // driplet isolation at MOD

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

	// match:

	if( fabs(dx) < sixcut && fabs(dy) < sixcut ) {

	  // average driplet and triplet at DUT:

	  double xa = 0.5 * ( xB + xc );
	  double ya = 0.5 * ( yB + yc );

	  if( driplets[jB].lk ) { // driplet linked to MOD
	    lsixlk = 1;
	    dddmin = driplets[jB].ttdmin;
	  }

	  // transform into DUT system:

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

	  if( ! rot90 ) { // straight, no Cu behind DUT
	    x4 = x8; // sixplet interpol
	    y4 = y8;
	  }

	} // six match

      } // driplets

      if( ! lsixlk ) continue; //only in-time six-tracks !

      if( dddmin < 0.3 ) continue; // require isolation at MOD

      if( fabs( x4 ) > 3.8 ) continue; // fiducial
      if( fabs( y4 ) > 3.9 ) continue;

      // reduce to 100x100 um region:

      double xmod = fmod( 9.000 + x4, 0.1 ); // [0,0.1] mm
      double ymod = fmod( 9.000 + y4, 0.1 ); // [0,0.1] mm
      if( chip0 == 142 || chip0 == 143 || chip0 == 144 )
	ymod = fmod( 9.050 + y4, 0.1 ); // [0,0.1] mm
      double xmod2 = fmod( 9.000 + x4, 0.2 ); // [0,0.2] mm
      double ymod2 = fmod( 9.000 + y4, 0.2 ); // [0,0.2] mm
      double xmod5 = fmod( 9.000 + x4, 0.05 );
      double ymod5 = fmod( 9.000 + y4, 0.05 );

      // from track x, y (at DUT) to sensor col, row:

      int kcol = x4 / ptchx[iDUT] + 0.5*nx[iDUT] - 1; // straight 50x50
      int krow = y4 / ptchy[iDUT] + 0.5*ny[iDUT];
      if( rot90 ) { // 25x100
	kcol = y4 / ptchx[iDUT] + 0.5*nx[iDUT]; // 0..398 even
	krow = x4 / ptchy[iDUT] + 0.5*ny[iDUT]; // 0..191 full
      }

      if( kcol < 0 ) kcol = 0;
      if( kcol >= nx[iDUT] ) kcol = nx[iDUT]-1;

      if( krow < 0 ) krow = 0;
      if( krow >= ny[iDUT] ) krow = ny[iDUT]-1;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT pixel clusters:

      int nm[99] = {0};

      for( vector<cluster>::iterator c = cl0[iDUT].begin(); c != cl0[iDUT].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;

	double Q0 = c->charge * norm; // cluster charge normalized to vertical incidence

	if( Q0 < 4 ) continue; // fresh

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

	  pcol[icol] += px->ph; // project cluster onto cols
	  prow[irow] += px->ph; // project cluster onto rows
	  qcol[icol] += q; // project cluster onto cols
	  qrow[irow] += q; // project cluster onto rows

	  sumph += px->ph;

	} // pix

	int ncol = colmax - colmin + 1;
	int nrow = rowmax - rowmin + 1;

	double P0 = sumph*norm;

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
	if( q12 > 1 ) eta = ( q1 - q2 ) / q12; // always negative
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

	// residuals:

	double cmsdx = cmsx - x4; // triplet extrapol
	double cmsdy = cmsy - y4;
	double cmsdu = cmsu - x4; // triplet extrapol
	double cmsdv = cmsv - y4;

	if( fabs( cmsdx ) < xcut ) {

	  if( fabs( cmsdy ) > ycut ) { // outlier

	    cout << "  " << iev
		 << ": ncl " << cl0[iDUT].size()
		 << ", dy " << cmsdy*1E3
		 << ", nrow " << nrow
		 << ", ncol " << ncol
		 << ", ph0 " << P0
		 << endl;

	    if( more && cl0[iDUT].size() == 1 ) {

	      // event display: pixels around track

	      TH2D hpxmap( "pxmap",
			   Form( "pixel map %i;col;row;PH [ADC]", iev ),
			   21, -10.5, 10.5, 21, -10.5, 10.5 );
	      hpxmap.GetXaxis()->SetNdivisions(211);
	      hpxmap.GetYaxis()->SetNdivisions(211);
	      hpxmap.SetMaximum(400);

 	      for( auto px : c->vpix ) { // C++17

		//hpxmap.Fill( px->col - kcol, px->row - krow, px->ph );
		hpxmap.Fill( px.col - kcol, px.row - krow, px.ph );

	      } // pix

	      hpxmap.Draw( "colz" );
	      c1.Update();

	      cout << "enter any key, q to stop" << endl;

	      while( !kbhit() )
		gSystem->ProcessEvents(); // ROOT

	      string any;
	      cin >> any;
	      if( any == Q )
		more = 0;

	    } // more

	  } // y cut

	} // x cut

      } // loop DUT clusters

      if( ldb ) cout << "    eff " << nm[49] << endl << flush;

    } // loop triplets iA

    gettimeofday( &tv, NULL );
    long s4 = tv.tv_sec; // seconds since 1.1.1970
    long u4 = tv.tv_usec; // microseconds
    zeit3 += s4 - s3 + ( u4 - u3 ) * 1e-6; // tracking

    if( ldb ) cout << "done ev " << iev << endl << flush;

    ++iev;

    if( syncmod ) { // shift all but MOD

      for( int ipl = 0; ipl < 6; ++ipl ) {
	cl0[ipl] = cl1[ipl]; // remember
	cl1[ipl] = cl[ipl]; // remember
      }

      cl0[iDUT] = cl1[iDUT]; // remember
      cl1[iDUT] = cl[iDUT]; // remember

    }

  } while( reader->NextEvent() && iev < lev );

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
  // done

  cout << endl;

  return 0;
}

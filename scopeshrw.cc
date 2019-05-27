
// Daniel Pitzl, DESY, Sep 2017
// telescope analysis with eudaq and ROC4Sens shallow

// row-wise readout

// make scopeshrw
// needs runs.dat
// needs align_24500.dat from tele

// scopesh 31110
// scopesh 31336  120 V
// scopesh 31342   20 V
// scopesh 31343   10 V

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
#include <cmath>

using namespace std;
using namespace eudaq;

class pixel {
 public:
  int col;
  int row;
  double ph;
  int ord;
  bool big;
  bool operator < (const pixel & pxObj ) const // col-wise ordering
  {
    return col < pxObj.col;
  }
};

struct cluster {
  vector <pixel> vpix; // Armin Burgmeier: list
  int size;
  int ncol, nrow;
  double col, row;
  double ph;
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
    double sumph = 0;
    c.big = 0;
    int minx = 9999;
    int maxx = 0;
    int miny = 9999;
    int maxy = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      double phpix = p->ph; // [ADC]
      sumph += phpix;
      c.col += (*p).col*phpix;
      c.row += (*p).row*phpix;
      if( p->big ) c.big = 1;
      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( sumph > 0 ) {
      c.col /= sumph;
      c.row /= sumph;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with non-positive ph" << endl;
    }

    c.ph = sumph;
    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;

    v.push_back(c); // add cluster to vector

    // look for a new seed = unused pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  // nothing left, return clusters

  delete [] gone;
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

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  string geoFileName( "geo.dat" );
  double DUTtilt0 = 19.3;
  double pbeam = 4.8;
  int chip0 = 110;

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

  } // alignFile

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

  } // hotFile

  ihotFile.close();

  for( int ipl = 0; ipl < 6; ++ipl )
    cout << ipl << ": hot " << hotset[ipl].size() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT:

  const double log10 = log(10);
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

  if( chip0 == 227 ) rot90 = 1; // KIT

  if( rot90 ) {
    cout << "DUT 90 degree rotated" << endl;
    return 90;
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

  if( chip0 == 118 ) upsigny = -1;

  if( run >= 36367 ) upsignx = -1; // May 2019

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

  if( DUTaligniteration == 0 )
    DUTtilt = DUTtilt0; // from runs.dat

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

  double xcut = 0.1;
  if( run == 31336 )
    xcut = 0.25; // skewed?
  if( chip0 == 118 )
    xcut = 0.05;

  double ycut = 0.2;
  if( fabs(DUTtilt) > 60 )
    ycut = 0.5;

  int rowcut = 0.3 * fabs( tan(DUTtilt*wt) ) * 0.15 / ptchy[iDUT];
  if( rowcut < 1 ) rowcut = 1;

  double dphcut = 12; // fresh

  if( chip0 >= 119 && chip0 <= 138 )
    dphcut = 25; // irrad1

  if( chip0 == 198 )
    dphcut = 20; // irrad dp2 8

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "scopeshrw" << run << ".root";

  TFile histoFile( rootFileName.str(  ).c_str(  ), "RECREATE" );

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
			    "ROC pixel hit map;ROC col;ROC row;ROC pixels",
			    155, 0, 155, 160, 0, 160 );
    }
    else {
      hcol[ipl] = TH1I( Form( "col%i", ipl ),
			Form( "%i col;col;plane %i pixels", ipl, ipl ), 
			max( 155, nx[ipl]/4 ), 0, nx[ipl] );
      hrow[ipl] = TH1I( Form( "row%i", ipl ),
			Form( "%i row;row;plane %i pixels", ipl, ipl ),
			max( 160, ny[ipl]/2 ), 0, ny[ipl] );
      hmap[ipl] = new TH2I( Form( "map%i", ipl ),
			    Form( "%i map;col;row;plane %i pixels", ipl, ipl ),
			    max( 155, nx[ipl]/4 ), 0, nx[ipl], max( 160, ny[ipl]/2 ), 0, ny[ipl] );
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

  TProfile dutphvsprev( "dutphvsprev", "Tsunami;previous PH [ADC];<PH> [ADC]", 120, -100, 500 );
  TProfile dutdphvsprev( "dutdphvsprev", "Tsunami;previous #DeltaPH [ADC];<#DeltaPH> [ADC]", 120, -100, 500 );
  TH1I dutphHisto( "dutph", "DUT PH;ADC-PED [ADC];pixels", 500, -100, 900 );
  TH1I dutdphHisto( "dutdph", "DUT #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );
  TH1I dutdph0Histo( "dutdph0", "DUT #DeltaPH;#DeltaPH [ADC];pixels", 500, -100, 900 );
  TH1I dutdpiHisto( "dutdpi", "DUT -#DeltaPH;-#DeltaPH [ADC];pixels", 1000, -200, 800 );
  TH1I dutdp2Histo( "dutdp2", "DUT #DeltaPH_{2};#DeltaPH_{2} [ADC];pixels", 1000, -200, 800 );

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

  TH1I dutadcHisto( "dutadc",
		    "DUT pixel ADC;pixel pulse height [ADC];pixels",
		    400, 0, 800 );
  TH1I dutcolHisto( "dutcol",
		    "DUT pixel column;pixel column;pixels",
		    nx[iDUT], -0.5, nx[iDUT]-0.5 );
  TH1I dutrowHisto( "dutrow",
		    "DUT pixel row;pixel row;pixels",
		    ny[iDUT], -0.5, ny[iDUT]-0.5 );

  TH1I dutph0Histo( "dutph0",
		   "normal fiducial cluster charge;normal fiducial cluster charge [ke];fiducial clusters",
		   240, 0, 1200 );

  TH1I dutnpxHisto( "dutnpx",
		     "DUT cluster size;cluster size [pixels];clusters",
		     52, 0.5, 52.5 );
  TH1I dutncolHisto( "dutncol",
		     "DUT cluster size;cluster size [columns];clusters",
		     52, 0.5, 52.5 );
  TH1I dutnrowHisto( "dutnrow",
		     "DUT cluster size;cluster size [rows];clusters",
		     80, 0.5, 80.5 );

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

  TH1I trixAHisto( "trixA", "triplets x at DUT;x [mm];triplets at DUT",
		   240, -12, 12 );
  TH1I triyAHisto( "triyA", "triplets y at DUT;y [mm];triplets at DUT",
		   120, -6, 6 );
  TH2I * trixyAHisto = new
    TH2I( "trixyA", "triplets x-y at DUT;x [mm];y [mm];triplets at DUT",
	  240, -12, 12, 120, -6, 6 );

  // DUT pixel vs triplets:

  TH1I z3Histo( "z3",
		"z3 should be zero;z3 [mm];triplets",
		100, -0.01, 0.01 );

  TH1I cmsdxpxHisto( "cmsdxpx",
		     "DUT pixel - Telescope dx;pixel - triplet #Deltax [mm];pixels",
		     200, -1, 1 );

  TH2I * cmsxvsx = new TH2I( "cmsxvsx",
			     "DUT vs Telescope x;track x [mm];DUT x [mm];track-cluster combinations",
			     160, -4, 4, 160, -4, 4 );
  TH2I * cmsyvsy = new TH2I( "cmsyvsy",
			     "DUT vs Telescope y;track y [mm];DUT y [mm];track-cluster combinations",
			     160, -4, 4, 160, -4, 4 );

  TH1I cmssxaHisto( "cmssxa",
		    "DUT + Telescope x;cluster + triplet #Sigmax [mm];clusters",
		    440, -11, 11 );
  TH1I cmsdxaHisto( "cmsdxa",
		    "DUT - Telescope x;cluster - triplet #Deltax [mm];clusters",
		    440, -11, 11 );
  TH1I cmsdxHisto( "cmsdx",
		   "DUT - Telescope x;cluster - triplet #Deltax [mm];clusters",
		   200, -0.5, 0.5 );
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
  TH2I * cmsdxvsev1 = new TH2I( "cmsdxvsev1",
				"DUT - Telescope x;event;#Deltax [mm]",
				100, 0, 50000, 100, -5, 5 );
  TH2I * cmsdxvsev2 = new TH2I( "cmsdxvsev2",
				"DUT - Telescope x;event;#Deltax [mm]",
				1000, 0, 1000*1000, 50, -5, 5 );

  TH1I cmsycHisto( "cmsyc",
		   "x-linked triplet at DUT y;triplet y at DUT [mm];x-linked triplets",
		   120, -6, 6 );
  TH1I cmssyaHisto( "cmssya",
		    "DUT + Telescope y;cluster + triplet #Sigmay [mm];clusters",
		    198, -99, 99 ); // shallow needs wide range
  TH1I cmsdyaHisto( "cmsdya",
		    "DUT - Telescope y;cluster - triplet #Deltay [mm];clusters",
		    198, -99, 99 ); // shallow needs wide range
  TH1I cmsdycHisto( "cmsdyc",
		    "#Deltay cut x;cluster - triplet #Deltay [mm];fiducial clusters",
		     200, -2, 2 );
  TProfile cmsdyvsx( "cmsdyvsx",
		     "DUT #Deltay vs x;x track [mm];<cluster - triplet #Deltay> [mm]",
		     50, -3.75, 3.75, -1, 1 );
  TProfile cmsdyvsy( "cmsdyvsy",
		     "DUT #Deltay vs y;y track [mm];<cluster - triplet #Deltay> [mm]",
		     80, -8, 8, -1, 1 );
  TProfile cmsdyvsty( "cmsdyvsty",
		      "DUT #Deltay vs #theta_{y};y track slope [rad];<cluster - triplet #Deltay> [mm]",
		      80, -0.002, 0.002, -1, 1 );

  TH1I trixAlkHisto( "trixAlk",
		    "linked triplet at DUT x;triplet x at DUT [mm];linked triplets",
		    240, -12, 12 );
  TH1I triyAlkHisto( "triyAlk",
		    "linked triplet at DUT y;triplet y at DUT [mm];linked triplets",
		    120, -6, 6 );
  TH2I * trixyAlkHisto = new
    TH2I( "trixyAlk", "triplets x-y at DUT;x [mm];y [mm];triplets at DUT",
	  240, -12, 12, 120, -6, 6 );

  TH1I cmscolHisto( "cmscol",
		    "DUT linked columns;DUT linked cluster column;linked clusters",
		    nx[iDUT], -0.5, nx[iDUT]-0.5 );
  TH1I cmsrowHisto( "cmsrow",
		    "DUT linked rows;DUT linked cluster row;linked clusters",
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

  TProfile cmsnrowvsy( "cmsnrowvsy",
		       "DUT cluster length vs y;y track [mm];<cluster length> [rows]",
		       100, -5, 5 );
  TProfile cmsnrowvsy0( "cmsnrowvsy0",
		       "DUT cluster length vs y;y track entry [mm];<cluster length> [rows]",
		       100, -5, 5 );
  TProfile cmsnrowvsy9( "cmsnrowvsy9",
		       "DUT cluster length vs y;y track exit [mm];<cluster length> [rows]",
		       100, -5, 5 );

  TProfile2D * cmsnpxvsxmym = new
    TProfile2D( "cmsnpxvsxmym",
	      "DUT cluster size vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];<cluster size> [pixels]",
		40, 0, 0.2, 40, 0, 0.2, 0, 20 );

  TH1I cmsph0Histo( "cmsph0",
		   "normal fiducial cluster charge;normal cluster charge [ADC];linked fiducial clusters",
		   240, 0, 1200 );

  TH1I cmspxphHisto( "cmspxph",
		    "DUT pixel charge linked;pixel charge [ADC];linked pixels",
		    200, 0, 400 );
  TH1I cmspxphlHisto( "cmspxphl",
		    "DUT pixel charge long linked;pixel charge [ADC];pixels in long linked clusters",
		    200, 0, 400 );

  TH1I npxrdHisto = TH1I( "npxrd",
			   "pixels in track road;pixels in track road;tracks",
			   120, 0.5, 120.5 );
  TH1I nrwrdHisto = TH1I( "nrwrd",
			   "rows in track road;rows in track road;tracks",
			   60, 0.5, 60.5 );

  TH1I dcol0Histo = TH1I( "dcol0",
			  "cluster-road 1st pixel;1st cluster-road [pixels];roads",
			  41, -20.5, 20.5 );
  TH1I dcol9Histo = TH1I( "dcol9",
			  "cluster-road lst pixel;lst cluster-road [pixels];roads",
			  41, -20.5, 20.5 );
  TH1I ncolrdHisto = TH1I( "ncolrd",
			   "columns in track road;pixel columns in track road;roads",
			   80, 0.5, 80.5 );

  TH1I dutpxphmidHisto =
    TH1I( "dutpxphmid",
	  "DUT pixel charge middle;pixel charge [ADC];middle pixels",
	  200, 0, 400 );
  TH1I dutpxphedgHisto =
    TH1I( "dutpxphedg",
	  "DUT pixel charge edge;pixel charge [ADC];edge pixels",
	  200, 0, 400 );

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

  TH1I npxrdphHisto = TH1I( "npxrdph",
			  "charged pixels in track road;charged pixels in track road;tracks",
			  80, 0.5, 80.5 );

  TH1I rowpHisto =
    TH1I( "rowp",
	  "DUT y charge;y charge [ADC];y bins in track road",
	  200, 0, 400 );
  TH1I rowp1Histo =
    TH1I( "rowp1",
	  "DUT 1st y charge;y charge [ADC];1st y bin in track road",
	  200, 0, 400 );
  TH1I rowp2Histo =
    TH1I( "rowp2",
	  "DUT 2nd y charge;y charge [ADC];2nd y bin in track road",
	  200, 0, 400 );
  TH1I rowp8Histo =
    TH1I( "rowp8",
	  "DUT last-1 y charge;y charge [ADC];lst-1 y bin in track road",
	  200, 0, 400 );
  TH1I rowp9Histo =
    TH1I( "rowp9",
	  "DUT last y charge;y charge [ADC];lst y bin in track road",
	  200, 0, 400 );
  TH1I rowpmHisto =
    TH1I( "rowpm",
	  "DUT mid y charge;y charge [ADC];mid y bin in track road",
	  200, 0, 400 );

  TH1I rowdyHisto =
    TH1I( "rowdy",
	  "DUT row dy;#Deltay [mm];rows",
	 160, -2, 2 );

  TH1I rowdHisto =
    TH1I( "rowd",
	  "track depth;track depth [#mum];pixels in track road",
	  100, -100, 100 );
  TH1I rowdphHisto =
    TH1I( "rowdph",
	  "track depth charged;track depth [#mum];charged pixels in track road",
	  100, -100, 100 );

  TProfile rowpvsdy =
    TProfile( "rowpvsdy",
	      "DUT charge vs track y;dy [#mum];<charge> [ADC]",
	      240, -1200, 1200 );
  TProfile rowpvsd =
    TProfile( "rowpvsd",
	      "DUT charge vs track depth;track depth [#mum];<pixel charge> [ADC]",
	      40, -100, 100 );
  TProfile rowpvsdxp =
    TProfile( "rowpvsdxp",
	      "DUT charge vs track depth, x > 0;track depth [#mum];<pixel charge> [ADC]",
	      40, -100, 100 );
  TProfile rowpvsdxm =
    TProfile( "rowpvsdxm",
	      "DUT charge vs track depth, x < 0;track depth [#mum];<pixel charge> [ADC]",
	      40, -100, 100 );
  TProfile rowp1vsd =
    TProfile( "rowp1vsd",
	      "DUT charge > 0 vs track depth;track depth [#mum];<pixel charge> [ADC]",
	      40, -100, 100 );

  TH2I * rowpdHisto = new
    TH2I( "rowpd",
	  "DUT charge vs track depth;track depth [#mum];pixel charge [ADC];entries",
	  40, -100, 100, 100, 0, 200 );

  TH1I rowpp60Histo =
    TH1I( "rowpp60",
	  "DUT pixel charge, d > 60 #mum;pixel charge [ADC];d > 60 #mum pixels",
	  200, 0, 400 );
  TH1I rowpp40Histo =
    TH1I( "rowpp40",
	  "DUT pixel charge, 40 < d < 60 #mum;pixel charge [ADC];40 < d 60 #mum pixels",
	  200, 0, 400 );
  TH1I rowpp20Histo =
    TH1I( "rowpp20",
	  "DUT pixel charge, 20 < d < 40 #mum;pixel charge [ADC];20 < d < 40 #mum pixels",
	  200, 0, 400 );
  TH1I rowpp00Histo =
    TH1I( "rowpp00",
	  "DUT pixel charge, 0 < d < 20 #mum;pixel charge [ADC];0 < d < 20 #mum pixels",
	  200, 0, 400 );
  TH1I rowpm20Histo =
    TH1I( "rowpm20",
	  "DUT pixel charge, -20 < d < 0 #mum;pixel charge [ADC];-20 < d < 0 #mum pixels",
	  200, 0, 400 );
  TH1I rowpm40Histo =
    TH1I( "rowpm40",
	  "DUT pixel charge, -40 < d < -20 #mum;pixel charge [ADC];-40 < d < -20 #mum pixels",
	  200, 0, 400 );
  TH1I rowpm60Histo =
    TH1I( "rowpm60",
	  "DUT pixel charge, -60 < d < -40 #mum;pixel charge [ADC];-60 < d < -40 #mum pixels",
	  200, 0, 400 );
  TH1I rowpm80Histo =
    TH1I( "rowpm80",
	  "DUT pixel charge, d < -60 #mum;pixel charge [ADC];d < -60 #mum pixels",
	  200, 0, 400 );

  TProfile2D * roipvsdxdy = new
    TProfile2D( "roipvsdxdy",
		"pixel signals at track;pixel - x track [#mum];pixel - y track [#mum];<pixel signal> [ADC]",
		50, -250, 250, 200, -1200, 800 );

  TProfile roirowpvsdy( "roirowpvsdy",
			"pixel signals at track;pixel - y track [#mum];<pixel signal> [ADC]",
			240, -1200, 1200 );

  TH1I roirowpHisto( "roirowp", "DUT row charge;row charge [ADC];rows", 200, -100, 700 );
  TH1I roinrowHisto( "roinrow", "DUT filled rows;filled rows;roads", 50, 0.5, 50.5 );
  TH1I roirow0Histo( "roirow0", "DUT 1st row;1st row;roads", 160, -0.5, 159.5 );
  TH1I roirow9Histo( "roirow9", "DUT lst row;lst row;roads", 160, -0.5, 159.5 );
  TH1I roirowszHisto( "roirowsz", "DUT row size;row size [cols];rows", 21, -0.5, 20.5 );

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

    if( chip0 == 114 && run == 31095 && iev == 145447 ) tludt = 0.001; // see cmsdxvsev2 and use ldbt
    if( chip0 == 114 && run == 31095 && iev == 146442 ) tludt = 0.001;
    if( chip0 == 114 && run == 31095 && iev == 147437 ) tludt = 0.001;

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

	hcol[ipl].Fill( ix+0.5 );
	hrow[ipl].Fill( iy+0.5 );
	hmap[ipl]->Fill( ix+0.5, iy+0.5 );

	// fill pixel block for clustering:

	pixel px;
	px.col = ix; // col
	px.row = iy; // row
	px.ph = ph;
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

    set <pixel> rowpx[160]; // per row = along track for straight PCB, sorted along col

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

	pixel px { col, row, ph, ipx, 0 };
	vpx.push_back(px);
	++ipx;

      } // roi px

      if( vpx.size() > 1199 ) {
	cout << "R4S vpx " << vpx.size() << " cleared in event " << iev << endl;
	vpx.clear();
      }

      // columns-wise common mode correction:

      set <pixel> compx[160]; // per row, sorted along col

      for( unsigned ipx = 0; ipx < vpx.size(); ++ipx ) {
	int row = vpx[ipx].row;
	pixel px { vpx[ipx].col, row, vpx[ipx].ph, vpx[ipx].ord, 0 };
	compx[row].insert(px); // sorted along row
      }

      for( unsigned row = 0; row < 160; ++row ) {

	if( compx[row].size() < 2 ) continue;

	auto px1 = compx[row].begin();
	auto px7 = compx[row].end(); --px7; // last

	int col1 = px1->col;
	int col7 = px7->col;
	double ph1 = px1->ph;
	double ph7 = px7->ph;

	auto px4 = px1;
	double phprev = px4->ph;
	double dphprev = 0;
	++px4;
	for( ; px4 != px7; ++px4 ) { // between 1 and 7, exclusively

	  int col4 = px4->col;
	  int row4 = px4->row;
	  double ph4 = px4->ph;

	  if( chip0 == -198 )
	    ph4 -= 0.13*phprev; // Tsunami or delta rays

	  double dph;
	  if( col4 - col1 < col7 - col4 )
	    dph = ph4 - ph1;
	  else
	    dph = ph4 - ph7;

	  dutphHisto.Fill( ph4 );
	  dutphvsprev.Fill( phprev, ph4 ); // Tsunami
	  dutdphHisto.Fill( dph );
	  dutdpiHisto.Fill( -dph );
	  dutdphvsprev.Fill( dphprev, dph );
	  if( col4 > 20 ) dutdph0Histo.Fill( dph ); // skip noisy regions

	  dph = ph4 - 0.5*(ph1+ph7); // 20.11.2018 less noise

	  dutdp2Histo.Fill( dph );

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
	  px.ord = pb.size(); // readout order
	  px.big = 0;

	  rowpx[row4].insert(px); // sorted along col

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

    // DUT clustering:

    hnpx[iDUT].Fill( pb.size() );

    cl[iDUT] = getClus(pb);

    hncl[iDUT].Fill( cl[iDUT].size() );

    for( unsigned icl = 0; icl < cl[iDUT].size(); ++ icl ) {

      hsiz[iDUT].Fill( cl[iDUT][icl].size );

    } // icl

    dutnpxvst2.Fill( evsec, pb.size() );
    dutnclvst2.Fill( evsec, cl[iDUT].size() );

    int DUTyld = 0;
    if( pb.size() ) DUTyld = 1; // no double counting: events with at least one px
    dutyldvst2.Fill( evsec, DUTyld );
    dutyldvst6.Fill( evsec/3600, DUTyld );

    if( ldb ) cout << "  DUT cl " << cl[iDUT].size() << endl << flush;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT:

    for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

      int colmin = 999;
      int colmax = -1;
      int rowmin = 999;
      int rowmax = -1;

      double phsum = 0;
      double sumcol = 0;
      double sumrow = 0;

      double phcol[nx[iDUT]];
      for( int icol = 0; icol < nx[iDUT]; ++icol ) phcol[icol] = 0;

      double phrow[ny[iDUT]];
      for( int irow = 0; irow < ny[iDUT]; ++irow ) phrow[irow] = 0;

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

	double ph = px->ph; // corrected
	if( ph < 0 ) continue;

	phsum += ph;
	phcol[icol] += ph; // project cluster onto cols
	phrow[irow] += ph; // project cluster onto rows

	sumcol += icol*ph;
	sumrow += irow*ph;

      } // pix

      int ncol = colmax - colmin + 1;
      int nrow = rowmax - rowmin + 1;

      if( colmin > 0 && colmax < nx[iDUT]-2 &&
	  rowmin > 0 && rowmax < ny[iDUT]-1 ) {
	dutph0Histo.Fill( phsum * norm );
	dutnpxHisto.Fill( c->size );
	dutncolHisto.Fill( ncol );
	dutnrowHisto.Fill( nrow );
      }

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
      double sxA = triplets[iA].sx;
      double syA = triplets[iA].sy;
      //double txy = sqrt( sxA*sxA + syA*syA ); // angle

      double zA = DUTz - zmA; // z DUT from mid of triplet
      double xA = xmA + sxA * zA; // triplet x at zDUT
      double yA = ymA + syA * zA;

      trixAHisto.Fill( xA );
      triyAHisto.Fill( yA );
      trixyAHisto->Fill( xA, yA );

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

      double zc = ( Nz*zA - Ny*ymA - Nx*xmA ) / ( Nx*sxA + Ny*syA + Nz ); // from zmA
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

      // define track road through pixels: (shallow tilt, no turn)
      // better: intersect at DUTz +- 0.5*DUTthickness

      double yrd9 = y4 + 0.5*DUTthickness * fabs( tan( DUTtilt*wt ) );
      double yrd0 = y4 - 0.5*DUTthickness * fabs( tan( DUTtilt*wt ) );

      double prd[155][160]; // charge
      double drd[155][160]; // depth

      for( int irow = 0; irow < 160; ++irow ) // x
	for( int icol = 0; icol < 155; ++icol ) // x
	  prd[icol][irow] = -1; // flag: pixel is outside track road

      int npxrd = 0;
      int rdrowmin = 159;
      int rdrowmax =  0;
      int rdcolmin = 154;
      int rdcolmax =  0;

      for( int icol = 0; icol < 155; ++icol ) { // x

	double xpix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT];

	if( fabs( x4 - xpix ) > 1.5 * ptchx[iDUT] ) continue; // 3-col road

	for( int irow = 0; irow < 160; ++irow ) { // y

	  double ypix = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // [mm]

	  drd[icol][irow] = ( y4 - ypix ) / fabs( tan( DUTtilt*wt ) ); // depth, 0 = mid

	  if( ypix > yrd9 + 1.5 * ptchy[iDUT] ) continue; // with smearing
	  if( ypix < yrd0 - 1.5 * ptchy[iDUT] ) continue;
	  if( ypix > yrd9 + 0.5 * ptchy[iDUT] ) continue; // tighter
	  if( ypix < yrd0 - 0.5 * ptchy[iDUT] ) continue;

	  ++npxrd;

	  prd[icol][irow] = 0; // pixel is inside track road

	  if( irow < rdrowmin ) rdrowmin = irow;
	  if( irow > rdrowmax ) rdrowmax = irow;
	  if( icol < rdcolmin ) rdcolmin = icol;
	  if( icol > rdcolmax ) rdcolmax = icol;

	} // row

      } // col

      npxrdHisto.Fill( npxrd );
      nrwrdHisto.Fill( rdrowmax-rdrowmin+1 );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT pixel clusters:

      int pxrowmin = 159;
      int pxrowmax =  0;
      int pxcolmin = 154;
      int pxcolmax =  0;

      for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;

	double P0 = c->ph * norm; // cluster charge normalized to vertical incidence

	int colmin = 999;
	int colmax = -1;
	int rowmin = 999;
	int rowmax = -1;

	double phcol[nx[iDUT]];
	for( int icol = 0; icol < nx[iDUT]; ++icol ) phcol[icol] = 0;

	double phrow[ny[iDUT]];
	for( int irow = 0; irow < ny[iDUT]; ++irow ) phrow[irow] = 0;

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

	  double ph = px->ph; // [ADC]
	  if( ph < 0 ) continue;

	  if( prd[icol][irow] > -0.5 )
	    prd[icol][irow] += ph; // road

	  phcol[icol] += ph; // project cluster onto cols
	  phrow[irow] += ph; // project cluster onto rows

	  // match to road in x:

	  double xpix = ( icol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT];
	  double dx = xpix - x4;
	  cmsdxpxHisto.Fill( dx );

	  if( fabs( dx ) > 1.5 * ptchx[iDUT] ) continue;

	  if( irow < pxrowmin ) pxrowmin = irow;
	  if( irow > pxrowmax ) pxrowmax = irow;
	  if( icol < pxcolmin ) pxcolmin = icol;
	  if( icol > pxcolmax ) pxcolmax = icol;

	} // pix

	int ncol = colmax - colmin + 1;
	int nrow = rowmax - rowmin + 1;
	int nlng = nrow;

	if( DUTaligniteration > 1 && nlng < rowcut ) continue;

	//crow = 0.5 * ( rowmin + rowmax ); // mid from head and tail, overwrite!

	// DUT - triplet:

	double cmsx = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	double cmsy = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm

	cmsxvsx->Fill( x4, cmsx ); // correlation
	cmsyvsy->Fill( y4, cmsy );

	// residuals for pre-alignment:

	cmssxaHisto.Fill( cmsx + x3 ); // rot, tilt and turn but no shift
	cmsdxaHisto.Fill( cmsx - x3 ); // peak

	double cmsdx = cmsx - x4; // triplet extrapol

	cmsdxHisto.Fill( cmsdx );

	double cmsdy = cmsy - y4; // with proper upsigny

	if( fabs(cmsdy) < ycut ) {

	  cmsdxcHisto.Fill( cmsdx ); // align: same sign
	  cmsdxvsx.Fill( x4, cmsdx ); // align: same sign
	  cmsdxvsy.Fill( y4, cmsdx ); // align: opposite sign
	  cmsdxvstx.Fill( sxA, cmsdx );
	  cmsdxvsev1->Fill( iev, cmsdx ); // sync stability
	  cmsdxvsev2->Fill( iev, cmsdx ); // sync stability

	  if( fabs( cmsdx ) < 0.02 )
	    if( ldbt )
	      cout << "\t\t dx " << cmsdx << endl;

	} // ycut

	if( fabs(cmsdx) < xcut ) {

	  cmsycHisto.Fill( yc );
	  cmssyaHisto.Fill( cmsy + y3 );
	  cmsdyaHisto.Fill( cmsy - y3 );
	  cmsdycHisto.Fill( cmsdy );
	  cmsdyvsx.Fill( x4, cmsdy );
	  cmsdyvsy.Fill( y4, cmsdy );
	  cmsdyvsty.Fill( syA, cmsdy );

	} // dx

	// xy cuts:

	if( fabs(cmsdx) < xcut  && fabs(cmsdy) < ycut ) {

	  trixAlkHisto.Fill( xA ); // telescope coordinates
	  triyAlkHisto.Fill( yA );
	  trixyAlkHisto->Fill( xA, yA );
	  cmscolHisto.Fill( ccol ); // map
	  cmsrowHisto.Fill( crow );

	  cmsrowminHisto.Fill( rowmin );
	  cmsrowmaxHisto.Fill( rowmax );
	  cmscolminHisto.Fill( colmin );
	  cmscolmaxHisto.Fill( colmax );
	  yrd0Histo.Fill( yrd0 );
	  yrd9Histo.Fill( yrd9 );

	  double cmsymax = ( rowmax + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm

 	  if( yrd0 > -3.8 && yrd9 < 3.8 ) {
	    cmsdyexHisto.Fill( cmsymax-yrd9 ); // should be zero
	    cmsdyinHisto.Fill( cmsymax-yrd0 ); // should be track length
	  }

	  cmsnpxHisto.Fill( c->size );
	  cmsncolHisto.Fill( ncol );
	  cmsnrowHisto.Fill( nrow );

	  cmsnrowvsy.Fill( y4, nrow );
	  cmsnrowvsy0.Fill( yrd0, nrow );
	  cmsnrowvsy9.Fill( yrd9, nrow );

	  cmsnpxvsxmym->Fill( xmod, ymod, c->size ); // cluster size map

	  cmsph0Histo.Fill( P0 ); // Landau

	  for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	    cmspxphHisto.Fill( px->ph );

	    if( nrow > rowcut )
	      cmspxphlHisto.Fill( px->ph );

	  } // pix

	} // cut xy

      } // loop DUT clusters

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // look at road:

      if( liso &&
 	  yrd0 > -3.8 && // fiducial contained road
 	  yrd9 <  3.8 ) {

	// telescope has out-of-time pile-up: require some match

	int npxrdph = 0;
	for( int irow = rdrowmin; irow <= rdrowmax; ++irow )
	  for( int icol = rdcolmin; icol <= rdcolmax; ++icol )
	    if( prd[icol][irow] > 0.1 )
	      ++npxrdph;

	if( rdcolmin > 0 && rdcolmax < 154 ) // fiducial road in x
	  npxrdphHisto.Fill( npxrdph );

	if( rdcolmin > 0 &&
	    rdcolmax < 154 && // fiducial road in x
	    npxrdph > rowcut ) {

	  bool ldbrd = 0;

	  for( int irow = rdrowmin; irow <= rdrowmax; ++irow ) { // y

	    double yrow = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // [mm]
	    double dy = y4 - yrow;

	    double P = 0;
	    double D = drd[rdcolmin][irow]; // depth

	    for( int icol = rdcolmin; icol <= rdcolmax; ++icol ) { // x

	      if( prd[icol][irow] < -0.5 ) continue; // pixel is outside track road

	      P += prd[icol][irow];
	      D  = drd[icol][irow];

	    } // cols in this row

	    rowpHisto.Fill( P );

	    if(      irow <= rdrowmin ) // 1st
	      rowp1Histo.Fill( P );
	    else if( irow <= rdrowmin+1 ) // 2nd
	      rowp2Histo.Fill( P );
	    else if( irow >= rdrowmax ) // lst
	      rowp9Histo.Fill( P );
	    else if( irow >= rdrowmax-1 )
	      rowp8Histo.Fill( P );
	    else {
	      rowpmHisto.Fill( P ); // mid
	      if( P < 0.1 )
		ldbrd = 0; // switch on for empty in-road debug
	    }

	    if( P > 0.1 )
	      rowdyHisto.Fill( dy );

	    rowdHisto.Fill( D*1E3 );
	    if( P > 0.1 )
	      rowdphHisto.Fill( D*1E3 );

	    rowpvsdy.Fill( dy*1E3, P );

	    rowpvsd.Fill( D*1E3, P );
	    if( P > 0.1 )
	      rowp1vsd.Fill( D*1E3, P );

	    rowpdHisto->Fill( D*1E3, P ); // 2D

	    if(      D > 0.060 )
	      rowpp60Histo.Fill( P );
	    else if( D > 0.040 )
	      rowpp40Histo.Fill( P );
	    else if( D > 0.020 )
	      rowpp20Histo.Fill( P );
	    else if( D > 0.000 )
	      rowpp00Histo.Fill( P );
	    else if( D >-0.020 )
	      rowpm20Histo.Fill( P );
	    else if( D >-0.040 )
	      rowpm40Histo.Fill( P );
	    else if( D >-0.060 )
	      rowpm60Histo.Fill( P );
	    else
	      rowpm80Histo.Fill( P );

	  } // loop irow

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
		cout << "  " << fixed << setprecision(1) << prd[icol][irow];
	      cout << setprecision(6) << endl;
	    }

	  } // ldb

	  // rowpix without dph cut:

	  int nrow = 0;
	  int row0 = 160;
	  int row9 = -1;

	  for( int row = 0; row < 160; ++row ) { // y

	    if( rowpx[row].size() < 5 ) continue; // expect roicol or zero

	    double rowp = 0;

	    double ypix = ( row + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // [mm]
	    double dy = ypix - y4;

	    for( auto px = rowpx[row].begin(); px != rowpx[row].end(); ++px ) {

	      double xpix = ( px->col + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT];

	      double dx = xpix - x4;

	      double P = px->ph;

	      roipvsdxdy->Fill( dx*1E3, dy*1E3, P );

	      if( fabs( x4 - xpix ) > 1.5 * ptchx[iDUT] ) continue;

	      rowp += P;

	    } // px

	    roirowpvsdy.Fill( dy*1E3, rowp );

	    if( ypix > yrd0 && // inside road
		ypix < yrd9 ) {

	      roirowpHisto.Fill( rowp );

	      if( rowp > 1 ) {
		++nrow;
		row9 = row;
		if( row < row0 ) row0 = row;
		roirowszHisto.Fill( rowpx[row].size() ); // 5 or more
	      }

	    } // track road

	  } // rows

	  roinrowHisto.Fill( nrow ); // length of ROI
	  roirow0Histo.Fill( row0 );
	  roirow9Histo.Fill( row9 );

	} // in-time

      } // fid road

    } // loop triplets iA

    if( ldb ) cout << "done ev " << iev << endl << flush;

    ++iev;

  } while( reader->NextEvent() && iev < lev );

  delete reader;

  cout << "done after " << iev << " events" << endl;
  cout << "resyncs " << nresync << endl;

  histoFile.Write();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT alignment:

  if( DUTaligniteration == 0 && cmsdxaHisto.GetEntries() > 999 ) {

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
      DUTalignx = fgp0->GetParameter(1);
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
      DUTalignx = fgp0->GetParameter(1);
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
      DUTaligny = fgp0->GetParameter(1);
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
      DUTaligny = fgp0->GetParameter(1);
    }

  } // cmsdxa

  // finer alignment:

  if( DUTaligniteration > 0 && cmsdxHisto.GetEntries() > 999 ) {

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
    DUTalignx += fgp0->GetParameter(1);

  } // iteration

    // dxvsy -> rot:

  if( DUTaligniteration > 1 && cmsdxvsy.GetEntries() > 999 ) {
    cmsdxvsy.Fit( "pol1", "q", "", -midx[iDUT]+0.2, midx[iDUT]-0.2 );
    TF1 * fdxvsy = cmsdxvsy.GetFunction( "pol1" );
    cout << endl << cmsdxvsy.GetTitle();
    cout << ": extra rot " << upsignx*upsigny*fdxvsy->GetParameter(1) << endl;
    DUTrot += upsignx*upsigny*fdxvsy->GetParameter(1);
  }

  // dxvstx -> dz:

  if( DUTaligniteration > 2 && cmsdxvstx.GetEntries() > 999 ) {
    cmsdxvstx.Fit( "pol1", "q", "", -0.002, 0.002 );
    TF1 * fdxvstx = cmsdxvstx.GetFunction( "pol1" );
    cout << endl << cmsdxvstx.GetTitle();
    cout << ": z shift " << upsignx*fdxvstx->GetParameter(1) << " mm" << endl;
    DUTz += upsignx*fdxvstx->GetParameter(1);
  }

  if( DUTaligniteration > 0 && cmsdycHisto.GetEntries() > 999 ) {

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
    DUTaligny += fgp0->GetParameter(1);

  } // iter

    // dyvsy -> tilt:

  if( DUTaligniteration > 0 && fabs(DUTtilt0) > 5 && cmsdyvsy.GetEntries() > 999 ) {
    cmsdyvsy.Fit( "pol1", "q", "", -3, 3 );
    TF1 * fdyvsy = cmsdyvsy.GetFunction( "pol1" );
    cout << endl << cmsdyvsy.GetTitle()
	 << ": slope " << fdyvsy->GetParameter(1)
	 << ", extra tilt " << fdyvsy->GetParameter(1)/wt*sa
	 << " deg"
	 << endl;
    DUTtilt += fdyvsy->GetParameter(1)/wt*sa;
  }

  // write new DUT alignment:
  /*
  ofstream DUTalignFile( DUTalignFileName.str() );

  DUTalignFile << "# DUT alignment for run " << run << endl;
  ++DUTaligniteration;
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
  */
  cout << endl;
  cout << "alignx " << DUTalignx << endl
       << "aligny " << DUTaligny << endl
       << "rot    " << DUTrot << endl
       << "tilt   " << DUTtilt << endl
       << "turn   " << DUTturn << endl
       << "dz     " << DUTz - zz[2] << endl
    ;

  cout
    << endl << "upsignx " << upsignx
    << endl << "upsigny " << upsigny
    << endl << "   xcut " << xcut
    << endl << "   ycut " << ycut
    << endl << " rowcut " << rowcut
    << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  cout << endl << histoFile.GetName() << endl;

  cout << endl;

  return 0;
}

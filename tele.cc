
// Daniel Pitzl, DESY, Feb 2016
// telescope analysis with eudaq
// triplet pre-alignment

// make tele
// tele -g geo_2015_01b.dat 13117  # empty telescope Jan 2015, narrow triplet
// tele -g geo_2015_04x.dat 19037  # tilted DUT 504
// tele -g geo_2015_07b.dat 20833  #  66k tilted DUT 504
// tele -g geo_2015_07b.dat 20842  # 521k tilted DUT 504
// tele -g geo_2016_03a.dat 23133
// tele -l 99999 -g geo_2016_03b.dat 23407

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TFile.h>
#include <TH1I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TSystem.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <map>
#include <set>
#include <cmath>

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
};

struct cluster {
  vector <pixel> vpix;
  int size;
  int ncol, nrow;
  double col, row;
  double mindxy;
};

struct triplet {
  double xm;
  double ym;
  double zm;
  double sx;
  double sy;
  vector <double> vx;
  vector <double> vy;
};

// globals:

pixel pb[999]; // global declaration: array of pixel hits
int fNHit; // global

//------------------------------------------------------------------------------
vector<cluster> getClus()
{
  // returns clusters with local coordinates
  // decodePixels should have been called before to fill pixel buffer pb 
  // simple clusterization
  // cluster search radius fCluCut ( allows fCluCut-1 empty pixels)

  const int fCluCut = 1; // clustering: 1 = no gap (15.7.2012)
  //const int fCluCut = 2;

  vector<cluster> v;
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
    do{
      growing = 0;
      for( int i = 0; i < fNHit; ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
	    //if( (   dr>=-fCluCut) && (dr<=fCluCut) 
	    //&& (dc>=-fCluCut) && (dc<=fCluCut) ) { // allow diagonal
            if( ( abs(dr) <= fCluCut && dc == 0 ) ||
		( abs(dc) <= fCluCut && dr == 0 ) ) { // only facing neighbours, same resolution
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
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      c.col += (*p).col;
      c.row += (*p).row;
      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    c.col /= c.size;
    c.row /= c.size;

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
  FileReader * reader;
  if( run < 100 )
    reader = new FileReader( runnum, "data/run0000$2R$X");
  else if( run < 1000 )
    reader = new FileReader( runnum, "data/run000$3R$X");
  else if( run < 10000 )
    reader = new FileReader( runnum, "data/run00$4R$X");
  else if( run < 100000 )
    reader = new FileReader( runnum, "data/run0$5R$X");
  else
    reader = new FileReader( runnum, "data/run$6R$X");

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // further arguments:

  int lev = 999222111; // last event
  string geoFileName( "geo.dat" );
  double mom = 4.8;

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

    if( !strcmp( argv[i], "-g" ) )
      geoFileName = argv[++i];

    if( !strcmp( argv[i], "-p" ) )
      mom = atof( argv[++i] ); // momentum

  } // argc

  const double ang = sqrt( 0.005*0.005 + pow( 0.002*4.8/mom, 2 ) );
  double f = 4.8/mom;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // geometry:

  int nx[6]; // x-pixels per plane
  int ny[6]; // y-pixels per plane
  double sizex[6]; // x size per plane
  double sizey[6]; // y size per plane
  double ptchx[6]; // x-pixel size
  double ptchy[6]; // y-pixel size
  double midx[6]; // x mid
  double midy[6]; // y mid

  double zz[6];

  for( int ipl = 0; ipl < 6; ++ipl )
    nx[ipl] = 0; // missing plane flag

  ifstream geoFile( geoFileName );

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

      if( ipl < 0 || ipl >= 6 ) {
	//cout << "wrong plane number " << ipl << endl;
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

    for( int ipl = 0; ipl < 6; ++ipl ) {
      if( nx[ipl] == 0 ) continue; // missing plane flag
      ptchx[ipl] = sizex[ipl] / nx[ipl]; // pixel size
      ptchy[ipl] = sizey[ipl] / ny[ipl];
      midx[ipl] = 0.5 * sizex[ipl]; // mid plane
      midy[ipl] = 0.5 * sizey[ipl]; // mid plane
    }

  } // geo scope

  geoFile.close();

  // for profile plots:
  //double drng = 0.1*f; // narrow spacing
  double drng = 0.2*f; // wide spacing

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // Create directory based on run number

  stringstream run_number;
  run_number << run;

  // Directory will be made in the current working directory (cwd)
  std::string outputDirectory = std::string(gSystem->WorkingDirectory())+"/run_"+run_number.str();

  // Check first that it does not already exit
  if(!gSystem->MakeDirectory(TString(outputDirectory)))
    gSystem->MakeDirectory(TString(outputDirectory));

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // hot pixels:

  ostringstream hotFileName; // output string stream

  hotFileName << "hot_" << run << ".dat";

  ifstream ihotFile( outputDirectory+"/"+hotFileName.str() );

  set <int> hotset[6];

  cout << endl;
  if( ihotFile.bad() || ! ihotFile.is_open() ) {
    cout << "no " << hotFileName.str() << ", will be created" << endl;
  }
  // can there be instructions between if and else ? no!
  else {

    cout << "read hot pixel list from " << hotFileName.str() << endl;

    string hash( "#" );
    string plane( "plane" );
    string pix( "pix" );

    int ipl = 0;

    while( ! ihotFile.eof() ) {

      string line;
      getline( ihotFile, line );
      cout << line << endl;

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
  // alignments:

  int aligniteration = 0;
  double alignx[6];
  double aligny[6];
  double alignz[6];
  double rotx[6];
  double roty[6];

  for( int ipl = 0; ipl < 6; ++ipl ) {

    alignx[ipl] = 0.000; // [mm] same sign as dxAB
    aligny[ipl] = 0.000; // [mm] same sign as dy
    alignz[ipl] = 0.000; // [mm]
    rotx[ipl] = 0.0000; // [rad] rot, same     sign dxvsy
    roty[ipl] = 0.0000; // [rad] rot, opposite sign dyvsx

  }

  ostringstream alignFileName; // output string stream

  alignFileName << "align_" << run << ".dat";

  ifstream ialignFile( outputDirectory+"/"+alignFileName.str() );

  cout << endl;
  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    cout << "no " << alignFileName.str() << ", will bootstrap" << endl;
    cout << endl;
  }
  // can there be instructions between if and else ? no!
  else {

    cout << "read alignment from " << alignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string plane( "plane" );
    string shiftx( "shiftx" );
    string shifty( "shifty" );
    string shiftz( "shiftz" );
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

      if( ipl < 0 || ipl >= 6 ) {
	//cout << "wrong plane number " << ipl << endl;
	continue;
      }

      double val;
      tokenizer >> val;
      if(      tag == shiftx )
	alignx[ipl] = val;
      else if( tag == shifty )
	aligny[ipl] = val;
      else if( tag == shiftz ) {
	alignz[ipl] = val;
	zz[ipl] += val;
      }
      else if( tag == rotxvsy )
	rotx[ipl] = val;
      else if( tag == rotyvsx )
	roty[ipl] = val;

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  ialignFile.close();

  // z position for triplet-driplet matching = DUT material:

  double DUTz = 70 + zz[2];

  if( run >= 13000 && run < 14000 ) // Jan 2015 empty Datura
    DUTz =  ( 1*zz[3] + 2*zz[2] ) / 3; // driplet more spacing

  if( run >= 23000 ) // 2016
    //DUTz = 50 + zz[2]; // Apr 2016, 24152 sixdxc 20.9 um
    //DUTz = 60 + zz[2]; // Apr 2016, 24152 sixdxc 18.9 um
    //DUTz = 70 + zz[2]; // Apr 2016, 24152 sixdxc 17.8 um
    DUTz = 80 + zz[2]; // Apr 2016, 24152 sixdxc 17.3 um
    //DUTz = 90 + zz[2]; // Apr 2016, 24152 sixdxc 17.7 um

  if( run >= 24290 ) // 603, tilt 180
    //DUTz = 50 + zz[2]; // Apr 2016, 24295 sixdxcsi 16.3 um
    DUTz = 60 + zz[2]; // Apr 2016, 24295 sixdxcsi 15.7 um
    //DUTz = 70 + zz[2]; // Apr 2016, 24295 sixdxcsi 16.2 um

  if( DUTz > zz[3] )
    DUTz = 0.5 * ( zz[2] + zz[3] );

  // DUT Cu window in x: from sixdtvsx

  double xminCu = -6.5;
  double xmaxCu =  6.5;
  if( run >= 20180 ) { // Jun 2015 504
    xminCu = -6.4;
    xmaxCu =  6.4;
  }
  if( run >= 24000 ) { // Apr 2016 504
    xminCu = -9.1;
    xmaxCu =  3.8;
  }
  if( run >= 24200 ) { // 603, check!
    xminCu = -9.1;
    xmaxCu =  3.8;
  }
  if( run >= 24397 ) { // 603 adjusted
    xminCu = -9.8;
    xmaxCu =  3.2;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "tele" << run << ".root";

  TFile* histoFile = new TFile( (outputDirectory+"/"+rootFileName.str(  )).c_str(  ), "RECREATE" );

  // book histos:

  TH1I hdtus = TH1I( "dtus", "time between events;time between events [us];events", 100, 0, 1000 );
  TH1I hdtms = TH1I( "dtms", "time between events;time between events [ms];events", 100, 0, 1000 );

  TH1I t1Histo = TH1I( "t1", "event time;event time [s];events", 100, 0, 1 );
  TH1I t2Histo = TH1I( "t2", "event time;event time [s];events", 500, 0, 500 );
  TH1I t3Histo = TH1I( "t3", "event time;event time [s];events", 150, 0, 1500 );
  TH1I t4Histo = TH1I( "t4", "event time;event time [s];events", 600, 0, 6000 );
  TH1I t5Histo = TH1I( "t5", "event time;event time [s];events", 600, 0, 60000 );

  TH1I hnpx[6];
  TH1I hcol0[6];
  TH1I hrow0[6];
  TH1I hcol[6];
  TH1I hrow[6];

  TH1I hncl[6];
  TH1I hsiz[6];
  TH1I hncol[6];
  TH1I hnrow[6];
  TH1I hmindxy[6];

  TH1I hdx[6];
  TH1I hdy[6];
  TH1I hdxc[6];
  TH1I hdyc[6];

  TProfile dxvsy[6];
  TProfile dyvsx[6];

  TProfile teleffx[6];

  for( int ipl = 0; ipl < 6; ++ipl ) {

    hnpx[ipl] = TH1I( Form( "npx%i", ipl ),
		      Form( "%i pixel per event;pixels;%i events", ipl, ipl ),
		      200, 0, 200 );

    hcol0[ipl] = TH1I( Form( "allcol%i", ipl ),
		      Form( "%i all col;col;%i all pixels", ipl, ipl ), 
		      max( 52, nx[ipl]/4 ), 0, nx[ipl] );
    hrow0[ipl] = TH1I( Form( "allrow%i", ipl ),
		      Form( "%i all row;row;%i all pixels", ipl, ipl ),
		      max( 80, ny[ipl]/2 ), 0, ny[ipl] );

    hcol[ipl] = TH1I( Form( "col%i", ipl ),
		      Form( "%i col;col;%i pixels", ipl, ipl ), 
		      max( 52, nx[ipl]/4 ), 0, nx[ipl] );
    hrow[ipl] = TH1I( Form( "row%i", ipl ),
		      Form( "%i row;row;%i pixels", ipl, ipl ),
		      max( 80, ny[ipl]/2 ), 0, ny[ipl] );

    hncl[ipl] = TH1I( Form( "ncl%i", ipl ),
		      Form( "plane %i cluster per event;cluster;plane %i events", ipl, ipl ),
		      51, -0.5, 50.5 );
    hsiz[ipl] = TH1I( Form( "clsz%i", ipl ),
		      Form( "%i cluster size;pixels/cluster;%i clusters", ipl, ipl ),
		      51, -0.5, 50.5 );
    hncol[ipl] = TH1I( Form( "ncol%i", ipl ), 
		       Form( "%i cluster size x;columns/cluster;%i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );
    hnrow[ipl] = TH1I( Form( "nrow%i", ipl ),
		       Form( "%i cluster size y;rows/cluster;%i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );

    hmindxy[ipl] = TH1I( Form( "mindxy%i", ipl ),
			 Form( "%i cluster isolation;distance to next cluster [mm];%i clusters",
			       ipl, ipl ),
			 100, 0, 3 );

    hdx[ipl] = TH1I( Form( "dx%im", ipl ),
		     Form( "%i-m dx;%i-m dx [mm];cluster pairs", ipl, ipl ),
		     100, -f, f );
    hdy[ipl] = TH1I( Form( "dy%im", ipl ),
		     Form( "%i-m dy;%i-m dy [mm];cluster pairs", ipl, ipl ),
		     100, -f, f );
    hdxc[ipl] = TH1I( Form( "dxc%im", ipl ),
		     Form( "%i-m dx;%i-m dx [mm];cluster pairs", ipl, ipl ),
		     100, -f, f );
    hdyc[ipl] = TH1I( Form( "dyc%im", ipl ),
		     Form( "%i-m dy;%i-m dy [mm];cluster pairs", ipl, ipl ),
		     100, -f, f );

    dxvsy[ipl] = TProfile( Form( "dx%imvsy", ipl ),
			   Form( "%i-m dx vs y;y [mm];<%i-m dx> [mm]", ipl, ipl ),
			   100, -midy[ipl], midy[ipl], -drng, drng );
    dyvsx[ipl] = TProfile( Form( "dy%imvsx", ipl ),
			   Form( "%i-m dy vs x;x [mm];<%i-m dy> [mm]", ipl, ipl ),
			   100, -midx[ipl], midx[ipl], -drng, drng );

    teleffx[ipl] = TProfile( Form( "eff%ivsx", ipl ),
			     Form( "plane %i efficiency;x [mm];plane %i efficiency", ipl, ipl ),
			     100, -midx[ipl], midx[ipl], -1, 2 );

  } // planes

  TH1I hdxCA[2];
  TH1I hdyCA[2];

  TH1I htridx[2];
  TH1I htridy[2];
  TProfile trimadxvstx[2];
  TProfile trimadyvsty[2];

  TH1I htridxc[2];
  TH1I htridyc[2];
  TH1I htridxci[2];
  TH1I htridyci[2];
  TH1I htridxct[2];
  TH1I htridyct[2];

  TH1I htridxs1[2];
  TH1I htridxs2[2];
  TH1I htridxs3[2];
  TH1I htridxs4[2];
  TH1I htridxs5[2];

  TH1I htridxc1[2];
  TH1I htridxc2[2];
  TH1I htridxc3[2];
  TH1I htridxc4[2];
  TH1I htridxc5[2];

  TProfile tridxvsx[2];
  TProfile tridxvsy[2];
  TProfile tridyvsx[2];
  TProfile tridyvsy[2];

  TProfile tridxvstx[2];
  TProfile tridyvsty[2];

  TH1I htrix[2];
  TH1I htriy[2];
  TH2I * htrixy[2];
  TH1I htritx[2];
  TH1I htrity[2];

  for( int itd = 0; itd < 2; ++itd ) {

    string tri("tri");
    string dri("dri");
    string std(tri);
    if( itd ) 
      std = dri;

    hdxCA[itd] = TH1I( Form( "%sdxCA", std.c_str() ),
		       Form( "%splet dx CA;%splet dx CA [mm];C-A pairs",
			     std.c_str(), std.c_str() ),
			100, -1, 1 );
    hdyCA[itd] = TH1I( Form( "%sdyCA", std.c_str() ),
		       Form( "%splet dy CA;%splet dy CA [mm];C-A pairs",
			     std.c_str(), std.c_str() ),
			100, -1, 1 );

    htridx[itd] = TH1I( Form( "%sdx", std.c_str() ),
			Form( "%splet dx;%splet dx [mm];%splets",
			      std.c_str(), std.c_str(), std.c_str() ),
			200, -0.1, 0.1 );

    htridy[itd] = TH1I( Form( "%sdy", std.c_str() ),
			Form( "%splet dy;%splet dy [mm];%splets",
			      std.c_str(), std.c_str(), std.c_str() ),
			200, -0.1, 0.1 );

    htridxc[itd] = TH1I( Form( "%sdxc", std.c_str() ),
			Form( "%splet dx;%splet dx [mm];%splets",
			      std.c_str(), std.c_str(), std.c_str() ),
			200, -0.1, 0.1 );
    htridyc[itd] = TH1I( Form( "%sdyc", std.c_str() ),
			Form( "%splet dy;%splet dy [mm];%splets",
			      std.c_str(), std.c_str(), std.c_str() ),
			200, -0.1, 0.1 );

    htridxci[itd] = TH1I( Form( "%sdxci", std.c_str() ),
			  Form( "isolated %splet dx;%splet dx [mm];isolated %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );
    htridyci[itd] = TH1I( Form( "%sdyci", std.c_str() ),
			  Form( "isolated %splet dy;%splet dy [mm];isolated %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );

    trimadxvstx[itd] =
      TProfile( Form( "%smadxvstx", std.c_str() ),
		Form( "%splet MAD(#Deltax) vs #theta_{x};%splet #theta_{x} [rad];%splet MAD(#Deltax) [mm]",
		      std.c_str(), std.c_str(), std.c_str() ),
		100, -0.005, 0.005, 0, 0.05 );
    trimadyvsty[itd] =
      TProfile( Form( "%smadyvsty", std.c_str() ),
		Form( "%splet MAD(#Deltay) vs #theta_{y};%splet #theta_{y} [rad];%splet MAD(#Deltay) [mm]",
		      std.c_str(), std.c_str(), std.c_str() ),
		100, -0.005, 0.005, 0, 0.05 );

    htridxct[itd] = TH1I( Form( "%sdxct", std.c_str() ),
			  Form( "%splet dx;%splet dx [mm];%splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );
    htridyct[itd] = TH1I( Form( "%sdyct", std.c_str() ),
			  Form( "%splet dy;%splet dy [mm];%splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );

    htridxs1[itd] = TH1I( Form( "%sdxs1", std.c_str() ),
			  Form( "%splet dx 1-px;1-px %splet dx [mm];1-px %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );
    htridxs2[itd] = TH1I( Form( "%sdxs2", std.c_str() ),
			  Form( "%splet dx 2-px;2-px %splet dx [mm];2-px %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );
    htridxs3[itd] = TH1I( Form( "%sdxs3", std.c_str() ),
			  Form( "%splet dx 3-px;3-px %splet dx [mm];3-px %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );
    htridxs4[itd] = TH1I( Form( "%sdxs4", std.c_str() ),
			  Form( "%splet dx 4-px;4-px %splet dx [mm];4-px %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );
    htridxs5[itd] = TH1I( Form( "%sdxs5", std.c_str() ),
			  Form( "%splet dx 5-px;5-px %splet dx [mm];5-px %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );

    htridxc1[itd] = TH1I( Form( "%sdxc1", std.c_str() ),
			  Form( "%splet dx 1-col;1-col %splet dx [mm];1-col %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );
    htridxc2[itd] = TH1I( Form( "%sdxc2", std.c_str() ),
			  Form( "%splet dx 2-col;2-col %splet dx [mm];2-col %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );
    htridxc3[itd] = TH1I( Form( "%sdxc3", std.c_str() ),
			  Form( "%splet dx 3-col;3-col %splet dx [mm];3-col %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );
    htridxc4[itd] = TH1I( Form( "%sdxc4", std.c_str() ),
			  Form( "%splet dx 4-col;4-col %splet dx [mm];4-col %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );
    htridxc5[itd] = TH1I( Form( "%sdxc5", std.c_str() ),
			  Form( "%splet dx 5-col;5-col %splet dx [mm];5-col %splets",
				std.c_str(), std.c_str(), std.c_str() ),
			  200, -0.1, 0.1 );

    tridxvsy[itd] = TProfile( Form( "%sdxvsy", std.c_str() ),
			      Form( "%splet dx vs y;%splet yB [mm];<%splets #Deltax> [mm]",
				    std.c_str(), std.c_str(), std.c_str() ),
			      100, -midy[itd], midy[itd], -0.05, 0.05 );
    tridyvsx[itd] = TProfile( Form( "%sdyvsx", std.c_str() ),
			      Form( "%splet dy vs x;%splet xB [mm];<%splets #Deltay> [mm]",
				    std.c_str(), std.c_str(), std.c_str() ),
			      100, -midx[itd], midx[itd], -0.05, 0.05 );

    tridxvsx[itd] = TProfile( Form( "%sdxvsx", std.c_str() ),
			      Form( "%splet dx vs x;%splet xB [mm];<%splets #Deltax> [mm]",
				    std.c_str(), std.c_str(), std.c_str() ),
			      100, -midx[itd], midx[itd], -0.05, 0.05 );
    tridyvsy[itd] = TProfile( Form( "%sdyvsy", std.c_str() ),
			      Form( "%splet dy vs y;%splet yB [mm];<%splets #Deltay> [mm]",
				    std.c_str(), std.c_str(), std.c_str() ),
			      100, -midy[itd], midy[itd], -0.05, 0.05 );

    tridxvstx[itd] =
      TProfile( Form( "%sdxvstx", std.c_str() ),
		Form( "%splet dx vs tx;%splet slope x [rad];<%splets #Deltax> [mm]",
		      std.c_str(), std.c_str(), std.c_str() ),
		80, -0.002, 0.002, -0.05, 0.05 );
    tridyvsty[itd] =
      TProfile( Form( "%sdyvsty", std.c_str() ),
		Form( "%splet dy vs ty;%splet slope y [rad];<%splets #Deltay> [mm]",
		      std.c_str(), std.c_str(), std.c_str() ),
		80, -0.002, 0.002, -0.05, 0.05 );

    htrix[itd] = TH1I( Form( "%sx", std.c_str() ),
		       Form( "%splet x;%splet x_{mid} [mm];%splets",
			     std.c_str(), std.c_str(), std.c_str() ),
		       240, -12, 12 );
    htriy[itd] = TH1I( Form( "%sy", std.c_str() ),
		       Form( "%splet y;%splet y_{mid} [mm];%splets",
			     std.c_str(), std.c_str(), std.c_str() ),
		       120, -6, 6 );
    htrixy[itd] = new
      TH2I( Form( "%sxy", std.c_str() ),
	    Form( "%splet x-y;%splet x_{mid} [mm];%splet y_{mid} [mm];%splets",
		  std.c_str(), std.c_str(), std.c_str(), std.c_str() ),
	    240, -12, 12, 120, -6, 6 );

    htritx[itd] = TH1I( Form( "%stx", std.c_str() ),
			Form( "%splet #theta_{x};%splet #theta_{x} [rad];%splets",
			      std.c_str(), std.c_str(), std.c_str() ),
			100, -ang, ang );
    htrity[itd] = TH1I( Form( "%sty", std.c_str() ),
			Form( "%splet #theta_{y};%splet #theta_{y} [rad];%splets",
			      std.c_str(), std.c_str(), std.c_str() ),
			100, -ang, ang );

  } // triplets and driplets

  TH1I hntri = TH1I( "ntri", "triplets per event;triplets;events", 51, -0.5, 50.5 );
  TH1I hndri = TH1I( "ndri", "driplets per event;driplets;events", 51, -0.5, 50.5 );
  TProfile ntrivst2 =
    TProfile( "ntrivst2",
	      "triplets per event vs time;time [s];<triplets/event>",
	      500, 0, 500, -0.5, 99.5 );

  TH1I hexdx[6];
  TH1I hexdy[6];

  TH1I hexdxc[6];
  TH1I hexdyc[6];

  TProfile exdxvsy[6];
  TProfile exdyvsx[6];

  TProfile exdxvstx[6];
  TProfile exdyvsty[6];

  TProfile exmadxvstx[6];
  TProfile exmadyvsty[6];

  for( int ipl = 0; ipl < 6; ++ipl ) {

    hexdx[ipl] = TH1I( Form( "exdx%i", ipl ),
		       Form( "ex dx %i;dx tri - plane %i [mm];triplet - cluster pairs", ipl, ipl ),
		       200, -0.2*f, 0.2*f );
    hexdy[ipl] = TH1I( Form( "exdy%i", ipl ),
		       Form( "ex dy %i;dy tri - plane %i [mm];triplet - cluster pairs", ipl, ipl ),
		       200, -0.2*f, 0.2*f );
    hexdxc[ipl] = TH1I( Form( "exdxc%i", ipl ),
			Form( "ex dx %i;dx tri - plane %i [mm];triplet - cluster pairs", ipl, ipl ),
			200, -0.2*f, 0.2*f );
    hexdyc[ipl] = TH1I( Form( "exdyc%i", ipl ),
			Form( "ex dy %i;dy tri - plane %i [mm];triplet - cluster pairs", ipl, ipl ),
			200, -0.2*f, 0.2*f );

    exdxvsy[ipl] = TProfile( Form( "exdxvsy%i", ipl ),
			     Form( "ex dx vs y %i;y at %i [mm];<#Deltax> [mm]", ipl, ipl ),
			     100, -midy[ipl], midy[ipl], -0.2*f, 0.2*f );
    exdyvsx[ipl] = TProfile( Form( "exdyvsx%i", ipl ),
			     Form( "ex dy vs x %i;x at %i [mm];<#Deltay> [mm]", ipl, ipl ),
			     100, -midx[ipl], midx[ipl], -0.2*f, 0.2*f );

    exdxvstx[ipl] =
      TProfile( Form( "exdxvstx%i", ipl ),
		Form( "dx vs tx at %i;slope x [rad];<#Deltax> at %i [mm]", ipl, ipl ),
		80, -0.002, 0.002, -0.1*f, 0.1*f );
    exdyvsty[ipl] =
      TProfile( Form( "exdyvsty%i", ipl ),
		Form( "dy vs ty at %i;slope y [rad];<#Deltay> at %i [mm]", ipl, ipl ),
		80, -0.002, 0.002, -0.1*f, 0.1*f );

    exmadxvstx[ipl] =
      TProfile( Form( "exmadxvstx%i", ipl ),
		Form( "MAD x vs tx at %i;slope x [rad];MAD #Deltax at %i [mm]", ipl, ipl ),
		80, -0.002*f, 0.002*f, 0, 0.1*f );
    exmadyvsty[ipl] =
      TProfile( Form( "exmadyvsty%i", ipl ),
		Form( "MAD y vs ty at %i;slope y [rad];MAD #Deltay at %i [mm]", ipl, ipl ),
		80, -0.002*f, 0.002*f, 0, 0.1*f );

  }  // ipl

  // driplets - triplets

  TH1I hsixdx = TH1I( "sixdx", "six dx;#Deltax [mm];triplet-driplet pairs", 100, -1*f, 1*f );
  TH1I hsixdy = TH1I( "sixdy", "six dy;#Deltay [mm];triplet-driplet pairs", 100, -1*f, 1*f );
  TH1I hsixdxc = TH1I( "sixdxc", "six dx;#Deltax [mm];triplet-driplet pairs", 200, -0.4*f, 0.4*f );
  TH1I hsixdxcsi = TH1I( "sixdxcsi", "six dx Si;#Deltax [mm];triplet-driplet pairs in Si", 200, -0.4*f, 0.4*f );
  TH1I hsixdxccu = TH1I( "sixdxccu", "six dx Cu;#Deltax [mm];triplet-driplet pairs in Cu", 200, -0.4*f, 0.4*f );
  TH1I hsixdxcsid = TH1I( "sixdxcsid", "six dx Si;#Deltax [mm];triplet-driplet pairs in Si", 200, -0.4*f, 0.4*f );
  TH1I hsixdyc = TH1I( "sixdyc", "six dy;#Deltay [mm];triplet-driplet pairs", 200, -0.4*f, 0.4*f );
  TH1I hsixdycsi = TH1I( "sixdycsi", "six dy Si;#Deltay [mm];triplet-driplet pairs in Si", 200, -0.4*f, 0.4*f );
  TH1I hsixdyccu = TH1I( "sixdyccu", "six dy Cu;#Deltay [mm];triplet-driplet pairs in Cu", 200, -0.4*f, 0.4*f );

  TProfile sixdxvsx =
    TProfile( "sixdxvsx",
	      "six #Deltax vs x;xB [mm];<driplet - triplet #Deltax [mm]",
	      220, -11, 11, -0.1*f, 0.1*f );
  TProfile sixmadxvsx =
    TProfile( "sixmadxvsx",
	      "six MAD x vs x;xB [mm];driplet - triplet MAD #Deltax [mm]",
	      220, -11, 11, 0, 0.1*f );
  TProfile sixmadxvsy =
    TProfile( "sixmadxvsy",
	      "six MAD x vs y;yB [mm];driplet - triplet MAD #Deltax [mm]",
	      110, -5.5, 5.5, 0, 0.1*f );
  TProfile sixmadxvstx =
    TProfile( "sixmadxvstx",
	      "six MAD x vs x;triplet #theta_{x} [rad];driplet - triplet MAD #Deltax [mm]",
	      80, -0.002, 0.002, 0, 0.1*f );
  TProfile sixmadxvsdtx =
    TProfile( "sixmadxvsdtx",
	      "six MAD x vs x;driplet-triplet #Delta#theta_{x} [rad];driplet - triplet MAD #Deltax [mm]",
	      80, -0.002, 0.002, 0, 0.1*f );
  TProfile sixdxvsy =
    TProfile( "sixdxvsy",
	      "six #Deltax vs y;yB [mm];<driplet - triplet #Deltax [mm]",
	      110, -5.5, 5.5, -0.5, 0.5 ); // for align
  TProfile sixdxvstx =
    TProfile( "sixdxvstx",
	      "six #Deltax vs slope x;slope x [rad];<driplet - triplet #Deltax> [mm]",
	      80, -0.002, 0.002, -0.1*f, 0.1*f );

  TProfile sixdyvsx =
    TProfile( "sixdyvsx",
	      "six #Deltay vs x;xB [mm];<driplet - triplet #Deltay [mm]",
	      220, -11, 11, -0.5, 0.5 ); // for align
  TProfile sixmadyvsx =
    TProfile( "sixmadyvsx",
	      "six MAD y vs x;xB [mm];driplet - triplet MAD #Deltay [mm]",
	      220, -11, 11, 0, 0.1*f );

  TProfile sixdyvsy =
    TProfile( "sixdyvsy",
	      "six #Deltay vs y;yB [mm];<driplet - triplet #Deltay [mm]",
	      110, -5.5, 5.5, -0.1*f, 0.1*f );
  TProfile sixdyvsty =
    TProfile( "sixdyvsty",
	      "six #Deltay vs slope y;slope y [rad];<driplet - triplet #Deltay> [mm]",
	      80, -0.002, 0.002, -0.1*f, 0.1*f );
  TProfile sixmadyvsy =
    TProfile( "sixmadyvsy",
	      "six MAD y vs y;yB [mm];driplet - triplet MAD #Deltay [mm]",
	      110, -5.5, 5.5, 0, 0.1*f );
  TProfile sixmadyvsty =
    TProfile( "sixmadyvsty",
	      "six MAD y vs #theta_{y};triplet #theta_{y} [rad];driplet - triplet MAD #Deltay [mm]",
	      80, -0.002, 0.002, 0, 0.1*f );
  TProfile sixmadyvsdty =
    TProfile( "sixmadyvsdty",
	      "six MAD y vs #Delta#theta_{y};driplet-triplet #Delta#theta_{y} [rad];driplet - triplet MAD #Deltay [mm]",
	      80, -0.002, 0.002, 0, 0.1*f );

  TProfile2D * sixdxyvsxy = new
    TProfile2D( "sixdxyvsxy",
		"driplet - triplet #Delta_{xy} vs x-y;x_{mid} [mm];y_{mid} [mm];<sqrt(#Deltax^{2}+#Deltay^{2})> [rad]",
		110, -11, 11, 55, -5.5, 5.5, 0, 0.1*f );

  TH1I hsixdtx =
    TH1I( "sixdtx",
	  "driplet slope x - triplet slope x;driplet slope x - triplet slope x;driplet-triplet pairs",
	  100, -0.0025*f, 0.0025*f );
  TH1I hsixdty =
    TH1I( "sixdty",
	  "driplet slope y - triplet slope y;driplet slope y - triplet slope y;driplet-triplet pairs",
	  100, -0.0025*f, 0.0025*f );     

  TH1I hsixdtxsi =
    TH1I( "sixdtxsi",
	  "driplet triplet #Delta#theta_{x} Si;driplet - triplet #Delta#theta_{x} [rad];driplet-triplet pairs in Si",
	  100, -0.0025*f, 0.0025*f );
  TH1I hsixdtxcu =
    TH1I( "sixdtxcu",
	  "driplet triplet #Delta#theta_{x} Cu;driplet - triplet #Delta#theta_{x} [rad];driplet-triplet pairs in Cu",
	  100, -0.010, 0.010 );

  TH1I hsixdtysi =
    TH1I( "sixdtysi",
	  "driplet triplet #Delta#theta_{y} Si;driplet - triplet #Delta#theta_{y} [rad];driplet-triplet pairs in Si",
	  100, -0.0025*f, 0.0025*f );
  TH1I hsixdtycu =
    TH1I( "sixdtycu",
	  "driplet triplet #Delta#theta_{y} Cu;driplet - triplet #Delta#theta_{y} [rad];driplet-triplet pairs in Cu",
	  100, -0.010*f, 0.010*f );

  TProfile sixdtvsx =
    TProfile( "sixdtvsx",
	      "driplet - triplet slope_{xy} vs x;x_{mid} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
	      110, -11, 11, 0, 0.1*f );
  TProfile2D * sixdtvsxy = new
    TProfile2D( "sixdtvsxy",
		"driplet - triplet slope_{xy} vs x-y;x_{mid} [mm];y_{mid} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [rad]",
		110, -11, 11, 55, -5.5, 5.5, 0, 0.1*f );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  int event_nr = 0;
  uint64_t evTLU0 = 0;
  const double fTLU = 384E6; // 384 MHz TLU clock
  uint64_t prevTLU = 0;

  map < int, int > pxmap[6];

  do {
    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() )
      eudaq::PluginManager::Initialize(evt);

    bool ldbg = 0;

    if( event_nr < 1  )
      ldbg = 1;

    if( lev < 100 )
      ldbg = 1;

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns
    if( event_nr < 2  ) // BORE has older time
      evTLU0 = evTLU;
    double evsec = (evTLU - evTLU0) / fTLU;
    t1Histo.Fill( evsec );
    t2Histo.Fill( evsec );
    t3Histo.Fill( evsec );
    t4Histo.Fill( evsec );
    t5Histo.Fill( evsec );

    double evdt = (evTLU - prevTLU) / fTLU;
    hdtus.Fill( evdt * 1E6 ); // [us]
    hdtms.Fill( evdt * 1E3 ); // [ms]
    prevTLU = evTLU;

    if( event_nr < 10 )
      cout << "tele processing  " << run << "." << event_nr << "  taken " << evsec << endl;
    else if( event_nr < 100 && event_nr%10 == 0 )
      cout << "tele processing  " << run << "." << event_nr << "  taken " << evsec << endl;
    else if( event_nr < 1000 && event_nr%100 == 0 )
      cout << "tele processing  " << run << "." << event_nr << "  taken " << evsec << endl;
    else if( event_nr%1000 == 0 )
      cout << "tele processing  " << run << "." << event_nr << "  taken " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    vector <cluster> cl[6];

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldbg ) std::cout << "PLANE " << plane.ID() << ": ";

      // /home/pitzl/eudaq/main/include/eudaq/CMSPixelHelper.hh

      int ipl = plane.ID();
      if ( ipl < 7 && run > 28000 )
	ipl = ipl - 1;
      if( ipl >= 6 ) continue; // only telescope here

      int npx = 0;

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldbg ) 
	  std::cout << plane.GetX(ipix)
		    << " " << plane.GetY(ipix)
		    << " " << plane.GetPixel(ipix) << " ";

	int ix = plane.GetX(ipix); // col pixel index
	int iy = plane.GetY(ipix); // row pixel index
	
	hcol0[ipl].Fill( ix );
	hrow0[ipl].Fill( iy );

	int ipx = ix*ny[ipl] + iy;
	if( pxmap[ipl].count(ipx) )
	  ++pxmap[ipl][ipx];
	else
	  pxmap[ipl].insert( make_pair( ipx, 1 ) );

	// fill pixel block for clustering:

	if( hotset[ipl].count(ipx) ) continue; // skip hot

	hcol[ipl].Fill( ix );
	hrow[ipl].Fill( iy );

	pb[npx].col = ix;
	pb[npx].row = iy;
	++npx;

	if( npx == 999 ) {
	  cout << "pixel buffer overflow in plane " << ipl
	       << ", event " << event_nr
	       << endl;
	  break;
	}

      } // pix
      
      hnpx[ipl].Fill(npx);

      if( ldbg ) std::cout << std::endl;

      // clustering:

      fNHit = npx; // for cluster search

      cl[ipl] = getClus();

      if( ldbg ) cout << "A clusters " << cl[ipl].size() << endl;

      hncl[ipl].Fill( cl[ipl].size() );

      for( vector<cluster>::iterator cA = cl[ipl].begin(); cA != cl[ipl].end(); ++cA ) {

	hsiz[ipl].Fill( cA->size );
	hncol[ipl].Fill( cA->ncol );
	hnrow[ipl].Fill( cA->nrow );

	// cluster isolation:

	double mindxy = 99.9;

	for( vector<cluster>::iterator cD = cl[ipl].begin(); cD != cl[ipl].end(); ++cD ) {
	  if( cD == cA ) continue; // self
	  double dx = (cD->col - cA->col)*ptchx[ipl];
	  double dy = (cD->row - cA->row)*ptchy[ipl];
	  double dxy = sqrt( dx*dx + dy*dy );
	  if( dxy < mindxy ) mindxy = dxy;
	}
	cA->mindxy = mindxy;
	hmindxy[ipl].Fill( mindxy );

      } // cl A

    } // planes

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // cluster pair correlations:

    for( int itd = 0; itd < 2; ++itd ) { // triplets 0-1-2 and driplets 3-4-5

      int im = 1; // mid plane triplet
      int ibeg = 0;
      int iend = 2;
      if( itd == 1 ) {
	im = 4; // mid plane driplet
	ibeg = 3;
	iend = 5;
      }

      for( vector<cluster>::iterator cA = cl[im].begin(); cA != cl[im].end(); ++cA ) {

	double xA = cA->col*ptchx[im] - alignx[im];
	double yA = cA->row*ptchy[im] - aligny[im];
	double xmid = xA - midx[im];
	double ymid = yA - midy[im];
	xA = xmid - ymid*rotx[im];
	yA = ymid + xmid*roty[im];

	for( int ipl = ibeg; ipl <= iend; ++ipl ) {

	  if( ipl == im ) continue;

	  for( vector<cluster>::iterator cB = cl[ipl].begin(); cB != cl[ipl].end(); ++cB ) {

	    double xB = cB->col*ptchx[ipl] - alignx[ipl];
	    double yB = cB->row*ptchy[ipl] - aligny[ipl];
	    double xmid = xB - midx[ipl];
	    double ymid = yB - midy[ipl];
	    xB = xmid - ymid*rotx[ipl];
	    yB = ymid + xmid*roty[ipl];

	    double dx = xB - xA;
	    double dy = yB - yA;
	    double dz = zz[ipl] - zz[im]; // signed
	    hdx[ipl].Fill( dx );
	    hdy[ipl].Fill( dy );
	    if( fabs(dy) < 0.003 * fabs(dz) ) { // beam divergence
	      hdxc[ipl].Fill( dx );
	      dxvsy[ipl].Fill( yB, dx );
	    }
	    if( fabs(dy) < 0.003 * fabs(dz) ) { // beam divergence
	      hdyc[ipl].Fill( dy );
	      dyvsx[ipl].Fill( xB, dy );
	    }

	  } // clusters

	} // ipl

      } // cl mid

    } // upstream and downstream internal correlations

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // telescope plane efficiency:

    double isoCut = 0.30; // [mm]
    double triCut = 0.05; // [mm]
    double effCut = 0.25; // [mm]

    for( int ipl = 0; ipl < 6; ++ipl ) {

      int ib = 1;
      int im = 2; // mid plane triplet
      int ie = 3;
      if(      ipl == 1 ) {
	ib = 0;
	im = 2;
	ie = 3;
      }
      else if( ipl == 2 ) {
	ib = 0;
	im = 1;
	ie = 3;
      }
      else if( ipl == 3 ) {
	ib = 1;
	im = 2;
	ie = 4;
      }
      else if( ipl == 4 ) {
	ib = 2;
	im = 3;
	ie = 5;
      }
      else if( ipl == 5 ) {
	ib = 2;
	im = 3;
	ie = 4;
      }

      double zD = zz[ipl];

      for( vector<cluster>::iterator cA = cl[ib].begin(); cA != cl[ib].end(); ++cA ) {

	if( cA->mindxy < isoCut ) continue;

	double xA = cA->col*ptchx[ib] - alignx[ib];
	double yA = cA->row*ptchy[ib] - aligny[ib];
	double xmid = xA - midx[ib];
	double ymid = yA - midy[ib];
	xA = xmid - ymid*rotx[ib];
	yA = ymid + xmid*roty[ib];
	double zA = zz[ib];

	for( vector<cluster>::iterator cC = cl[ie].begin(); cC != cl[ie].end(); ++cC ) {

	  if( cC->mindxy < isoCut ) continue;

	  double xC = cC->col*ptchx[ie] - alignx[ie];
	  double yC = cC->row*ptchy[ie] - aligny[ie];
	  double xmid = xC - midx[ie];
	  double ymid = yC - midy[ie];
	  xC = xmid - ymid*rotx[ie];
	  yC = ymid + xmid*roty[ie];
	  double zC = zz[ie];

	  double dx2 = xC - xA;
	  double dy2 = yC - yA;
	  double dz02 = zC - zA;

	  if( fabs( dx2 ) > ang * dz02 ) continue; // angle cut
	  if( fabs( dy2 ) > ang * dz02 ) continue; // angle cut

	  double xavg2 = 0.5*(xA + xC);
	  double yavg2 = 0.5*(yA + yC);
	  double zavg2 = 0.5*(zA + zC);
	  double slpx = ( xC - xA ) / dz02; // slope x
	  double slpy = ( yC - yA ) / dz02; // slope y

	  for( vector<cluster>::iterator cB = cl[im].begin(); cB != cl[im].end(); ++cB ) {

	    double xB = cB->col*ptchx[im] - alignx[im]; // stretch and shift
	    double yB = cB->row*ptchy[im] - aligny[im];
	    double xmid = xB - midx[im];
	    double ymid = yB - midy[im];
	    xB = xmid - ymid*rotx[im];
	    yB = ymid + xmid*roty[im];
	    double zB = zz[im];

	    // interpolate track to B:

	    double dz = zB - zavg2;
	    double xk = xavg2 + slpx * dz; // triplet at B
	    double yk = yavg2 + slpy * dz;

	    double dx3 = xB - xk;
	    double dy3 = yB - yk;

	    if( fabs( dx3 ) > triCut ) continue;
	    if( fabs( dy3 ) > triCut ) continue;

	    // interpolate track to D:

	    double da = zD - zavg2;
	    double xi = xavg2 + slpx * da; // triplet at D
	    double yi = yavg2 + slpy * da;

	    if( fabs( xi ) > 10.4 ) continue; // fiducial
	    if( fabs( yi ) >  5.2 ) continue; // fiducial

	    // eff pl:

	    int nm = 0;

	    for( vector<cluster>::iterator cD = cl[ipl].begin(); cD != cl[ipl].end(); ++cD ) {

	      double xD = cD->col*ptchx[ipl] - alignx[ipl];
	      double yD = cD->row*ptchy[ipl] - aligny[ipl];
	      double xmid = xD - midx[ipl];
	      double ymid = yD - midy[ipl];
	      xD = xmid - ymid*rotx[ipl];
	      yD = ymid + xmid*roty[ipl];

	      double dx4 = xD - xi;
	      double dy4 = yD - yi;

	      if( fabs( dx4 ) > effCut ) continue;
	      if( fabs( dy4 ) > effCut ) continue;

	      ++nm;

	      if( nm > 0 ) break; // one link is enough

	    } // cl D

	    teleffx[ipl].Fill( xi, nm );

	  } // cl B

	} // cl C

      } // cl A

    } // eff planes

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets 1 vs 2-0:
    // driplets 4 vs 5-3:

    vector <triplet> triplets;
    vector <triplet> driplets;

    double tricut = 0.05; // [mm]

    for( int itd = 0; itd < 2; ++itd ) { // triplets 0-1-2 and driplets 3-4-5

      int im = 1; // mid plane triplet
      int ib = 0;
      int ie = 2;
      if( itd == 1 ) {
	im = 4; // mid plane driplet
	ib = 3;
	ie = 5;
      }

      for( vector<cluster>::iterator cA = cl[ib].begin(); cA != cl[ib].end(); ++cA ) {

	double xA = cA->col*ptchx[ib] - alignx[ib];
	double yA = cA->row*ptchy[ib] - aligny[ib];
	double xmid = xA - midx[ib];
	double ymid = yA - midy[ib];
	xA = xmid - ymid*rotx[ib];
	yA = ymid + xmid*roty[ib];
	double zA = zz[ib];

	for( vector<cluster>::iterator cC = cl[ie].begin(); cC != cl[ie].end(); ++cC ) {

	  double xC = cC->col*ptchx[ie] - alignx[ie];
	  double yC = cC->row*ptchy[ie] - aligny[ie];
	  double xmid = xC - midx[ie];
	  double ymid = yC - midy[ie];
	  xC = xmid - ymid*rotx[ie];
	  yC = ymid + xmid*roty[ie];
	  double zC = zz[ie];

	  double dx2 = xC - xA;
	  double dy2 = yC - yA;
	  double dz02 = zC - zA;
	  hdxCA[itd].Fill( dx2 );
	  hdyCA[itd].Fill( dy2 );

	  if( fabs( dx2 ) > ang * dz02 ) continue; // angle cut
	  if( fabs( dy2 ) > ang * dz02 ) continue; // angle cut

	  double xavg2 = 0.5*(xA + xC);
	  double yavg2 = 0.5*(yA + yC);
	  double zavg2 = 0.5*(zA + zC);
	  double slpx = ( xC - xA ) / dz02; // slope x
	  double slpy = ( yC - yA ) / dz02; // slope y

	  for( vector<cluster>::iterator cB = cl[im].begin(); cB != cl[im].end(); ++cB ) {

	    double xB = cB->col*ptchx[im] - alignx[im]; // stretch and shift
	    double yB = cB->row*ptchy[im] - aligny[im];
	    double xmid = xB - midx[im];
	    double ymid = yB - midy[im];
	    xB = xmid - ymid*rotx[im];
	    yB = ymid + xmid*roty[im];
	    double zB = zz[im];

	    // interpolate track to B:

	    double dz = zB - zavg2;
	    double xk = xavg2 + slpx * dz; // triplet at B
	    double yk = yavg2 + slpy * dz;

	    double dx3 = xB - xk;
	    double dy3 = yB - yk;

	    htridx[itd].Fill( dx3 );
	    htridy[itd].Fill( dy3 );

	    bool iso = 1;
	    if( cA->mindxy < isoCut ) iso = 0;
	    if( cC->mindxy < isoCut ) iso = 0;
	    if( cB->mindxy < isoCut ) iso = 0;

	    if( fabs( dy3 ) < 0.02 ) {

	      htridxc[itd].Fill( dx3 );
	      if( iso ) htridxci[itd].Fill( dx3 );
	      trimadxvstx[itd].Fill( slpx, fabs(dx3) );
	      if( fabs( slpx < 0.001 ) )
		htridxct[itd].Fill( dx3 );

	      if(      cB->size == 1 )
		htridxs1[itd].Fill( dx3 ); // 4.2 um
	      else if( cB->size == 2 )
		htridxs2[itd].Fill( dx3 ); // 4.0 um
	      else if( cB->size == 3 )
		htridxs3[itd].Fill( dx3 ); // 3.8 um
	      else if( cB->size == 4 )
		htridxs4[itd].Fill( dx3 ); // 4.3 um
	      else
		htridxs5[itd].Fill( dx3 ); // 3.6 um

	      if(      cB->ncol == 1 )
		htridxc1[itd].Fill( dx3 ); // 4.0 um
	      else if( cB->ncol == 2 )
		htridxc2[itd].Fill( dx3 ); // 4.1 um
	      else if( cB->ncol == 3 )
		htridxc3[itd].Fill( dx3 ); // 3.6 um
	      else if( cB->ncol == 4 )
		htridxc4[itd].Fill( dx3 ); // 3.5 um
	      else
		htridxc5[itd].Fill( dx3 ); // 4.1 um

	      tridxvsx[itd].Fill( xk, dx3 );
	      tridxvsy[itd].Fill( yk, dx3 );

	    } // dy

	    if( fabs( dx3 ) < 0.02 ) {

	      htridyc[itd].Fill( dy3 );
	      if( iso ) htridyci[itd].Fill( dy3 );
	      tridyvsx[itd].Fill( xk, dy3 );
	      tridyvsy[itd].Fill( yk, dy3 );
	      trimadyvsty[itd].Fill( slpy, fabs(dy3) );
	      if( fabs( slpy < 0.001 ) )
		htridyct[itd].Fill( dy3 );

	    }

	    // store triplets:

	    if( fabs( dx3 ) < tricut && fabs( dy3 ) < tricut ) {

	      triplet tri;
	      tri.xm = xavg2;
	      tri.ym = yavg2;
	      tri.zm = zavg2;
	      tri.sx = slpx;
	      tri.sy = slpy;

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

	      if( itd )
		driplets.push_back(tri);
	      else
		triplets.push_back(tri);

	      htrix[itd].Fill( xavg2 );
	      htriy[itd].Fill( yavg2 );
	      htrixy[itd]->Fill( xavg2, yavg2 );
	      htritx[itd].Fill( slpx );
	      htrity[itd].Fill( slpy );

	    } // triplet

	    // check z spacing: A-B as baseline

	    if( fabs( dx3 ) < tricut && fabs( dy3 ) < tricut ) {
	      double dzAB = zB - zA;
	      double tx = ( xB - xA ) / dzAB; // slope x
	      double ty = ( yB - yA ) / dzAB; // slope y
	      double dz = zC - zB;
	      double xk = xB + tx * dz; // at C
	      double yk = yB + ty * dz; // at C
	      double dx = xC - xk;
	      double dy = yC - yk;
	      tridxvstx[itd].Fill( tx, dx ); // adjust zpos, same sign
	      tridyvsty[itd].Fill( ty, dy );
	    }

	  } // cl B

	} // cl C

      } // cl A

    } // triplets and driplets

    hntri.Fill( triplets.size() );
    ntrivst2.Fill( evsec, triplets.size() );
    hndri.Fill( driplets.size() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // extrapolate triplets to each downstream plane
    // dy vs ty: dz

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // i = A = upstream

      double avxA = triplets[iA].xm;
      double avyA = triplets[iA].ym;
      double avzA = triplets[iA].zm;
      double slxA = triplets[iA].sx;
      double slyA = triplets[iA].sy;

      for( int ipl = 3; ipl <= 5; ++ipl ) {

	// triplet at plane:

	double zA = zz[ipl] - avzA; // z from mid of triplet to plane
	double xA = avxA + slxA * zA; // triplet at mid
	double yA = avyA + slyA * zA;

	for( vector<cluster>::iterator cC = cl[ipl].begin(); cC != cl[ipl].end(); ++cC ) {

	  double xC = cC->col*ptchx[ipl] - alignx[ipl];
	  double yC = cC->row*ptchy[ipl] - aligny[ipl];
	  double xmid = xC - midx[ipl];
	  double ymid = yC - midy[ipl];
	  xC = xmid - ymid*rotx[ipl];
	  yC = ymid + xmid*roty[ipl];

	  double dx = xC - xA;
	  double dy = yC - yA;
	  hexdx[ipl].Fill( dx );
	  hexdy[ipl].Fill( dy );
	  if( fabs( dy ) < 0.5 ) {
	    hexdxc[ipl].Fill( dx );
	    exdxvsy[ipl].Fill( yC, dx );
	    exdxvstx[ipl].Fill( slxA, dx );
	    exmadxvstx[ipl].Fill( slxA, fabs(dx) );
	  }
	  if( fabs( dx ) < 0.5 ) {
	    hexdyc[ipl].Fill( dy );
	    exdxvsy[ipl].Fill( xC, dy );
	    exdyvsty[ipl].Fill( slyA, dy );
	    exmadyvsty[ipl].Fill( slyA, fabs(dy) );
	  }

	} // clus

      } // planes

    } // triplets

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // extrapolate driplets to each upstream plane

    for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // j = B = downstream

      double avxB = driplets[jB].xm;
      double avyB = driplets[jB].ym;
      double avzB = driplets[jB].zm;
      double slxB = driplets[jB].sx;
      double slyB = driplets[jB].sy;

      for( int ipl = 0; ipl <= 2; ++ipl ) {

	// driplet at plane:

	double zB = zz[ipl] - avzB; // z from mid of driplet to plane
	double xB = avxB + slxB * zB; // driplet at mid
	double yB = avyB + slyB * zB;

	for( vector<cluster>::iterator cC = cl[ipl].begin(); cC != cl[ipl].end(); ++cC ) {

	  double xC = cC->col*ptchx[ipl] - alignx[ipl];
	  double yC = cC->row*ptchy[ipl] - aligny[ipl];
	  double xmid = xC - midx[ipl];
	  double ymid = yC - midy[ipl];
	  xC = xmid - ymid*rotx[ipl];
	  yC = ymid + xmid*roty[ipl];

	  double dx = xC - xB;
	  double dy = yC - yB;
	  hexdx[ipl].Fill( dx );
	  hexdy[ipl].Fill( dy );
	  if( fabs( dy ) < 0.5 ) {
	    hexdxc[ipl].Fill( dx );
	    exdxvsy[ipl].Fill( yC, dx );
	    exdxvstx[ipl].Fill( slxB, dx );
	    exmadxvstx[ipl].Fill( slxB, fabs(dx) );
	  }
	  if( fabs( dx ) < 0.5 ) {
	    hexdyc[ipl].Fill( dy );
	    exdxvsy[ipl].Fill( xC, dy );
	    exdyvsty[ipl].Fill( slyB, dy );
	    exmadyvsty[ipl].Fill( slyB, fabs(dy) );
	  }

	} // clus

      } // planes

    } // driplets

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // match triplets and driplets, measure offset

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // i = A = upstream

      double avxA = triplets[iA].xm;
      double avyA = triplets[iA].ym;
      double avzA = triplets[iA].zm;
      double slxA = triplets[iA].sx;
      double slyA = triplets[iA].sy;

      // triplet at DUT:

      double zA = DUTz - avzA; // z from mid of triplet to mid driplet
      double xA = avxA + slxA * zA; // triplet at mid
      double yA = avyA + slyA * zA;

      for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // j = B = downstream

	double avxB = driplets[jB].xm;
	double avyB = driplets[jB].ym;
	double avzB = driplets[jB].zm;
	double slxB = driplets[jB].sx;
	double slyB = driplets[jB].sy;

	// driplet at DUT:

	double zB = DUTz - avzB; // z from mid of triplet to mid
	double xB = avxB + slxB * zB; // triplet at mid
	double yB = avyB + slyB * zB;

	// driplet - triplet:

	double dx = xB - xA;
	double dy = yB - yA;
	double dxy = sqrt( dx*dx + dy*dy );
	double dtx = slxB - slxA;
	double dty = slyB - slyA;
	double dtxy = sqrt( dtx*dtx + dty*dty );

	hsixdx.Fill( dx ); // for align fit
	hsixdy.Fill( dy ); // for align fit

	if( fabs(dy) < 0.1 ) {

	  hsixdxc.Fill( dx );
	  if( xA > xminCu && xA < xmaxCu ) // no Cu
	    hsixdxcsi.Fill( dx );
	  else
	    hsixdxccu.Fill( dx );

	  sixdxvsx.Fill( xA, dx );
	  sixmadxvsx.Fill( xA, fabs(dx) );
	  if( xA > xminCu && xA < xmaxCu ) { // no Cu
	    sixdxvsy.Fill( yA, dx );
	    sixdxvstx.Fill( slxA, dx );
	    sixmadxvsy.Fill( yA, fabs(dx) );
	    sixmadxvstx.Fill( slxA, fabs(dx) );
	    sixmadxvsdtx.Fill( dtx, fabs(dx) ); // U-shape
	    if( fabs( dtx ) < 0.0005 )
	      hsixdxcsid.Fill( dx );
	  } // Si
	} // dy

	if( fabs(dx) < 0.1 ) {

	  hsixdyc.Fill( dy );
	  if( xA > xminCu && xA < xmaxCu ) // no Cu
	    hsixdycsi.Fill( dy );
	  else
	    hsixdyccu.Fill( dy );

	  sixdyvsx.Fill( xA, dy );
	  sixmadyvsx.Fill( xA, fabs(dy) );
	  if( xA > xminCu && xA < xmaxCu ) { // no Cu
	    sixdyvsy.Fill( yA, dy );
	    sixdyvsty.Fill( slyA, dy );
	    sixmadyvsy.Fill( yA, fabs(dy) );
	    sixmadyvsty.Fill( slyA, fabs(dy) );
	    sixmadyvsdty.Fill( dty, fabs(dy) ); // U-shape
	  }
	}

	// compare slopes:

	if( fabs(dy) < 0.1 && fabs(dx) < 0.1 ) {

	  sixdxyvsxy->Fill( xA, yA, dxy );

	  hsixdtx.Fill( dtx );
	  if( xA > xminCu && xA < xmaxCu ) { // no Cu
	    hsixdtxsi.Fill( dtx );
	    hsixdtysi.Fill( dty );
	  }
	  else {
	    hsixdtxcu.Fill( dtx );
	    hsixdtycu.Fill( dty );
	  }
	  hsixdty.Fill( dty ); // width: 0.3 mrad
	  sixdtvsx.Fill( xA, dtxy );
	  sixdtvsxy->Fill( xA, yA, dtxy );

	} // match

	// matched tri-dri => GBL

	if( fabs(dy) < 0.5 &&
	    fabs(dx) < 0.5 &&
	    fabs( dtx ) < 0.001 &&
	    fabs( dty ) < 0.001 ) {

	} // match

      } // driplets

    } // triplets

    ++event_nr;

  } while( reader->NextEvent() && event_nr < lev );

  //delete reader;

  cout << "done after " << event_nr << " events" << endl;
  histoFile->Write();
  histoFile->Close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // hot pixels:

  ofstream hotFile( outputDirectory+"/"+hotFileName.str() );

  hotFile << "# telescope hot pixel ist for run " << run << endl;

  for( int ipl = 0; ipl < 6; ++ipl ) {
    hotFile << endl;
    hotFile << "plane " << ipl << endl;
    int nmax = 0;
    int ntot = 0;
    int nhot = 0;
    for( map < int, int >::iterator jpx = pxmap[ipl].begin(); jpx != pxmap[ipl].end(); ++ jpx ) {
      int nhit = jpx->second;
      ntot += nhit;
      if( nhit > nmax ) nmax = nhit;
      if( nhit > event_nr/128 ) {
	++nhot;
	int ipx = jpx->first;
	int ix = ipx/ny[ipl];
	int iy = ipx%ny[ipl];
	hotFile << "pix " << setw(4) << ix << setw(5) << iy << endl;
      }
    } // jpx
    cout << ipl
	 << ": active " << pxmap[ipl].size()
	 << ", sum " << ntot
	 << ", max " << nmax
	 << ", hot " << nhot
	 << endl;
  } // ipl

  cout << "hot pixel list written to " << hotFileName.str() << endl;

  hotFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignment fits:

  for( int ipl = 0; ipl < 6; ++ipl ) {

    double nb = hdx[ipl].GetNbinsX();
    double ne = hdx[ipl].GetSumOfWeights();
    double nm = hdx[ipl].GetMaximum();

    if( nm < 99 ) continue;

    cout << endl << hdx[ipl].GetTitle() << endl;
    cout << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl;
    cout << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;
    cout << "  at " << hdx[ipl].GetBinCenter( hdx[ipl].GetMaximumBin() )
	 << endl;

    TF1 * fgp0x = new TF1( "fgp0x", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0x->SetParameter( 0, nm ); // amplitude
    fgp0x->SetParameter( 1, hdx[ipl].GetBinCenter( hdx[ipl].GetMaximumBin() ) );
    fgp0x->SetParameter( 2, 0.05 ); // sigma
    fgp0x->SetParameter( 3, hdx[ipl].GetBinContent(1) ); // BG
    hdx[ipl].Fit( "fgp0x", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0x->GetParameter(0)
	 << endl << "mid " << fgp0x->GetParameter(1)
	 << endl << "sig " << fgp0x->GetParameter(2)
	 << endl << " BG " << fgp0x->GetParameter(3)
	 << endl;

    // dy:

    nb = hdy[ipl].GetNbinsX();
    ne = hdy[ipl].GetSumOfWeights();
    nm = hdy[ipl].GetMaximum();
    cout << endl << hdy[ipl].GetTitle() << endl;
    cout << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl;
    cout << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;
    cout << "  at " << hdy[ipl].GetBinCenter( hdy[ipl].GetMaximumBin() )
	 << endl;

    TF1 * fgp0y = new TF1( "fgp0y", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0y->SetParameter( 0, nm ); // amplitude
    fgp0y->SetParameter( 1, hdy[ipl].GetBinCenter( hdy[ipl].GetMaximumBin() ) );
    fgp0y->SetParameter( 2, 0.05 ); // sigma
    fgp0y->SetParameter( 3, hdy[ipl].GetBinContent(1) ); // BG
    hdy[ipl].Fit( "fgp0y", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0y->GetParameter(0)
	 << endl << "mid " << fgp0y->GetParameter(1)
	 << endl << "sig " << fgp0y->GetParameter(2)
	 << endl << " BG " << fgp0y->GetParameter(3)
	 << endl;

    alignx[ipl] += fgp0x->GetParameter(1);
    aligny[ipl] += fgp0y->GetParameter(1);

    // x-y rotation from profiles:

    if( aligniteration ) {

      dxvsy[ipl].Fit( "pol1", "q", "", -midy[ipl], midy[ipl] );
      TF1 * fdxvsy = dxvsy[ipl].GetFunction( "pol1" );
      cout << endl << dxvsy[ipl].GetTitle() << " slope " << fdxvsy->GetParameter(1) << endl;

      dyvsx[ipl].Fit( "pol1", "q", "", -midx[ipl], midx[ipl] );
      TF1 * fdyvsx = dyvsx[ipl].GetFunction( "pol1" );
      cout << endl << dyvsx[ipl].GetTitle() << " slope " << fdyvsx->GetParameter(1) << endl;

      rotx[ipl] += fdxvsy->GetParameter(1);
      roty[ipl] -= fdyvsx->GetParameter(1); // sign

    }

  } // ipl

  // z-shift of last planes:

  if( aligniteration ) {

    for( int itd = 0; itd < 2; ++itd ) {
      cout << endl;
      int ipl = 2+3*itd;
      tridxvstx[itd].Fit( "pol1", "q", "", -0.002, 0.002 );
      TF1 * f1 = tridxvstx[itd].GetFunction( "pol1" );
      cout << tridxvstx[itd].GetTitle()
	   << " dz " << f1->GetParameter(1)
	   << " plane " << ipl
	   << " new zpos " << zz[2+3*itd] + f1->GetParameter(1)
	   << endl;
      alignz[ipl] += f1->GetParameter(1);
    }
  }

  // driplet vs triplet:

  if( aligniteration && hsixdx.GetMaximum() > 99 ) {

    double nb = hsixdx.GetNbinsX();
    double ne = hsixdx.GetSumOfWeights();
    double nm = hsixdx.GetMaximum();

    cout << endl << hsixdx.GetTitle() << endl;
    cout << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl;
    cout << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;
    cout << "  at " << hsixdx.GetBinCenter( hsixdx.GetMaximumBin() )
	 << endl;

    TF1 * fgp0x = new TF1( "fgp0x", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0x->SetParameter( 0, nm ); // amplitude
    fgp0x->SetParameter( 1, hsixdx.GetBinCenter( hsixdx.GetMaximumBin() ) );
    fgp0x->SetParameter( 2, 0.05 ); // sigma
    fgp0x->SetParameter( 3, hsixdx.GetBinContent(1) ); // BG
    hsixdx.Fit( "fgp0x", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0x->GetParameter(0)
	 << endl << "mid " << fgp0x->GetParameter(1)
	 << endl << "sig " << fgp0x->GetParameter(2)
	 << endl << " BG " << fgp0x->GetParameter(3)
	 << endl;

    // dy:

    nb = hsixdy.GetNbinsX();
    ne = hsixdy.GetSumOfWeights();
    nm = hsixdy.GetMaximum();
    cout << endl << hsixdy.GetTitle() << endl;
    cout << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl;
    cout << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;
    cout << "  at " << hsixdy.GetBinCenter( hsixdy.GetMaximumBin() )
	 << endl;

    TF1 * fgp0y = new TF1( "fgp0y", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0y->SetParameter( 0, nm ); // amplitude
    fgp0y->SetParameter( 1, hsixdy.GetBinCenter( hsixdy.GetMaximumBin() ) );
    fgp0y->SetParameter( 2, 0.05 ); // sigma
    fgp0y->SetParameter( 3, hsixdy.GetBinContent(1) ); // BG
    hsixdy.Fit( "fgp0y", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0y->GetParameter(0)
	 << endl << "mid " << fgp0y->GetParameter(1)
	 << endl << "sig " << fgp0y->GetParameter(2)
	 << endl << " BG " << fgp0y->GetParameter(3)
	 << endl;

    // update driplet planes:

    for( int ipl = 3; ipl < 6; ++ipl ) {
      alignx[ipl] += fgp0x->GetParameter(1);
      aligny[ipl] += fgp0y->GetParameter(1);
    }

    // x-y rotation from profiles:

    sixdxvsy.Fit( "pol1", "q", "", -midy[2], midy[2] );
    TF1 * fdxvsy = sixdxvsy.GetFunction( "pol1" );
    cout << endl << sixdxvsy.GetTitle() << " slope " << fdxvsy->GetParameter(1) << endl;

    sixdyvsx.Fit( "pol1", "q", "", -midx[2], midx[2] );
    TF1 * fdyvsx = sixdyvsx.GetFunction( "pol1" );
    cout << endl << sixdyvsx.GetTitle() << " slope " << fdyvsx->GetParameter(1) << endl;

    for( int ipl = 3; ipl < 6; ++ipl ) {
      rotx[ipl] += fdxvsy->GetParameter(1);
      roty[ipl] -= fdyvsx->GetParameter(1); // sign
    }

    // dz from dy vs ty:

  } // aligniteration > 0

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write alignment to file:

  ofstream alignFile( outputDirectory+"/"+alignFileName.str() );

  alignFile << "# telescope alignment for run " << run << endl;
  ++aligniteration;
  alignFile << "iteration " << aligniteration << endl;

  for( int ipl = 0; ipl < 6; ++ipl ) {
    alignFile << endl;
    alignFile << "plane " << ipl << endl;
    alignFile << "shiftx " << alignx[ipl] << endl;
    alignFile << "shifty " << aligny[ipl] << endl;
    alignFile << "shiftz " << alignz[ipl] << endl;
    alignFile << "rotxvsy " << rotx[ipl] << endl;
    alignFile << "rotyvsx " << roty[ipl] << endl;
  } // ipl

  alignFile.close();

  cout << endl
       << "wrote telescope alignment iteration " << aligniteration
       << " to " << alignFileName.str()
       << endl;
  if( aligniteration == 1 )
    cout << "need one more align iteration: please run again!" << endl;

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  return 0;
}

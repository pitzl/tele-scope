
// Daniel Pitzl, DESY, Feb 2016, Jan 2018, May 2019
// telescope analysis with eudaq and openMP
// updated from eudeq53
// produces align_run.dat

// make tele
// tele -g geo_2018_06r.dat -p 5.6 -l 99999 33095

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TFile.h>
#include <TH1I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <list>
#include <map>
#include <set>
#include <cmath>
#include <time.h> // clock_gettime

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int nn;
};

struct cluster {
  vector <pixel> vpix;
  //int size; int ncol, nrow;
  unsigned scr; // compressed
  float col, row;
  float mindxy;
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

bool ldbg = 0; // global

//------------------------------------------------------------------------------
vector<cluster> getClus( vector <pixel> pb, int fCluCut = 1 ) // 1 = no gap
{
  // returns clusters with pixel coordinates
  // next-neighbour topological clustering (allows fCluCut-1 empty pixels)

  vector <cluster> vc;
  if( pb.size() == 0 ) return vc;

  int * gone = new int[pb.size()];
  for( unsigned i = 0; i < pb.size(); ++i )
    gone[i] = 0;

  unsigned seed = 0;

  while( seed < pb.size() ) {

    // start a new cluster:

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do{
      growing = 0;
      for( unsigned i = 0; i < pb.size(); ++i ) {
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

    // count pixel neighbours:

    for( vector <pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {
      vector <pixel>::iterator q = p;
      ++q;
      for( ; q != c.vpix.end(); ++q )
	if( abs( p->col - q->col ) <= 1 &&abs( p->row - q->row ) <= 1 ) {
	  ++p->nn;
	  ++q->nn;
	}
    }

    c.col = 0;
    c.row = 0;
    int sumnn = 0;
    int minx = 9999;
    int maxx = 0;
    int miny = 9999;
    int maxy = 0;

    for( vector <pixel>::iterator p = c.vpix.begin(); p != c.vpix.end(); ++p ) {

      int nn = max( 1, p->nn ); // neighbours
      sumnn += nn;
      c.col += p->col * nn;
      c.row += p->row * nn;
      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;

    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    //c.col /= c.vpix.size();
    //c.row /= c.vpix.size();
    c.col /= sumnn; // weighted cluster center
    c.row /= sumnn;

    //c.size = c.vpix.size();
    //c.ncol = maxx-minx+1;
    //c.nrow = maxy-miny+1;
    c.scr = c.vpix.size() + 1024 * ( (maxx-minx+1) + 1024*(maxy-miny+1) ); // compressed size, ncol nrow
    c.mindxy = 999;

    c.vpix.clear(); // save space

    vc.push_back(c); // add cluster to vector

    // look for a new seed = unused pixel:

    while( ( ++seed < pb.size() ) && gone[seed] );

  } // while over seeds

  delete [] gone;

  return vc; // vector of clusters

} // getClus

//------------------------------------------------------------------------------
list < vector <cluster> > oneplane( unsigned ipl, list < vector <pixel> > pxlist )
{
  list < vector < cluster > > clist;

  for( auto ev = pxlist.begin(); ev != pxlist.end(); ++ev ) {

    vector <pixel> pb = *ev;

    vector <cluster> vcl = getClus(pb);

    if( ldbg ) cout << ipl << " clusters " << vcl.size() << endl;

    for( vector<cluster>::iterator cA = vcl.begin(); cA != vcl.end(); ++cA ) {

      // cluster isolation:

      vector<cluster>::iterator cD = cA;
      ++cD;
      for( ; cD != vcl.end(); ++cD ) {
	double dx = cD->col - cA->col;
	double dy = cD->row - cA->row;
	double dxy = sqrt( dx*dx + dy*dy );
	if( dxy < cA->mindxy ) cA->mindxy = dxy;
	if( dxy < cD->mindxy ) cD->mindxy = dxy;
      }

    } // cl A

    clist.push_back(vcl);

  }

  return clist;

} // oneplane

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc < 4 ) {
    cout << "format: tele -g geo_year_mon.dat run" << endl;
    return 1;
  }

  // run number = last arg

  string runnum( argv[argc-1] );
  int run = atoi( argv[argc-1] );

  cout << "run " << run << endl;
  FileReader * reader;
  if( run < 100 )
    reader = new FileReader( runnum.c_str(), "data/run0000$2R$X" );
  else if( run < 1000 )
    reader = new FileReader( runnum.c_str(), "data/run000$3R$X" );
  else if( run < 10000 )
    reader = new FileReader( runnum.c_str(), "data/run00$4R$X" );
  else if( run < 100000 )
    reader = new FileReader( runnum.c_str(), "data/run0$5R$X" );
  else
    reader = new FileReader( runnum.c_str(), "data/run$6R$X" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // further arguments:

  int lev = 999222111; // last event
  string geoFileName{ "geo.dat" };
  double mom = 4.8;

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

    if( !strcmp( argv[i], "-g" ) )
      geoFileName = argv[++i];

    if( !strcmp( argv[i], "-p" ) )
      mom = atof( argv[++i] ); // momentum

  } // argc

  double f = 4.8/mom;
  const double ang = sqrt( 0.005*0.005 + pow( 0.002*f, 2 ) ); // [rad]

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

  if( geoFile.bad() || ! geoFile.is_open() ) {
    cout << "Error opening " << geoFileName << endl;
    return 1;
  }

  cout << "read geometry from " << geoFileName << endl;

  { // open local scope

    string hash{ "#" };
    string plane{ "plane" };
    string type{ "type" };
    string sizexs{ "sizex" };
    string sizeys{ "sizey" };
    string npixelx{ "npixelx" };
    string npixely{ "npixely" };
    string zpos{ "zpos" };

    int ipl = -1;
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

      if( ipl < 0 || ipl > 5 ) { // Mimosa 1.6
	cout << "geo wrong plane number " << ipl << endl;
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
      ptchx[ipl] = sizex[ipl] / nx[ipl]; // pixel size 21.2/1152=18.4
      ptchy[ipl] = sizey[ipl] / ny[ipl];
      midx[ipl] = 0.5 * sizex[ipl]; // mid plane
      midy[ipl] = 0.5 * sizey[ipl]; // mid plane
    }

  } // geo scope

  geoFile.close();

  // for profile plots:
  //double drng = 0.1*f; // narrow spacing
  double drng = 0.2*f; // wide spacing [mm]

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // hot pixels:

  ostringstream hotFileName; // output string stream

  hotFileName << "hot_" << run << ".dat";

  ifstream ihotFile( hotFileName.str() );

  set <int> hotset[9];

  cout << endl;
  if( ihotFile.bad() || ! ihotFile.is_open() ) {
    cout << "no " << hotFileName.str() << ", will be created" << endl;
  }
  // can there be instructions between if and else ? no!
  else {

    cout << "read hot pixel list from " << hotFileName.str() << endl;

    string hash{ "#" };
    string plane{ "plane" };
    string pix{ "pix" };

    int ipl = -1;

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

      if( ipl < 0 || ipl > 5 ) {
	cout << "hot pixel wrong plane number " << ipl << endl;
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

  for( int ipl = 0; ipl < 9; ++ipl )
    cout << ipl << ": hot " << hotset[ipl].size() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignments:

  int aligniteration = 0;
  double alignx[9];
  double aligny[9];
  double alignz[9];
  double rotx[9];
  double roty[9];

  for( int ipl = 0; ipl < 9; ++ipl ) {

    alignx[ipl] = 0.000; // [mm] same sign as dxAB
    aligny[ipl] = 0.000; // [mm] same sign as dy
    alignz[ipl] = 0.000; // [mm]
    rotx[ipl] = 0.0000; // [rad] rot, same     sign dxvsy
    roty[ipl] = 0.0000; // [rad] rot, opposite sign dyvsx

  }

  ostringstream alignFileName; // output string stream

  alignFileName << "align_" << run << ".dat";

  ifstream ialignFile( alignFileName.str() );

  cout << endl;
  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    cout << "no " << alignFileName.str() << ", will bootstrap" << endl;
    cout << endl;
  }
  // can there be instructions between if and else ? no!
  else {

    cout << "read alignment from " << alignFileName.str() << endl;

    string hash{ "#" };
    string iteration{ "iteration" };
    string plane{ "plane" };
    string shiftx{ "shiftx" };
    string shifty{ "shifty" };
    string shiftz{ "shiftz" };
    string rotxvsy{ "rotxvsy" };
    string rotyvsx{ "rotyvsx" };

    int ipl = -1;

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

      if( tag == iteration ) {
	tokenizer >> aligniteration;
	continue;
      }

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 0 || ipl > 5 ) {
	cout << "align wrong plane number " << ipl << endl;
	continue;
      }

      double val;
      tokenizer >> val;
      if(      tag == shiftx )
	alignx[ipl] = val;
      else if( tag == shifty )
	aligny[ipl] = val;
      else if( tag == shiftz )
	alignz[ipl] = val;
      else if( tag == rotxvsy )
	rotx[ipl] = val;
      else if( tag == rotyvsx )
	roty[ipl] = val;

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  ialignFile.close();

  if( aligniteration == 0 ) f *= 3; // wider binning

  // z position for triplet-driplet matching = DUT material:

  double DUTz = 0.5 * ( zz[3] + zz[4] ); // 0.5*(304+560) = 432

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "tele" << run << ".root";

  TFile* histoFile = new TFile( rootFileName.str(  ).c_str(  ), "RECREATE" );

  // book histos:

  TH1I hdtus = TH1I( "dtus", "time between events;time between events [us];events", 100, 0, 1000 );
  TH1I hdtms = TH1I( "dtms", "time between events;time between events [ms];events", 100, 0, 1000 );

  TH1I t1Histo = TH1I( "t1", "event time;event time [s];events / 10 ms", 100, 0, 1 );
  TH1I t2Histo = TH1I( "t2", "event time;event time [s];events / s", 300, 0, 300 );
  TH1I t3Histo = TH1I( "t3", "event time;event time [s];events / 10 s", 150, 0, 1500 );
  TH1I t4Histo = TH1I( "t4", "event time;event time [s];events / 10 s", 600, 0, 6000 );
  TH1I t5Histo = TH1I( "t5", "event time;event time [s];events / min", 1000, 0, 60000 );

  TH1I hnpx0[9];
  TH1I hpivot[9];
  TH1I hcol0[9];
  TH1I hrow0[9];
  TH2I * hmap0[9];

  TH1I hcol[9];
  TH1I hrow[9];
  TH1I hbool[9];
  TH1I hnpx[9];

  TH1I hncl[9];
  TH1I hccol[9];
  TH1I hcrow[9];
  TH1I hnpix[9];
  TH1I hncol[9];
  TH1I hnrow[9];
  TH1I hmindxy[9];

  TH2I * hxx[9];
  TH1I hdx[9];
  TH1I hdy[9];
  TProfile dxvsx[9];
  TProfile dxvsy[9];
  TProfile dyvsx[9];
  TProfile dyvsy[9];

  TH1I hdx4[9];
  TH1I hdy4[9];
  TProfile dx4vsy[9];
  TProfile dy4vsx[9];

  TProfile effvsx[9];
  TProfile effvsxm[9];
  TProfile effvsym[9];
  TProfile2D * effvsxmym[9];

  for( int ipl = 0; ipl <= 5; ++ipl ) {

    hpivot[ipl] = TH1I( Form( "pivot%i", ipl ),
			Form( "plane %i pivot;pivot;plane %i events", ipl, ipl ),
			800, 0, 8000 );
    hnpx0[ipl] = TH1I( Form( "allnpx%i", ipl ),
		      Form( "plane %i all pixels per event;all pixels / event;plane %i events", ipl, ipl ),
		      200, 0, 200 );
    hcol0[ipl] = TH1I( Form( "allcol%i", ipl ),
		       Form( "plane %i all pix col;col;all plane %i pixels", ipl, ipl ), 
		       nx[ipl]/4, 0, nx[ipl] );
    hrow0[ipl] = TH1I( Form( "allrow%i", ipl ),
		       Form( "plane %i all pix row;row;all plane %i pixels", ipl, ipl ),
		       ny[ipl]/2, 0, ny[ipl] );
    hmap0[ipl] = new TH2I( Form( "map%i", ipl ),
			   Form( "plane %i map all;col;row;all plane %i pixels", ipl, ipl ),
			   nx[ipl]/4, 0, nx[ipl], ny[ipl]/2, 0, ny[ipl] );

    hnpx[ipl] = TH1I( Form( "npx%i", ipl ),
		      Form( "plane %i cool pixel per event;cool pixels / event;plane %i events", ipl, ipl ),
		      200, 0, 200 );
    hcol[ipl] = TH1I( Form( "col%i", ipl ),
		      Form( "plane %i cool pix col;col;cool plane %i pixels", ipl, ipl ), 
		      nx[ipl]/4, 0, nx[ipl] );
    hrow[ipl] = TH1I( Form( "row%i", ipl ),
		      Form( "plane %i cool pix row;row;cool plane %i pixels", ipl, ipl ),
		      ny[ipl]/2, 0, ny[ipl] );

    hncl[ipl] = TH1I( Form( "ncl%i", ipl ),
		      Form( "plane %i cluster per event;cluster;plane %i events", ipl, ipl ),
		      51, -0.5, 50.5 );
    hccol[ipl] = TH1I( Form( "ccol%i", ipl ),
		       Form( "plane %i cluster col;col;%i clusters", ipl, ipl ), 
		       nx[ipl]/4, 0, nx[ipl] );
    hcrow[ipl] = TH1I( Form( "crow%i", ipl ),
		       Form( "plane %i cluster row;row;%i clusters", ipl, ipl ),
		       ny[ipl]/2, 0, ny[ipl] );
    hnpix[ipl] = TH1I( Form( "npix%i", ipl ),
		       Form( "plane %i cluster size;pixels/cluster;%i clusters", ipl, ipl ),
		       80, 0.5, 80.5 );
    hncol[ipl] = TH1I( Form( "ncol%i", ipl ), 
		       Form( "plane %i cluster size x;columns/cluster;%i clusters", ipl, ipl ),
		       50, 0.5, 50.5 );
    hnrow[ipl] = TH1I( Form( "nrow%i", ipl ),
		       Form( "plane %i cluster size y;rows/cluster;%i clusters", ipl, ipl ),
		       50, 0.5, 50.5 );

    hmindxy[ipl] = TH1I( Form( "mindxy%i", ipl ),
			 Form( "plane %i cluster isolation;distance to next cluster [pixels];%i clusters",
			       ipl, ipl ),
			 100, 0, 10 );

    hxx[ipl] = new TH2I( Form( "xx%im", ipl ),
			 Form( "plane %i-m x-x;mid  plane x [mm];plane %i x [mm];cluster pairs", ipl, ipl ),
			 nx[ipl]/4, -midx[ipl], midx[ipl], nx[ipl]/4, -midx[ipl], midx[ipl] );

    hdx[ipl] = TH1I( Form( "dx%im", ipl ),
		     Form( "plane %i-m dx;%i-m dx [mm];cluster pairs", ipl, ipl ),
		     200, -2, 2 );
    hdy[ipl] = TH1I( Form( "dy%im", ipl ),
		     Form( "plane %i-m dy;%i-m dy [mm];cluster pairs", ipl, ipl ),
		     200, -2, 2 );

    dxvsx[ipl] = TProfile( Form( "dx%imvsx", ipl ),
			   Form( "plane %i-m dx vs x;x [mm];<%i-m dx> [mm]", ipl, ipl ),
			   100, -midx[ipl], midx[ipl], -drng, drng );
    dxvsy[ipl] = TProfile( Form( "dx%imvsy", ipl ),
			   Form( "plane %i-m dx vs y;y [mm];<%i-m dx> [mm]", ipl, ipl ),
			   100, -midy[ipl], midy[ipl], -drng, drng );
    dyvsx[ipl] = TProfile( Form( "dy%imvsx", ipl ),
			   Form( "plane %i-m dy vs x;x [mm];<%i-m dy> [mm]", ipl, ipl ),
			   100, -midx[ipl], midx[ipl], -drng, drng );
    dyvsy[ipl] = TProfile( Form( "dy%imvsy", ipl ),
			   Form( "plane %i-m dy vs y;y [mm];<%i-m dy> [mm]", ipl, ipl ),
			   100, -midy[ipl], midy[ipl], -drng, drng );

    hdx4[ipl] = TH1I( Form( "dx4%i", ipl ),
		      Form( "plane %i dx4;plane %i dx4 [#mum];track-cluster pairs", ipl, ipl ),
		      200, -100, 100 );
    hdy4[ipl] = TH1I( Form( "dy4%i", ipl ),
		      Form( "plane %i dy4;plane %i dy4 [#mum];track-cluster pairs", ipl, ipl ),
		      200, -100, 100 );
    dx4vsy[ipl] = TProfile( Form( "dx4%ivsy", ipl ),
			    Form( "plane %i dx4 vs y;y [mm];plane %i <dx4> [#mum]", ipl, ipl ),
			    100, -midy[ipl], midy[ipl], -100, 100 );
    dy4vsx[ipl] = TProfile( Form( "dy4%ivsx", ipl ),
			    Form( "plane %i dy4 vs x;x [mm];plane %i <dy4> [#mum]", ipl, ipl ),
			    100, -midx[ipl], midx[ipl], -100, 100 );

    effvsx[ipl] = TProfile( Form( "eff%ivsx", ipl ),
			     Form( "plane %i efficiency;x [mm];plane %i efficiency", ipl, ipl ),
			     100, -midx[ipl], midx[ipl], -1, 2 );
    effvsxm[ipl] = TProfile( Form( "eff%ivsxm", ipl ),
			     Form( "plane %i efficiency;x mod 36.8 [#mum];plane %i efficiency", ipl, ipl ),
			     74, 0, 2*ptchx[ipl]*1E3, -1, 2 );
    effvsym[ipl] = TProfile( Form( "eff%ivsym", ipl ),
			     Form( "plane %i efficiency;y mod 36.8 [#mum];plane %i efficiency", ipl, ipl ),
			     74, 0, 2*ptchy[ipl]*1E3, -1, 2 );
    effvsxmym[ipl] = new
      TProfile2D( Form( "eff%ivsxmym", ipl ),
		  Form( "plane %i efficiency;x mod 36.8 [#mum];y mod 36.8 [#mum];plane %i efficiency", ipl, ipl ),
		  74, 0, 2*ptchx[ipl]*1E3, 74, 0, 2*ptchy[ipl]*1E3, -1, 2 );

  } // planes

  TH1I hdxCA[2];
  TH1I hdyCA[2];
  TProfile dxCAvsx[2];
  TProfile dyCAvsy[2];

  TH1I htridx[2];
  TH1I htridy[2];

  TH1I htridxc[2];
  TH1I htridyc[2];
  TH1I htridxci[2];
  TH1I htridyci[2];
  TH1I htridxct[2];
  TH1I htridyct[2];
  TH1I htridyc1[2];
  TH1I htridyc2[2];
  TH1I htridyc3[2];
  TH1I htridyc4[2];
  TH1I htridyc5[2];
  TH1I htridyc6[2];
  TH1I htridycg[2];

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
  TH1I htridxc6[2];
  TH1I htridxcg[2];

  TProfile tridxvsx[2];
  TProfile tridxvsy[2];
  TProfile tridxvsxm[2];
  TH2I * tridxvsdy[2];
  TProfile2D * tridxvsxmym[2];
  TProfile2D * tridyvsxmym[2];
  TProfile2D * tridxycvsxmym[2];
  TProfile2D * tridxyvsxmym[2];
  TProfile2D * tridxy2vsxmym[2];
  TProfile trimadxvsxm[2];
  TProfile2D * trimadxvsxmym[2];
  TProfile tridxvstx[2];
  TProfile trimadxvstx[2];

  TProfile tridyvsx[2];
  TProfile tridyvsy[2];
  TProfile tridyvsym[2];
  TProfile trimadyvsym[2];
  TProfile2D * trimadyvsxmym[2];
  TProfile tridyvsty[2];
  TProfile trimadyvsty[2];

  TProfile tridxCvsx[2];
  TProfile tridxCvsy[2];
  TProfile tridyCvsx[2];
  TProfile tridyCvsy[2];
  TProfile tridxCvsax[2];
  TProfile tridyCvsay[2];

  TH1I htrix[2];
  TH1I htriy[2];
  TH2I * htrixy[2];
  TH1I htritx[2];
  TH1I htrity[2];

  TH1I htrincol[2];
  TH1I htrinrow[2];
  TH1I htrinpix[2];
  TH1I htrixmod[2];
  TH1I htrixmod1[2];
  TH1I htrixmod2[2];
  TH1I htrixmod3[2];
  TH1I htrixmod4[2];
  TH1I htrixmod5[2];
  TH1I htrixmod6[2];
  TProfile trincolvsxm[2];
  TProfile trinrowvsym[2];
  TProfile2D * trinpixvsxmym[2];
  TProfile2D * trinpixgvsxmym[2];

  for( int itd = 0; itd < 2; ++itd ) {

    string tri{"tri"};
    string dri{"dri"};
    string tds{tri};
    if( itd ) 
      tds = dri;
    int ipl = 1+3*itd; // 1 or 4

    hdxCA[itd] = TH1I( Form( "%sdxCA", tds.c_str() ),
		       Form( "%splet dx CA;%splet dx CA [mm];C-A pairs",
			     tds.c_str(), tds.c_str() ),
			100, -1, 1 );
    hdyCA[itd] = TH1I( Form( "%sdyCA", tds.c_str() ),
		       Form( "%splet dy CA;%splet dy CA [mm];C-A pairs",
			     tds.c_str(), tds.c_str() ),
			100, -1, 1 );

    dxCAvsx[itd] = TProfile( Form( "%sdxCAvsx", tds.c_str() ),
			     Form( "%splet dx CA vs x;%splet xC [mm];<%splets #Deltax CA> [mm]",
				   tds.c_str(), tds.c_str(), tds.c_str() ),
			     100, -midx[ipl], midx[ipl], -1, 1 );
    dyCAvsy[itd] = TProfile( Form( "%sdyCAvsy", tds.c_str() ),
			     Form( "%splet dy CA vs y;%splet yC [mm];<%splets #Deltay CA> [mm]",
				   tds.c_str(), tds.c_str(), tds.c_str() ),
			     100, -midy[ipl], midy[ipl], -1, 1 );

    htridx[itd] = TH1I( Form( "%sdx", tds.c_str() ),
			Form( "%splet dx;%splet dx [#mum];%splets",
			      tds.c_str(), tds.c_str(), tds.c_str() ),
			200, -100, 100 );

    htridy[itd] = TH1I( Form( "%sdy", tds.c_str() ),
			Form( "%splet dy;%splet dy [#mum];%splets",
			      tds.c_str(), tds.c_str(), tds.c_str() ),
			200, -100, 100 );
    htridyc[itd] = TH1I( Form( "%sdyc", tds.c_str() ),
			Form( "%splet dy;%splet dy [#mum];%splets",
			      tds.c_str(), tds.c_str(), tds.c_str() ),
			200, -100, 100 );
    htridyci[itd] = TH1I( Form( "%sdyci", tds.c_str() ),
			  Form( "isolated %splet dy;%splet dy [#mum];isolated %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );

    htridyct[itd] = TH1I( Form( "%sdyct", tds.c_str() ),
			  Form( "%splet dy;%splet dy [#mum];%splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );

    htridyc1[itd] = TH1I( Form( "%sdyc1", tds.c_str() ),
			  Form( "%splet dy 1-row;1-row %splet dy [#mum];1-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridyc2[itd] = TH1I( Form( "%sdyc2", tds.c_str() ),
			  Form( "%splet dy 2-row;2-row %splet dy [#mum];2-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridyc3[itd] = TH1I( Form( "%sdyc3", tds.c_str() ),
			  Form( "%splet dy 3-row;3-row %splet dy [#mum];3-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridyc4[itd] = TH1I( Form( "%sdyc4", tds.c_str() ),
			  Form( "%splet dy 4-row;4-row %splet dy [#mum];4-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridyc5[itd] = TH1I( Form( "%sdyc5", tds.c_str() ),
			  Form( "%splet dy 5-row;5-row %splet dy [#mum];5-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridyc6[itd] = TH1I( Form( "%sdyc6", tds.c_str() ),
			  Form( "%splet dy 6-row;6-row %splet dy [#mum];6-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridycg[itd] = TH1I( Form( "%sdycg", tds.c_str() ),
			  Form( "%splet dy good row;good-row %splet dy [#mum];good-row %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );

    htridxc[itd] = TH1I( Form( "%sdxc", tds.c_str() ),
			Form( "%splet dx;%splet dx [#mum];%splets",
			      tds.c_str(), tds.c_str(), tds.c_str() ),
			200, -100, 100 );
    htridxci[itd] = TH1I( Form( "%sdxci", tds.c_str() ),
			  Form( "isolated %splet dx;%splet dx [#mum];isolated %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridxct[itd] = TH1I( Form( "%sdxct", tds.c_str() ),
			  Form( "%splet dx;%splet dx [#mum];%splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );

    htridxs1[itd] = TH1I( Form( "%sdxs1", tds.c_str() ),
			  Form( "%splet dx 1-px;1-px %splet dx [#mum];1-px %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridxs2[itd] = TH1I( Form( "%sdxs2", tds.c_str() ),
			  Form( "%splet dx 2-px;2-px %splet dx [#mum];2-px %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridxs3[itd] = TH1I( Form( "%sdxs3", tds.c_str() ),
			  Form( "%splet dx 3-px;3-px %splet dx [#mum];3-px %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridxs4[itd] = TH1I( Form( "%sdxs4", tds.c_str() ),
			  Form( "%splet dx 4-px;4-px %splet dx [#mum];4-px %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridxs5[itd] = TH1I( Form( "%sdxs5", tds.c_str() ),
			  Form( "%splet dx 5-px;5-px %splet dx [#mum];5-px %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );

    htridxc1[itd] = TH1I( Form( "%sdxc1", tds.c_str() ),
			  Form( "%splet dx 1-col;1-col %splet dx [#mum];1-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridxc2[itd] = TH1I( Form( "%sdxc2", tds.c_str() ),
			  Form( "%splet dx 2-col;2-col %splet dx [#mum];2-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridxc3[itd] = TH1I( Form( "%sdxc3", tds.c_str() ),
			  Form( "%splet dx 3-col;3-col %splet dx [#mum];3-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridxc4[itd] = TH1I( Form( "%sdxc4", tds.c_str() ),
			  Form( "%splet dx 4-col;4-col %splet dx [#mum];4-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridxc5[itd] = TH1I( Form( "%sdxc5", tds.c_str() ),
			  Form( "%splet dx 5-col;5-col %splet dx [#mum];5-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridxc6[itd] = TH1I( Form( "%sdxc6", tds.c_str() ),
			  Form( "%splet dx 6-col;6-col %splet dx [#mum];6-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );
    htridxcg[itd] = TH1I( Form( "%sdxcg", tds.c_str() ),
			  Form( "%splet dx good col;good-col %splet dx [#mum];good-col %splets",
				tds.c_str(), tds.c_str(), tds.c_str() ),
			  200, -100, 100 );

    tridxvsx[itd] = TProfile( Form( "%sdxvsx", tds.c_str() ),
			      Form( "%splet dx vs x;%splet xB [mm];<%splets #Deltax> [#mum]",
				    tds.c_str(), tds.c_str(), tds.c_str() ),
			      100, -midx[ipl], midx[ipl], -20, 20 );
    tridxvsy[itd] = TProfile( Form( "%sdxvsy", tds.c_str() ),
			      Form( "%splet dx vs y;%splet yB [mm];<%splets #Deltax> [#mum]",
				    tds.c_str(), tds.c_str(), tds.c_str() ),
			      100, -midy[ipl], midy[ipl], -20, 20 );
    tridxvsxm[itd] = TProfile( Form( "%sdxvsxm", tds.c_str() ),
			       Form( "%splet dx vs xmod;%splet xB mod 36.8 [#mum];<%splets #Deltax> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       74, 0, 2*ptchx[ipl]*1E3, -20, 20 );
    tridxvsdy[itd] = new
      TH2I( Form( "%sdxvsdy", tds.c_str() ),
	    Form( "%splet #Deltay vs #Deltax;%splet #Deltax ;%splet #Deltay [#mum];%splets",
		  tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
	    80, -20, 20, 80, -20, 20 );
    tridxvsxmym[itd] = new
      TProfile2D( Form( "%sdxvsxmym", tds.c_str() ),
		  Form( "%splet #Deltax map;%splet xB mod 36.8 [#mum];%splet yB mod 36.8 [#mum];%splet <#Deltax> [#mum]",
			tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
		  74, 0, 2*ptchx[ipl]*1E3, 74, 0, 2*ptchy[ipl]*1E3, -20, 20 );
    tridyvsxmym[itd] = new
      TProfile2D( Form( "%sdyvsxmym", tds.c_str() ),
		  Form( "%splet #Deltay map;%splet xB mod 36.8 [#mum];%splet yB mod 36.8 [#mum];%splet <#Deltay> [#mum]",
			tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
		  74, 0, 2*ptchx[ipl]*1E3, 74, 0, 2*ptchy[ipl]*1E3, -20, 20 );
    tridxycvsxmym[itd] = new
      TProfile2D( Form( "%sdxycvsxmym", tds.c_str() ),
		  Form( "%splet #Deltax*y map;%splet xB mod 36.8 [#mum];%splet yB mod 36.8 [#mum];%splet <#Deltax#Deltay> [#mum^{2}]",
			tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
		  74, 0, 2*ptchx[ipl]*1E3, 74, 0, 2*ptchy[ipl]*1E3, -400, 400 );
    tridxyvsxmym[itd] = new
      TProfile2D( Form( "%sdxyvsxmym", tds.c_str() ),
		Form( "%splet #Deltax+y map;%splet xB mod 36.8 [#mum];%splet yB mod 36.8 [#mum];%splet <#Deltax+y> [#mum^{2}]",
		      tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
		74, 0, 2*ptchx[ipl]*1E3, 74, 0, 2*ptchy[ipl]*1E3, -20, 20 );
    tridxy2vsxmym[itd] = new
      TProfile2D( Form( "%sdxy2vsxmym", tds.c_str() ),
		Form( "%splet #Deltaxy map;%splet xB mod 36.8 [#mum];%splet yB mod 36.8 [#mum];%splet <#Deltaxy> [#mum^{2}]",
		      tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
		74, 0, 2*ptchx[ipl]*1E3, 74, 0, 2*ptchy[ipl]*1E3, 0, 2500 );
    trimadxvsxm[itd] =
      TProfile( Form( "%smadxvsxm", tds.c_str() ),
		Form( "%splet MAD(#Deltax) map;%splet xB mod 36.8 [#mum];%splet MAD(#Deltax) [#mum]",
		      tds.c_str(), tds.c_str(), tds.c_str() ),
		74, 0, 2*ptchx[ipl]*1E3, 0, 50 );

    trimadxvsxmym[itd] = new
      TProfile2D( Form( "%smadxvsxmym", tds.c_str() ),
		Form( "%splet MAD(#Deltax) vs xmod;%splet xB mod 36.8 [#mum];%splet yB mod 36.8 [#mum];%splet MAD(#Deltax) [#mum]",
		      tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
		74, 0, 2*ptchx[ipl]*1E3, 74, 0, 2*ptchy[ipl]*1E3, 0, 50 );

    tridxvstx[itd] = TProfile( Form( "%sdxvstx", tds.c_str() ),
			       Form( "%splet dx vs tx;%splet slope x [mrad];<%splets #Deltax> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       80, -2, 2, -20, 20 );
    trimadxvstx[itd] =
      TProfile( Form( "%smadxvstx", tds.c_str() ),
		Form( "%splet MAD(#Deltax) vs #theta_{x};%splet #theta_{x} [mrad];%splet MAD(#Deltax) [#mum]",
		      tds.c_str(), tds.c_str(), tds.c_str() ),
		100, -5, 5, 0, 50 );

    tridyvsx[itd] = TProfile( Form( "%sdyvsx", tds.c_str() ),
			      Form( "%splet dy vs x;%splet xB [mm];<%splets #Deltay> [#mum]",
				    tds.c_str(), tds.c_str(), tds.c_str() ),
			      100, -midx[ipl], midx[ipl], -20, 20 );
    tridyvsy[itd] = TProfile( Form( "%sdyvsy", tds.c_str() ),
			      Form( "%splet dy vs y;%splet yB [mm];<%splets #Deltay> [#mum]",
				    tds.c_str(), tds.c_str(), tds.c_str() ),
			      100, -midy[ipl], midy[ipl], -20, 20 );
    tridyvsym[itd] = TProfile( Form( "%sdyvsym", tds.c_str() ),
			       Form( "%splet dy vs ymod;%splet yB mod 36.8 [#mum];<%splets #Deltay> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       74, 0, 2*ptchy[ipl]*1E3, -20, 20 );
    trimadyvsym[itd] =
      TProfile( Form( "%smadyvsym", tds.c_str() ),
		Form( "%splet MAD(#Deltay) vs ymod;%splet yB mod 36.8 [#mum];%splet MAD(#Deltay) [#mum]",
		      tds.c_str(), tds.c_str(), tds.c_str() ),
		74, 0, 2*ptchy[ipl]*1E3, 0, 50 );

    trimadyvsxmym[itd] = new
      TProfile2D( Form( "%smadyvsxmym", tds.c_str() ),
		Form( "%splet MAD(#Deltay) vs xmod ymod;%splet xB mod 36.8 [#mum];%splet yB mod 36.8 [#mum];%splet MAD(#Deltay) [#mum]",
		      tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
		74, 0, 2*ptchx[ipl]*1E3, 74, 0, 2*ptchy[ipl]*1E3, 0, 50 );

    tridyvsty[itd] = TProfile( Form( "%sdyvsty", tds.c_str() ),
			       Form( "%splet dy vs ty;%splet slope y [mrad];<%splets #Deltay> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       80, -2, 2, -20, 20 );
    trimadyvsty[itd] =
      TProfile( Form( "%smadyvsty", tds.c_str() ),
		Form( "%splet MAD(#Deltay) vs #theta_{y};%splet #theta_{y} [mrad];%splet MAD(#Deltay) [#mum]",
		      tds.c_str(), tds.c_str(), tds.c_str() ),
		100, -5, 5, 0, 50 );

    tridxCvsx[itd] = TProfile( Form( "%sdxCvsx", tds.c_str() ),
			       Form( "%splet dxC vs x;%splet x at C [mm];<%splets #Deltax C> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       100, -midx[ipl], midx[ipl], -20, 20 );
    tridxCvsy[itd] = TProfile( Form( "%sdxCvsy", tds.c_str() ),
			       Form( "%splet dxC vs y;%splet y at C [mm];<%splets #Deltax C> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       100, -midy[ipl], midy[ipl], -20, 20 );
    tridyCvsx[itd] = TProfile( Form( "%sdyCvsx", tds.c_str() ),
			       Form( "%splet dyC vs x;%splet x at C [mm];<%splets #Deltay C> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       100, -midx[ipl], midx[ipl], -20, 20 );
    tridyCvsy[itd] = TProfile( Form( "%sdyCvsy", tds.c_str() ),
			       Form( "%splet dy vs y;%splet y at C [mm];<%splets #Deltay C> [#mum]",
				     tds.c_str(), tds.c_str(), tds.c_str() ),
			       100, -midy[ipl], midy[ipl], -20, 20 );
    tridxCvsax[itd] = TProfile( Form( "%sdxCvsax", tds.c_str() ),
				Form( "%splet dxC vs axAB;%splet AB slope x [mrad];<%splets #Deltax C> [#mum]",
				      tds.c_str(), tds.c_str(), tds.c_str() ),
				80, -2, 2, -20, 20 );
    tridyCvsay[itd] = TProfile( Form( "%sdyCvsay", tds.c_str() ),
				Form( "%splet dyC vs ayAB;%splet AB slope y [mrad];<%splets #Deltay C> [#mum]",
				      tds.c_str(), tds.c_str(), tds.c_str() ),
				80, -2, 2, -20, 20 );

    htrix[itd] = TH1I( Form( "%sx", tds.c_str() ),
		       Form( "%splet x;%splet x_{mid} [mm];%splets",
			     tds.c_str(), tds.c_str(), tds.c_str() ),
		       240, -12, 12 );
    htriy[itd] = TH1I( Form( "%sy", tds.c_str() ),
		       Form( "%splet y;%splet y_{mid} [mm];%splets",
			     tds.c_str(), tds.c_str(), tds.c_str() ),
		       120, -6, 6 );
    htrixy[itd] = new
      TH2I( Form( "%sxy", tds.c_str() ),
	    Form( "%splet x-y;%splet x_{mid} [mm];%splet y_{mid} [mm];%splets",
		  tds.c_str(), tds.c_str(), tds.c_str(), tds.c_str() ),
	    240, -12, 12, 120, -6, 6 );

    htritx[itd] = TH1I( Form( "%stx", tds.c_str() ),
			Form( "%splet #theta_{x};%splet #theta_{x} [mrad];%splets",
			      tds.c_str(), tds.c_str(), tds.c_str() ),
			200, -ang*1E3, ang*1E3 );
    htrity[itd] = TH1I( Form( "%sty", tds.c_str() ),
			Form( "%splet #theta_{y};%splet #theta_{y} [mrad];%splets",
			      tds.c_str(), tds.c_str(), tds.c_str() ),
			200, -ang*1E3, ang*1E3 );

    htrincol[itd] = TH1I( Form( "%sncol", tds.c_str() ),
			  Form( "%s cluster size x;columns/cluster;%s clusters on tracks",
				tds.c_str(), tds.c_str() ),
			  50, 0.5, 50.5 );
    htrinrow[itd] = TH1I( Form( "%snrow", tds.c_str() ),
			  Form( "%s cluster size y;rows/cluster;%s clusters on tracks",
				tds.c_str(), tds.c_str() ),
			  50, 0.5, 50.5 );
    htrinpix[itd] = TH1I( Form( "%snpix", tds.c_str() ),
			  Form( "%s cluster size;pixel/cluster;%s clusters on tracks",
				tds.c_str(), tds.c_str() ),
			  50, 0.5, 50.5 );

    htrixmod[itd] = TH1I( Form( "%sxmod", tds.c_str() ),
			   Form( "%s xmod;track x mod 38.8 #mum;%s clusters",
				 tds.c_str(), tds.c_str() ),
			   74, 0, 2*ptchx[ipl]*1E3 );
    htrixmod1[itd] = TH1I( Form( "%sxmod1", tds.c_str() ),
			   Form( "%s xmod for col 1;track x mod 38.8 #mum;%s 1-column clusters",
				 tds.c_str(), tds.c_str() ),
			   74, 0, 2*ptchx[ipl]*1E3 );
    htrixmod2[itd] = TH1I( Form( "%sxmod2", tds.c_str() ),
			   Form( "%s xmod for col 2;track x mod 38.8 #mum;%s 2-column clusters",
				 tds.c_str(), tds.c_str() ),
			   74, 0, 2*ptchx[ipl]*1E3 );
    htrixmod3[itd] = TH1I( Form( "%sxmod3", tds.c_str() ),
			   Form( "%s xmod for col 3;track x mod 38.8 #mum;%s 3-column clusters",
				 tds.c_str(), tds.c_str() ),
			   74, 0, 2*ptchx[ipl]*1E3 );
    htrixmod4[itd] = TH1I( Form( "%sxmod4", tds.c_str() ),
			   Form( "%s xmod for col 4;track x mod 38.8 #mum;%s 4-column clusters",
				 tds.c_str(), tds.c_str() ),
			   74, 0, 2*ptchx[ipl]*1E3 );
    htrixmod5[itd] = TH1I( Form( "%sxmod5", tds.c_str() ),
			   Form( "%s xmod for col 5;track x mod 38.8 #mum;%s 5-column clusters",
				 tds.c_str(), tds.c_str() ),
			   74, 0, 2*ptchx[ipl]*1E3 );
    htrixmod6[itd] = TH1I( Form( "%sxmod6", tds.c_str() ),
			   Form( "%s xmod for col 6;track x mod 38.8 #mum;%s 6-column clusters",
				 tds.c_str(), tds.c_str() ),
			   74, 0, 2*ptchx[ipl]*1E3 );

    trincolvsxm[itd] = TProfile( Form( "%sncolvsxm", tds.c_str() ),
				 Form( "%s cluster size;x mod 73.6 [#mum];%s <cluster size> [columns]",
				       tds.c_str(), tds.c_str() ),
				 74, 0, 4*ptchx[ipl]*1E3, 0, 2.5 );
    trinrowvsym[itd] = TProfile( Form( "%snrowvsym", tds.c_str() ),
				 Form( "%s cluster size;y mod 73.6 [#mum];%s <cluster size> [rows]",
				       tds.c_str(), tds.c_str() ),
				 74, 0, 4*ptchy[ipl]*1E3, 0, 2.5 );
    trinpixvsxmym[itd] = new
      TProfile2D( Form( "%snpixvsxmym", tds.c_str() ),
		  Form( "%s cluster size;x mod 73.6 [#mum];y mod 73.6 [#mum];%s <cluster size> [pixels]",
				tds.c_str(), tds.c_str() ),
		  74, 0, 4*ptchx[ipl]*1E3, 74, 0, 4*ptchy[ipl]*1E3, 0, 4.5 );
    trinpixgvsxmym[itd] = new
      TProfile2D( Form( "%snpixgvsxmym", tds.c_str() ),
		  Form( "%s cluster size;x mod 73.6 [#mum];y mod 73.6 [#mum];%s <cluster size> [pixels]",
				tds.c_str(), tds.c_str() ),
		  74, 0, 4*ptchx[ipl]*1E3, 74, 0, 4*ptchy[ipl]*1E3, 0, 4.5 );

  } // triplets and driplets

  TH1I hntri = TH1I( "ntri", "triplets per event;triplets;events", 51, -0.5, 50.5 );
  TH1I hndri = TH1I( "ndri", "driplets per event;driplets;events", 51, -0.5, 50.5 );
  TProfile ntrivsev =
    TProfile( "ntrivsev", "triplets per event vs time;time [events];<triplets/event>",
	      250, 0, 250*1000, -0.5, 99.5 );

  TH1I hexdx[9];
  TH1I hexdy[9];

  TH1I hexdxc[9];
  TH1I hexdyc[9];

  TProfile exdxvsy[9];
  TProfile exdyvsx[9];

  TProfile exdxvstx[9];
  TProfile exdyvsty[9];

  TProfile exmadxvstx[9];
  TProfile exmadyvsty[9];

  for( int ipl = 0; ipl <= 5; ++ipl ) {

    hexdx[ipl] = TH1I( Form( "exdx%i", ipl ),
		       Form( "ex dx %i;dx tri - plane %i [#mum];triplet - cluster pairs", ipl, ipl ),
		       200, -200*f, 200*f );
    hexdy[ipl] = TH1I( Form( "exdy%i", ipl ),
		       Form( "ex dy %i;dy tri - plane %i [#mum];triplet - cluster pairs", ipl, ipl ),
		       200, -200*f, 200*f );
    hexdxc[ipl] = TH1I( Form( "exdxc%i", ipl ),
			Form( "ex dx %i;dx tri - plane %i [#mum];triplet - cluster pairs", ipl, ipl ),
			200, -200*f, 200*f );
    hexdyc[ipl] = TH1I( Form( "exdyc%i", ipl ),
			Form( "ex dy %i;dy tri - plane %i [#mum];triplet - cluster pairs", ipl, ipl ),
			200, -200*f, 200*f );

    exdxvsy[ipl] = TProfile( Form( "exdxvsy%i", ipl ),
			     Form( "ex dx vs y %i;y at %i [mm];<#Deltax> [#mum]", ipl, ipl ),
			     100, -midy[ipl], midy[ipl], -200*f, 200*f );
    exdyvsx[ipl] = TProfile( Form( "exdyvsx%i", ipl ),
			     Form( "ex dy vs x %i;x at %i [mm];<#Deltay> [#mum]", ipl, ipl ),
			     100, -midx[ipl], midx[ipl], -200*f, 200*f );

    exdxvstx[ipl] =
      TProfile( Form( "exdxvstx%i", ipl ),
		Form( "dx vs tx at %i;slope x [mrad];<#Deltax> at %i [#mum]", ipl, ipl ),
		80, -2, 2, -100*f, 100*f );
    exdyvsty[ipl] =
      TProfile( Form( "exdyvsty%i", ipl ),
		Form( "dy vs ty at %i;slope y [mrad];<#Deltay> at %i [#mum]", ipl, ipl ),
		80, -2, 2, -100*f, 100*f );

    exmadxvstx[ipl] =
      TProfile( Form( "exmadxvstx%i", ipl ),
		Form( "MAD x vs tx at %i;slope x [mrad];MAD #Deltax at %i [#mum]", ipl, ipl ),
		80, -2*f, 2*f, 0, 100*f );
    exmadyvsty[ipl] =
      TProfile( Form( "exmadyvsty%i", ipl ),
		Form( "MAD y vs ty at %i;slope y [mrad];MAD #Deltay at %i [#mum]", ipl, ipl ),
		80, -2*f, 2*f, 0, 100*f );

  }  // ipl

  // driplets - triplets

  TH1I hsixdx = TH1I( "sixdx", "six dx;#Deltax [mm];triplet-driplet pairs", 800, -2*f, 2*f );
  TH1I hsixdy = TH1I( "sixdy", "six dy;#Deltay [mm];triplet-driplet pairs", 400, -1*f, 1*f );

  TH1I hsixdxc = TH1I( "sixdxc", "six dx;#Deltax [#mum];triplet-driplet pairs", 400, -200*f, 200*f );
  TH1I hsixdxcsid = TH1I( "sixdxcsid", "six dx Si;#Deltax [#mum];triplet-driplet pairs in Si", 400, -200*f, 200*f );
  TH1I hsixdyc = TH1I( "sixdyc", "six dy;#Deltay [#mum];triplet-driplet pairs", 400, -200*f, 200*f );

  TH2I * hsixxy = new
    TH2I( "sixxy", "sixplet x-y;sixplet x_{mid} [mm];sixplet y_{mid} [mm];sixplets",
	  240, -12, 12, 120, -6, 6 );

  TProfile sixdxvsx =
    TProfile( "sixdxvsx",
	      "six #Deltax vs x;xB [mm];<driplet - triplet #Deltax [#mum]",
	      220, -11, 11, -100*f, 100*f );
  TProfile sixmadxvsx =
    TProfile( "sixmadxvsx",
	      "six MAD x vs x;xB [mm];driplet - triplet MAD #Deltax [#mum]",
	      220, -11, 11, 0, 100*f );
  TProfile sixmadxvsy =
    TProfile( "sixmadxvsy",
	      "six MAD x vs y;yB [mm];driplet - triplet MAD #Deltax [#mum]",
	      110, -5.5, 5.5, 0, 100*f );
  TProfile sixmadxvstx =
    TProfile( "sixmadxvstx",
	      "six MAD x vs x;triplet #theta_{x} [mrad];driplet - triplet MAD #Deltax [#mum]",
	      80, -2, 2, 0, 100*f );
  TProfile sixmadxvsdtx =
    TProfile( "sixmadxvsdtx",
	      "six MAD x vs x;driplet-triplet #Delta#theta_{x} [mrad];driplet - triplet MAD #Deltax [#mum]",
	      80, -2, 2, 0, 100*f );
  TProfile sixdxvsy =
    TProfile( "sixdxvsy",
	      "six #Deltax vs y;yB [mm];<driplet - triplet #Deltax [#mum]",
	      110, -5.5, 5.5, -100*f, 100*f );
  TProfile sixdxvstx =
    TProfile( "sixdxvstx",
	      "six #Deltax vs slope x;slope x [mrad];<driplet - triplet #Deltax> [#mum]",
	      80, -2, 2, -100*f, 100*f );

  TProfile sixdyvsx =
    TProfile( "sixdyvsx",
	      "six #Deltay vs x;xB [mm];<driplet - triplet #Deltay [#mum]",
	      220, -11, 11, -100*f, 100*f );
  TProfile sixmadyvsx =
    TProfile( "sixmadyvsx",
	      "six MAD y vs x;xB [mm];driplet - triplet MAD #Deltay [#mum]",
	      220, -11, 11, 0, 100*f );

  TProfile sixdyvsy =
    TProfile( "sixdyvsy",
	      "six #Deltay vs y;yB [mm];<driplet - triplet #Deltay [#mum]",
	      110, -5.5, 5.5, -100*f, 100*f );
  TProfile sixdyvsty =
    TProfile( "sixdyvsty",
	      "six #Deltay vs slope y;slope y [mrad];<driplet - triplet #Deltay> [#mum]",
	      80, -2, 2, -100*f, 100*f );
  TProfile sixmadyvsy =
    TProfile( "sixmadyvsy",
	      "six MAD y vs y;yB [mm];driplet - triplet MAD #Deltay [#mum]",
	      110, -5.5, 5.5, 0, 100*f );
  TProfile sixmadyvsty =
    TProfile( "sixmadyvsty",
	      "six MAD y vs #theta_{y};triplet #theta_{y} [mrad];driplet - triplet MAD #Deltay [#mum]",
	      80, -2, 2, 0, 100*f );
  TProfile sixmadyvsdty =
    TProfile( "sixmadyvsdty",
	      "six MAD y vs #Delta#theta_{y};driplet-triplet #Delta#theta_{y} [mrad];driplet - triplet MAD #Deltay [#mum]",
	      80, -2, 2, 0, 100*f );

  TProfile2D * sixdxyvsxy = new
    TProfile2D( "sixdxyvsxy",
		"driplet - triplet #Delta_{xy} vs x-y;x_{mid} [mm];y_{mid} [mm];<sqrt(#Deltax^{2}+#Deltay^{2})> [mrad]",
		110, -11, 11, 55, -5.5, 5.5, 0, 100*f );

  TH1I hsixdtx( "sixdtx",
		"driplet slope x - triplet slope x;driplet slope x - triplet slope x [mrad];driplet-triplet pairs",
		200, -5*f, 5*f );
  TH1I hsixdty( "sixdty",
		"driplet slope y - triplet slope y;driplet slope y - triplet slope y [mrad];driplet-triplet pairs",
		200, -5*f, 5*f );

  TProfile sixdtvsx =
    TProfile( "sixdtvsx",
	      "driplet - triplet slope_{xy} vs x;x_{mid} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [mrad]",
	      220, -11, 11, 0, 100*f );
  TProfile2D * sixdtvsxy = new
    TProfile2D( "sixdtvsxy",
		"driplet - triplet slope_{xy} vs x-y;x_{mid} [mm];y_{mid} [mm];<sqrt(#Delta#theta_{x}^{2}+#Delta#theta_{y}^{2})> [mrad]",
		110, -11, 11, 55, -5.5, 5.5, 0, 100*f );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  timespec ts;
  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s0 = ts.tv_sec; // seconds since 1.1.1970
  long f0 = ts.tv_nsec; // nanoseconds
  double zeit1 = 0; // read
  double zeit2 = 0; // clus
  double zeit3 = 0; // track

  int iev = 0;
  uint64_t evTLU0 = 0;
  const double fTLU = 384E6; // 384 MHz TLU clock
  uint64_t prevTLU = 0;

  map < int, int > pxmap[9]; // for hot pixels

  list < vector < pixel > > pxlist[9];

  do {

    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() ) {
      cout << "Begin Of Run Event" << endl << flush;
      eudaq::PluginManager::Initialize(evt);
    }

    if( iev < 0  )
      ldbg = 1;

    if( lev < 100 )
      ldbg = 1;

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns
    if( iev < 2  ) // BORE has older time
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

    if( iev < 10 || ldbg )
      cout << "tele reading  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 100 && iev%10 == 0 )
      cout << "tele reading  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 1000 && iev%100 == 0 )
      cout << "tele reading  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev < 10000 && iev%1000 == 0 )
      cout << "tele reading  " << run << "." << iev << "  taken " << evsec << endl;
    else if( iev%10000 == 0 )
      cout << "tele reading  " << run << "." << iev << "  taken " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard( evt );

    string MIM{"MIMOSA26"};
    int mpl = 0; // Mimosa planes start at 0 in eudaq1.6

    if( ldbg ) cout << "planes " << sevt.NumPlanes() << endl << flush;

    for( size_t ipl = 0; ipl < sevt.NumPlanes(); ++ipl ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(ipl);

      if( ldbg )
	cout
	  << "  " << mpl
	  << ": plane " << plane.ID()
	  << " " << plane.Type() // NI
	  << " " << plane.Sensor() // MIMOSA26
	  << " frames " << plane.NumFrames() // 2
	  << " pivot " << plane.PivotPixel() // 6830
	  << " total " << plane.TotalPixels() // 663552
	  << " hits " << plane.HitPixels() // 5
	  << flush;

      //0: plane 1 NI MIMOSA26 frames 2 pivot 6830 total 663552 hits 5: 486 296 1: 635 68 1: 509 307 1: 510 307 1: 509 308 1

      if( plane.Sensor() != MIM ) {
	if( ldbg ) cout << endl << flush;
	continue;
      }

      hpivot[mpl].Fill( plane.PivotPixel() );

      vector <double> pxl = plane.GetPixels<double>();

      hnpx0[mpl].Fill( pxl.size() );

      vector <pixel> pb; // for clustering

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldbg ) 
	  cout << " :"
	       << " " << plane.GetX(ipix) // col
	       << " " << plane.GetY(ipix) // row
	       << " " << plane.GetPixel(ipix) // charge
	       << flush;

	int ix = plane.GetX(ipix); // col pixel index
	int iy = plane.GetY(ipix); // row pixel index

	hcol0[mpl].Fill( ix );
	hrow0[mpl].Fill( iy );
	hmap0[mpl]->Fill( ix, iy );

	int ipx = ix * ny[mpl] + iy;

	if( ldbg )
	  cout << " " << ipx << flush;

	if( pxmap[mpl].count(ipx) )
	  ++pxmap[mpl][ipx];
	else
	  pxmap[mpl].insert( make_pair( ipx, 1 ) ); // slow

	if( hotset[mpl].count(ipx) ) {
	  if( ldbg )
	    cout << " hot" << flush;
	  continue; // skip hot
	}

	// fill pixel block for clustering:

	hcol[mpl].Fill( ix );
	hrow[mpl].Fill( iy );

	pixel px;
	px.col = ix;
	px.row = iy;
	px.nn = 1; // init neighbours for best resolution
	pb.push_back(px);

	if( ldbg )
	  cout << " " << pb.size() << flush;

	if( pb.size() == 999 ) {
	  cout << "pixel buffer overflow in plane " << mpl
	       << ", event " << iev
	       << endl;
	  break;
	}

      } // pix

      if( ldbg )
	cout << " done" << endl << flush;

      hnpx[mpl].Fill( pb.size() );

      // for clustering:

      pxlist[mpl].push_back(pb);

      ++mpl;
      if( mpl > 5 ) break; // skip others

    } // planes

    ++iev;

  } while( reader->NextEvent() && iev < lev ); // event loop

  delete reader;

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s1 = ts.tv_sec; // seconds since 1.1.1970
  long f1 = ts.tv_nsec; // nanoseconds
  zeit1 += s1 - s0 + ( f1 - f0 ) * 1e-9; // read

  cout << "read " << iev << " events"
	 << " in " << s1 - s0 + ( f1 - f0 ) * 1e-9 << " s"
       << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // hot pixels:

  cout << endl << "Mimosa hot pixel list for run " << run << endl;

  ofstream hotFile( hotFileName.str() );

  hotFile << "# telescope hot pixel list for run " << run
	  << " with " << iev << " events"
	  << endl;

  for( int ipl = 0; ipl <= 5; ++ipl ) {
    hotFile << endl;
    hotFile << "plane " << ipl << endl;
    int nmax = 0;
    int ntot = 0;
    int nhot = 0;
    for( map < int, int >::iterator jpx = pxmap[ipl].begin(); jpx != pxmap[ipl].end(); ++ jpx ) {
      int nhit = jpx->second;
      ntot += nhit;
      if( nhit > nmax ) nmax = nhit;
      if( nhit > iev/128 ) {
	++nhot;
	int ipx = jpx->first;
	int ix = ipx/ny[ipl];
	int iy = ipx%ny[ipl];
	hotFile << "pix "
		<< setw(4) << ix
		<< setw(5) << iy
		<< "  " << nhit
		<< endl;
      }
    } // jpx
    cout
      << "  " << ipl
      << ": active " << pxmap[ipl].size()
      << ", sum " << ntot
      << ", max " << nmax
      << ", hot " << nhot
      << endl;
  } // ipl

  cout << "hot pixel list written to " << hotFileName.str() << endl;

  hotFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // make clusters:

  cout << endl << "parallel clustering" << flush;

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s2 = ts.tv_sec; // seconds since 1.1.1970
  long f2 = ts.tv_nsec; // nanoseconds

  list < vector <cluster> > clist[9];

  //#pragma omp sections // test, not parallel
#pragma omp parallel sections
  {
#pragma omp section
    {
      clist[0] = oneplane( 0, pxlist[0] );
    }
#pragma omp section
    {
      clist[1] = oneplane( 1, pxlist[1] );
    }
#pragma omp section
    {
      clist[2] = oneplane( 2, pxlist[2] );
    }
#pragma omp section
    {
      clist[3] = oneplane( 3, pxlist[3] );
    }
#pragma omp section
    {
      clist[4] = oneplane( 4, pxlist[4] );
    }
#pragma omp section
    {
      clist[5] = oneplane( 5, pxlist[5] );
    }

  } // parallel

  clock_gettime( CLOCK_REALTIME, &ts );
  time_t s3 = ts.tv_sec; // seconds since 1.1.1970
  long f3 = ts.tv_nsec; // nanoseconds
  zeit2 += s3 - s2 + ( f3 - f2 ) * 1e-9; // cluster

  cout << " in " << s3 - s2 + ( f3 - f2 ) * 1e-9 << " s" << endl;

  for( unsigned ipl = 0; ipl <= 5; ++ipl )
    pxlist[ipl].clear(); // memory

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignment iterations:

  int maxiter = aligniteration + 1;
  if( maxiter < 8 ) maxiter = 8;

  for( ; aligniteration < maxiter; ++aligniteration ) {

    cout << endl << "alignment iteration " << aligniteration << endl;

    for( unsigned ipl = 0; ipl <= 5; ++ipl ) {
      hxx[ipl]->Reset();
      hdx[ipl].Reset();
      hdy[ipl].Reset();
      dxvsx[ipl].Reset();
      dxvsy[ipl].Reset();
      dyvsx[ipl].Reset();
      dyvsy[ipl].Reset();
      hdx4[ipl].Reset();
      hdy4[ipl].Reset();
      dx4vsy[ipl].Reset();
      dy4vsx[ipl].Reset();
      effvsx[ipl].Reset();
      effvsxm[ipl].Reset();
      effvsym[ipl].Reset();
      effvsxmym[ipl]->Reset();
    } // pl

    for( int itd = 0; itd < 2; ++itd ) { // triplets 0-1-2 and driplets 3-4-5

      hdxCA[itd].Reset();
      hdyCA[itd].Reset();
      dxCAvsx[itd].Reset();
      dyCAvsy[itd].Reset();

      htridx[itd].Reset();
      htridy[itd].Reset();

      htridxc[itd].Reset();
      htridxci[itd].Reset();
      htridxct[itd].Reset();

      tridxvsx[itd].Reset();
      tridxvsy[itd].Reset();
      tridxvsxm[itd].Reset();
      tridxvsdy[itd]->Reset();
      tridxvsxmym[itd]->Reset();
      tridyvsxmym[itd]->Reset();
      tridxycvsxmym[itd]->Reset();
      tridxyvsxmym[itd]->Reset();
      tridxy2vsxmym[itd]->Reset();
      trimadxvsxm[itd].Reset();
      trimadxvsxmym[itd]->Reset();
      tridxvstx[itd].Reset();
      trimadxvstx[itd].Reset();

      htridxs1[itd].Reset();
      htridxs2[itd].Reset();
      htridxs3[itd].Reset();
      htridxs4[itd].Reset();
      htridxs5[itd].Reset();

      htridxc1[itd].Reset();
      htridxc2[itd].Reset();
      htridxc3[itd].Reset();
      htridxc4[itd].Reset();
      htridxc5[itd].Reset();
      htridxc6[itd].Reset();
      htridxcg[itd].Reset();

      htridyc[itd].Reset();
      htridyci[itd].Reset();
      htridyct[itd].Reset();
      htridyc1[itd].Reset();
      htridyc2[itd].Reset();
      htridyc3[itd].Reset();
      htridyc4[itd].Reset();
      htridyc5[itd].Reset();
      htridyc6[itd].Reset();
      htridycg[itd].Reset();

      tridyvsx[itd].Reset();
      tridyvsy[itd].Reset();
      tridyvsym[itd].Reset();
      trimadyvsym[itd].Reset();
      trimadyvsxmym[itd]->Reset();
      tridyvsty[itd].Reset();
      trimadyvsty[itd].Reset();

      htrix[itd].Reset();
      htriy[itd].Reset();
      htrixy[itd]->Reset();
      htritx[itd].Reset();
      htrity[itd].Reset();

      htrincol[itd].Reset();
      htrinrow[itd].Reset();
      htrinpix[itd].Reset();
      htrixmod1[itd].Reset();
      htrixmod2[itd].Reset();
      htrixmod3[itd].Reset();
      htrixmod4[itd].Reset();
      htrixmod5[itd].Reset();
      htrixmod6[itd].Reset();
      trincolvsxm[itd].Reset();
      trinrowvsym[itd].Reset();
      trinpixvsxmym[itd]->Reset();
      trinpixgvsxmym[itd]->Reset();

      tridxCvsx[itd].Reset();
      tridxCvsy[itd].Reset();
      tridyCvsx[itd].Reset();
      tridyCvsy[itd].Reset();
      tridxCvsax[itd].Reset();
      tridyCvsay[itd].Reset();

      hntri.Reset();
      ntrivsev.Reset();
      hndri.Reset();

    }

    for( unsigned ipl = 0; ipl <= 5; ++ipl ) {

      hexdx[ipl].Reset();
      hexdy[ipl].Reset();
      hexdxc[ipl].Reset();
      exdxvsy[ipl].Reset();
      exdxvstx[ipl].Reset();
      exmadxvstx[ipl].Reset();
      hexdyc[ipl].Reset();
      exdyvsx[ipl].Reset();
      exdyvsty[ipl].Reset();
      exmadyvsty[ipl].Reset();

    }

    hsixdx.Reset();
    hsixdy.Reset();
    hsixdxc.Reset();
    sixdxvsx.Reset();
    sixmadxvsx.Reset();
    sixdxvsy.Reset();
    sixdxvstx.Reset();
    sixmadxvsy.Reset();
    sixmadxvstx.Reset();
    sixmadxvsdtx.Reset();
    hsixdxcsid.Reset();

    hsixdyc.Reset();
    sixdyvsx.Reset();
    sixmadyvsx.Reset();
    sixdyvsy.Reset();
    sixdyvsty.Reset();
    sixmadyvsy.Reset();
    sixmadyvsty.Reset();
    sixmadyvsdty.Reset();

    hsixxy->Reset();
    sixdxyvsxy->Reset();
    hsixdtx.Reset();
    hsixdty.Reset();
    sixdtvsx.Reset();
    sixdtvsxy->Reset();

    // loop over events, correlate planes:

    list < vector <cluster> >::iterator evi[9];
    for( unsigned ipl = 0; ipl <= 5; ++ipl )
      evi[ipl] = clist[ipl].begin();

    unsigned nev = 0;

    cout << "tracking ev";

    for( ; evi[0] != clist[0].end();
	 ++evi[0], ++evi[1], ++evi[2], ++evi[3], ++evi[4], ++evi[5] ) {

      vector <cluster> cl[9]; // Mimosa planes
      for( unsigned ipl = 0; ipl <= 5; ++ipl )
	cl[ipl] = *evi[ipl];

      ++nev;
      if( nev%10000 == 0 )
	cout << " " << nev << flush;

      // final cluster plots:

      if( aligniteration == maxiter-1 ) {

	for( unsigned ipl = 0; ipl <= 5; ++ipl ) {

	  hncl[ipl].Fill( cl[ipl].size() );

	  for( vector<cluster>::iterator cA = cl[ipl].begin(); cA != cl[ipl].end(); ++cA ) {

	    hccol[ipl].Fill( cA->col );
	    hcrow[ipl].Fill( cA->row );
	    unsigned nrowA = cA->scr/(1024*1024);
	    unsigned ncolA = (cA->scr - nrowA*1024*1024)/1024;
	    unsigned npixA = cA->scr % 1024;
	    hnpix[ipl].Fill( npixA );
	    hncol[ipl].Fill( ncolA );
	    hnrow[ipl].Fill( nrowA );
	    hmindxy[ipl].Fill( cA->mindxy );

	  } // clus

	} // planes

      } // laster iter

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

	// A = mid plane:

	for( vector<cluster>::iterator cA = cl[im].begin(); cA != cl[im].end(); ++cA ) {

	  double xA = cA->col*ptchx[im] - alignx[im];
	  double yA = cA->row*ptchy[im] - aligny[im];
	  double xmid = xA - midx[im];
	  double ymid = yA - midy[im];
	  xA = xmid - ymid*rotx[im];
	  yA = ymid + xmid*roty[im];

	  for( int ipl = ibeg; ipl <= iend; ++ipl ) {

	    if( ipl == im ) continue;

	    double sign = ipl - im; // along track: -1 or 1

	    // B = A +- 1

	    for( vector<cluster>::iterator cB = cl[ipl].begin(); cB != cl[ipl].end(); ++cB ) {

	      double xB = cB->col*ptchx[ipl] - alignx[ipl];
	      double yB = cB->row*ptchy[ipl] - aligny[ipl];
	      double xmid = xB - midx[ipl];
	      double ymid = yB - midy[ipl];
	      xB = xmid - ymid*rotx[ipl];
	      yB = ymid + xmid*roty[ipl];

	      double dx = xB - xA;
	      double dy = yB - yA;
	      hxx[ipl]->Fill( xA, xB );
	      hdx[ipl].Fill( dx ); // for shift: fixed sign
	      hdy[ipl].Fill( dy );
	      dxvsx[ipl].Fill( xB, dx*sign ); // for turn angle: sign along track
	      dxvsy[ipl].Fill( yB, dx      );
	      dyvsx[ipl].Fill( xB, dy      );
	      dyvsy[ipl].Fill( yB, dy*sign ); // for tilt angle: sign along track

	    } // clusters

	  } // ipl

	} // cl mid

      } // upstream and downstream internal correlations

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // telescope plane efficiency:

      double isoCut = 0.30/ptchx[1]; // [px]
      double triCut = 0.05; // [mm]
      double effCut = 0.25; // [mm]

      for( int ipl = 0; ipl <= 5; ++ipl ) {

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

	double zD = zz[ipl] + alignz[ipl];

	for( vector<cluster>::iterator cA = cl[ib].begin(); cA != cl[ib].end(); ++cA ) {

	  if( cA->mindxy < isoCut ) continue;

	  double xA = cA->col*ptchx[ib] - alignx[ib];
	  double yA = cA->row*ptchy[ib] - aligny[ib];
	  double zA = zz[ib] + alignz[ib];
	  double xmid = xA - midx[ib];
	  double ymid = yA - midy[ib];
	  xA = xmid - ymid*rotx[ib];
	  yA = ymid + xmid*roty[ib];

	  for( vector<cluster>::iterator cC = cl[ie].begin(); cC != cl[ie].end(); ++cC ) {

	    if( cC->mindxy < isoCut ) continue;

	    double xC = cC->col*ptchx[ie] - alignx[ie];
	    double yC = cC->row*ptchy[ie] - aligny[ie];
	    double zC = zz[ie] + alignz[ie];
	    double xmid = xC - midx[ie];
	    double ymid = yC - midy[ie];
	    xC = xmid - ymid*rotx[ie];
	    yC = ymid + xmid*roty[ie];

	    double dx2 = xC - xA;
	    double dy2 = yC - yA;
	    double dzCA = zC - zA;

	    if( fabs( dx2 ) > ang * dzCA ) continue; // angle cut
	    if( fabs( dy2 ) > ang * dzCA ) continue; // angle cut

	    double xavg2 = 0.5*(xA + xC);
	    double yavg2 = 0.5*(yA + yC);
	    double zavg2 = 0.5*(zA + zC);

	    double slpx = ( xC - xA ) / dzCA; // slope x
	    double slpy = ( yC - yA ) / dzCA; // slope y

	    for( vector<cluster>::iterator cB = cl[im].begin(); cB != cl[im].end(); ++cB ) {

	      double xB = cB->col*ptchx[im] - alignx[im]; // stretch and shift
	      double yB = cB->row*ptchy[im] - aligny[im];
	      double zB = zz[im] + alignz[im];
	      double xmid = xB - midx[im];
	      double ymid = yB - midy[im];
	      xB = xmid - ymid*rotx[im];
	      yB = ymid + xmid*roty[im];

	      // interpolate track to B:

	      double dz = zB - zavg2;
	      double xm = xavg2 + slpx * dz; // triplet at B
	      double ym = yavg2 + slpy * dz;

	      double dxm = xB - xm;
	      double dym = yB - ym;

	      if( fabs( dxm ) > triCut ) continue;
	      if( fabs( dym ) > triCut ) continue;

	      // inter/extrapolate track to D:

	      double da = zD - zavg2;
	      double xi = xavg2 + slpx * da; // triplet at D
	      double yi = yavg2 + slpy * da;

	      // transform into local frame:

	      double xr = xi + yi*rotx[ipl] + alignx[ipl];
	      double yr = yi - xi*roty[ipl] + aligny[ipl];

	      if( fabs( xr ) > 10.4 ) continue; // fiducial
	      if( fabs( yr ) >  5.2 ) continue; // fiducial

	      // eff pl:

	      int nm = 0;

	      for( vector<cluster>::iterator cD = cl[ipl].begin(); cD != cl[ipl].end(); ++cD ) {

		double xD = cD->col*ptchx[ipl] - midx[ipl];
		double yD = cD->row*ptchy[ipl] - midy[ipl];

		double dx4 = xD - xr;
		double dy4 = yD - yr;

		if( fabs( dy4 ) < effCut ) {
		  hdx4[ipl].Fill( dx4 );
		  dx4vsy[ipl].Fill( yr, dx4 );
		}
		if( fabs( dx4 ) < effCut ) {
		  hdy4[ipl].Fill( dy4 );
		  dy4vsx[ipl].Fill( xr, dy4 );
		}

		if( fabs( dx4 ) > effCut ) continue;
		if( fabs( dy4 ) > effCut ) continue;

		++nm;

		if( nm > 0 ) break; // one link is enough

	      } // cl D

	      effvsx[ipl].Fill( xr, nm );

	      double xmod2 = fmod( xr + sizex[ipl] + 0.5*ptchx[ipl], 2*ptchx[ipl] );
	      double ymod2 = fmod( yr + sizey[ipl] + 0.5*ptchy[ipl], 2*ptchy[ipl] );
	      effvsxm[ipl].Fill( xmod2*1E3, nm );
	      effvsym[ipl].Fill( ymod2*1E3, nm );
	      effvsxmym[ipl]->Fill( xmod2*1E3, ymod2*1E3, nm );

	    } // cl B

	  } // cl C

	} // cl A

      } // eff planes

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // triplets 1 vs 2-0:
      // driplets 4 vs 5-3:

      vector <triplet> triplets;
      vector <triplet> driplets;

      double tricut = 0.050; // [mm]

      for( int itd = 0; itd < 2; ++itd ) { // triplets 0-1-2 and driplets 3-4-5

	int ib = 0;
	int ie = 2;
	int im = 1; // mid plane triplet
	if( itd == 1 ) {
	  ib = 3;
	  ie = 5;
	  im = 4; // mid plane driplet
	}
	double zA = zz[ib] + alignz[ib];
	double zC = zz[ie] + alignz[ie];
	double zB = zz[im] + alignz[im];
	double dzCA = zC - zA;

	for( vector<cluster>::iterator cA = cl[ib].begin(); cA != cl[ib].end(); ++cA ) {

	  double xA = cA->col*ptchx[ib] - alignx[ib];
	  double yA = cA->row*ptchy[ib] - aligny[ib];
	  double xmid = xA - midx[ib];
	  double ymid = yA - midy[ib];
	  xA = xmid - ymid*rotx[ib];
	  yA = ymid + xmid*roty[ib];

	  unsigned nrowA = cA->scr/(1024*1024);
	  unsigned ncolA = (cA->scr - nrowA*1024*1024)/1024;
	  //unsigned npixA = cA->scr % 1024;
	  bool goodncolA = 1;
	  if( ncolA == 2 || ncolA > 4 ) goodncolA = 0;
	  bool goodnrowA = 1;
	  if( nrowA == 2 || nrowA > 4 ) goodnrowA = 0;

	  for( vector<cluster>::iterator cC = cl[ie].begin(); cC != cl[ie].end(); ++cC ) {

	    double xC = cC->col*ptchx[ie] - alignx[ie];
	    double yC = cC->row*ptchy[ie] - aligny[ie];
	    double xmid = xC - midx[ie];
	    double ymid = yC - midy[ie];
	    xC = xmid - ymid*rotx[ie];
	    yC = ymid + xmid*roty[ie];

	    double dx2 = xC - xA;
	    double dy2 = yC - yA;
	    hdxCA[itd].Fill( dx2 );
	    hdyCA[itd].Fill( dy2 );
	    if( fabs( dy2 ) < 0.001 * dzCA )
	      dxCAvsx[itd].Fill( xC, dx2 );
	    if( fabs( dx2 ) < 0.001 * dzCA )
	      dyCAvsy[itd].Fill( yC, dy2 );

	    if( fabs( dx2 ) > ang * dzCA ) continue; // angle cut
	    if( fabs( dy2 ) > ang * dzCA ) continue; // angle cut

	    double xavg2 = 0.5*(xA + xC);
	    double yavg2 = 0.5*(yA + yC);
	    double zavg2 = 0.5*(zA + zC);

	    double slpx = ( xC - xA ) / dzCA; // slope x
	    double slpy = ( yC - yA ) / dzCA; // slope y

	    // interpolate track to B:

	    double dz = zB - zavg2;
	    double xm = xavg2 + slpx * dz; // triplet at B
	    double ym = yavg2 + slpy * dz;

	    // transform into local frame:

	    double xr = xm + ym*rotx[im] + alignx[im];
	    double yr = ym - xm*roty[im] + aligny[im];

	    double xmod2 = fmod( xr + sizex[im] + 0.5*ptchx[im], 2*ptchx[im] );
	    double ymod2 = fmod( yr + sizey[im] + 0.5*ptchy[im], 2*ptchy[im] );
	    double xmod4 = fmod( xr + sizex[im] + 0.5*ptchx[im], 4*ptchx[im] );
	    double ymod4 = fmod( yr + sizey[im] + 0.5*ptchy[im], 4*ptchy[im] );

	    unsigned nrowC = cC->scr/(1024*1024);
	    unsigned ncolC = (cC->scr - nrowC*1024*1024)/1024;
	    //unsigned npixC = cC->scr % 1024;
	    bool goodncolC = 1;
	    if( ncolC == 2 || ncolC > 4 ) goodncolC = 0;
	    bool goodnrowC = 1;
	    if( nrowC == 2 || nrowC > 4 ) goodnrowC = 0;

	    for( vector<cluster>::iterator cB = cl[im].begin(); cB != cl[im].end(); ++cB ) {

	      double xB = cB->col*ptchx[im] - alignx[im]; // stretch and shift
	      double yB = cB->row*ptchy[im] - aligny[im];
	      double xmid = xB - midx[im];
	      double ymid = yB - midy[im];
	      xB = xmid - ymid*rotx[im];
	      yB = ymid + xmid*roty[im];

	      double dxm = xB - xm;
	      double dym = yB - ym;

	      htridx[itd].Fill( dxm*1E3 );
	      htridy[itd].Fill( dym*1E3 );

	      bool iso = 1;
	      if( cA->mindxy < isoCut ) iso = 0;
	      if( cC->mindxy < isoCut ) iso = 0;
	      if( cB->mindxy < isoCut ) iso = 0;

	      unsigned nrowB = cB->scr/(1024*1024);
	      unsigned ncolB = (cB->scr - nrowB*1024*1024)/1024;
	      unsigned npixB = cB->scr % 1024;

	      if( ncolB > 99 )
		cout << "scrB " << cB->scr
		     << ", nrow " << nrowB
		     << ", ncol " << ncolB
		     << ", npix " << npixB
		     << endl;

	      bool goodncolB = 1;
	      if( ncolB == 2 || ncolB > 4 ) goodncolB = 0;
	      bool goodnrowB = 1;
	      if( nrowB == 2 || nrowB > 4 ) goodnrowB = 0;

	      if( fabs( dym ) < 0.02 ) {

		htridxc[itd].Fill( dxm*1E3 );
		if( iso ) htridxci[itd].Fill( dxm*1E3 );
		tridxvsx[itd].Fill( xm, dxm*1E3 );
		tridxvsy[itd].Fill( ym, dxm*1E3 );
		tridxvsxm[itd].Fill( xmod2*1E3, dxm*1E3 );
		trimadxvsxm[itd].Fill( xmod2*1E3, fabs(dxm)*1E3 );
		trimadxvsxmym[itd]->Fill( xmod2*1E3, ymod2*1E3, fabs(dxm)*1E3 );
		tridxvstx[itd].Fill( slpx*1E3, dxm*1E3 ); // adjust zpos, same sign
		trimadxvstx[itd].Fill( slpx*1E3, fabs(dxm)*1E3 ); // U-shape
		if( fabs( slpx < 0.001 ) )
		  htridxct[itd].Fill( dxm*1E3 );

		if(      npixB == 1 )
		  htridxs1[itd].Fill( dxm*1E3 ); // 4.2 um
		else if( npixB == 2 )
		  htridxs2[itd].Fill( dxm*1E3 ); // 4.0 um
		else if( npixB == 3 )
		  htridxs3[itd].Fill( dxm*1E3 ); // 3.8 um
		else if( npixB == 4 )
		  htridxs4[itd].Fill( dxm*1E3 ); // 4.3 um
		else
		  htridxs5[itd].Fill( dxm*1E3 ); // 3.6 um

		if(      ncolB == 1 )
		  htridxc1[itd].Fill( dxm*1E3 ); // 3.5 um
		else if( ncolB == 2 )
		  htridxc2[itd].Fill( dxm*1E3 ); // 4.2 um
		else if( ncolB == 3 )
		  htridxc3[itd].Fill( dxm*1E3 ); // 3.5 um
		else if( ncolB == 4 )
		  htridxc4[itd].Fill( dxm*1E3 ); // 3.5 um
		else if( ncolB == 5 )
		  htridxc5[itd].Fill( dxm*1E3 ); // 
		else
		  htridxc6[itd].Fill( dxm*1E3 ); // 

		if( goodncolA && goodncolB && goodncolC )
		  htridxcg[itd].Fill( dxm*1E3 ); // 

	      } // dy

	      if( fabs( dxm ) < 0.02 ) {

		htridyc[itd].Fill( dym*1E3 );
		if( iso ) htridyci[itd].Fill( dym*1E3 );
		tridyvsx[itd].Fill( xm, dym*1E3 );
		tridyvsy[itd].Fill( ym, dym*1E3 );
		tridyvsym[itd].Fill( ymod2*1E3, dym*1E3 );
		trimadyvsym[itd].Fill( ymod2*1E3, fabs(dym)*1E3 );
		trimadyvsxmym[itd]->Fill( xmod2*1E3, ymod2*1E3, fabs(dym)*1E3 );
		tridyvsty[itd].Fill( slpy*1E3, dym*1E3 );
		trimadyvsty[itd].Fill( slpy*1E3, fabs(dym)*1E3 ); // U-shape
		if( fabs( slpy < 0.001 ) )
		  htridyct[itd].Fill( dym*1E3 );

		if(      nrowB == 1 )
		  htridyc1[itd].Fill( dym*1E3 ); // 
		else if( nrowB == 2 )
		  htridyc2[itd].Fill( dym*1E3 ); // 
		else if( nrowB == 3 )
		  htridyc3[itd].Fill( dym*1E3 ); // 
		else if( nrowB == 4 )
		  htridyc4[itd].Fill( dym*1E3 ); // 
		else if( nrowB == 5 )
		  htridyc5[itd].Fill( dym*1E3 ); // 
		else
		  htridyc6[itd].Fill( dym*1E3 ); // 

		if( goodnrowA && goodnrowB && goodnrowC )
		  htridycg[itd].Fill( dym*1E3 ); // 

	      }

	      if( fabs( dxm ) < 0.020 && fabs( dym ) < 0.020 ) { // for tight z spacing

		tridxvsxmym[itd]->Fill( xmod2*1E3, ymod2*1E3, dxm * 1E3 );
		tridyvsxmym[itd]->Fill( xmod2*1E3, ymod2*1E3, dym * 1E3 );
		tridxvsdy[itd]->Fill( dxm*1E3, dym*1E3 ); // correlated residuals? no
		tridxycvsxmym[itd]->Fill( xmod2*1E3, ymod2*1E3, dxm * dym * 1E6 ); // correlation map
		tridxyvsxmym[itd]->Fill( xmod2*1E3, ymod2*1E3, ( dxm + dym ) * 1E3 ); // shift
		tridxy2vsxmym[itd]->Fill( xmod2*1E3, ymod2*1E3, ( dxm*dxm + dym*dym ) * 1E6 ); // resolution

	      }

	      // store triplets:

	      if( fabs( dxm ) < tricut && fabs( dym ) < tricut ) {

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
		htritx[itd].Fill( slpx*1E3 );
		htrity[itd].Fill( slpy*1E3 );

		htrincol[itd].Fill( ncolB );
		htrinrow[itd].Fill( nrowB );
		htrinpix[itd].Fill( npixB );

		htrixmod[itd].Fill( xmod2*1E3 );
		if(      ncolB == 1 )
		  htrixmod1[itd].Fill( xmod2*1E3 );
		else if( ncolB == 2 )
		  htrixmod2[itd].Fill( xmod2*1E3 );
		else if( ncolB == 3 )
		  htrixmod3[itd].Fill( xmod2*1E3 );
		else if( ncolB == 4 )
		  htrixmod4[itd].Fill( xmod2*1E3 );
		else if( ncolB == 5 )
		  htrixmod5[itd].Fill( xmod2*1E3 );
		else
		  htrixmod6[itd].Fill( xmod2*1E3 );

		trincolvsxm[itd].Fill( xmod4*1E3, ncolB );
		trinrowvsym[itd].Fill( ymod4*1E3, nrowB );
		trinpixvsxmym[itd]->Fill( xmod4*1E3, ymod4*1E3, npixB );
		if( goodnrowA && goodnrowC && goodncolA && goodncolC ) // better resolution, lower stat
		  trinpixgvsxmym[itd]->Fill( xmod4*1E3, ymod4*1E3, npixB );

		// check z spacing: A-B as baseline

		double dzAB = zB - zA;
		double ax = ( xB - xA ) / dzAB; // slope x
		double ay = ( yB - yA ) / dzAB; // slope y
		double dz = zC - zB;
		double xk = xB + ax * dz; // at C
		double yk = yB + ay * dz; // at C
		double dx = xC - xk;
		double dy = yC - yk;
		tridxCvsx[itd].Fill( xk, dx*1E3 );
		tridxCvsy[itd].Fill( yk, dx*1E3 );
		tridyCvsx[itd].Fill( xk, dy*1E3 );
		tridyCvsy[itd].Fill( yk, dy*1E3 );
		tridxCvsax[itd].Fill( ax*1E3, dx*1E3 ); // adjust zpos, same sign
		tridyCvsay[itd].Fill( ay*1E3, dy*1E3 );

	      } // triplet

	    } // cl B

	  } // cl C

	} // cl A

      } // triplets and driplets

      hntri.Fill( triplets.size() );
      ntrivsev.Fill( nev, triplets.size() );
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

	  double zA = zz[ipl] + alignz[ipl] - avzA; // z from mid of triplet to plane
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
	    hexdx[ipl].Fill( dx*1E3 );
	    hexdy[ipl].Fill( dy*1E3 );
	    if( fabs( dy ) < 0.5 ) {
	      hexdxc[ipl].Fill( dx*1E3 );
	      exdxvsy[ipl].Fill( yC, dx*1E3 );
	      exdxvstx[ipl].Fill( slxA*1E3, dx*1E3 );
	      exmadxvstx[ipl].Fill( slxA*1E3, fabs(dx)*1E3 );
	    }
	    if( fabs( dx ) < 0.5 ) {
	      hexdyc[ipl].Fill( dy*1E3 );
	      exdyvsx[ipl].Fill( xC, dy*1E3 );
	      exdyvsty[ipl].Fill( slyA*1E3, dy*1E3 );
	      exmadyvsty[ipl].Fill( slyA*1E3, fabs(dy)*1E3 );
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

	  double zB = zz[ipl] + alignz[ipl] - avzB; // z from mid of driplet to plane
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
	    hexdx[ipl].Fill( dx*1E3 );
	    hexdy[ipl].Fill( dy*1E3 );
	    if( fabs( dy ) < 0.5 ) {
	      hexdxc[ipl].Fill( dx*1E3 );
	      exdxvsy[ipl].Fill( yC, dx*1E3 );
	      exdxvstx[ipl].Fill( slxB*1E3, dx*1E3 );
	      exmadxvstx[ipl].Fill( slxB*1E3, fabs(dx)*1E3 );
	    }
	    if( fabs( dx ) < 0.5 ) {
	      hexdyc[ipl].Fill( dy*1E3 );
	      exdyvsx[ipl].Fill( xC, dy*1E3 );
	      exdyvsty[ipl].Fill( slyB*1E3, dy*1E3 );
	      exmadyvsty[ipl].Fill( slyB*1E3, fabs(dy)*1E3 );
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

	  if( fabs(dy) < 0.100 ) {

	    hsixdxc.Fill( dx*1E3 );

	    sixdxvsx.Fill( xA, dx*1E3 );
	    sixmadxvsx.Fill( xA, fabs(dx)*1E3 );
	    sixdxvsy.Fill( yA, dx*1E3 );
	    sixdxvstx.Fill( slxA*1E3, dx*1E3 );
	    sixmadxvsy.Fill( yA, fabs(dx)*1E3 );
	    sixmadxvstx.Fill( slxA*1E3, fabs(dx)*1E3 );
	    sixmadxvsdtx.Fill( dtx*1E3, fabs(dx)*1E3 ); // U-shape
	    if( fabs( dtx ) < 0.0005 )
	      hsixdxcsid.Fill( dx*1E3 );

	  } // dy

	  if( fabs(dx) < 0.100 ) {

	    hsixdyc.Fill( dy*1E3 );

	    sixdyvsx.Fill( xA, dy*1E3 );
	    sixmadyvsx.Fill( xA, fabs(dy)*1E3 );
	    sixdyvsy.Fill( yA, dy*1E3 );
	    sixdyvsty.Fill( slyA*1E3, dy*1E3 );
	    sixmadyvsy.Fill( yA, fabs(dy)*1E3 );
	    sixmadyvsty.Fill( slyA*1E3, fabs(dy)*1E3 );
	    sixmadyvsdty.Fill( dty*1E3, fabs(dy)*1E3 ); // U-shape

	  }

	  // compare slopes:

	  if( fabs(dy) < 0.100 && fabs(dx) < 0.100 ) {

	    hsixxy->Fill( xA, yA );
	    sixdxyvsxy->Fill( xA, yA, dxy );

	    hsixdtx.Fill( dtx*1E3 );
	    hsixdty.Fill( dty*1E3 );
	    sixdtvsx.Fill( xA, dtxy );
	    sixdtvsxy->Fill( xA, yA, dtxy );

	  } // match

	} // driplets

      } // triplets

    } // events

    cout << endl;

    clock_gettime( CLOCK_REALTIME, &ts );
    time_t s3 = ts.tv_sec; // seconds since 1.1.1970
    long f3 = ts.tv_nsec; // nanoseconds
    zeit3 += s3 - s2 + ( f3 - f2 ) * 1e-9; // track

    clock_gettime( CLOCK_REALTIME, &ts );
    time_t s9 = ts.tv_sec; // seconds since 1.1.1970
    long f9 = ts.tv_nsec; // nanoseconds

    cout << endl
	 << "done after " << iev << " events"
	 << " in " << s9 - s0 + ( f9 - f0 ) * 1e-9 << " s"
	 << " (read " << zeit1
	 << " s, cluster " << zeit2
	 << " s, tracking " << zeit3 << " s)"
	 << endl;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // alignment fits:

    if( maxiter == aligniteration + 1 ) {      
      histoFile->Write(); // before fitting
      histoFile->Close();
    }

    cout << endl << "alignment fits:" << endl;

    for( int ipl = 0; ipl <= 5; ++ipl ) {

      double nb = hdx[ipl].GetNbinsX();
      double ne = hdx[ipl].GetSumOfWeights();
      double nm = hdx[ipl].GetMaximum();

      cout << endl << hdx[ipl].GetTitle() << endl;
      if( nm < 99 ) {
	cout << "  peak " << nm << " not enough" << endl;
	continue;
      }

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

      alignx[ipl] += fgp0x->GetParameter(1);

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

      aligny[ipl] += fgp0y->GetParameter(1);

      // x-y rotation:

      if( aligniteration >= 1 && dyvsx[ipl].GetEntries() > 999 ) {

	dyvsx[ipl].Fit( "pol1", "q", "", -midx[ipl], midx[ipl] );
	TF1 * fdyvsx = dyvsx[ipl].GetFunction( "pol1" );
	cout << endl << dyvsx[ipl].GetTitle()
	     << " slope " << fdyvsx->GetParameter(1)
	     << " = rot"
	     << endl;
	roty[ipl] -= fdyvsx->GetParameter(1); // sign

      }

      if( aligniteration >= 1 && dxvsy[ipl].GetEntries() > 999 ) {

	dxvsy[ipl].Fit( "pol1", "q", "", -midy[ipl], midy[ipl] );
	TF1 * fdxvsy = dxvsy[ipl].GetFunction( "pol1" );
	cout << endl << dxvsy[ipl].GetTitle()
	     << " slope " << fdxvsy->GetParameter(1)
	     << " = rot"
	     << endl;
	rotx[ipl] += fdxvsy->GetParameter(1);

      }

    } // ipl

    // z-shift:

    if( aligniteration >= 3 ) {

      for( int itd = 0; itd < 2; ++itd ) { // upstream and downstream

	cout << endl;

	int ipl = 2+3*itd; // 2 or 5
	/*
	if( tridxCvsax[itd].GetEntries() > 999 ) {

	  tridxCvsax[itd].Fit( "pol1", "q", "", -1, 1 ); // um vs mrad
	  TF1 * f1 = tridxCvsax[itd].GetFunction( "pol1" );
	  alignz[ipl-2] += 0.5*f1->GetParameter(1);
	  alignz[ipl]   += 0.5*f1->GetParameter(1);
	  cout << tridxCvsax[itd].GetTitle()
	       << " dz " << f1->GetParameter(1)
	       << ", plane " << ipl-2
	       << " new zpos " << zz[ipl-2] + alignz[ipl-2]
	       << ", plane " << ipl
	       << " new zpos " << zz[ipl] + alignz[ipl]
	       << endl;

	}
	*/
	if( tridxvstx[itd].GetEntries() > 999 ) {

	  tridxvstx[itd].Fit( "pol1", "q", "", -1, 1 ); // um vs mrad
	  TF1 * f1 = tridxvstx[itd].GetFunction( "pol1" );
	  alignz[ipl-2] -= 0.5*f1->GetParameter(1);
	  alignz[ipl]   -= 0.5*f1->GetParameter(1);
	  cout << tridxvstx[itd].GetTitle()
	       << " dz " << f1->GetParameter(1)
	       << ", plane " << ipl-2
	       << " new zpos " << zz[ipl-2] + alignz[ipl-2]
	       << ", plane " << ipl
	       << " new zpos " << zz[ipl] + alignz[ipl]
	       << endl;

	}

      } // itd

    } // aligniteration

    // driplet vs triplet:

    if( aligniteration >= 4 && hsixdx.GetMaximum() > 99 ) {

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

      for( int ipl = 3; ipl <= 5; ++ipl ) {
	alignx[ipl] += fgp0x->GetParameter(1);
	aligny[ipl] += fgp0y->GetParameter(1);
      }

    } // aligniteration

    if( aligniteration >= 5 && hsixdx.GetMaximum() > 99 ) {

      // x-y rotation from profiles:

      if( sixdxvsy.GetEntries() > 999 ) {
	sixdxvsy.Fit( "pol1", "q", "", -midy[2], midy[2] );
	TF1 * fdxvsy = sixdxvsy.GetFunction( "pol1" );
	cout << endl << sixdxvsy.GetTitle()
	     << " slope " << fdxvsy->GetParameter(1)
	     << " mu/mm" << endl;
	for( int ipl = 3; ipl <= 5; ++ipl )
	  rotx[ipl] += fdxvsy->GetParameter(1)*1E-3;
      }

      if( sixdyvsx.GetEntries() > 999 ) {
	sixdyvsx.Fit( "pol1", "q", "", -midx[2], midx[2] );
	TF1 * fdyvsx = sixdyvsx.GetFunction( "pol1" );
	cout << endl << sixdyvsx.GetTitle()
	     << " slope " << fdyvsx->GetParameter(1)
	     << " mu/mm" << endl;
	for( int ipl = 3; ipl <= 5; ++ipl )
	  roty[ipl] -= fdyvsx->GetParameter(1)*1E-3; // sign
      }

    } // aligniteration

      // dz from dx vs tx:

    if( aligniteration >= 6 && sixdxvstx.GetEntries() > 999 ) {

      sixdxvstx.Fit( "pol1", "q", "", -1, 1 ); // [mu/mrad]
      TF1 * f1 = sixdxvstx.GetFunction( "pol1" );
      cout << endl << sixdxvstx.GetTitle()
	   << "  dz " << f1->GetParameter(1)
	   << endl;

      for( int ipl = 3; ipl <= 5; ++ipl )
	alignz[ipl] += f1->GetParameter(1);

    } // aligniteration

  } // aligniterations

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write alignment to file:

  ofstream alignFile( alignFileName.str() );

  alignFile << "# telescope alignment for run " << run << endl;

  alignFile << "iteration " << aligniteration << endl;

  for( int ipl = 0; ipl <= 5; ++ipl ) {
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
  if( aligniteration <= 7 )
    cout << endl << "need more align iterations: please run again!" << endl;

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  return 0;
}

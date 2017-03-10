
// Daniel Pitzl, DESY, Feb 2016
// telescope analysis with eudaq and DUT

// make scope
// scope -l 9999 20833  # tilted DUT 504
// scope 19037
// scope -s 23403 # re-sync
// scope -t 3.2 23494 # apply threshold [ke]

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TFile.h>
#include <TH1I.h> // counting
#include <TH1D.h> // weighted counts
#include <TH2I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>
#include <TSystem.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <set>
#include <cmath>

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int adc;
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
      c.charge += Qpix;
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

    if( ! ( c.charge == 0 ) ) {
      c.col /= sumQ;
      c.row /= sumQ;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with zero charge" << endl;
    }

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
  bool useMODasREF = false;
  cout << "use MOD as reference plane" << useMODasREF << endl;
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // further arguments:

  // Default runs file
  string runFileName( "runs.dat" );

  int lev = 999222111; // last event

  bool syncdut = 0; // re-sync required ?
  bool syncref = 0; // re-sync required ?
  bool syncmod = 0; // re-sync required ?

  double thr = 0; // offline pixel threshold [ke]

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

    if( !strcmp( argv[i], "-t" ) )
      thr = atof( argv[++i] ); // [ke]

    if( !strcmp( argv[i], "-s" ) ) {
      syncdut = 1;
      syncref = 1;
      syncmod = 1;
    }

    if( !strcmp( argv[i], "-d" ) )
      syncdut = 1;

    if( !strcmp( argv[i], "-r" ) )
      syncref = 1;

    if( !strcmp( argv[i], "-f" ) )
      runFileName = argv[++i];

  } // argc

  if( syncdut )
    cout << "re-sync DUT" << endl;
  if( syncref )
    cout << "re-sync REF" << endl;

  cout << "apply offline pixel threshold at " << thr << " ke" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  string geoFileName( "geo.dat" );
  double DUTtilt0 = 19.3;
  double MODtilt0 = 16.5;
  double MODturn0 = 27.8;
  double pbeam = 5.6;
  int chip0 = 504;
  int refchip0 = 501;
  string gainFileName( "gain.dat" );
  string refgainFileName( "gain.dat" );
  int weib = 3;

  ifstream runsFile( runFileName );

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
    string WEIB( "weib" );
    string GAIN( "gain" );
    string REFCHIP( "refchip" );
    string REFGAIN( "refgain" );
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

      if( tag == REFGAIN ) {
	tokenizer >> refgainFileName;
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

      if( tag == REFCHIP ) {
	tokenizer >> refchip0;
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
	<< "  Weibull version " << weib << endl
	<< "  REF chip " << refchip0 << endl
	<< "  REF gain file " << refgainFileName << endl
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
  // Create a directory based on run number
  stringstream run_number;
  run_number << run;

  // Directory created will be created in the current working directory (cwd)
  std::string outputDirectory = std::string(gSystem->WorkingDirectory())+"/run_"+run_number.str();

  // Check first if directory does not already exist
  if(!gSystem->MakeDirectory(TString(outputDirectory)))
    gSystem->MakeDirectory(TString(outputDirectory));

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignments:

  int aligniteration = 0;
  double alignx[9];
  double aligny[9];
  double rotx[9];
  double roty[9];

  ostringstream alignFileName; // output string stream

  alignFileName << "align_" << run << ".dat";

  ifstream ialignFile( outputDirectory+"/"+alignFileName.str() );

  cout << endl;

  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    cout << "Error opening " << alignFileName.str() << endl
	 << "  please do: tele -g " << geoFileName << " " << run << endl
	 << endl;
    return 1;
  }
  // can there be instructions between if and else ?
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

  ifstream ihotFile( outputDirectory+"/"+hotFileName.str() );

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
  // DUT:

  double wt = atan(1.0) / 45.0; // pi/180 deg

  double qwid = 3.5; // [ke] for Moyal

  bool rot90 = 0; // 504
  //bool rot90 = 1; // 506

  double upsign = 1; // 504

  double upsignx = upsign;
  double upsigny = upsign;
  if( chip0 == 603 && run > 24290 ) // tilt 180 deg
    upsigny = -upsign;
  
  if( rot90 ) {
    upsignx =  1;
    upsigny = -1;
  }

  int iDUT = 7;

  int DUTaligniteration = 0;
  double DUTalignx = 0.0;
  double DUTaligny = 0.0;
  double DUTrot = 0.0;
  double DUTturn = 0;
  double DUTtilt = DUTtilt0; // [deg]
  double DUTz = 40 + zz[2];

  ostringstream DUTalignFileName; // output string stream

  DUTalignFileName << "alignDUT_" << run << ".dat";

  ifstream iDUTalignFile( outputDirectory+"/"+DUTalignFileName.str() );

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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // REF:

  int iREF = 8;

  int REFaligniteration = 0;
  double REFalignx = 0.0;
  double REFaligny = 0.0;
  double REFrot = 0.0;
  double REFz = 55 + zz[5];

  ostringstream REFalignFileName; // output string stream

  REFalignFileName << "alignREF_" << run << ".dat";

  ifstream iREFalignFile( outputDirectory+"/"+REFalignFileName.str() );

  cout << endl;

  if( iREFalignFile.bad() || ! iREFalignFile.is_open() ) {
    cout << "no " << REFalignFileName.str() << ", will bootstrap" << endl;
  }
  else {

    cout << "read REFalignment from " << REFalignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string alignx( "alignx" );
    string aligny( "aligny" );
    string rot( "rot" );
    string dz( "dz" );

    while( ! iREFalignFile.eof() ) {

      string line;
      getline( iREFalignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) 
	tokenizer >> REFaligniteration;

      double val;
      tokenizer >> val;
      if(      tag == alignx )
	REFalignx = val;
      else if( tag == aligny )
	REFaligny = val;
      else if( tag == rot )
	REFrot = val;
      else if( tag == dz )
	REFz = val + zz[5];

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  iREFalignFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOD:

  int iMOD = 6;

  int MODaligniteration = 0;
  double MODalignx = 0.0;
  double MODaligny = 0.0;
  double MODrot = 0.0;
  double MODtilt = MODtilt0; // [deg]
  double MODturn = MODturn0; // [deg]
  double MODz = 47 + zz[4];

  ostringstream MODalignFileName; // output string stream

  MODalignFileName << "alignMOD_" << run << ".dat";

  ifstream iMODalignFile( outputDirectory+"/"+MODalignFileName.str() );

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
  // DUT gain:

  double p0[52][80]; // Weibull
  double p1[52][80];
  double p2[52][80];
  double p3[52][80];
  double p4[52][80];
  double p5[52][80];

  ifstream gainFile( gainFileName );

  if(! gainFile ) {
    cout << "gain file " << gainFileName << " not found" << endl;
    return 1;
  }
  else {
    cout << endl << "using DUT pixel gain file " << gainFileName << endl;
    char ih[99];
    int col;
    int row;
    while( gainFile >> ih ) {
      gainFile >> col;
      gainFile >> row;
      if( col < 0 || col > 51 || row < 0 || row > 79 ) {
	cout << "invalid pixel in gain file " << col << " " << row << endl;
	continue;
      }
      gainFile >> p0[col][row];
      gainFile >> p1[col][row];
      gainFile >> p2[col][row];
      gainFile >> p3[col][row];
      gainFile >> p4[col][row];
      gainFile >> p5[col][row]; // gain ratio
    } // while

  } // gainFile

  double ke = 0.250; // to get q0f peak at 22 ke

  if( chip0 == 500 )
    ke = 0.305; // 14393 to get q0f peak at 22 ke no eps in Q

  if( chip0 == 504 ) {
    ke = 0.254; // 14614 to get q0f peak at 22 ke
    if( run >= 19019 ) // Apr 2015
      ke = 0.243; // 19037 to get q0f peak at 22 ke
    if( run >= 20785 ) // Jul 2015
      ke = 0.235;
    if( run >= 20823 ) // 12.7.2015 chiller 17
      ke = 0.244; // to get q0f peak at 22 ke
    if( run >= 23000 ) // Mar 2016 chiller 17
      ke = 0.259; // to get q0f peak at 22 ke
    if( run >= 23339 ) // Mar 2016 6 GeV trim 36
      ke = 0.249; // to get q0f peak at 22 ke
    if( run >= 23350 ) // Mar 2016 chiller 17
      ke = 0.259; // to get q0f peak at 22 ke
    if( run >= 23518 ) // Mar 2016 chiller 17
      ke = 0.254; // to get q0f peak at 22 ke
  }

  if( chip0 == 506 ) {
    ke = 0.290; // 14654 to get q0f peak at 22 ke
    if( run >= 19441 ) // chiller off
      ke = 0.263;
    if( run >= 19582 ) // chiller off
      ke = 0.268;
    if( run >= 19702 ) // chiller on
      ke = 0.257;
    if( run >= 20744 ) // chiller off
      ke = 0.268;
  }

  if( chip0 == 603 )
    ke = 0.288; // 24434

  // tsunami correction:

  double eps = 0.0;
  if( chip0 == 504 )
    eps = 0.05; // 2016 sy 7.66
  //eps = 0.06; // 2015 sy 7.74
  //eps = 0.04; // 2016 sy 7.72

  if( chip0 == 603 )
    eps = 0.11;

  // correct trend in q vs y:

  double yco = 0.0;
  if( chip0 == 504 )
    yco = 0.0009;

  // skew correction, depends on tilt:

  double skwmid = 0.0;
  double skwwid = 0.05;
  double skwrng = 0.0; // [um]
  double skwoff = 0.0; // [um]
  double skwslp = 0.0; // [um]

  if( rot90 ) {

    if(      DUTtilt <  4.0 ) { // 14705 (1.3)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff = 0.1; // [um]
      skwslp =-429; // [um]
    }
    else if( DUTtilt <  8.0 ) { // 14703 (5.9)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff = 1.5; // [um]
      skwslp =-335; // [um]
    }
    else if( DUTtilt < 12.0 ) { // 14701 (10.0)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff = 0.1; // [um]
      skwslp =-291; // [um]
    }
    else if( DUTtilt < 14.5 ) { // 14699 (14.1)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.1; // [um]
      skwslp =-233.4; // [um]
    }
    else if( DUTtilt < 15.5 ) { // 20178 (15.1)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.7; // [um]
      skwslp =-182.2; // [um]
    }
    else if( DUTtilt < 16.5 ) { // 14697 (16.0)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.4; // [um]
      skwslp =-199.6; // [um]
    }
    else if( DUTtilt < 17.5 ) { // 20176 (17.0)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.6; // [um]
      skwslp =-157.6; // [um]
    }
    else if( DUTtilt < 18.5 ) { // 20174 (18.0)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.6; // [um]
      skwslp =-139.7; // [um]
    }
    else if( DUTtilt < 19.5 ) { // 20171 (19.0)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.4; // [um]
      skwslp =-122.6; // [um]
    }
    else if( DUTtilt < 20.5 ) { // 20169 (20.0)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.7; // [um]
      skwslp =-100.8; // [um]
    }
    else if( DUTtilt < 21.5 ) { // 20167 (21.0)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.8; // [um]
      skwslp =-86.0; // [um]
    }
    else if( DUTtilt < 22.5 ) { // 20165 (22.0)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.7; // [um]
      skwslp =-70.8; // [um]
    }
    else if( DUTtilt < 23.8 ) { // 19942 (23.4) nskw
      skwmid = 0.054;
      skwwid = 0.042;
      skwrng = 0.4; // [um]
      skwoff =-1.7; // [um]
      skwslp =-38.6; // [um/skw]
    }
    else if( DUTtilt < 24.8 ) { // 19733 (24.4)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.4; // [um]
      skwslp =-39.5; // [um/skw]
    }
    else if( DUTtilt < 25.8 ) { // 19938 (25.2) nskw
      skwmid = 0.092;
      skwwid = 0.043;
      skwrng = 0.4; // [um]
      skwoff =-0.1; // [um]
      skwslp =-6.2; // [um/skw]
    }
    else if( DUTtilt < 26.0 ) { // 20028 (25.95) nskw lin
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // [um]
      skwoff = 0.5; // [um]
      skwslp = 8.4; // [um/skw]
    }
    else if( DUTtilt < 26.3 ) { // 19731 (26.05)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.4; // [um]
      skwslp =-9.7; // [um/skw]
    }
    else if( DUTtilt < 26.8 ) { // 19947 (26.5) nskw
      skwmid = 0.078;
      skwwid = 0.2285;
      skwrng = 8.2; // [um]
      skwoff = 1.6; // [um]
      skwslp =-18.3; // [um/skw]
    }
    else if( DUTtilt < 27.3 ) { // 19975 (27.0), 19805 (27.2)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.6; // [um]
      skwslp = 21; // [um/skw]
    }
    else if( DUTtilt < 27.6 ) { // 19951 (27.4) nskw
      skwmid = 0.0985;
      skwwid = 0.0363;
      skwrng = 0.5; // [um]
      skwoff = 3.7; // [um]
      skwslp = 33.1; // [um/skw]
    }
    else if( DUTtilt < 28.1 ) { // 19729 (27.8), 19737 (28.0)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.8; // [um]
      skwslp = 24; // [um/skw]
    }
    else if( DUTtilt < 28.8 ) { // 19953 (28.3) nskw
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.9; // [um]
      skwslp = 39.3; // [um/skw]
    }
    else if( DUTtilt < 29.4 ) { // 19955 (29.2) nskw
      skwmid = 0.006;
      skwwid = 0.022;
      skwrng =-0.6; // [um]
      skwoff = 0.9; // [um]
      skwslp = 59.2; // [um/skw]
    }
    else if( DUTtilt < 29.9 ) { // 19726 (29.5)
      // .x fitskwlin.C("cmsdy0vsskwcol", -0.15, 0.15 )
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.8; // [um]
      skwslp = 57; // [um/skw]
    }
    else if( DUTtilt < 30.4 ) { // 19957 (30.1) nskw
      skwmid =-0.062;
      skwwid = 0.087;
      skwrng =-0.7; // [um]
      skwoff =-4.3; // [um]
      skwslp = 81.9; // [um/skw]
    }
    else if( DUTtilt < 31.4 ) { // 19724 (31.3)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0; // lin
      skwoff =-0.8; // [um]
      skwslp = 86; // [um/skw]
    }
    else if( DUTtilt < 32.1 ) { // 19867 (31.5)
      skwmid = 0.0;
      skwwid = 0.0117;
      skwrng = 0.8; // [um]
      skwoff =-1.3; // [um]
      skwslp = 83; // [um/skw]
    }
    else if( DUTtilt < 33.9 ) { // 19722 (33.0)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0; // lin
      skwoff =-0.8; // [um]
      skwslp = 107; // [um/skw]
    }
    else if( DUTtilt < 35.0 ) { // 19720 (34.8)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0; // lin
      skwoff =-0.6; // [um]
      skwslp = 128; // [um/skw]
    }
    else if( DUTtilt < 35.7 ) { // 19859 (35.1)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0; // lin
      skwoff =-1.0; // [um]
      skwslp = 142; // [um/skw]
    }
    else if( DUTtilt < 37.4 ) { // 19718 (36.5)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0; // lin
      skwoff =-0.8; // [um]
      skwslp = 134; // [um/skw]
    }
    else if( DUTtilt < 39.2 ) { // 19716 (38.3)
      skwmid = 0.0;
      skwwid = 0.0117;
      skwrng =-3.3; // [um]
      skwoff =-0.9; // [um]
      skwslp = 105; // [um/skw]
    }
    else if( DUTtilt < 40.05 ) { // 19714 (40.0)
      skwmid = 0.0026;
      skwwid =-0.0123;
      skwrng = 3.6; // [um]
      skwoff =-1.3; // [um]
      skwslp = 72.7; // [um/skw]
    }
    else if( DUTtilt < 40.3 ) { // 19442 (40.1) nskw
      skwmid =-0.0356;
      skwwid = 0.0116;
      skwrng = 3.2; // [um]
      skwoff =-3.0; // [um]
      skwslp = 104.3; // [um/skw]
    }
    else if( DUTtilt < 40.7 ) { // 20750 (40.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0130;
      skwrng =-3.5; // [um]
      skwoff = 1.6; // [um]
      skwslp = 65.3; // [um/skw]
    }
    else if( DUTtilt < 42.3 ) { // 19443 (42.1) skw/ncol
      skwmid =-0.0551;
      skwwid = 0.0256;
      skwrng = 2.6; // [um]
      skwoff =-3.1; // [um]
      skwslp = 50.5; // [um/skw]
    }
    else if( DUTtilt < 43.5 ) { // 20751 (42.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.05;
      skwrng = 0.0; // [um]
      skwoff = 0.0; // [um]
      skwslp = 0.0; // [um/skw]
    }
    else if( DUTtilt < 44.3 ) { // 19444 (44.1) skw/ncol
      skwmid =-0.054;
      skwwid = 0.041;
      skwrng = 4.9; // [um]
      skwoff =-0.2; // [um]
      skwslp =-60.3; // [um/skw]
    }
    else if( DUTtilt < 45.5 ) { // 20752 (44.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0085;
      skwrng = 1.4; // [um]
      skwoff = 0.2; // [um]
      skwslp =-104; // [um/skw]
    }
    else if( DUTtilt < 46.3 ) { // 19445 (46.1) skw/ncol
      skwmid =-0.0012;
      skwwid = 0.014;
      skwrng = 3.7; // [um]
      skwoff = 0.4; // [um]
      skwslp =-86.9; // [um/skw]
    }
    else if( DUTtilt < 47.5 ) { // 20753 (46.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0150;
      skwrng = 4.2; // [um]
      skwoff = 0.5; // [um]
      skwslp =-133; // [um/skw]
    }
    else if( DUTtilt < 48.3 ) { // 19446 (48.1) skw/ncol
      skwmid =-0.0008;
      skwwid = 0.018;
      skwrng = 6.4; // [um]
      skwoff = 0.5; // [um]
      skwslp =-71.2; // [um/skw]
    }
    else if( DUTtilt < 49.5 ) { // 20754 (48.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0180;
      skwrng = 6.4; // [um]
      skwoff = 0.7; // [um]
      skwslp =-75; // [um/skw]
    }
    else if( DUTtilt < 50.3 ) { // 19447 (50.1) nskw
      skwmid =-0.0019;
      skwwid = 0.0211;
      skwrng = 8.2; // [um]
      skwoff = 0.7; // [um]
      skwslp =-27.5; // [um/skw]
    }
    else if( DUTtilt < 51.5 ) { // 20755 (50.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0213;
      skwrng = 7.8; // [um]
      skwoff = 1.0; // [um]
      skwslp =-62.6; // [um/skw]
    }
    else if( DUTtilt < 52.3 ) { // 19448 (52.1) nskw
      skwmid =-0.005;
      skwwid = 0.025;
      skwrng = 8.7; // [um]
      skwoff = 0.3; // [um]
      skwslp =-62.4; // [um/skw]
    }
    else if( DUTtilt < 52.9 ) { // 20756 (52.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.043;
      skwrng = 20.6; // [um]
      skwoff = 0.9; // [um]
      skwslp =-332; // [um/skw]
    }
    else if( DUTtilt < 53.6 ) { // 19809 (53.1)
      skwmid = 0.0;
      skwwid = 0.112;
      skwrng = 453; // [um]
      skwoff =-1.3; // [um]
      skwslp =-3878; // [um/skw]
    }
    else if( DUTtilt < 54.05 ) { // 19845 (54.0)
      skwmid = 0.0;
      skwwid = 0.138;
      skwrng = 478; // [um]
      skwoff =-1.9; // [um]
      skwslp =-3426; // [um/skw]
    }
    else if( DUTtilt < 54.3 ) { // 19449 (54.1) nskw
      skwmid =-0.0148;
      skwwid = 0.0229;
      skwrng = 5.2; // [um]
      skwoff = 0.3; // [um]
      skwslp =-70.1; // [um/skw]
    }
    else if( DUTtilt < 54.7 ) { // 20757 (54.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.602;
      skwrng = 9910; // [um]
      skwoff = 0.2; // [um]
      skwslp =-16544; // [um/skw]
    }
    else if( DUTtilt < 55.2 ) { // 19811 (54.8) nskw
      skwmid = 0.0;
      skwwid = 0.163;
      skwrng = 329; // [um]
      skwoff =-1.9; // [um]
      skwslp =-2095; // [um/skw]
    }
    else if( DUTtilt < 56.3 ) { // 19847 (55.8)
      skwmid = 0.0;
      skwwid = 0.00985;
      skwrng =-1.1; // [um]
      skwoff =-1.8; // [um]
      skwslp =-197; // [um/skw]
    }
    else if( DUTtilt < 56.6 ) { // 20758 (56.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.01965;
      skwrng =-1.85; // [um]
      skwoff =-0.35; // [um]
      skwslp =-239.7; // [um/skw]
    }
    else if( DUTtilt < 57.2 ) { // 19813 (56.7) nskw
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff = 1.9; // [um]
      skwslp =-297; // [um/skw]
    }
    else if( DUTtilt < 58.1 ) { // 19850 (57.6)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0; // lin
      skwoff =-1.9; // [um]
      skwslp =-334; // [um/skw]
    }
    else if( DUTtilt < 59.0 ) { // 19815 (58.5)
      skwmid = 0.0;
      skwwid = 0.0131;
      skwrng = 3.8; // [um]
      skwoff =-1.9; // [um]
      skwslp =-334; // [um/skw]
    }
    else if( DUTtilt < 59.45 ) { // 19852 (59.4)
      skwmid = 0.0;
      skwwid = 0.0145;
      skwrng = 5.1; // [um]
      skwoff = 16.6; // [um]
      skwslp =-278; // [um/skw]
    }
    else if( DUTtilt < 59.9 ) { // 20779 (59.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0147;
      skwrng = 4.1; // [um]
      skwoff = 19.1; // [um]
      skwslp =-288; // [um/skw]
    }
    else if( DUTtilt < 60.4 ) { // 19817 (60.3)
      skwmid = 0.0;
      skwwid = 0.0154;
      skwrng = 5.7; // [um]
      skwoff =-1.8; // [um]
      skwslp =-228; // [um/skw]
    }
    else if( DUTtilt < 60.9 ) { // 20760 (60.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0173;
      skwrng = 4.8; // [um]
      skwoff = 0.4; // [um]
      skwslp =-271; // [um/skw]
    }
    else if( DUTtilt < 61.3 ) { // 19854 (61.2)
      skwmid = 0.0;
      skwwid = 0.0254;
      skwrng = 9.2; // [um]
      skwoff =-0.5; // [um]
      skwslp =-399; // [um/skw]
    }
    else if( DUTtilt < 61.7 ) { // 20778 (61.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.027;
      skwrng = 4.0; // [um]
      skwoff = 0.6; // [um]
      skwslp =-320; // [um/skw]
    }
    else if( DUTtilt < 62.2 ) { // 19819 (62.0)
      skwmid = 0.0;
      skwwid = 0.0324;
      skwrng = 8.1; // [um]
      skwoff =-1.8; // [um]
      skwslp =-412; // [um/skw]
    }
    else if( DUTtilt < 62.8 ) { // 20761 (62.5) nskw lin
      skwmid = 0.0; // fix
      skwwid = 0.02;
      skwrng = 0.0; // [um]
      skwoff = 0.1; // [um]
      skwslp =-310; // [um/skw]
    }
    else if( DUTtilt < 63.2 ) { // 19856 (63.0)
      skwmid = 0.0;
      skwwid =-0.0076;
      skwrng = 1.0; // [um]
      skwoff =-1.1; // [um]
      skwslp =-371; // [um/skw]
    }
    else if( DUTtilt < 63.7 ) { // 20777 (63.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.028;
      skwrng =-2.9; // [um]
      skwoff = 1.0; // [um]
      skwslp =-385; // [um/skw]
    }
    else if( DUTtilt < 64.2 ) { // 19821 (63.9)
      skwmid = 0.0;
      skwwid = 0.038;
      skwrng =-7.9; // [um]
      skwoff =-1.4; // [um]
      skwslp =-321; // [um/skw]
    }
    else if( DUTtilt < 64.8 ) { // 20762 (64.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.047;
      skwrng =-7.5; // [um]
      skwoff = 0.4; // [um]
      skwslp =-374; // [um/skw]
    }
    else if( DUTtilt < 65.55 ) { // 20776 (65.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.153;
      skwrng = 95.8; // [um]
      skwoff = 0.5; // [um]
      skwslp =-1074; // [um/skw]
    }
    else if( DUTtilt < 66.1 ) { // 19823 (65.6)
      skwmid = 0.0;
      skwwid = 0.218;
      skwrng = 3382; // [um]
      skwoff =-1.5; // [um]
      skwslp =-15887; // [um/skw]
    }
    else if( DUTtilt < 66.9 ) { // 20763 (66.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.019;
      skwrng = 1.7; // [um]
      skwoff = 0.2; // [um]
      skwslp =-464; // [um/skw]
    }
    else if( DUTtilt < 68.0 ) { // 19826 (67.5)
      skwmid = 0.0;
      skwwid = 0.213;
      skwrng = 755; // [um]
      skwoff =-1.7; // [um]
      skwslp =-4005; // [um/skw]
    }
    else if( DUTtilt < 68.9 ) { // 20764 (68.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0266;
      skwrng =-2.6; // [um]
      skwoff =-0.5; // [um]
      skwslp =-532; // [um/skw]
    }
    else if( DUTtilt < 69.7 ) { // 19828 (69.2)
      skwmid = 0.0;
      skwwid = 0.129;
      skwrng =-310; // [um]
      skwoff =-1.3; // [um]
      skwslp = 1622; // [um/skw]
    }
    else if( DUTtilt < 70.9 ) { // 20765 (70.5) nskw
      skwmid = 0.0; // fix
      skwwid =-0.0104;
      skwrng =-1.2; // [um]
      skwoff = 0.6; // [um]
      skwslp =-610; // [um/skw]
    }
    else if( DUTtilt < 71.6 ) { // 19830 (71.1)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-1.2; // [um]
      skwslp =-679; // [um/skw]
    }
    else if( DUTtilt < 72.7 ) { // 20766 (72.5) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0325;
      skwrng =-6.1; // [um]
      skwoff = 0.4; // [um]
      skwslp =-658; // [um/skw]
    }
    else if( DUTtilt < 73.3 ) { // 19832 (72.8)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-1.2; // [um]
      skwslp =-888; // [um/skw]
    }
    else if( DUTtilt < 74.6 ) { // 20767 (74.4) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0285;
      skwrng =-5.9; // [um]
      skwoff = 0.6; // [um]
      skwslp =-720; // [um/skw]
    }
    else if( DUTtilt < 75.2 ) { // 19834 (74.7)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.9; // [um]
      skwslp =-1007; // [um/skw]
    }
    else if( DUTtilt < 76.7 ) { // 19836 (76.4), 20768 (76.4)
      skwmid = 0.0; // fix
      skwwid = 0.06;
      skwrng = 0.0; // lin
      skwoff = 0.2; // [um]
      skwslp =-1156; // [um/skw]
    }
    else if( DUTtilt < 77.1 ) { // 20160 (76.9) nskw
      skwmid = 0.0; // fix
      skwwid = 0.163;
      skwrng =-253.5; // [um]
      skwoff =-6.8; // [um]
      skwslp = 735; // [um/skw]
    }
    else if( DUTtilt < 78.0 ) { // 20150 (77.9) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0143;
      skwrng =-2.7; // [um]
      skwoff =-2.9; // [um]
      skwslp =-874; // [um/skw]
    }
    else if( DUTtilt < 78.3 ) { // 19838 (78.2)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff = 0.4; // [um]
      skwslp =-1370; // [um/skw]
    }
    else if( DUTtilt < 78.7 ) { // 20769 (78.4) nskw
      skwmid = 0.0; // fix
      skwwid =-0.0804;
      skwrng = 61.1; // [um]
      skwoff = 2.2; // [um]
      skwslp =-347; // [um/skw]
    }
    else if( DUTtilt < 79.2 ) { // 20158 (78.9) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0179;
      skwrng =-2.2; // [um]
      skwoff =-6.2; // [um] misaligned
      skwslp =-877; // [um/skw]
    }
    else if( DUTtilt < 79.7 ) { // 20770 (79.4) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0059;
      skwrng =-0.5; // [um]
      skwoff = 0.0; // [um]
      skwslp =-1069; // [um/skw]
    }
    else if( DUTtilt < 79.94 ) { // 20151 (79.9) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0162;
      skwrng =-2.5; // [um]
      skwoff =-0.7; // [um]
      skwslp =-993; // [um/skw]
    }
    else if( DUTtilt < 80.2 ) { // 19840 (80.0)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff = 0.6; // [um]
      skwslp =-1607; // [um/skw]
    }
    else if( DUTtilt < 80.7 ) { // 20771 (80.3) nskw
      skwmid = 0.0; // fix
      skwwid =-0.0278;
      skwrng =-3.9; // [um]
      skwoff = 2.8; // [um]
      skwslp =-1130; // [um/skw]
    }
    else if( DUTtilt < 81.4 ) { // 20157 (80.9) nskw
      skwmid = 0.0; // fix
      skwwid = 0.1204;
      skwrng =-400; // [um]
      skwoff =-0.4; // [um]
      skwslp = 2241; // [um/skw]
    }
    else if( DUTtilt < 81.6 ) { // 20772 (81.3) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0194;
      skwrng =-6.2; // [um]
      skwoff = 0.9; // [um]
      skwslp =-1195; // [um/skw]
    }
    else if( DUTtilt < 81.825 ) { // 19843 (81.81)
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // lin
      skwoff =-0.8; // [um]
      skwslp =-1903; // [um/skw]
    }
    else if( DUTtilt < 82.1 ) { // 20152 (81.84) nskw
      skwmid = 0.0; // fix
      skwwid = 0.1075;
      skwrng =-433; // [um]
      skwoff = 1.0; // [um]
      skwslp = 2704; // [um/skw]
    }
    else if( DUTtilt < 82.6 ) { // 20773 (82.2) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0520;
      skwrng =-39; // [um]
      skwoff = 3.1; // [um]
      skwslp =-818; // [um/skw]
    }
    else if( DUTtilt < 83.3 ) { // 20155 (82.9) nskw
      skwmid = 0.0; // fix
      skwwid = 0.024;
      skwrng =-15.0; // [um]
      skwoff = 1.0; // [um]
      skwslp =-977; // [um/skw]
    }
    else if( DUTtilt < 83.6 ) { // 20774 (83.2) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0197;
      skwrng =-16; // [um]
      skwoff = 1.6; // [um]
      skwslp =-1260; // [um/skw]
    }
    else if( DUTtilt < 84.3 ) { // 20154 (83.8) nskw
      skwmid = 0.0; // fix
      skwwid = 0.016;
      skwrng =-5.2; // [um]
      skwoff =-0.6; // [um]
      skwslp =-1312; // [um/skw]
    }
    else if( DUTtilt < 84.6 ) { // 20775 (84.2) nskw
      skwmid = 0.0; // fix
      skwwid = 0.076;
      skwrng =-180; // [um]
      skwoff = 2.8; // [um]
      skwslp = 417; // [um/skw]
    }
    else
      skwrng =   0.0;
  } // rot90
  else {
    if( DUTtilt <  3.5 ) { // 20245 (3.0) nskw
      skwmid = 0.0; // fix
      skwwid = 0.1088;
      skwrng =-296; // [um]
      skwoff =-1.9; // [um]
      skwslp = 2249; // [um]
    }
    else if( DUTtilt <  5.5 ) { // 20243 (5.0) nskw
      skwmid = 0.0; // fix
      skwwid = 0.136;
      skwrng =-267; // [um]
      skwoff =-1.5; // [um]
      skwslp = 1613; // [um]
    }
    else if( DUTtilt <  7.5 ) { // 20241 (7.0) nskw
      skwmid = 0.0; // fix
      skwwid = 0.1346;
      skwrng =-238; // [um]
      skwoff =-0.9; // [um]
      skwslp = 1467; // [um]
    }
    else if( DUTtilt <  9.5 ) { // 20239 (9.0) nskw
      skwmid = 0.0; // fix
      skwwid = 0.1408;
      skwrng =-170; // [um]
      skwoff =-0.7; // [um]
      skwslp = 967.5; // [um]
    }
    else if( DUTtilt < 11.5 ) { // 20237 (11.0) nskw
      skwmid = 0.0; // fix
      skwwid = 0.1524;
      skwrng =-131; // [um]
      skwoff =-0.3; // [um]
      skwslp = 669; // [um]
    }
    else if( DUTtilt < 13.5 ) { // 20235 (13.0) nskw
      skwmid = 0.0; // fix
      skwwid = 0.1509;
      skwrng =-97.1; // [um]
      skwoff =-0.3; // [um]
      skwslp = 498; // [um]
    }
    else if( DUTtilt < 14.5 ) { // 20230 (14.0) nskw
      skwmid = 0.0; // fix
      skwwid = 0.1432;
      skwrng =-103; // [um]
      skwoff =-0.1; // [um]
      skwslp = 584; // [um]
    }
    else if( DUTtilt < 15.5 ) { // 20225 (15.0) nskw
      skwmid = 0.0; //
      skwwid = 0.1576;
      skwrng =-114; // linear
      skwoff =-0.2; // [um]
      skwslp = 615; // [um]
    }
    else if( DUTtilt < 16.5 ) { // 20223 (16.0) nskw
      skwmid = 0.0; //
      skwwid = 0.2195;
      skwrng =-194; // [um]
      skwoff =-0.1; // [um]
      skwslp = 800; // [um]
    }
    else if( DUTtilt < 17.5 ) { // 20221 (17.0) nskw
      skwmid = 0.0; //
      skwwid = 0.2551;
      skwrng =-328; // [um]
      skwoff =-0.1; // [um]
      skwslp = 1219; // [um]
    }
    else if( DUTtilt < 18.5 ) { // 20219 (18.0) nskw
      skwmid = 0.0; //
      skwwid = 0.384;
      skwrng =-870; // [um]
      skwoff =-0.1; // [um]
      skwslp = 2225; // [um]
    }
    else if( DUTtilt < 19.2 ) { // 20182 (19.0) nskw lin
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // linear
      skwoff = 0.0; // [um]
      skwslp =-5.5; // [um]
    }
    else if( DUTtilt < 19.7 ) { // 20194 (19.5) nskw lin
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // linear
      skwoff = 0.0; // [um]
      skwslp = 1.6; // [um]
    }
    else if( DUTtilt < 20.5 ) { // 20184 (20.0) nskw lin
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // linear
      skwoff = 0.1; // [um]
      skwslp = 10.9; // [um]
    }
    else if( DUTtilt < 21.5 ) { // 20186 (21.0) nskw
      skwmid = 0.0;
      skwwid = 0.1701;
      skwrng = 32.7; // [um]
      skwoff = 0.1; // [um]
      skwslp =-157; // [um]
    }
    else if( DUTtilt < 22.5 ) { // 20188 (22.0) nskw lin
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // linear
      skwoff = 0.0; // [um]
      skwslp = 38.3; // [um]
    }
    else if( DUTtilt < 23.5 ) { // 20190 (23.0) nskw lin
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // linear
      skwoff = 0.1; // [um]
      skwslp = 47.4; // [um]
    }
    else if( DUTtilt < 24.5 ) { // 20192 (24.0) nskw lin
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // linear
      skwoff = 1.6; // [um]
      skwslp = 53.7; // [um]
    }
    else if( DUTtilt < 26.6 ) { // 20202 (26.0) nskw
      skwmid = 0.0;
      skwwid = 0.1463;
      skwrng = 60.6; // [um]
      skwoff =-0.5; // [um]
      skwslp =-334; // [um]
    }
    else if( DUTtilt < 27.0 ) { // 20785 (26.6) nskw lin
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // [um]
      skwoff =-1.0; // [um]
      skwslp = 76.5; // [um]
    }
    else if( DUTtilt < 28.0 ) { // 20799 (27.6) nskw lin
      skwmid = 0.0;
      skwwid = 0.056;
      skwrng = 0; // [um]
      skwoff =-0.3; // [um]
      skwslp = 62; // [um]
    }
    else if( DUTtilt < 28.3 ) { // 20204 (28.0) nskw
      skwmid = 0.0;
      skwwid = 0.2044;
      skwrng = 196; // 
      skwoff =-0.4; // [um]
      skwslp =-898; // [um]
    }
    else if( DUTtilt < 29.0 ) { // 20786 (28.6) nskw
      skwmid = 0.0;
      skwwid =-0.01176;
      skwrng = 2.7; // [um]
      skwoff =-0.2; // [um]
      skwslp = 46.1; // [um]
    }
    else if( DUTtilt < 30.0 ) { // 20799 (29.6) nskw lin
      skwmid = 0.0;
      skwwid = 0.056;
      skwrng = 0; // [um]
      skwoff =-0.2; // [um]
      skwslp = 21; // [um]
    }
    else if( DUTtilt < 30.3 ) { // 20206 (30.0) nskw lin
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // linear
      skwoff =-0.4; // [um]
      skwslp =-5.4; // [um]
    }
    else if( DUTtilt < 31.0 ) { // 20787 (30.6) nskw lin
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0.0; // linear
      skwoff = 0.4; // [um]
      skwslp =-12.8; // [um]
    }
    else if( DUTtilt < 32.0 ) { // 20797 (31.6) nskw lin
      skwmid = 0.0;
      skwwid = 0.05;
      skwrng = 0; // [um]
      skwoff =-0.3; // [um]
      skwslp =-59; // [um]
    }
    else if( DUTtilt < 32.3 ) { // 20208 (32.0) nskw
      skwmid = 0.0; // fix
      skwwid = 0.01908;
      skwrng =-1.9; // [um]
      skwoff =-0.4; // [um]
      skwslp =-79.6; // [um]
    }
    else if( DUTtilt < 33.0 ) { // 20788 (32.6) nskw
      skwmid = 0.0;
      skwwid = 0.195;
      skwrng =-188.5; // [um]
      skwoff =-0.3; // [um]
      skwslp = 865; // [um]
    }
    else if( DUTtilt < 34.0 ) { // 20796 (33.6) nskw
      skwmid = 0.0;
      skwwid = 0.011;
      skwrng = 1.9; // [um]
      skwoff =-0.2; // [um]
      skwslp =-120; // [um]
    }
    else if( DUTtilt < 34.3 ) { // 20210 (34.0) nskw
      skwmid = 0.0;
      skwwid = 0.0110;
      skwrng = 1.5; // [um]
      skwoff =-0.3; // [um]
      skwslp =-150; // [um]
    }
    else if( DUTtilt < 35.0 ) { // 20789 (34.6) nskw
      skwmid = 0.0;
      skwwid = 0.013;
      skwrng = 2.7; // [um]
      skwoff =-0.2; // [um]
      skwslp =-124; // [um]
    }
    else if( DUTtilt < 36.0 ) { // 20795 (35.6) nskw
      skwmid = 0.0;
      skwwid = 0.015;
      skwrng = 3.6; // [um]
      skwoff =-0.4; // [um]
      skwslp =-99; // [um]
    }
    else if( DUTtilt < 36.3 ) { // 20212 (36.0) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0148;
      skwrng = 3.3; // [um]
      skwoff =-0.4; // [um]
      skwslp =-120; // [um]
    }
    else if( DUTtilt < 37.0 ) { // 20790 (36.6) nskw
      skwmid = 0.0;
      skwwid = 0.0159;
      skwrng = 4.1; // [um]
      skwoff =-0.4; // [um]
      skwslp =-67.2; // [um]
    }
    else if( DUTtilt < 38.0 ) { // 20794 (37.6) nskw
      skwmid = 0.0;
      skwwid = 0.018;
      skwrng = 5.1; // [um]
      skwoff =-0.3; // [um]
      skwslp =-50; // [um]
    }
    else if( DUTtilt < 38.3 ) { // 20214 (38.0) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0182;
      skwrng = 4.6; // [um]
      skwoff =-0.4; // [um]
      skwslp =-77.3; // [um]
    }
    else if( DUTtilt < 39.0 ) { // 20791 (38.6) nskw
      skwmid = 0.0;
      skwwid = 0.0194;
      skwrng = 5.2; // [um]
      skwoff =-0.4; // [um]
      skwslp =-45.3; // [um]
    }
    else if( DUTtilt < 40.0 ) { // 20793 (39.6) nskw
      skwmid = 0.0;
      skwwid = 0.025;
      skwrng = 6.3; // [um]
      skwoff =-0.6; // [um]
      skwslp =-96; // [um]
    }
    else if( DUTtilt < 40.3 ) { // 20216 (40.0) nskw
      skwmid = 0.0; // fix
      skwwid = 0.0246;
      skwrng = 5.4; // [um]
      skwoff =-0.6; // [um]
      skwslp =-116; // [um]
    }
    else if( DUTtilt < 41.0 ) { // 20792 (40.6) nskw
      skwmid = 0.0;
      skwwid = 0.053;
      skwrng = 31; // [um]
      skwoff =-0.5; // [um]
      skwslp =-493; // [um]
    }

  } // not rot90

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // REF gain:

  double r0[52][80]; // Weibull
  double r1[52][80];
  double r2[52][80];
  double r3[52][80];
  double r4[52][80];
  double r5[52][80];

  ifstream refgainFile( refgainFileName );

  if(! refgainFile ) {
    cout << "REF gain file " << refgainFileName << " not found" << endl;
    return 1;
  }
  else {
    cout << endl << "using REF pixel gain file " << refgainFileName << endl;
    char ih[99];
    int col;
    int row;
    while( refgainFile >> ih ) {
      refgainFile >> col;
      refgainFile >> row;
      if( col < 0 || col > 51 || row < 0 || row > 79 ) {
	cout << "invalid pixel in gain file " << col << " " << row << endl;
	continue;
      }
      refgainFile >> r0[col][row];
      refgainFile >> r1[col][row];
      refgainFile >> r2[col][row];
      refgainFile >> r3[col][row];
      refgainFile >> r4[col][row];
      refgainFile >> r5[col][row]; // gain ratio
    } // while

  } // gainFile

  double refke = 0.250; // to get q0f peak at 22 ke

  if(      refchip0 == 500 )
    refke = 0.336;

  else if( refchip0 == 501 )
    refke = 0.305;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "scope" << run << ".root";

  TFile* histoFile = new TFile( (outputDirectory+"/"+rootFileName.str(  )).c_str(  ), "RECREATE" );

  // book histos:
 
  TH1I hdtus = TH1I( "dtus", "time between events;time between events [us];events", 100, 0, 1000 );
  TH1I hdtms = TH1I( "dtms", "time between events;time between events [ms];events", 100, 0, 1000 );

  TH1I t1Histo = TH1I( "t1", "event time;event time [s];events", 100, 0, 1 );
  TH1I t2Histo = TH1I( "t2", "event time;event time [s];events", 300, 0, 300 );
  TH1I t3Histo = TH1I( "t3", "event time;event time [s];events", 150, 0, 1500 );
  TH1I t4Histo = TH1I( "t4", "event time;event time [s];events", 100, 0, 10000 );
  TH1I t5Histo = TH1I( "t5", "event time;event time [s];events", 600, 0, 60000 );

  TH1I hcol[9];
  TH1I hrow[9];
  TH1I hnpx[9];
  TH2I * hmap[9];

  TH1I hncl[9];
  TH1I hsiz[9];
  TH1I hncol[9];
  TH1I hnrow[9];

  for( int ipl = 0; ipl < 9; ++ipl ) {

    hcol[ipl] = TH1I( Form( "col%i", ipl ),
		      Form( "%i col;col;%i pixels", ipl, ipl ), 
		      max( 52, nx[ipl]/4 ), 0, nx[ipl] );
    hrow[ipl] = TH1I( Form( "row%i", ipl ),
		      Form( "%i row;row;%i pixels", ipl, ipl ),
		      max( 80, ny[ipl]/2 ), 0, ny[ipl] );
    hmap[ipl] = new TH2I( Form( "map%i", ipl ),
			  Form( "%i map;col;row;%i pixels", ipl, ipl ),
			  max( 52, nx[ipl]/4 ), 0, nx[ipl], max( 80, ny[ipl]/2 ), 0, ny[ipl] );

    hnpx[ipl] = TH1I( Form( "npx%i", ipl ),
		      Form( "%i pixel per event;pixels;%i events", ipl, ipl ),
		      200, 0, 200 );

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

  } // planes

  TProfile dutnpxvst2 =
    TProfile( "dutnpxvst2",
	      "DUT pixels vs time;time [s];DUT pixels per event",
	      150, 0, 1500, -0.5, 99.5 );
  TProfile dutnclvst2 =
    TProfile( "dutnclvst2",
	      "DUT clusters vs time;time [s];DUT clusters per event with pixels",
	      150, 0, 1500, 0.5, 99.5 );
  TProfile dutyldvst2 =
    TProfile( "dutyldvst2",
	      "DUT yield vs time;time [s];DUT events with pixels",
	      150, 0, 1500, -0.5, 1.5 );

  // driplets:

  TH1I hdx35 = TH1I( "dx35", "3-5 dx;3-5 dx [mm];cluster pairs", 100, -1, 1 );
  TH1I hdy35 = TH1I( "dy35", "3-5 dy;3-5 dy [mm];cluster pairs", 100, -1, 1 );

  TH1I hdridx = TH1I( "dridx", "driplet dx;driplet dx [mm];driplets", 100, -0.1, 0.1 );
  TH1I hdridy = TH1I( "dridy", "driplet dy;driplet dy [mm];driplets", 100, -0.1, 0.1 );

  TH1I hdridxc = TH1I( "dridxc", "driplet dx;driplet dx [mm];driplets", 100, -0.1, 0.1 );
  TH1I hdridyc = TH1I( "dridyc", "driplet dy;driplet dy [mm];driplets", 100, -0.1, 0.1 );

  TProfile dridxvsy =
    TProfile( "dridxvsy",
	      "driplet dx vs y;driplet yB [mm];<driplets #Deltax> [mm]",
	      110, -5.5, 5.5, -0.05, 0.05 );
  TProfile dridyvsx =
    TProfile( "dridyvsx",
	      "driplet dy vs x;driplet xB [mm];<driplets #Deltay> [mm]",
	      110, -11, 11, -0.05, 0.05 );

  TProfile dridxvstx =
    TProfile( "dridxstx",
	      "driplet dx vs slope x;driplet slope x [rad];<driplets #Deltax> [mm]",
	      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile dridyvsty =
    TProfile( "dridysty",
	      "driplet dy vs slope y;driplet slope y [rad];<driplets #Deltay> [mm]",
	      60, -0.003, 0.003, -0.05, 0.05 );

  TH1I ndriHisto = TH1I( "ndri", "driplets;driplets;events", 51, -0.5, 50.5 );

  TH1I dddmin1Histo = TH1I( "dddmin1",
			    "telescope driplets isolation at ref;driplets min #Delta_{xy} [mm];driplet pairs",
			    100, 0, 1 );
  TH1I dddmin2Histo = TH1I( "dddmin2",
			    "telescope driplets isolation at ref;driplets min #Delta_{xy} [mm];driplet pairs",
			    150, 0, 15 );

  TH1I dddmin1modHisto = TH1I( "dddmin1mod",
			    "telescope driplets isolation at mod;driplets min #Delta_{xy} [mm];driplet pairs",
			    100, 0, 1 );
  TH1I dddmin2modHisto = TH1I( "dddmin2mod",
			    "telescope driplets isolation at mod;driplets min #Delta_{xy} [mm];driplet pairs",
			    150, 0, 15 );

  // REF vs driplets:

  TH1I refsxaHisto = TH1I( "refsxa",
			   "REF + driplet x;REF cluster + driplet #Sigmax [mm];REF clusters",
			   440, -11, 11 );
  TH1I refdxaHisto = TH1I( "refdxa",
			   "REF - driplet x;REF cluster - driplet #Deltax [mm];REF clusters",
			   440, -11, 11 );

  TH1I refsyaHisto = TH1I( "refsya",
			   "REF + driplet y;REF cluster + driplet #Sigmay [mm];REF clusters",
			   220, -5.5, 5.5 );
  TH1I refdyaHisto = TH1I( "refdya",
			   "REF - driplet y;REF cluster - driplet #Deltay [mm];REF clusters",
			   220, -5.5, 5.5 );

  TH1I refdxHisto = TH1I( "refdx",
			   "REF - driplet x;REF cluster - driplet #Deltax [mm];REF clusters",
			   200, -1, 1 );
  TH1I refdyHisto = TH1I( "refdy",
			   "REF - driplet y;REF cluster - driplet #Deltay [mm];REF clusters",
			   200, -1, 1 );

  TProfile refdyvsx =
    TProfile( "refdyvsx",
	      "REF #Deltay vs x;x track [mm];<cluster - driplet #Deltay> [mm]",
	      50, -3.75, 3.75, -0.5, 0.5 );
  TProfile refdyvsy =
    TProfile( "refdyvsy",
	      "REF #Deltay vs y;y track [mm];<cluster - driplet #Deltay> [mm]",
	      76, -3.8, 3.8, -0.5, 0.5 );
  TProfile refdyvsty =
    TProfile( "refdyvsty",
	      "REF #Deltay vs #theta_{y};y track slope [rad];<cluster - driplet #Deltay> [mm]",
	      80, -0.002, 0.002, -0.5, 0.5 );

  TH1I refnpxHisto =
    TH1I( "refnpx",
	  "REF linked clusters;REF cluster size [pixels];linked REF cluster",
	  20, 0.5, 20.5 );

  TH1I refqHisto = TH1I( "refq",
			 "REF linked clusters;REF cluster charge [ke];linked REF cluster",
			 80, 0, 80 );

  TProfile2D * refnpxvsxmym = new
    TProfile2D( "refnpxvsxmym",
		"REF cluster size vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];REF <cluster size> [pixels]",
		60, 0, 0.3, 40, 0, 0.2, 0, 20 );

  TProfile2D * refqxvsxmym = new
    TProfile2D( "refqxvsxmym",
	      "REF cluster charge vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];REF <cluster charge> [ke]",
		60, 0, 0.3, 40, 0, 0.2, 0, 0.1 );

  TH1I reflkxBHisto = TH1I( "reflkxb",
			   "linked driplet at REF x;driplet x at REF [mm];linked driplets",
			   220, -11, 11 );
  TH1I reflkyBHisto = TH1I( "reflkyb",
			   "linked driplet at REF y;driplet y at REF [mm];linked driplets",
			   120, -6, 6 );
  TH1I reflkxHisto = TH1I( "reflkx",
			   "linked driplet at REF x;driplet x at REF [mm];linked driplets",
			   60, -4.5, 4.5 );
  TH1I reflkyHisto = TH1I( "reflky",
			   "linked driplet at REF y;driplet y at REF [mm];linked driplets",
			   90, -4.5, 4.5 );

  TH1I reflkcolHisto = TH1I( "reflkcol",
			     "REF linked col;REF linked col;linked REF cluster",
			     52, 0, 52 );
  TH1I reflkrowHisto = TH1I( "reflkrow",
			     "REF linked row;REF linked row;linked REF cluster",
			     80, 0, 80 );
  TProfile reflkvst1 =
    TProfile( "reflkvst1",
	      "driplet-REF links vs time;time [s];driplets with REF links",
	      300, 0, 300, -0.5, 1.5 );
  TProfile reflkvst2 =
    TProfile( "reflkvst2",
	      "driplet-REF links vs time;time [s];driplets with REF links",
	      150, 0, 1500, -0.5, 1.5 );
  TProfile reflkvst3 =
    TProfile( "reflkvst3",
	      "driplet-REF links vs time;time [s];driplets with REF links",
	      300, 0, 30000, -0.5, 1.5 );

  TH1I ndrilkHisto = TH1I( "ndrilk", "driplet - REF links;driplet - REF links;events",
			    11, -0.5, 10.5 );

  // MOD vs driplets:

  TH1I modsxaHisto = TH1I( "modsxa",
			   "MOD + driplet x;MOD cluster + driplet #Sigmax [mm];MOD clusters",
			   320, -32, 32 );
  TH1I moddxaHisto = TH1I( "moddxa",
			   "MOD - driplet x;MOD cluster - driplet #Deltax [mm];MOD clusters",
			   320, -32, 32 );

  TH1I modsyaHisto = TH1I( "modsya",
			   "MOD + driplet y;MOD cluster + driplet #Sigmay [mm];MOD clusters",
			   110, -5.5, 5.5 );
  TH1I moddyaHisto = TH1I( "moddya",
			   "MOD - driplet y;MOD cluster - driplet #Deltay [mm];MOD clusters",
			   110, -5.5, 5.5 );

  TH1I moddxHisto = TH1I( "moddx",
			   "MOD - driplet x;MOD cluster - driplet #Deltax [mm];MOD clusters",
			   200, -1, 1 );
  TProfile moddxvsx =
    TProfile( "moddxvsx",
	      "MOD #Deltax vs x;x track [mm];<cluster - driplet #Deltax> [mm]",
	      432, -32.4, 32.4, -0.5, 0.5 );
  TProfile moddxvsy =
    TProfile( "moddxvsy",
	      "MOD #Deltax vs y;y track [mm];<cluster - driplet #Deltax> [mm]",
	      76, -3.8, 3.8, -0.5, 0.5 );
  TProfile moddxvstx =
    TProfile( "moddxvstx",
	      "MOD #Deltax vs #theta_{x};x track slope [rad];<cluster - driplet #Deltax> [mm]",
	      80, -0.002, 0.002, -0.5, 0.5 );

  TH1I moddyHisto = TH1I( "moddy",
			   "MOD - driplet y;MOD cluster - driplet #Deltay [mm];MOD clusters",
			   200, -1, 1 );
  TProfile moddyvsx =
    TProfile( "moddyvsx",
	      "MOD #Deltay vs x;x track [mm];<cluster - driplet #Deltay> [mm]",
	      432, -32.4, 32.4, -0.5, 0.5 );
  TProfile moddyvsy =
    TProfile( "moddyvsy",
	      "MOD #Deltay vs y;y track [mm];<cluster - driplet #Deltay> [mm]",
	      76, -3.8, 3.8, -0.5, 0.5 );
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

  TProfile2D * modnpxvsxmym = new
    TProfile2D( "modnpxvsxmym",
		"MOD cluster size vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];MOD <cluster size> [pixels]",
		60, 0, 0.3, 40, 0, 0.2, 0, 20 );

  TProfile2D * modqxvsxmym = new
    TProfile2D( "modqxvsxmym",
	      "MOD cluster charge vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];MOD <cluster charge> [ke]",
		60, 0, 0.3, 40, 0, 0.2, 0, 0.1 );

  TH1I modlkxBHisto = TH1I( "modlkxb",
			   "linked driplet at MOD x;driplet x at MOD [mm];linked driplets",
			   220, -11, 11 );
  TH1I modlkyBHisto = TH1I( "modlkyb",
			   "linked driplet at MOD y;driplet y at MOD [mm];linked driplets",
			   120, -6, 6 );
  TH1I modlkxHisto = TH1I( "modlkx",
			   "linked driplet at MOD x;driplet x at MOD [mm];linked driplets",
			   480, -4.5, 36 );
  TH1I modlkyHisto = TH1I( "modlky",
			   "linked driplet at MOD y;driplet y at MOD [mm];linked driplets",
			   180, -9, 9 );
  TH2I modlkxvsyHisto = TH2I( "modlkxvsy",
			      "linked driplet at MOD x vs y;driplet x at MOD [mm];driplet y at MOD [mm]",
			      480, -4.5, 36 , 180, -9, 9);

  TH1I modlkcolHisto = TH1I( "modlkcol",
			     "MOD linked col;MOD linked col;linked MOD cluster",
			     52, 0, 52 );
  TH1I modlkrowHisto = TH1I( "modlkrow",
			     "MOD linked row;MOD linked row;linked MOD cluster",
			     80, 0, 80 );

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
	      300, 0, 30000, -0.5, 1.5 );

  // DUT:

  TH1I dutpxq1stHisto =
    TH1I( "dutpxq1st",
	  "DUT pixel charge 1st;1st pixel charge [ke];1st pixels",
	  100, 0, 25 );

  TH1I dutpxq2ndHisto =
    TH1I( "dutpxq2nd",
	  "DUT pixel charge 2nd;2nd pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd02Histo =
    TH1I( "dutpxq2nd02",
	  "DUT pixel charge 2nd - 0.02 q1;2nd - 0.02 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd03Histo =
    TH1I( "dutpxq2nd03",
	  "DUT pixel charge 2nd - 0.03 q1;2nd - 0.03 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd04Histo =
    TH1I( "dutpxq2nd04",
	  "DUT pixel charge 2nd - 0.04 q1;2nd - 0.04 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd05Histo =
    TH1I( "dutpxq2nd05",
	  "DUT pixel charge 2nd - 0.05 q1;2nd - 0.05 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd06Histo =
    TH1I( "dutpxq2nd06",
	  "DUT pixel charge 2nd - 0.06 q1;2nd - 0.06 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd07Histo =
    TH1I( "dutpxq2nd07",
	  "DUT pixel charge 2nd - 0.07 q1;2nd - 0.07 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd08Histo =
    TH1I( "dutpxq2nd08",
	  "DUT pixel charge 2nd - 0.08 q1;2nd - 0.08 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd09Histo =
    TH1I( "dutpxq2nd09",
	  "DUT pixel charge 2nd - 0.09 q1;2nd - 0.09 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd10Histo =
    TH1I( "dutpxq2nd10",
	  "DUT pixel charge 2nd - 0.10 q1;2nd - 0.10 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd11Histo =
    TH1I( "dutpxq2nd11",
	  "DUT pixel charge 2nd - 0.11 q1;2nd - 0.11 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd12Histo =
    TH1I( "dutpxq2nd12",
	  "DUT pixel charge 2nd - 0.12 q1;2nd - 0.12 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd13Histo =
    TH1I( "dutpxq2nd13",
	  "DUT pixel charge 2nd - 0.13 q1;2nd - 0.13 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd14Histo =
    TH1I( "dutpxq2nd14",
	  "DUT pixel charge 2nd - 0.14 q1;2nd - 0.14 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd15Histo =
    TH1I( "dutpxq2nd15",
	  "DUT pixel charge 2nd - 0.15 q1;2nd - 0.15 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );

  TH1I dutcolHisto =
    TH1I( "dutcol",
	  "DUT pixel column;pixel column;pixels",
	  52, -0.5, 51.5 );
  TH1I dutrowHisto =
    TH1I( "dutrow",
	  "DUT pixel row;pixel row;pixels",
	  80, -0.5, 79.5 );

  TH1I dutq0fHisto =
    TH1I( "dutq0f",
	  "normal fiducial cluster charge;normal fiducial cluster charge [ke];fiducial clusters",
	  160, 0, 80 );

  TH1I dutncolHisto =
    TH1I( "dutncol",
	  "DUT cluster size;cluster size [columns];clusters",
	  52, 0.5, 52.5 );
  TH1I dutnrowHisto =
    TH1I( "dutnrow",
	  "DUT cluster size;cluster size [rows];clusters",
	  80, 0.5, 80.5 );
  TH1I dutcolminHisto =
    TH1I( "dutcolmin",
	  "DUT first cluster column;first cluster column;clusters",
	  52, -0.5, 51.5 );
  TH1I dutcolmaxHisto =
    TH1I( "dutcolmax",
	  "DUT last cluster column;last cluster column;clusters",
	  52, -0.5, 51.5 );
  TH1I dutcol0qHisto =
    TH1I( "dutcol0q",
	  "DUT first column charge;first column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol0oddqHisto =
    TH1I( "dutcol0oddq",
	  "DUT odd first column charge;odd first column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol0eveqHisto =
    TH1I( "dutcol0eveq",
	  "DUT eve first column charge;eve first column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol9qHisto =
    TH1I( "dutcol9q",
	  "DUT last column charge;last column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol1qHisto =
    TH1I( "dutcol1q",
	  "DUT 2nd column charge;2nd column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol2qHisto =
    TH1I( "dutcol2q",
	  "DUT 3rd column charge;3rd column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol3qHisto =
    TH1I( "dutcol3q",
	  "DUT 4th column charge;4th column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol4qHisto =
    TH1I( "dutcol4q",
	  "DUT 5th column charge;5th column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol5qHisto =
    TH1I( "dutcol5q",
	  "DUT 6th column charge;6th column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol6qHisto =
    TH1I( "dutcol6q",
	  "DUT 7th column charge;7th column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol7qHisto =
    TH1I( "dutcol7q",
	  "DUT 8th column charge;8th column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol8qHisto =
    TH1I( "dutcol8q",
	  "DUT 9th column charge;9th column charge [ke];clusters",
	  100, 0, 50 );

  // triplets:

  TH1I hdx02 = TH1I( "dx02", "0-2 dx;0-2 dx [mm];cluster pairs", 100, -1, 1 );
  TH1I hdy02 = TH1I( "dy02", "0-2 dy;0-2 dy [mm];cluster pairs", 100, -1, 1 );

  TH1I htridx = TH1I( "tridx", "triplet dx;triplet dx [mm];triplets", 100, -0.1, 0.1 );
  TH1I htridy = TH1I( "tridy", "triplet dy;triplet dy [mm];triplets", 100, -0.1, 0.1 );

  TH1I htridxc = TH1I( "tridxc", "triplet dx;triplet dx [mm];triplets", 100, -0.05, 0.05 );
  TH1I htridyc = TH1I( "tridyc", "triplet dy;triplet dy [mm];triplets", 100, -0.05, 0.05 );

  TH1I htridxc1 = TH1I( "tridxc1", "triplet dx 1-col;1-col triplet dx [mm];1-col triplets",
			100, -0.05, 0.05 );
  TH1I htridxc2 = TH1I( "tridxc2", "triplet dx 2-col;2-col triplet dx [mm];2-col triplets",
			100, -0.05, 0.05 );
  TH1I htridxc3 = TH1I( "tridxc3", "triplet dx 3-col;3-col triplet dx [mm];3-col triplets",
			100, -0.05, 0.05 );
  TH1I htridxc4 = TH1I( "tridxc4", "triplet dx 4-col;4-col triplet dx [mm];4-col triplets",
			100, -0.05, 0.05 );
  TH1I htridxc5 = TH1I( "tridxc5", "triplet dx 5-col;5-col triplet dx [mm];5-col triplets",
			100, -0.05, 0.05 );

  TH1I htridxs1 = TH1I( "tridxs1", "triplet dx 1-px;1-px triplet dx [mm];1-px triplets",
			100, -0.05, 0.05 );
  TH1I htridxs2 = TH1I( "tridxs2", "triplet dx 2-px;2-px triplet dx [mm];2-px triplets",
			100, -0.05, 0.05 );
  TH1I htridxs3 = TH1I( "tridxs3", "triplet dx 3-px;3-px triplet dx [mm];3-px triplets",
			100, -0.05, 0.05 );
  TH1I htridxs4 = TH1I( "tridxs4", "triplet dx 4-px;4-px triplet dx [mm];4-px triplets",
			100, -0.05, 0.05 );
  TH1I htridxs5 = TH1I( "tridxs5", "triplet dx 5-px;5-px triplet dx [mm];5-px triplets",
			100, -0.05, 0.05 );
  TProfile tridxvsy =
    TProfile( "tridxvsy",
	      "triplet dx vs y;triplet yB [mm];<triplets #Deltax> [mm]",
	      110, -5.5, 5.5, -0.05, 0.05 );
  TProfile tridyvsx =
    TProfile( "tridyvsx",
	      "triplet dy vs x;triplet xB [mm];<triplets #Deltay> [mm]",
	      110, -11, 11, -0.05, 0.05 );
  TProfile tridxvstx =
    TProfile( "tridxvstx",
	      "triplet dx vs slope x;triplet slope x [rad];<triplets #Deltax> [mm]",
	      60, -0.003, 0.003, -0.05, 0.05 );
  TProfile tridyvsty =
    TProfile( "tridyvsty",
	      "triplet dy vs slope y;triplet slope y [rad];<triplets #Deltay> [mm]",
	      60, -0.003, 0.003, -0.05, 0.05 );

  TH1I ntriHisto = TH1I( "ntri", "triplets;triplets;events", 51, -0.5, 50.5 );

  TH1I trixHisto = TH1I( "trix", "triplets at z DUT;x [mm];triplets",
			  240, -12, 12 );
  TH1I triyHisto = TH1I( "triy", "triplets at z DUT;y [mm];triplets",
			  120, -6, 6 );
  TH2I * trixyHisto = new
    TH2I( "trixy", "triplets at z DUT;x [mm];y [mm];triplets",
	  240, -12, 12, 120, -6, 6 );
  TH1I tritxHisto = TH1I( "tritx", "triplet slope x;slope x [rad];triplets",
			    100, -0.005, 0.005 );
  TH1I trityHisto = TH1I( "trity", "triplet slope y;slope y [rad];triplets",
			    100, -0.005, 0.005 );

  TH1I ttdxHisto = TH1I( "ttdx", "telescope triplets;triplets #Deltax [mm];triplet pairs",
			 100, -5, 5 );
  TH1I ttdx1Histo = TH1I( "ttdx1", "telescope triplets;triplets #Deltax [mm];triplet pairs",
			  100, -0.5, 0.5 );
  TH1I ttdmin1Histo = TH1I( "ttdmin1",
			    "telescope triplets isolation;triplets min #Delta_{xy} [mm];triplet pairs",
			    100, 0, 1 );
  TH1I ttdmin2Histo = TH1I( "ttdmin2",
			    "telescope triplets isolation;triplets min #Delta_{xy} [mm];triplet pairs",
			    150, 0, 15 );

  // dripets - triplets:

  TH1I hsixdx = TH1I( "sixdx", "six dx;dx [mm];triplet-driplet pairs", 100, -1, 1 );
  TH1I hsixdy = TH1I( "sixdy", "six dy;dy [mm];triplet-driplet pairs", 100, -1, 1 );
  TH1I hsixdxc = TH1I( "sixdxc", "six dx;dx [mm];triplet-driplet pairs", 100, -0.5, 0.5 );
  TH1I hsixdyc = TH1I( "sixdyc", "six dy;dy [mm];triplet-driplet pairs", 100, -0.5, 0.5 );

  TProfile sixdxvsy =
    TProfile( "sixdxvsy",
	      "six #Deltax vs y;yB [mm];<driplet - triplet #Deltax> [mm]",
	      100, -5, 5, -0.5, 0.5 );
  TProfile sixdyvsx =
    TProfile( "sixdyvsx",
	      "six #Deltay vs x;xB [mm];<driplet - triplet #Deltay> [mm]",
	      200, -10, 10, -0.5, 0.5 );

  TProfile sixdxvstx =
    TProfile( "sixdxvstx",
	      "six #Deltax vs slope x;slope x [rad];<driplet - triplet #Deltax> [mm]",
	      100, -0.002, 0.002, -0.5, 0.5 );
  TProfile sixdyvsty =
    TProfile( "sixdyvsty",
	      "six #Deltay vs slope y;slope y [rad];<driplet - triplet #Deltay> [mm]",
	      100, -0.002, 0.002, -0.5, 0.5 );

  TH1I hsixdslpx =
    TH1I( "sixdslpx",
	  "driplet slope x - triplet slope x;driplet slope x - triplet slope x;driplet-triplet pairs",
	  100, -0.005, 0.005 );     
  TH1I hsixdslpy =
    TH1I( "sixdslpy",
	  "driplet slope y - triplet slope y;driplet slope y - triplet slope y;driplet-triplet pairs",
	  100, -0.005, 0.005 );     

  TH2I * sixxyHisto = new
    TH2I( "sixxy", "sixplets at z DUT;x [mm];y [mm];sixplets",
	  240, -12, 12, 120, -6, 6 );

  // DUT pixel vs triplets:

  TH1I z3Histo = TH1I( "z3",
		       "z3 should be zero;z3 [mm];triplets",
		       100, -0.01, 0.01 );

  TH1I cmssxaHisto = TH1I( "cmssxa",
			   "DUT + Telescope x;cluster + triplet #Sigmax [mm];clusters",
			   440, -11, 11 );
  TH1I cmsdxaHisto = TH1I( "cmsdxa",
			   "DUT - Telescope x;cluster - triplet #Deltax [mm];clusters",
			   440, -11, 11 );

  TH1I cmssyaHisto = TH1I( "cmssya",
			   "DUT + Telescope y;cluster + triplet #Sigmay [mm];clusters",
			   440, -11, 11 ); // shallow needs wide range
  TH1I cmsdyaHisto = TH1I( "cmsdya",
			   "DUT - Telescope y;cluster - triplet #Deltay [mm];clusters",
			   440, -11, 11 );

  TH1I cmsdxHisto = TH1I( "cmsdx",
			   "DUT - Telescope x;cluster - triplet #Deltax [mm];clusters",
			   200, -0.5, 0.5 );
  TH1I cmsdyHisto = TH1I( "cmsdy",
			   "DUT - Telescope y;cluster - triplet #Deltay [mm];clusters",
			   500, -0.5, 0.5 );

  TH1I cmslkxHisto = TH1I( "cmslkx",
			   "linked triplet at DUT x;triplet x at DUT [mm];linked triplets",
			   220, -11, 11 );
  TH1I cmslkyHisto = TH1I( "cmslky",
			   "linked triplet at DUT y;triplet y at DUT [mm];linked triplets",
			   110, -5.5, 5.5 );

  TH1I cmscolHisto = TH1I( "cmscol",
			   "DUT linked columns;DUT linked cluster column;linked clusters",
			   52, 0, 52 );
  TH1I cmsrowHisto = TH1I( "cmsrow",
			   "DUT linked rows;DUT linked cluster row;linked clusters",
			   80, 0, 80 );

  TH1I cmsdxfHisto =
    TH1I( "cmsdxf",
	  "fiducial #Deltax ;cluster - triplet #Deltax [mm];fiducial clusters",
	  200, -1, 1 );
  TH1I cmsdx0fHisto =
    TH1I( "cmsdx0f",
	  "fiducial #Deltax0 ;cluster - triplet #Deltax0 [mm];fiducial clusters",
	  200, -1, 1 );
  TH1I cmsdxfq9Histo =
    TH1I( "cmsdxfq9",
	  "fiducial #Deltax q9;cluster - triplet #Deltax [mm];fiducial q9 clusters",
	  200, -1, 1 );
  TH1I cmsdxfcHisto =
    TH1I( "cmsdxfc",
	  "fiducial #Deltax cut y;cluster - triplet #Deltax [mm];fiducial clusters",
	  200, -1, 1 );

  TProfile cmsdxvsx =
    TProfile( "cmsdxvsx",
	      "#Deltax vs x;x track [mm];<cluster - triplet #Deltax> [mm]",
	      50, -3.75, 3.75, -0.5, 0.5 );
  TProfile cmsdxvsy =
    TProfile( "cmsdxvsy",
	      "#Deltax vs y;y track [mm];<cluster - triplet #Deltax> [mm]",
	      76, -3.8, 3.8, -0.5, 0.5 );
  TProfile cmsdxvstx =
    TProfile( "cmsdxvstx",
	      "#Deltax vs #theta_{x};x track slope [rad];<cluster - triplet #Deltax> [mm]",
	      80, -0.004, 0.004, -0.5, 0.5 );

  TH1I cmsdyfq0Histo =
    TH1I( "cmsdyfq0",
	  "fiducial #Deltay q0;cluster - triplet #Deltay [mm];fiducial q0 clusters",
	  200, -1, 1 );
  TH1I cmsdyfq9Histo =
    TH1I( "cmsdyfq9",
	  "fiducial #Deltay q9;cluster - triplet #Deltay [mm];fiducial q9 clusters",
	  200, -1, 1 );
  TH1I cmsdyfcHisto =
    TH1I( "cmsdyfc",
	  "fiducial #Deltay cut x;cluster - triplet #Deltay [mm];fiducial clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq0Histo =
    TH1I( "cmsdyfcq0",
	  "fiducial #Deltay cut x q0;cluster - triplet #Deltay [mm];fiducial q0 clusters",
	  200, -1, 1 );
  TH1I cmsdyfcq9Histo =
    TH1I( "cmsdyfcq9",
	  "fiducial #Deltay cut x q9;cluster - triplet #Deltay [mm];fiducial q9 clusters",
	  200, -1, 1 );
  TH1I cmsdyfcq8Histo =
    TH1I( "cmsdyfcq8",
	  "fiducial #Deltay cut x q8;cluster - triplet #Deltay [mm];fiducial q8 clusters",
	  200, -1, 1 );
  TH1I cmsdyfcq1Histo =
    TH1I( "cmsdyfcq1",
	  "fiducial #Deltay cut x q1;cluster - triplet #Deltay [mm];fiducial q1 clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq2Histo =
    TH1I( "cmsdyfcq2",
	  "fiducial #Deltay cut x q2;cluster - triplet #Deltay [mm];fiducial q2 clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3Histo =
    TH1I( "cmsdyfcq3",
	  "fiducial #Deltay cut x q3;cluster - triplet #Deltay [mm];fiducial q3 clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdy0fcq3Histo =
    TH1I( "cmsdy0fcq3",
	  "fiducial #Deltay0 cut x q3;cluster - triplet #Deltay0 [mm];fiducial q3 clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3dotHisto =
    TH1I( "cmsdyfcq3dot",
	  "fiducial #Deltay cut x q3 dot;cluster - triplet #Deltay [mm];fiducial q3 dot clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3nodHisto =
    TH1I( "cmsdyfcq3nod",
	  "fiducial #Deltay cut x q3 no dot;cluster - triplet #Deltay [mm];fiducial q3 no dot clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3oddHisto =
    TH1I( "cmsdyfcq3odd",
	  "fiducial #Deltay cut x q3 odd col;cluster - triplet #Deltay [mm];fiducial q3 odd clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3eveHisto =
    TH1I( "cmsdyfcq3eve",
	  "fiducial #Deltay cut x q3 eve col;cluster - triplet #Deltay [mm];fiducial q3 eve clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3midHisto =
    TH1I( "cmsdyfcq3mid",
	  "fiducial #Deltay cut x q3 mid;cluster - triplet #Deltay [mm];fiducial q3 mid clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3rimHisto =
    TH1I( "cmsdyfcq3rim",
	  "fiducial #Deltay cut x q3 rim;cluster - triplet #Deltay [mm];fiducial q3 rim clusters",
	  500, -0.5, 0.5 );

  TH1I cmsdyfcq4Histo =
    TH1I( "cmsdyfcq4",
	  "fiducial #Deltay cut x q4;cluster - triplet #Deltay [mm];fiducial q4 clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq4nodHisto =
    TH1I( "cmsdyfcq4nod",
	  "fiducial #Deltay cut x q4 no dot;cluster - triplet #Deltay [mm];fiducial q4 no dot clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq4oddHisto =
    TH1I( "cmsdyfcq4odd",
	  "fiducial #Deltay cut x q4 odd col;cluster - triplet #Deltay [mm];fiducial q4 odd clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq4eveHisto =
    TH1I( "cmsdyfcq4eve",
	  "fiducial #Deltay cut x q4 eve col;cluster - triplet #Deltay [mm];fiducial q4 eve clusters",
	  500, -0.5, 0.5 );

  TProfile cmsrmsyvsq =
    TProfile( "cmsrmsyvsq",
	      "DUT #Deltay vs Q0;normal cluster charge [ke];MAD(#Deltay) [mm]",
	      150, 0, 150, 0, 0.2 );

  TProfile cmsdyvsx =
    TProfile( "cmsdyvsx",
	      "DUT #Deltay vs x;x track [mm];<cluster - triplet #Deltay> [mm]",
	      50, -3.75, 3.75, -0.5, 0.5 );
  TProfile cmsdyvsy =
    TProfile( "cmsdyvsy",
	      "DUT #Deltay vs y;y track [mm];<cluster - triplet #Deltay> [mm]",
	      76, -3.8, 3.8, -0.2, 0.2 );
  TProfile cmsdyvsty =
    TProfile( "cmsdyvsty",
	      "DUT #Deltay vs #theta_{y};y track slope [rad];<cluster - triplet #Deltay> [mm]",
	      80, -0.002, 0.002, -0.2, 0.2 );

  TProfile cmsqrowvsd =
    TProfile( "cmsqrowvsd",
	      "DUT charge vs depth;track depth [mm];<row charge> [ke]",
	      40, -0.2, 0.2, 0, 0.65 ); // 0.65 = cut off 1.5 ke
  TProfile cmsqrowvsdc =
    TProfile( "cmsqrowvsdc",
	      "DUT charge vs depth fiducial;track depth [mm];<fiducial row charge> [ke]",
	      40, -0.2, 0.2, 0, 0.65 );
  TH1I cmsqrownegHisto =
    TH1I( "cmsqrowneg",
	  "row charge negative depth;row charge [ke];linked fiducial deep rows",
	  160, 0, 80 );
  TH1I cmsqrowposHisto =
    TH1I( "cmsqrowpos",
	  "row charge positive depth;row charge [ke];linked fiducial shallow rows",
	  160, 0, 80 );

  TH1I cmsq0fHisto =
    TH1I( "cmsq0f",
	  "normal fiducial cluster charge;normal cluster charge [ke];linked fiducial clusters",
	  160, 0, 80 );
  TH1I cmsq0fdotHisto =
    TH1I( "cmsq0fdot",
	  "normal fiducial cluster charge dot;normal cluster charge dot [ke];linked fiducial dot clusters",
	  160, 0, 80 );
  TH1I cmsq0fnodHisto =
    TH1I( "cmsq0fnod",
	  "normal fiducial cluster charge no dot;normal cluster charge no dot [ke];linked fiducial no dot clusters",
	  160, 0, 80 );
  TH1I cmsq0foddHisto =
    TH1I( "cmsq0fodd",
	  "normal fiducial cluster charge odd col;normal cluster charge odd col [ke];linked fiducial odd col clusters",
	  160, 0, 80 );
  TH1I cmsq0feveHisto =
    TH1I( "cmsq0feve",
	  "normal fiducial cluster charge eve col;normal cluster charge eve col [ke];linked fiducial eve col clusters",
	  160, 0, 80 );
  TH1I cmsq0fcHisto =
    TH1I( "cmsq0fc",
	  "normal fiducial cluster charge core;normal cluster charge core [ke];linked fiducial core clusters",
	  160, 0, 80 );
  TH1I cmsq0fpHisto =
    TH1I( "cmsq0fp",
	  "normal fiducial cluster charge peri;normal cluster charge peri [ke];linked fiducial peri clusters",
	  160, 0, 80 );

  TProfile cmsqxvst1 =
    TProfile( "cmsqxvst1",
	      "DUT cluster charge vs time;time [s];<cluster charge> [ke]",
	      300, 0, 300, 0, 0.1 ); // 0.1 = cut off 8 ke
  TProfile cmsqxvst2 =
    TProfile( "cmsqxvst2",
	      "DUT cluster charge vs time;time [s];<cluster charge> [ke]",
	      150, 0, 1500, 0, 0.1 );
  TProfile cmsqxvst3 =
    TProfile( "cmsqxvst3",
	      "DUT cluster charge vs time;time [s];<cluster charge> [ke]",
	      300, 0, 30000, 0, 0.1 );
  TProfile cmsqxvsx =
    TProfile( "cmsqxvsx",
	      "DUT cluster charge vs x;x track [mm];<cluster charge> [ke]",
	      50, -3.75, 3.75, 0, 0.1 );
  TProfile cmsqxvsy =
    TProfile( "cmsqxvsy",
	      "DUT cluster charge vs y;y track [mm];<cluster charge> [ke]",
	      76, -3.8, 3.8, 0, 0.1 );
  TProfile cmsqxvsxm =
    TProfile( "cmsqxvsxm",
	      "DUT cluster charge vs xmod;x track mod 0.3 [mm];<cluster charge> [ke]",
	      60, 0, 0.3, 0, 0.1 ); // 0.1 = cut off 8 ke
  TProfile cmsqxvsym =
    TProfile( "cmsqxvsym",
	      "DUT cluster charge vs ymod;y track mod 0.2 [mm];<cluster charge> [ke]",
	      40, 0, 0.2, 0, 0.1 );
  TProfile cmsqxvsymn =
    TProfile( "cmsqxvsymn",
	      "DUT cluster charge vs ymod no dot;y track mod 0.2 [mm];<no dot cluster charge> [ke]",
	      40, 0, 0.2, 0, 0.1 );
  TProfile cmsqxvsymd =
    TProfile( "cmsqxvsymd",
	      "DUT cluster charge vs ymod dot;y track mod 0.2 [mm];<dot cluster charge> [ke]",
	      40, 0, 0.2, 0, 0.1 );
  TProfile cmsqxvsxmd =
    TProfile( "cmsqxvsxmd",
	      "DUT cluster charge vs xmod dot;x track mod 0.3 [mm];<dot cluster charge> [ke]",
	      60, 0, 0.3, 0, 0.1 );
  TProfile2D * cmsqxvsxmym = new
    TProfile2D( "cmsqxvsxmym",
	      "DUT cluster charge vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];<cluster charge> [ke]",
		120, 0, 0.3, 80, 0, 0.2, 0, 0.1 );

  TH1I cmsnpxHisto =
    TH1I( "cmsnpx",
	  "linked DUT cluster size;cluster size [pixels];linked fiducial clusters",
	  80, 0.5, 80.5 );
  TH1I cmsncolHisto =
    TH1I( "cmsncol",
	  "linked DUT cluster size;cluster size [columns];linked fiducial clusters",
	  52, 0.5, 52.5 );
  TH1I cmsnrowHisto =
    TH1I( "cmsnrow",
	  "linked DUT cluster size;cluster size [rows];linked fiducial clusters",
	  80, 0.5, 80.5 );
  TProfile cmsnpxvsq =
    TProfile( "cmsnpxvsq",
	      "DUT cluster size vs Q0;normal cluster charge [ke];<cluster size> [pixels]",
	      150, 0, 150, 0, 20 );

  TH1I cmsrmsrowHisto =
    TH1I( "cmsrmsrow",
	  "DUT row rms;DUT cluster row rms;DUT linked clusters",
	  100, 0, 1 );
  TH1I cmsskwrowHisto =
    TH1I( "cmsskwrow",
	  "DUT row skw;DUT cluster row skw;DUT linked clusters",
	  120, -0.3, 0.3 );
  TProfile cmsrmsrowvsq0 =
    TProfile( "cmsrmsrowvsq0",
	      "DUT cluster row rms vs Q0;normal cluster charge [ke];<cluster row rms>",
	      140, 10, 150, 0, 1 );
  TProfile cmsrmsrowvsym =
    TProfile( "cmsrmsrowvsym",
	      "DUT cluster row rms vs ymod;y track mod 0.2 [mm];<cluster row rms>",
	      40, 0, 0.2, 0, 1 );
  TProfile cmsskwrowvsym =
    TProfile( "cmsskwrowvsym",
	      "DUT cluster row skw vs ymod;y track mod 0.2 [mm];<cluster row skw>",
	      40, 0, 0.2, -0.3, 0.3 );

  TProfile cmsdyvsrmsrow =
    TProfile( "cmsdyvsrmsrow",
	      "DUT y resid vs row rms;cluster row rms;<#Deltay> [mm]",
	      40, 0, 1, -0.15, 0.15 );
  TProfile cmsrmsyvsrmsrow =
    TProfile( "cmsrmsyvsrmsrow",
	      "DUT y MAD vs row rms;cluster row rms;MAD(#Deltay) [mm]",
	      40, 0, 1, 0, 0.1 );

  TProfile cmsdyvsskwrow =
    TProfile( "cmsdyvsskwrow",
	      "DUT y resid vs row skw;cluster row skw;<#Deltay> [mm]",
	      80, -0.2, 0.2, -0.15, 0.15 );
  TProfile cmsdy0vsskwrow =
    TProfile( "cmsdy0vsskwrow",
	      "DUT y resid0 vs row skw;cluster row skw;<#Deltay0> [mm]",
	      80, -0.2, 0.2, -0.15, 0.15 );
  TProfile cmsskwrowvsdy =
    TProfile( "cmsskwrowvsdy",
	      "DUT row skw vs y resid;#Deltay [#mum];<cluster row skw>",
	      50, -50, 50, -0.3, 0.3 );
  TProfile cmsskwrowvsdy0 =
    TProfile( "cmsskwrowvsdy0",
	      "DUT row skw vs y0 resid;#Deltay0 [#mum];<cluster row skw>",
	      50, -50, 50, -0.3, 0.3 );

  TProfile cmsrmsyvsskwrow =
    TProfile( "cmsrmsyvsskwrow",
	      "DUT y MAD vs row skw;cluster row skw;MAD(#Deltay) [mm]",
	      80, -0.2, 0.2, 0, 0.1 );

  TH1I cmspxqHisto =
    TH1I( "cmspxq",
	  "DUT pixel charge linked;pixel charge [ke];linked pixels",
	  100, 0, 25 );

  TH1I cmspxqoddHisto =
    TH1I( "cmspxqodd",
	  "DUT pixel charge linked odd col;pixel charge [ke];linked odd col pixels",
	  100, 0, 25 );
  TH1I cmspxqeveHisto =
    TH1I( "cmspxqeve",
	  "DUT pixel charge linked even col;pixel charge [ke];linked even col pixels",
	  100, 0, 25 );
  TH1I cmspxq1Histo =
    TH1I( "cmspxq1",
	  "DUT pixel charge linked 1-px;pixel charge [ke];linked 1-px pixels",
	  100, 0, 25 );
  TH1I cmspxq2Histo =
    TH1I( "cmspxq2",
	  "DUT pixel charge linked 2-px;pixel charge [ke];linked 2-px pixels",
	  100, 0, 25 );
  TH1I cmspxq3Histo =
    TH1I( "cmspxq3",
	  "DUT pixel charge linked 3-px;pixel charge [ke];linked 3-px pixels",
	  100, 0, 25 );
  TProfile cmspxqvsq =
    TProfile( "cmspxqvsq",
	      "DUT pixel charge vs Q0;normal cluster charge [ke];<pixel charge> [ke",
	      150, 0, 150, 0, 50 );
  TProfile cmspxqvsxm =
    TProfile( "cmspxqvsxm",
	      "DUT pixel charge vs xmod;x track mod 0.3 [mm];<pixel charge> [ke]",
	      60, 0, 0.3, 0, 50 );
  TProfile cmspxqvsym =
    TProfile( "cmspxqvsym",
	      "DUT pixel charge vs ymod;y track mod 0.2 [mm];<pixel charge> [ke]",
	      40, 0, 0.2, 0, 50 );

  TProfile cmsdxvsxm =
    TProfile( "cmsdxvsxm",
	      "DUT x residual vs xmod;x track mod 0.3 [mm];<x residual> [mm]",
	      60, 0, 0.3, -0.2, 0.2 );
  TProfile cmsdxvsym =
    TProfile( "cmsdxvsym",
	      "DUT x residual vs ymod;y track mod 0.2 [mm];<x residual> [mm]",
	      40, 0, 0.2, -0.2, 0.2 );
  TProfile cmsdyvsxm =
    TProfile( "cmsdyvsxm",
	      "DUT y residual vs xmod;x track mod 0.3 [mm];<y residual> [mm]",
	      60, 0, 0.3, -0.2, 0.2 );
  TProfile cmsdyvsym =
    TProfile( "cmsdyvsym",
	      "DUT y residual vs ymod;y track mod 0.2 [mm];<y residual> [mm]",
	      40, 0, 0.2, -0.2, 0.2 );
  TProfile cmsdyvsymodd =
    TProfile( "cmsdyvsymodd",
	      "DUT y residual vs ymod odd col;y track mod 0.2 [mm];<y residual> odd col [mm]",
	      40, 0, 0.2, -0.2, 0.2 );
  TProfile cmsdyvsymeve =
    TProfile( "cmsdyvsymeve",
	      "DUT y residual vs ymod even col;y track mod 0.2 [mm];<y residual> even col [mm]",
	      40, 0, 0.2, -0.2, 0.2 );

  TProfile cmsrmsxvsx =
    TProfile( "cmsrmsxvsx",
	      "DUT x resolution vs x;x track [mm];MAD(#Deltax) [mm]",
	      50, -3.75, 3.75, 0, 0.2 );
  TProfile cmsrmsyvsx =
    TProfile( "cmsrmsyvsx",
	      "DUT y resolution vs x;x track [mm];MAD(#Deltay) [mm]",
	      50, -3.75, 3.75, 0, 0.2 );
  TProfile cmsrmsxvsy =
    TProfile( "cmsrmsxvsy",
	      "DUT x resolution vs y;y track [mm];MAD(#Deltax) [mm]",
	      76, -3.8, 3.8, 0, 0.2 );
  TProfile cmsrmsyvsy =
    TProfile( "cmsrmsyvsy",
	      "DUT y resolution vs y;y track [mm];MAD(#Deltay) [mm]",
	      76, -3.8, 3.8, 0, 0.2 );

  TProfile cmsrmsyvstx =
    TProfile( "cmsrmsyvstx",
	      "DUT y resolution vs slope x;track slope x [rad];MAD(#Deltay) [mm]",
	      60, -0.003, 0.003, 0, 0.2 );
  TProfile cmsrmsyvsty =
    TProfile( "cmsrmsyvsty",
	      "DUT y resolution vs slope y;track slope y [rad];MAD(#Deltay) [mm]",
	      60, -0.003, 0.003, 0, 0.2 );

  TProfile cmsrmsxvsxm =
    TProfile( "cmsrmsxvsxm",
	      "DUT x resolution vs xmod;x track mod 0.3 [mm];MAD(#Deltax) [mm]",
	      60, 0, 0.3, 0, 0.2 );
  TProfile cmsrmsyvsxm =
    TProfile( "cmsrmsyvsxm",
	      "DUT y resolution vs xmod;x track mod 0.3 [mm];MAD(#Deltay) [mm]",
	      60, 0, 0.3, 0, 0.2 );
  TProfile cmsrmsxvsym =
    TProfile( "cmsrmsxvsym",
	      "DUT x resolution vs ymod;y track mod 0.2 [mm];MAD(#Deltax) [mm]",
	      40, 0, 0.2, 0, 0.2 );
  TProfile cmsrmsyvsym =
    TProfile( "cmsrmsyvsym",
	      "DUT y resolution vs ymod;y track mod 0.2 [mm];MAD(#Deltay) [mm]",
	      40, 0, 0.2, 0, 0.2 );

  TProfile2D * cmsrmsxvsxmym = new
    TProfile2D( "cmsrmsxvsxmym",
		"DUT x resolution vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];MAD(#Deltax) [mm]",
		60, 0, 0.3, 40, 0, 0.2, 0, 0.2 );
  TProfile2D * cmsrmsyvsxmym = new
    TProfile2D( "cmsrmsyvsxmym",
		"DUT y resolution vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];MAD(#Deltay) [mm]",
		60, 0, 0.3, 40, 0, 0.2, 0, 0.2 );

  TProfile cmsdyvst1 =
    TProfile( "cmsdyvst1",
	      "DUT #Deltay vs time;time [s];<cluster - triplet #Deltay> [mm]",
	      300, 0, 300, -0.1, 0.1 );
  TProfile cmsdyvst2 =
    TProfile( "cmsdyvst2",
	      "DUT #Deltay vs time;time [s];<cluster - triplet #Deltay> [mm]",
	      100, 0, 1000, -0.1, 0.1 );

  TProfile cmsrmsyvst1 =
    TProfile( "cmsrmsyvst1",
	      "DUT y resolution vs time;time [s];MAD(cluster - triplet #Deltay) [mm]",
	      300, 0, 300, -0.1, 0.1 );
  TProfile cmsrmsyvst2 =
    TProfile( "cmsrmsyvst2",
	      "DUT y resolution vs time;time [s];MAD(cluster - triplet #Deltay) [mm]",
	      100, 0, 1000, -0.1, 0.1 );

  TProfile cmsncolvsx =
    TProfile( "cmsncolvsx",
	      "DUT cluster size vs x;x track [mm];<cluster size> [columns]",
	      50, -3.75, 3.75, 0, 20 );
  TProfile cmsncolvsy =
    TProfile( "cmsncolvsy",
	      "DUT cluster size vs y;y track [mm];<cluster size> [columns]",
	      76, -3.8, 3.8, 0, 20 );

  TProfile cmsncolvsxm =
    TProfile( "cmsncolvsxm",
	      "DUT cluster size vs xmod;x track mod 0.3 [mm];<cluster size> [columns]",
	      60, 0, 0.3, 0, 20 );
  TProfile cmsnrowvsxm =
    TProfile( "cmsnrowvsxm",
	      "DUT cluster size vs xmod;x track mod 0.3 [mm];<cluster size> [rows]",
	      60, 0, 0.3, 0, 20 );

  TProfile cmsncolvsym =
    TProfile( "cmsncolvsym",
	      "DUT cluster size vs ymod;y track mod 0.2 [mm];<cluster size> [columns]",
	      40, 0, 0.2, 0, 20 );
  TProfile cmsnrowvsym =
    TProfile( "cmsnrowvsym",
	      "DUT cluster size vs ymod;y track mod 0.2 [mm];<cluster size> [rows]",
	      40, 0, 0.2, 0, 20 );
  TProfile2D * cmsnpxvsxmym = new
    TProfile2D( "cmsnpxvsxmym",
	      "DUT cluster size vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];<cluster size> [pixels]",
		120, 0, 0.3, 80, 0, 0.2, 0, 20 );

  TH2I * trixylkHisto = new
    TH2I( "trixylk",
	  "triplets with DUT;x [mm];y [mm];DUT-linked triplets",
	  240, -12, 12, 120, -6, 6 );
  TH2I * sixxylkHisto = new
    TH2I( "sixxylk",
	  "REF-linked sixplets at z DUT;x [mm];y [mm];REF-linked sixplets",
	  240, -12, 12, 120, -6, 6 );
  TH2I * sixxyeffHisto = new
    TH2I( "sixxyeff",
	  "REF-linked sixplets with DUT;x [mm];y [mm];DUT-REF-linked sixplets",
	  240, -12, 12, 120, -6, 6 );

  TH1I effdminHisto = TH1I( "effdmin",
			 "min DUT - triplet distance;min DUT - triplet #Delta_{xy} [mm];REF linked fiducial iso triplets",
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
	  150, 0, 150 );
  TH1I effclq1Histo =
    TH1I( "effclq1",
	  "nearest cluster charge dxy > 0.1;DUT cluster charge [ke];dxy > 0.1 mm nearest cluster",
	  150, 0, 150 );
  TH1I effclq2Histo =
    TH1I( "effclq2",
	  "nearest cluster charge dxy > 0.2;DUT cluster charge [ke];dxy > 0.2 mm nearest cluster",
	  150, 0, 150 );
  TH1I effclq3Histo =
    TH1I( "effclq3",
	  "nearest cluster charge dxy > 0.3;DUT cluster charge [ke];dxy > 0.3 mm nearest cluster",
	  150, 0, 150 );
  TH1I effclq4Histo =
    TH1I( "effclq4",
	  "nearest cluster charge dxy > 0.4;DUT cluster charge [ke];dxy > 0.4 mm nearest cluster",
	  150, 0, 150 );
  TH1I effclq5Histo =
    TH1I( "effclq5",
	  "nearest cluster charge dxy > 0.5;DUT cluster charge [ke];dxy > 0.5 mm nearest cluster",
	  150, 0, 150 );
  TH1I effclq6Histo =
    TH1I( "effclq6",
	  "nearest cluster charge dxy > 0.6;DUT cluster charge [ke];dxy > 0.6 mm nearest cluster",
	  150, 0, 150 );
  TH1I effclq7Histo =
    TH1I( "effclq7",
	  "nearest cluster charge dxy > 0.7;DUT cluster charge [ke];dxy > 0.7 mm nearest cluster",
	  150, 0, 150 );
  TH1I effclq8Histo =
    TH1I( "effclq8",
	  "nearest cluster charge dxy > 0.8;DUT cluster charge [ke];dxy > 0.8 mm nearest cluster",
	  150, 0, 150 );
  TH1I effclq9Histo =
    TH1I( "effclq9",
	  "nearest cluster charge dxy > 0.9;DUT cluster charge [ke];dxy > 0.9 mm nearest cluster",
	  150, 0, 150 );

  TH1I effclqrHisto =
    TH1I( "effclqr",
	  "nearest cluster charge, 1 driplet-REF;DUT cluster charge [ke];mono nearest clusters",
	  150, 0, 150 );

  TProfile2D * effvsxy = new
    TProfile2D( "effvsxy",
		"DUT efficiency vs x;x track at DUT [mm];y track at DUT [mm];efficiency",
		60, -4.5, 4.5, 90, -4.5, 4.5, -1, 2 ); // bin = pix
  TProfile effvsx =
    TProfile( "effvsx",
	      "DUT efficiency vs x;x track at DUT [mm];efficiency",
	      54, -4.05, 4.05, -1, 2 ); // bin = col
  TProfile effvsy =
    TProfile( "effvsy",
	      "DUT efficiency vs y;y track at DUT [mm];efficiency",
	      80, -4, 4, -1, 2 ); // bin = row

  TProfile effvst1 =
    TProfile( "effvst1",
	      "DUT efficiency vs time;time [s];efficiency",
	      300, 0, 300, -1, 2 );
  TProfile effvst2 =
    TProfile( "effvst2",
	      "DUT efficiency vs time;time [s];efficiency",
	      150, 0, 1500, -1, 2 );
  TProfile effvst3 =
    TProfile( "effvst3",
	      "DUT efficiency vs time;time [s];efficiency",
	      300, 0, 30000, -1, 2 );
  TProfile effvst4 =
    TProfile( "effvst4",
	      "DUT efficiency vs time;time [s];efficiency",
	      600, 0, 60000, -1, 2 );

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
		"DUT efficiency vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];efficiency",
		60, 0, 0.3, 40, 0, 0.2, -1, 2 );
  TProfile effvsxm =
    TProfile( "effvsxm",
	      "DUT efficiency vs xmod;x track mod 0.3 [mm];efficiency",
	      60, 0, 0.3, -1, 2 );

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

  int event_nr = 0;
  uint64_t evTLU0 = 0;
  const double fTLU = 384E6; // 384 MHz TLU clock
  uint64_t prevTLU = 0;

  vector <cluster> cl0[9]; // remember from previous event

  do {
    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() )
      eudaq::PluginManager::Initialize(evt);

    bool ldbg = 0;

    if( event_nr <  0 )
      ldbg = 1;

    if( lev < 10 )
      ldbg = 1;

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns
    if( event_nr < 2  )
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
      cout << "scope processing  " << run << "." << event_nr << "  taken " << evsec << endl;
    else if( event_nr < 100 && event_nr%10 == 0 )
      cout << "scope processing  " << run << "." << event_nr << "  taken " << evsec << endl;
    else if( event_nr < 1000 && event_nr%100 == 0 )
      cout << "scope processing  " << run << "." << event_nr << "  taken " << evsec << endl;
    else if( event_nr%1000 == 0 )
      cout << "scope processing  " << run << "." << event_nr << "  taken " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    vector <cluster> cl[9];

    int DUTyld = 0;

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldbg ) std::cout << "PLANE " << plane.ID() << ": ";

      // /home/pitzl/eudaq/main/include/eudaq/CMSPixelHelper.hh

      int ipl = plane.ID();

      if ( ipl < 7 && run > 28000 )
	ipl = ipl - 1;

      int npx = 0;

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldbg ) 
	  std::cout << plane.GetX(ipix)
		    << " " << plane.GetY(ipix)
		    << " " << plane.GetPixel(ipix) << " ";

	int ix = plane.GetX(ipix); // global column 0..415
	int iy = plane.GetY(ipix); // global row 0..159
	int adc = plane.GetPixel(ipix); // ADC 0..255

	// skip hot pixels:

	int ipx = ix*ny[ipl] + iy;
	if( hotset[ipl].count(ipx) ) continue;

	double q = adc;

	if( ipl == iDUT &&
	    adc > 0 &&
	    ix >= 0 && ix < 52 &&
	    iy >= 0 && iy < 80 ) {

	  double Ared = adc - p4[ix][iy]; // p4 is asymptotic maximum

	  if( Ared >= 0 )
	    Ared = -0.1; // avoid overflow

	  double a3 = p3[ix][iy]; // positive
	  if( weib == 3 )
	    q = p1[ix][iy] *
	      ( pow( -log( -Ared / a3 ), 1/p2[ix][iy] ) - p0[ix][iy] ) * ke;
	  // q = ( (-ln(-(A-p4)/p3))^1/p2 - p0 )*p1

	  if( q < thr ) continue; // offline threshold [ke]

	}  // DUT

	if( ipl == iREF &&
	    adc > 0 &&
	    ix >= 0 && ix < 52 &&
	    iy >= 0 && iy < 80 ) {

	  double Ared = adc - r4[ix][iy]; // r4 is asymptotic maximum

	  if( Ared >= 0 )
	    Ared = -0.1; // avoid overflow

	  double a3 = r3[ix][iy]; // positive
	  if( weib == 3 )
	    q = r1[ix][iy] *
	      ( pow( -log( -Ared / a3 ), 1/r2[ix][iy] ) - r0[ix][iy] ) * refke;
	  // q = ( (-ln(-(A-r4)/r3))^1/r2 - r0 )*r1

	}  // REF

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

	  // flip for upper ROCs into local addresses:

	  if( iy > 79 ) {
	    roc = 15 - roc; // 15..8
	    col = 51 - col; // 51..0
	    row = 159 - iy; // 79..0
	  }

	  ix = x;
	  iy = y;

	} // MOD

	hcol[ipl].Fill( ix );
	hrow[ipl].Fill( iy );
	hmap[ipl]->Fill( ix, iy );

	// fill pixel block for clustering:

	pb[npx].col = ix;
	pb[npx].row = iy;
	pb[npx].adc = adc;
	pb[npx].q = q;
	pb[npx].ord = npx; // readout order
	pb[npx].big = 0;
	++npx;

	if( ipl == iMOD ) {

	  // double big pixels:
	  // 0+1
	  // 2..51
	  // 52+53

	  int col = xm % 52; // 0..51

	  if( col == 0 ) {
	    pb[npx].col = ix-1; // double
	    pb[npx].row = iy;
	    pb[npx-1].adc *= 0.5;
	    pb[npx-1].q *= 0.5;
	    pb[npx].adc = 0.5*adc;
	    pb[npx].q = 0.5*q;
	    pb[npx].big = 1;
	    ++npx;
	  }

	  if( col == 51 ) {
	    pb[npx].col = ix+1; // double
	    pb[npx].row = iy;
	    pb[npx-1].adc *= 0.5;
	    pb[npx-1].q *= 0.5;
	    pb[npx].adc = 0.5*adc;
	    pb[npx].q = 0.5*q;
	    pb[npx].big = 1;
	    ++npx;
	  }

	  if( ym == 79 ) {
	    pb[npx].col = ix; // double
	    pb[npx].row = 80;
	    pb[npx-1].adc *= 0.5;
	    pb[npx-1].q *= 0.5;
	    pb[npx].adc = 0.5*adc;
	    pb[npx].q = 0.5*q;
	    pb[npx].big = 1;
	    ++npx;
	  }

	  if( ym == 80 ) {
	    pb[npx].col = ix; // double
	    pb[npx].row = 81;
	    pb[npx-1].adc *= 0.5;
	    pb[npx-1].q *= 0.5;
	    pb[npx].adc = 0.5*adc;
	    pb[npx].q = 0.5*q;
	    pb[npx].big = 1;
	    ++npx;
	  }

	} // MOD

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

      if( ldbg ) cout << "clusters " << cl[ipl].size() << endl;

      hncl[ipl].Fill( cl[ipl].size() );

      for( vector<cluster>::iterator c = cl[ipl].begin(); c != cl[ipl].end(); ++c ) {

	hsiz[ipl].Fill( c->size );
	hncol[ipl].Fill( c->ncol );
	hnrow[ipl].Fill( c->nrow );

      } // cl

      if( ipl == iDUT ) {
	if( npx ) DUTyld = 1;
	dutnpxvst2.Fill( evsec, npx );
	dutnclvst2.Fill( evsec, cl[ipl].size() );
      }

    } // planes

    dutyldvst2.Fill( evsec, DUTyld );

    if( ! syncdut )
      cl0[iDUT] = cl[iDUT];
    if( ! syncref )
      cl0[iREF] = cl[iREF];
    if( ! syncmod )
      cl0[iMOD] = cl[iMOD];

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make driplets 3+5-4:

    vector <triplet> driplets;
    vector <triplet> dripletsmod;

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

	if( fabs( dx2 ) > 0.005 * dz35 ) continue; // angle cut
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

    ndriHisto.Fill( driplets.size() );
    dripletsmod = driplets;
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // driplets vs REF:

    double xcutREF = 0.12;
    double ycutREF = 0.08;
    int nm = 0;
    int nmmod = 0;
    int ndrilk = 0;

    const double cr = cos( REFrot );
    const double sr = sin( REFrot );

    double xcutMOD = 0.18;
    double ycutMOD = 0.12;

    for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // jB = downstream

      double avx = driplets[jB].xm;
      double avy = driplets[jB].ym;
      double avz = driplets[jB].zm;
      double slx = driplets[jB].sx;
      double sly = driplets[jB].sy;

      double zB = REFz - avz; // z REF from mid of driplet
      double xB = avx + slx * zB; // driplet impact point on REF
      double yB = avy + sly * zB;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // dri vs dri: isolation at REF

      double dddmin = 99.9;

      for( unsigned int jj = 0; jj < driplets.size(); ++jj ) {

	if( jj == jB ) continue;

	double zj = REFz - driplets[jj].zm;
	double xj = driplets[jj].xm + driplets[jj].sx * zj; // triplet impact point on DUT
	double yj = driplets[jj].ym + driplets[jj].sy * zj;

	double dx = xB - xj;
	double dy = yB - yj;
	double dd = sqrt( dx*dx + dy*dy );
	if( dd < dddmin ) dddmin = dd;

      } // jj

      dddmin1Histo.Fill( dddmin );
      dddmin2Histo.Fill( dddmin );
      driplets[jB].ttdmin = dddmin;


      // transform into REF system:

      double x3 = xB*cr + yB*sr; // rotate
      double y3 = yB*cr - xB*sr;

      double x4 = -x3 + REFalignx; // invert and shift to mid
      double y4 =  y3 + REFaligny; // shift to mid

      double xmod = fmod( 9.000 + x4, 0.3 ); // [0,0.3] mm, 2 pixel wide
      double ymod = fmod( 9.000 + y4, 0.2 ); // [0,0.2] mm

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // driplets vs REF clusters:

      for( vector<cluster>::iterator c = cl0[iREF].begin(); c != cl0[iREF].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;
	double refx = ( ccol + 0.5 - nx[iREF]/2 ) * ptchx[iREF]; // -3.9..3.9 mm
	double refy = ( crow + 0.5 - ny[iREF]/2 ) * ptchy[iREF]; // -4..4 mm
	double q = c->charge;
	double qx = exp( -q/qwid);

	int npx = c->size;

	// residuals for pre-alignment:

	refsxaHisto.Fill( refx + x3 ); // peak
	refdxaHisto.Fill( refx - x3 ); // 

	refsyaHisto.Fill( refy + y3 ); // 
	refdyaHisto.Fill( refy - y3 ); // peak

	double refdx = refx - x4;
	double refdy = refy - y4;

	if( fabs( refdx ) < xcutREF ) {
	  refdyHisto.Fill( refdy );
	  refdyvsx.Fill( x4, refdy );
	  refdyvsy.Fill( y4, refdy );
	  refdyvsty.Fill( sly, refdy );
	}

	if( fabs( refdy ) < ycutREF ) {
	  refdxHisto.Fill( refdx );
	}

	if( fabs( refdx ) < xcutREF &&
	    fabs( refdy ) < ycutREF &&
	    fabs( x4 ) < 3.9 && // fiducial at REF
	    fabs( y4 ) < 3.9 ) {
	  refnpxHisto.Fill( npx );
	  refqHisto.Fill( q );
	  refnpxvsxmym->Fill( xmod, ymod, npx );
	  refqxvsxmym->Fill( xmod, ymod, qx );
	}

	if( fabs( refdx ) < xcutREF && fabs( refdy ) < ycutREF ) { // 99.93  99.97

	  reflkxBHisto.Fill( xB );
	  reflkyBHisto.Fill( yB );
	  reflkxHisto.Fill( x4 );
	  reflkyHisto.Fill( y4 );
	  reflkcolHisto.Fill( ccol );
	  reflkrowHisto.Fill( crow );

	  if( fabs( x4 ) < 3.9 && // fiducial at REF
	      fabs( y4 ) < 3.9 ) { // 3.8 same eff

	    driplets[jB].lk = 1;
	    nm = 1; // we have a REF-driplet match in this event
	    ++ndrilk;

	  } // REF fiducial

	} // REF link x and y

      } // REF

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // driplet at MOD:

      zB = MODz - avz; // z MOD from mid of driplet

      // intersect inclined track with tilted MOD plane:

      double zc = (Nzm*zB - Nym*avy - Nxm*avx) / (Nxm*slx + Nym*sly + Nzm); // from avz
      double yc = avy + sly * zc;
      double xc = avx + slx * zc;

      double dzc = zc + avz - MODz; // from MOD z0 [-8,8] mm
      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // dri vs dri: isolation at MOD
      
      dddmin = 99.9;
      
      for( unsigned int jj = 0; jj < dripletsmod.size(); ++jj ) {
	
	if( jj == jB ) continue;
	
	double zj = MODz - driplets[jj].zm;
	double xj = dripletsmod[jj].xm + dripletsmod[jj].sx * zj; // triplet impact point on DUT
	double yj = dripletsmod[jj].ym + dripletsmod[jj].sy * zj;
	
	double dx = xc - xj;
	double dy = yc - yj;
	double dd = sqrt( dx*dx + dy*dy );
	if( dd < dddmin ) dddmin = dd;
	
      } // jj
      
      dddmin1modHisto.Fill( dddmin );
      dddmin2modHisto.Fill( dddmin );
      dripletsmod[jB].ttdmin = dddmin;


      // transform into MOD system: (passive).
      // large rotations don't commute: careful with order

      double x1 = com*xc - som*dzc; // turn o
      double y1 = yc;
      double z1 = som*xc + com*dzc;

      double x2 = x1;
      double y2 = cam*y1 + sam*z1; // tilt a

      x3 = cfm*x2 + sfm*y2; // rot
      y3 =-sfm*x2 + cfm*y2;

      x4 =-x3 + MODalignx; // shift to mid
      y4 = y3 + MODaligny; // invert y, shift to mid

      xmod = fmod( 9.000 + x4, 0.3 ); // [0,0.3] mm, 2 pixel wide
      ymod = fmod( 9.000 + y4, 0.2 ); // [0,0.2] mm

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // driplets vs MOD clusters:

      for( vector<cluster>::iterator c = cl0[iMOD].begin(); c != cl0[iMOD].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;
	double modx = ( ccol + 0.5 - nx[iMOD]/2 ) * ptchx[iMOD]; // -3.9..3.9 mm
	double mody = ( crow + 0.5 - ny[iMOD]/2 ) * ptchy[iMOD]; // -4..4 mm
	double q = c->charge;
	double qx = exp( -q/qwid);

	int npx = c->size;

	// residuals for pre-alignment:

	modsxaHisto.Fill( modx + x3 ); // peak
	moddxaHisto.Fill( modx - x3 ); // 

	modsyaHisto.Fill( mody + y3 ); // 
	moddyaHisto.Fill( mody - y3 ); // peak

	double moddx = modx - x4;
	double moddy = mody - y4;

	if( fabs( moddx ) < xcutMOD &&
	    c->big == 0 ) {
	  moddyHisto.Fill( moddy );
	  moddyvsx.Fill( x4, moddy );
	  moddyvsy.Fill( y4, moddy );
	  moddyvsty.Fill( sly, moddy );
	}

	if( fabs( moddy ) < ycutMOD &&
	    c->big == 0 ) {
	  moddxHisto.Fill( moddx );
	  moddxvsx.Fill( x4, moddx );
	  moddxvsy.Fill( y4, moddx );
	  moddxvstx.Fill( slx, moddx );
	}

	if( fabs( moddx ) < xcutMOD &&
	    fabs( moddy ) < ycutMOD &&
	    c->big == 0 ) {
	  modnpxHisto.Fill( npx );
	  modqHisto.Fill( q );
	  modnpxvsxmym->Fill( xmod, ymod, npx );
	  modqxvsxmym->Fill( xmod, ymod, qx );
	  dripletsmod[jB].lk = 1; //store driplets that have match
	  nmmod = 1; // we have a MOD-driplet match in this event
	}

	if( fabs( moddx ) < xcutMOD &&
	    fabs( moddy ) < ycutMOD ) {

	  modlkxBHisto.Fill( xB );
	  modlkyBHisto.Fill( yB );
	  modlkxHisto.Fill( x4 );
	  modlkyHisto.Fill( y4 );
	  modlkcolHisto.Fill( ccol );
	  modlkrowHisto.Fill( crow );
	  modlkxvsyHisto.Fill( x4, y4 );

	} // MOD link x and y

      } // MOD

    } // driplets

    reflkvst1.Fill( evsec, nm ); // REF yield vs time
    reflkvst2.Fill( evsec, nm );
    reflkvst3.Fill( evsec, nm );
    ndrilkHisto.Fill( ndrilk );

    modlkvst1.Fill( evsec, nmmod ); // Mod yield vs time
    modlkvst2.Fill( evsec, nmmod );
    modlkvst3.Fill( evsec, nmmod );

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
      dutpxq2ndHisto.Fill( q2 );

      dutpxq2nd02Histo.Fill( q2 - 0.02*q1 ); // 
      dutpxq2nd03Histo.Fill( q2 - 0.03*q1 ); // good for 506
      dutpxq2nd04Histo.Fill( q2 - 0.04*q1 ); // too little
      dutpxq2nd05Histo.Fill( q2 - 0.05*q1 ); // good
      dutpxq2nd06Histo.Fill( q2 - 0.06*q1 ); // best for digV2.1
      dutpxq2nd07Histo.Fill( q2 - 0.07*q1 ); // too much
      dutpxq2nd08Histo.Fill( q2 - 0.08*q1 ); // 
      dutpxq2nd09Histo.Fill( q2 - 0.09*q1 ); // 
      dutpxq2nd10Histo.Fill( q2 - 0.10*q1 ); // 
      dutpxq2nd11Histo.Fill( q2 - 0.11*q1 ); // 
      dutpxq2nd12Histo.Fill( q2 - 0.12*q1 ); // 
      dutpxq2nd13Histo.Fill( q2 - 0.13*q1 ); // 
      dutpxq2nd14Histo.Fill( q2 - 0.14*q1 ); // 
      dutpxq2nd15Histo.Fill( q2 - 0.15*q1 ); // 

    } // DUT cl

    // apply:

    for( vector<cluster>::iterator c = cl0[iDUT].begin(); c != cl0[iDUT].end(); ++c ) {

      // pix in clus:

      for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	if( px->ord == 0 ) continue;

	// look for previous, in all clusters:

	double qprv = 0;

	for( vector<cluster>::iterator d = cl0[iDUT].begin(); d != cl0[iDUT].end(); ++d )
	  for( vector<pixel>::iterator qx = d->vpix.begin(); qx != d->vpix.end(); ++qx )
	    if( qx->ord == px->ord - 1 ) 
	      qprv = qx->q;

	// apply tsunami correction:

	px->q -= eps*qprv; // overwrite!
	//px->q -= 1.4; // [ke] worse

      } // pix

    } // DUT clusters

    // correct row-trend:

    for( vector<cluster>::iterator c = cl0[iDUT].begin(); c != cl0[iDUT].end(); ++c )
      for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px )
	px->q *= ( 1 + yco * ( px->row - 40 ) );

    // tsunami corrected clusters:

    for( vector<cluster>::iterator c = cl0[iDUT].begin(); c != cl0[iDUT].end(); ++c ) {

      int colmin = 99;
      int colmax = -1;
      int rowmin = 99;
      int rowmax = -1;

      double qsum = 0;
      double sumcol = 0;
      double sumrow = 0;

      double qcol[52];
      for( int icol = 0; icol < 52; ++icol ) qcol[icol] = 0;

      double qrow[80];
      for( int irow = 0; irow < 80; ++irow ) qrow[irow] = 0;

      for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	int icol = px->col;
	int irow = px->row;

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

      // cluster charge after tsunami correction:

      c->charge = qsum; // overwrite ! creates shift in cmsqvsxmym
      c->col = sumcol / qsum; // overwrite !
      c->row = sumrow / qsum; // overwrite !

      if( colmin > 0 && colmax < 51 && rowmin > 0 && rowmax <79 ) {
	dutq0fHisto.Fill( qsum * norm );
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

	if( fabs( dx2 ) > 0.005 * dz02 ) continue; // angle cut
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

	} // cl B

      } // cl C

    } // cl A

    ntriHisto.Fill( triplets.size() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets:

    double xcut = 0.4;
    double ycut = 0.2;

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

      double avx = triplets[iA].xm;
      double avy = triplets[iA].ym;
      double avz = triplets[iA].zm;
      double slx = triplets[iA].sx;
      double sly = triplets[iA].sy;
      double txy = sqrt( slx*slx + sly*sly );

      double zA = DUTz - avz; // z DUT from mid of triplet
      double xA = avx + slx * zA; // triplet impact point on DUT
      double yA = avy + sly * zA;

      trixHisto.Fill( xA );
      triyHisto.Fill( yA );
      trixyHisto->Fill( xA, yA );
      tritxHisto.Fill( slx );
      trityHisto.Fill( sly );

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
      // match triplet and driplet for efficiency:

      bool lsixlk = 0;
      double sixcut = 0.1; // [mm]
      double dddmin = 99.9; // driplet isolation at REF
      double sixdslp = 0.099; // [rad]

      vector <triplet> dripletsForEff;
      dripletsForEff = driplets;
      if (useMODasREF){
	dripletsForEff = dripletsmod;
      }

      for( unsigned int jB = 0; jB < dripletsForEff.size(); ++jB ) { // j = B = downstream

	double avxB = dripletsForEff[jB].xm;
	double avyB = dripletsForEff[jB].ym;
	double avzB = dripletsForEff[jB].zm;
	double slxB = dripletsForEff[jB].sx;
	double slyB = dripletsForEff[jB].sy;

	// driplet at mid:

	double zB = DUTz - avzB; // z from mid of triplet to mid
	double xB = avxB + slxB * zB; // triplet at mid
	double yB = avyB + slyB * zB;

	// driplet - triplet:

	double dx = xB - xA;
	double dy = yB - yA;

	hsixdx.Fill( dx ); // for align fit
	hsixdy.Fill( dy ); // for align fit
	if( fabs(dy) < sixcut ) {
	  hsixdxc.Fill( dx );
	  sixdxvsy.Fill( yB, dx );
	  sixdxvstx.Fill( slx, dx );
	}
	if( fabs(dx) < sixcut ) {
	  hsixdyc.Fill( dy );
	  sixdyvsx.Fill( xB, dy );
	  sixdyvsty.Fill( sly, dy );
	}

	// compare slopes:

	if( fabs(dy) < 0.1 && fabs(dx) < 0.1 ) {
	  hsixdslpx.Fill( slxB - slx );
	  hsixdslpy.Fill( slyB - sly ); // width: 0.3 mrad
	}

	if( fabs(dx) < sixcut && fabs(dy) < sixcut ) {
	  sixxyHisto->Fill( xA, yA );
	  if( dripletsForEff[jB].lk ) {
	    lsixlk = 1;
	    dddmin = dripletsForEff[jB].ttdmin;
	    sixdslp = sqrt( pow( slxB - slx, 2 ) + pow( slyB - sly, 2 ) );
	  }
	}

      } // dripletsForEff

      if( lsixlk && dddmin < 0.6 ) liso = 0; // require isolation at REF

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // intersect inclined track with tilted DUT plane:

      double zc = (Nz*zA - Ny*avy - Nx*avx) / (Nx*slx + Ny*sly + Nz); // from avz
      double yc = avy + sly * zc;
      double xc = avx + slx * zc;

      double dzc = zc + avz - DUTz; // from DUT z0 [-8,8] mm

      // transform into DUT system: (passive).
      // large rotations don't commute: careful with order

      double x1 = co*xc - so*dzc; // turn o
      double y1 = yc;
      double z1 = so*xc + co*dzc;

      double x2 = x1;
      double y2 = ca*y1 + sa*z1; // tilt a
      double z2 =-sa*y1 + ca*z1; // should be zero (in DUT plane). is zero

      double x3 = cf*x2 + sf*y2; // rot
      double y3 =-sf*x2 + cf*y2;
      double z3 = z2; // should be zero (in DUT plane). is zero

      z3Histo.Fill( z3 ); // is zero

      double x4 = upsignx*x3 + DUTalignx; // shift to mid
      double y4 =-upsigny*y3 + DUTaligny; // invert y, shift to mid

      bool fiducial = 1;
      if( fabs( x4 ) > 3.9 ) fiducial = 0; // skip big col
      if( fabs( y4 ) > 3.9 ) fiducial = 0; // skip 1st and last row

      // reduce to 2x2 pixel region:

      double xmod = fmod( 9.000 + x4, 0.3 ); // [0,0.3] mm, 2 pixel wide
      double ymod = fmod( 9.000 + y4, 0.2 ); // [0,0.2] mm
      if( rot90 ) { // x = col = yt, y = row = xt
	xmod = fmod( 9.000 + y4, 0.3 ); // [0,0.3] mm, 2 pixel wide
	ymod = fmod( 9.000 + x4, 0.2 ); // [0,0.2] mm
      }

      // bias dot, from cmsqvsxmym:

      bool ldot = 1;
      if( xmod < 0.105 ) ldot = 0; // dot at x = 125
      if( xmod > 0.195 ) ldot = 0; // and at x = 175
      if( DUTtilt < 6 ) {
	if( ymod < 0.055 ) ldot = 0; // dot at y =  75
	if( ymod > 0.195 ) ldot = 0; // dot at y = 175
	if( ymod > 0.095 && ymod < 0.155 ) ldot = 0; // band between dots
      }

      bool ymid = 0;
      if( ymod > 0.040 && ymod < 0.060 ) ymid = 1;
      if( ymod > 0.140 && ymod < 0.160 ) ymid = 1;

      // pixel core, 2x2 pixel region:

      bool lcore = 1;
      if( xmod < 0.015 ) lcore = 0; // outer edge, see cmsncolvsxm
      if( xmod > 0.285 ) lcore = 0; // outer edge
      if( xmod > 0.105 && xmod < 0.195 ) lcore = 0; // cmsqxvsxmym
      if( ymod > 0.030 && ymod < 0.070 ) lcore = 0; // inner edge
      if( ymod > 0.130 && ymod < 0.170 ) lcore = 0; // inner edge

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT pixel clusters:

      int nm = 0; 
      double dmin = 19.9; // [mm]
      double dxmin = 9;
      double dymin = 9;
      double clQ0 = 0;

      for( vector<cluster>::iterator c = cl0[iDUT].begin(); c != cl0[iDUT].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;
	double ccol0 = ccol; // before correction
	double crow0 = crow; // before correction

	int oddcol = int(ccol) % 2;

	int colmin = 99;
	int colmax = -1;
	int rowmin = 99;
	int rowmax = -1;

	double qcol[52];
	for( int icol = 0; icol < 52; ++icol ) qcol[icol] = 0;

	double qrow[80];
	for( int irow = 0; irow < 80; ++irow ) qrow[irow] = 0;

	for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	  int icol = px->col;
	  if( icol <  0 ) continue;
	  if( icol > 51 ) continue;

	  int irow = px->row;
	  if( irow <  0 ) continue;
	  if( irow > 79 ) continue;

	  if( icol < colmin ) colmin = icol;
	  if( icol > colmax ) colmax = icol;
	  if( irow < rowmin ) rowmin = irow;
	  if( irow > rowmax ) rowmax = irow;

	  double q = px->q; // [ke] corrected
	  if( q < 0 ) continue;

	  qcol[icol] += q; // project cluster onto cols
	  qrow[irow] += q; // project cluster onto rows

	} // pix

	int ncol = colmax - colmin + 1;
	int nrow = rowmax - rowmin + 1;

	bool fiducialc = fiducial;
	if( colmin ==  0 ) fiducialc = 0;
	if( colmax == 51 ) fiducialc = 0;
	if( rowmin ==  0 ) fiducialc = 0;
	if( rowmax == 79 ) fiducialc = 0;

	// eta-algo in rows:

	double q1 = 0; // highest charge
	double q2 = 0; // 2nd highest
	int i1 = 99;
	int i2 = 99;
	double sumq = 0;
	double sumrow = 0;
	double sumrow2 = 0;
	double sumrow3 = 0;

	for( int irow = rowmin; irow <= rowmax; ++irow ) {

	  double q = qrow[irow];
	  sumq += q;
	  sumrow += irow*q;
	  double drow = irow - crow;
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

	double rmsrow = 0;
	double skwrow = 0;
	if( nrow > 1 ) {
	  rmsrow = sqrt( sumrow2/sumq );
	  skwrow = sumrow3/sumq*8/nrow/nrow/nrow; // normalized 3rd moment
	  // from fitskw:
	  crow -= // overwrite!
	    ( skwoff + skwslp * (skwrow-skwmid) + // Sun 5.7.2015
	      skwrng * sin( (skwrow-skwmid) / skwwid ) ) * 1E-3 / ptchy[iDUT];
	}

	double q12 = q1 + q2;
	double eta = 0;
	if( q12 > 1 ) eta = ( q1 - q2 ) / q12;
	if( i2 > i1 ) eta = -eta;

	if( nrow == -2 ) { // inactive, need gain cal

	  // correct even/odd PH effect in readout direction:

	  double A1 = q1; // larger
	  double A2 = q2; // smaller

	  //double A12 = A1 + A2;

	  // cross talk for 2-row clusters:
	  // A1 = (1-cx)U1 + cx U2
	  // A2 = (1-cx)U2 + cx U1

	  // .x fitgp0.C("cmsdyfctq3d",-22,22) (no bias dot)
	  // minimize bias in cmsdyvsym.Draw()

	  //double cx =-0.02; // chip 205, run 10891, sy 7.07
	  double cx = 0.00; // chip 205, run 10891, sy 6.98+-0.05, run 6701 6.98
	  //double cx = 0.02; // chip 205, run 10891, sy 7.10, run 6701 7.06
	  //double cx = 0.04; // chip 205, run 10891, sy 7.21

	  //if( chip0 == 203 ) cx = 0.12; // cmsrmsyvseta worse
	  //if( chip0 == 203 ) cx =-0.12; // cmsrmsyvseta better but tilted sy 11.0
	  //if( chip0 == 203 ) cx =-0.25; // cmsrmsyvseta flat and hump at 0.5, sy 12.0

	  // analog chip 113, run 11553 .x fittp0.C("cmsdyfctq3d",-22,22)
	  // cx -0.02  sy 7.29
	  // cx  0.00  sy 7.22
	  // cx  0.02  sy 7.30 +- 0.035

	  //if( chip0 > 399 ) cx = 0.05; // digV2.1 sy 8.47 
	  //if( chip0 > 399 ) cx = 0.08; // digV2.1 sy 8.11 run 12521
	  if( chip0 > 399 ) cx = 0.10; // digV2.1 sy 8.07
	  //if( chip0 > 399 ) cx = 0.12; // digV2.1 sy 8.14

	  if( chip0 > 499 ) cx = 0.00; // chip 500  sy 8.28
	  //if( chip0 > 499 ) cx =  0.02; // chip 500  sy 8.52
	  //if( chip0 > 499 ) cx = -0.02; // chip 500  sy 8.32

	  // unfold = invert cross talk matrix:

	  double det = 1-2*cx; // Determinante

	  double U1 = ( (1-cx)*A1 - cx*A2 ) / det;
	  double U2 = ( (1-cx)*A2 - cx*A1 ) / det; // U1+U2 = A12

	  // threshold:

	  if( U1 < 1.5 ) U1 = 1.5; // [ke] threshold
	  if( U2 < 1.5 ) U2 = 1.5; // This alone improves resolution !

	  double U12 = U1 + U2;
	  if( fabs(cx) > 0.01 ) // for chips 400
	    crow = ( i1*U1 + i2*U2 ) / U12; // overwrite !

	  if( U12 > 1 ) eta = ( U1 - U2 ) / U12;
	  if( i2 > i1 ) eta = -eta;

	} // 2-row

	// column cluster:

	sumq = 0;
	double sumcol = 0;
	double sumcol2 = 0;
	double sumcol3 = 0;
	double p1 = 0; // highest charge
	double p2 = 0; // 2nd highest
	int j1 = 99;
	int j2 = 99;

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

	//double rmscol = 0;
	double skwcol = 0;
	if( ncol > 1 ) {
	  //rmscol = sqrt( sumcol2/sumq );
	  skwcol = sumcol3/sumq*8/ncol/ncol/ncol; // normalized 3rd moment
	  // from fitskw:
	  ccol -= // overwrite!
	    ( skwoff + skwslp * (skwcol-skwmid) + // Fri 3.7.2015
	      skwrng * sin( (skwcol-skwmid) / skwwid ) ) * 1E-3 / ptchx[iDUT];
	}

	double p12 = p1 + p2;
	double uta = 0;
	if( p12 > 1 ) uta = ( p1 - p2 ) / p12;
	if( j2 > j1 ) uta = -uta;

	if( rot90 )
	  eta = uta;

	double Q0 = c->charge * norm; // cluster charge normalized to vertical incidence
	double Qx = exp(-Q0/qwid);

	bool lq = 0;
	if( Q0 > 17 &&  Q0 < 30 ) lq = 1; // for Q0 Landau peak at 22

	// DUT - triplet:

	double cmsx = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	double cmsy = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm

	double cmsx0 = ( ccol0 + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // before corr
	double cmsy0 = ( crow0 + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // 

	if( rot90 ) {
	  cmsx = ( crow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	  cmsy = ( ccol + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	  cmsx0 = ( crow0 + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	  cmsy0 = ( ccol0 + 0.5 - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	}

	// residuals for pre-alignment:

	cmssxaHisto.Fill( cmsx + x3 ); // rot, tilt and turn but no shift
	cmsdxaHisto.Fill( cmsx - x3 ); // peak

	cmssyaHisto.Fill( cmsy + y3 ); // peak
	cmsdyaHisto.Fill( cmsy - y3 );

	// residuals:

	double cmsdx = cmsx - x4; // new style
	double cmsdy = cmsy - y4;
	double cmsdx0 = cmsx0 - x4; // before skw correction
	double cmsdy0 = cmsy0 - y4;

	cmsdxHisto.Fill( cmsdx );
	cmsdyHisto.Fill( cmsdy );

	if( chip0 == 500 ) { // 28.3.2015
	  if( oddcol )
	    cmsdy -= 1.2E-3;
	  else
	    cmsdy += 1.2E-3;
	}
	if( chip0 == 504 ) {
	  if( oddcol )
	    cmsdy -= 0.0E-3; // 2016
	  else
	    cmsdy += 0.5E-3; // 2016
	}

	if( fabs(cmsdx) < xcut  && fabs(cmsdy) < ycut ) {

	  cmslkxHisto.Fill( xc );
	  cmslkyHisto.Fill( yc );
	  cmscolHisto.Fill( ccol ); // map
	  cmsrowHisto.Fill( crow );

	}

	// cuts:

	if( fiducialc &&
	    ( cl0[iDUT].size() == 1 || chip0 > 600 ) && // single cluster = clean events
	    liso && lsixlk ) { // require REF link (cleaner)

	  // for dx:

	  cmsdxfHisto.Fill( cmsdx );
	  cmsdx0fHisto.Fill( cmsdx0 );
	  if( Q0 > 77 )  // 1% Landau tail
	    cmsdxfq9Histo.Fill( cmsdx );

	  if( fabs(cmsdy) < ycut ) {

	    cmsdxfcHisto.Fill( cmsdx );

	    // profiles for alignment:

	    if( lq ) {
	      cmsdxvsx.Fill( x4, cmsdx ); //cmsdxvsx->Fit("pol1","","",-3.5,3.5)
	      cmsdxvsy.Fill( y4, cmsdx );
	      cmsdxvstx.Fill( slx, cmsdx );
	    } // lq

	  } // dy

	  // for dy:

	  if(      Q0 < 15 ) // broken cluster
	    cmsdyfq0Histo.Fill( cmsdy );
	  else if( Q0 > 77 ) // 1% Landau tail
	    cmsdyfq9Histo.Fill( cmsdy );

	  if( fabs(cmsdx) < xcut ) {

	    cmsdyfcHisto.Fill( cmsdy );
	    cmsrmsyvsq.Fill( Q0, fabs(cmsdy) ); //resolution vs charge

	    if(      Q0 < 15 ) // broken cluster
	      cmsdyfcq0Histo.Fill( cmsdy );

	    else if( Q0 < 17 )
	      ;

	    else if( Q0 > 77 ) // 1% Landau tail
	      cmsdyfcq9Histo.Fill( cmsdy ); // 82 um

	    else if( Q0 > 35 ) // Landau tail
	      cmsdyfcq8Histo.Fill( cmsdy );

	    else { // 17..35

	      cmsdyfcq1Histo.Fill( cmsdy );

	      if( Q0 < 30 ) {
		cmsdyfcq2Histo.Fill( cmsdy );
	      }

	      if( Q0 < 25 ) {

		cmsdyfcq3Histo.Fill( cmsdy );
		cmsdy0fcq3Histo.Fill( cmsdy0 );

		if( ldot )
		  cmsdyfcq3dotHisto.Fill( cmsdy );
		else {
		  cmsdyfcq3nodHisto.Fill( cmsdy );
		  if( oddcol )
		    cmsdyfcq3oddHisto.Fill( cmsdy );
		  else
		    cmsdyfcq3eveHisto.Fill( cmsdy );
		}

		if( ymid ) {
		  cmsdyfcq3midHisto.Fill( cmsdy );
		}
		else
		  cmsdyfcq3rimHisto.Fill( cmsdy );

	      }

	      if( Q0 > 19 && Q0 < 23 ) {

		cmsdyfcq4Histo.Fill( cmsdy );

		if( !ldot) {

		  cmsdyfcq4nodHisto.Fill( cmsdy );

		  if( oddcol )
		    cmsdyfcq4oddHisto.Fill( cmsdy );
		  else
		    cmsdyfcq4eveHisto.Fill( cmsdy );

		}

	      }

	    } // Q0 > 17

	    if( lq ) {
	      cmsdyvsx.Fill( x4, cmsdy );
	      cmsdyvsy.Fill( y4, cmsdy );
	      cmsdyvsty.Fill( sly, cmsdy );
	    } // lq

	    // depth charge profile:

	    for( int irow = rowmin; irow <= rowmax; ++irow ) {

	      double q = qrow[irow];
	      double qx = exp(-q/qwid);

	      double yrow = ( irow + 0.5 - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	      double depth = ( y4 - yrow ) / tan(DUTtilt*wt); // [-thick/2,thick/2]

	      cmsqrowvsd.Fill( depth, qx );
	      if( irow > rowmin && irow < rowmax ) {
		cmsqrowvsdc.Fill( depth, qx );
		if( depth < 0 )
		  cmsqrownegHisto.Fill( q );
		else
		  cmsqrowposHisto.Fill( q );
	      }

	    } // rows

	  } // dx

	  // cut in x and y:

	  if( fabs(cmsdx) < xcut  && fabs(cmsdy) < ycut ) {

	    cmsq0fHisto.Fill( Q0 ); // Landau

	    if( ldot ) { // sensor bias dot
	      cmsq0fdotHisto.Fill( Q0 );
	    }
	    else { // no dot

	      cmsq0fnodHisto.Fill( Q0 );

	      if( oddcol )
		cmsq0foddHisto.Fill( Q0 ); // 21.84
	      else
		cmsq0feveHisto.Fill( Q0 ); // 21.51

	      if( lcore ) { // pixel core
		cmsq0fcHisto.Fill( Q0 );
	      }
	      else {
		cmsq0fpHisto.Fill( Q0 );
	      } // not core

	    } // no dot

	    // cluster charge profiles, exponential weighting Qx:

	    cmsqxvst1.Fill( evsec, Qx );
	    cmsqxvst2.Fill( evsec, Qx );
	    cmsqxvst3.Fill( evsec, Qx );
	    cmsqxvsx.Fill( x4, Qx );
	    cmsqxvsy.Fill( y4, Qx );
	    cmsqxvsxm.Fill( xmod, Qx ); // Q within pixel
	    cmsqxvsym.Fill( ymod, Qx ); // Q within pixel

	    if( ! ldot )
	      cmsqxvsymn.Fill( ymod, Qx ); // q within pixel

	    if( ( xmod > 120 && xmod < 140 ) ||
		( xmod > 160 && xmod < 180 ) ) // dot
	      cmsqxvsymd.Fill( ymod, Qx ); // q within pixel

	    if( ( ymod >  60 && ymod <  80 ) ||
		( ymod > 160 && ymod < 180 ) ) // dot
	      cmsqxvsxmd.Fill( xmod, Qx ); // q within pixel

	    cmsqxvsxmym->Fill( xmod, ymod, Qx ); // cluster charge profile

	    if( lq ) {
	      cmsnpxHisto.Fill( c->size );
	      cmsncolHisto.Fill( ncol );
	      cmsnrowHisto.Fill( nrow );
	    }
	    cmsnpxvsq.Fill( Q0, c->size ); // rising

	    if( nrow > 1 ) {
	      cmsrmsrowHisto.Fill( rmsrow );
	      cmsskwrowHisto.Fill( skwrow );
	      cmsrmsrowvsq0.Fill( Q0, rmsrow );
	      cmsrmsrowvsym.Fill( ymod, rmsrow/nrow ); // within pixel
	      cmsskwrowvsym.Fill( ymod, skwrow ); // within pixel
	    }

	    for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	      cmspxqHisto.Fill( px->q );
	      if( oddcol )
		cmspxqoddHisto.Fill( px->q );
	      else
		cmspxqeveHisto.Fill( px->q );
	      /*
	      if( px->ord == 0 )
		cmspxq1stHisto.Fill( px->q );
	      else
		cmspxq2ndHisto.Fill( px->q );
	      */
	      if( c->size == 1 ) cmspxq1Histo.Fill( px->q );
	      if( c->size == 2 ) cmspxq2Histo.Fill( px->q ); // flat
	      if( c->size >= 3 ) cmspxq3Histo.Fill( px->q ); // hump at low q
	      cmspxqvsq.Fill( Q0, px->q );
	      cmspxqvsxm.Fill( xmod, px->q );
	      cmspxqvsym.Fill( ymod, px->q );

	    } // pix

	    if( lq ) { // Landau peak

	      cmsdxvsxm.Fill( xmod, cmsdx );
	      cmsdxvsym.Fill( ymod, cmsdx );
	      cmsdyvsxm.Fill( xmod, cmsdy );
	      if( ! ldot ) {
		cmsdyvsym.Fill( ymod, cmsdy );
		if( oddcol ) 
		  cmsdyvsymodd.Fill( ymod, cmsdy );
		else
		  cmsdyvsymeve.Fill( ymod, cmsdy );
	      }

	      if( nrow > 1 ) {
		cmsdyvsrmsrow.Fill( rmsrow/ncol, cmsdy );
		cmsrmsyvsrmsrow.Fill( rmsrow/ncol, fabs(cmsdy) );
		cmsdyvsskwrow.Fill( skwrow, cmsdy );
		cmsdy0vsskwrow.Fill( skwrow, cmsdy0 ); // before correction
		cmsskwrowvsdy.Fill( cmsdy*1E3, skwrow );
		cmsskwrowvsdy0.Fill( cmsdy0*1E3, skwrow ); // before correction
		cmsrmsyvsskwrow.Fill( skwrow, fabs(cmsdy) );
	      }

	      // resolution profiles:

	      cmsrmsxvsx.Fill( x4, fabs(cmsdx) ); // resolution across cols
	      cmsrmsyvsx.Fill( x4, fabs(cmsdy) ); // resolution across cols
	      cmsrmsxvsy.Fill( y4, fabs(cmsdx) ); // resolution across rows
	      cmsrmsyvsy.Fill( y4, fabs(cmsdy) ); // resolution across rows

	      cmsrmsyvstx.Fill( slx, fabs(cmsdy) );
	      cmsrmsyvsty.Fill( sly, fabs(cmsdy) );

	      cmsrmsxvsxm.Fill( xmod, fabs(cmsdx) ); // resolution within pixel
	      cmsrmsyvsxm.Fill( xmod, fabs(cmsdy) ); // resolution within pixel
	      if( !ldot ) {
		cmsrmsxvsym.Fill( ymod, fabs(cmsdx) ); // resolution within pixel
		cmsrmsyvsym.Fill( ymod, fabs(cmsdy) ); // resolution within pixel
	      }
	      cmsrmsyvsxmym->Fill( xmod, ymod, fabs(cmsdy) ); // resolution within pixel
	      cmsrmsxvsxmym->Fill( xmod, ymod, fabs(cmsdx) ); // resolution within pixel

	      // pump oscillation ?

	      cmsdyvst1.Fill( evsec, cmsdy );
	      cmsdyvst2.Fill( evsec, cmsdy );
	      cmsrmsyvst1.Fill( evsec, fabs(cmsdy) );
	      cmsrmsyvst2.Fill( evsec, fabs(cmsdy) );

	      // cluster size profiles:

	      cmsncolvsx.Fill( x4, ncol ); // no trend
	      cmsncolvsy.Fill( y4, ncol ); // 

	      cmsncolvsxm.Fill( xmod, ncol );
	      cmsnrowvsxm.Fill( xmod, nrow );
	      if( ( xmod > 0.020 && xmod < 0.130 ) ||
		  ( xmod > 0.170 && xmod < 0.280 ) ) {
		cmsncolvsym.Fill( ymod, ncol ); // within pixel
		cmsnrowvsym.Fill( ymod, nrow ); // within pixel
	      }
	      cmsnpxvsxmym->Fill( xmod, ymod, c->size ); // cluster size map

	    } // q Landau peak

	  } // dx and dy

	} // fiducial and iso

	// wide linking cut in x and y for efficiency:

	double dxy = sqrt( cmsdx*cmsdx + cmsdy*cmsdy );
	//if( fabs( cmsdx ) < xcut && fabs( cmsdy ) < ycut ) {
	//if( dxy < 0.6 ) {
	if( dxy < 0.9 ) {
	  nm = 1;
	  trixylkHisto->Fill( xA, yA );
	}
	if( dxy < dmin ) {
	  dmin = dxy;
	  clQ0 = Q0;
	  dxmin = cmsdx;
	  dymin = cmsdy;
	}

      } // loop DUT clusters

      // DUT efficiency vs isolated REF-linked fiducial tracks:

      if( lsixlk
	  //&& cl0[iDUT].size() < 2 // empty or single cluster, same eff
	  ) {

	double fidx0 =-3.9;
	double fidx9 = 3.9;
	double fidy0 =-3.9;
	double fidy9 = 3.9;
	if( run == 23311 ) { // 504, 34M
	  fidx0 =-3.4;
	  fidy9 = 3.8;
	}
	if( run >= 24421 && run <= 24429 ) { // 603 Vdig 2
	  fidy0 =-3.0;
	}

	if( x4 > fidx0 && x4 < fidx9 &&
	    y4 > fidy0 && y4 < fidy9 ) { // fiducial
	  effvsdmin.Fill( dddmin, nm ); // at REF, small effect
	  if( dddmin > 0.4 )
	    effvstmin.Fill( ttdmin, nm ); // at DUT, flat
	}

	//if( dddmin > 0.4 && ttdmin > 1.4 ) { // stronger iso [mm] 99.94
	//if( dddmin > 0.4 && ttdmin > 0.6 ) { // iso [mm] 99.93
	if( dddmin > 0.6 ) { // iso [mm] 99.94

	  sixxylkHisto->Fill( xA, yA );
	  if( nm ) sixxyeffHisto->Fill( xA, yA );

	  effvsxy->Fill( x4, y4, nm );

	  if( y4 > fidy0 && y4 < fidy9 )
	    effvsx.Fill( x4, nm );

	  if( x4 > fidx0 && x4 < fidx9 )
	    effvsy.Fill( y4, nm );

	  if( x4 > fidx0 && x4 < fidx9 &&
	      y4 > fidy0 && y4 < fidy9 ) { // fiducial

	    effdminHisto.Fill( dmin );
	    if( nm == 0 ) {
	      effdmin0Histo.Fill( dmin );
	      effrxmin0Histo.Fill( dxmin/dmin ); // 0
	      effrymin0Histo.Fill( dymin/dmin ); // +-1
	      effdxmin0Histo.Fill( dxmin );
	      effdymin0Histo.Fill( dymin );
	    }
	    effclq0Histo.Fill( clQ0 );
	    if( ndrilk == 1 ) effclqrHisto.Fill( clQ0 );
	    if( dmin > 0.1 ) effclq1Histo.Fill( clQ0 );
	    if( dmin > 0.2 ) effclq2Histo.Fill( clQ0 );
	    if( dmin > 0.3 ) effclq3Histo.Fill( clQ0 );
	    if( dmin > 0.4 ) effclq4Histo.Fill( clQ0 );
	    if( dmin > 0.5 ) effclq5Histo.Fill( clQ0 );
	    if( dmin > 0.6 ) effclq6Histo.Fill( clQ0 );
	    if( dmin > 0.7 ) effclq7Histo.Fill( clQ0 );
	    if( dmin > 0.8 ) effclq8Histo.Fill( clQ0 );
	    if( dmin > 0.9 ) effclq9Histo.Fill( clQ0 );
	    effvst1.Fill( evsec, nm );
	    effvst2.Fill( evsec, nm );
	    effvst3.Fill( evsec, nm );
	    effvst4.Fill( evsec, nm );
	    effvsxt->Fill( evsec, x4, nm );
	    effvsntri.Fill( triplets.size(), nm );
	    effvsndri.Fill( dripletsForEff.size(), nm );
	    effvsxmym->Fill( xmod, ymod, nm );
	    effvsxm.Fill( xmod, nm );
	    effvstx.Fill( slx, nm );
	    effvsty.Fill( sly, nm );
	    effvstxy.Fill( txy, nm ); // no effect
	    effvsdslp.Fill( sixdslp, nm ); // no effect

	  } // fiducial

	} // iso

      } // six

    } // loop triplets iA

    ++event_nr;

    if( syncdut )
      cl0[iDUT] = cl[iDUT]; // remember for re-sync
    if( syncref )
      cl0[iREF] = cl[iREF];
    if( syncmod )
      cl0[iMOD] = cl[iMOD];

  } while( reader->NextEvent() && event_nr < lev );

  cout << "done after " << event_nr << " events" << endl;
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
      fgp0->SetParameter( 0, moddxaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, moddxaHisto.GetBinCenter( moddxaHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, moddxaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddxaHisto.GetBinContent(1) ); // BG
      moddxaHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODalignx = fgp0->GetParameter(1);
    }
    else {
      cout << endl << modsxaHisto.GetTitle()
	   << " bin " << modsxaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, modsxaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, modsxaHisto.GetBinCenter( modsxaHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, modsxaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, modsxaHisto.GetBinContent(1) ); // BG
      modsxaHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODalignx = fgp0->GetParameter(1);
    }

    if( moddyaHisto.GetMaximum() > modsyaHisto.GetMaximum() ) {
      cout << endl << moddyaHisto.GetTitle()
	   << " bin " << moddyaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, moddyaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, moddyaHisto.GetBinCenter( moddyaHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, moddyaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddyaHisto.GetBinContent(1) ); // BG
      moddyaHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODaligny = fgp0->GetParameter(1);
    }
    else {
      cout << endl << modsyaHisto.GetTitle()
	   << " bin " << modsyaHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, modsyaHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, modsyaHisto.GetBinCenter( modsyaHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, modsyaHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, modsyaHisto.GetBinContent(1) ); // BG
      modsyaHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODaligny = fgp0->GetParameter(1);
    }

    // finer alignment:

    if( MODaligniteration > 0 && fabs( newMODalignx - MODalignx ) < 0.1 ) {

      cout << endl << moddxHisto.GetTitle()
	   << " bin " << moddxHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, moddxHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, moddxHisto.GetBinCenter( moddxHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 8*moddxHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddxHisto.GetBinContent(1) ); // BG
      moddxHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODalignx = MODalignx + fgp0->GetParameter(1);

    }

    if( MODaligniteration > 0 && fabs( newMODaligny - MODaligny ) < 0.1 ) {

      cout << endl << moddyHisto.GetTitle()
	   << " bin " << moddyHisto.GetBinWidth(1)
	   << endl;
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, moddyHisto.GetMaximum() ); // amplitude
      fgp0->SetParameter( 1, moddyHisto.GetBinCenter( moddyHisto.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 5*moddyHisto.GetBinWidth(1) ); // sigma
      fgp0->SetParameter( 3, moddyHisto.GetBinContent(1) ); // BG
      moddyHisto.Fit( "fgp0", "q" );
      cout << "Fit Gauss + BG:"
	   << endl << "  A " << fgp0->GetParameter(0)
	   << endl << "mid " << fgp0->GetParameter(1)
	   << endl << "sig " << fgp0->GetParameter(2)
	   << endl << " BG " << fgp0->GetParameter(3)
	   << endl;
      newMODaligny = MODaligny + fgp0->GetParameter(1);

      // dyvsx -> rot

      moddyvsx.Fit( "pol1", "q", "", -midx[iMOD]+0.2, midx[iMOD]-0.2 );
      TF1 * fdyvsx = moddyvsx.GetFunction( "pol1" );
      cout << endl << moddyvsx.GetTitle()
	   << ": extra rot " << fdyvsx->GetParameter(1) << endl;
      MODrot += fdyvsx->GetParameter(1);

      // dxvsx -> turn:

      moddxvsx.Fit( "pol1", "q", "", -midx[iMOD]+0.2, midx[iMOD]-0.2 );
      TF1 * fdxvsx = moddxvsx.GetFunction( "pol1" );
      cout << endl << moddxvsx.GetTitle()
	   << ": slope " << fdxvsx->GetParameter(1)
	   << ", extra turn " << fdxvsx->GetParameter(1)/wt/max(som,0.01)
	   << " deg"
	   << endl;
      MODturn += fdxvsx->GetParameter(1)/wt/max(som,0.01); // [deg] min 0.6 deg

      // dyvsy -> tilt:

      moddyvsy.Fit( "pol1", "q", "", -midy[iMOD]+0.2, midy[iMOD]-0.2 );
      TF1 * fdyvsy = moddyvsy.GetFunction( "pol1" );
      cout << endl << moddyvsy.GetTitle()
	   << ": slope " << fdyvsy->GetParameter(1)
	   << ", extra tilt " << fdyvsy->GetParameter(1)/wt/max(sam,0.01)
	   << " deg"
	   << endl;
      MODtilt += fdyvsy->GetParameter(1)/wt/max(sam,0.01); // [deg] min 0.6 deg

      // dyvsty -> dz:

      moddyvsty.Fit( "pol1", "q", "", -0.002, 0.002 );
      TF1 * fdyvsty = moddyvsty.GetFunction( "pol1" );
      cout << endl << moddyvsty.GetTitle()
	   << ": z shift " << fdyvsty->GetParameter(1)
	   << " mm (add to MODz)"
	   << endl;
      MODz += fdyvsty->GetParameter(1);
    }

    // write new MOD alignment:

    ofstream MODalignFile( outputDirectory+"/"+MODalignFileName.str() );
    
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
	 << " to " << MODalignFileName.str()
	 << endl;

  } // MOD

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // REF alignment:

  double newREFalignx = REFalignx;
  double newREFaligny = REFaligny;

  if( refdxaHisto.GetMaximum() > refsxaHisto.GetMaximum() ) {
    cout << endl << refdxaHisto.GetTitle()
	 << " bin " << refdxaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refdxaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refdxaHisto.GetBinCenter( refdxaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, refdxaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refdxaHisto.GetBinContent(1) ); // BG
    refdxaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFalignx = fgp0->GetParameter(1);
  }
  else {
    cout << endl << refsxaHisto.GetTitle()
	 << " bin " << refsxaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refsxaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refsxaHisto.GetBinCenter( refsxaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, refsxaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refsxaHisto.GetBinContent(1) ); // BG
    refsxaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFalignx = fgp0->GetParameter(1);
  }

  if( refdyaHisto.GetMaximum() > refsyaHisto.GetMaximum() ) {
    cout << endl << refdyaHisto.GetTitle()
	 << " bin " << refdyaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refdyaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refdyaHisto.GetBinCenter( refdyaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, refdyaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refdyaHisto.GetBinContent(1) ); // BG
    refdyaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFaligny = fgp0->GetParameter(1);
  }
  else {
    cout << endl << refsyaHisto.GetTitle()
	 << " bin " << refsyaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refsyaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refsyaHisto.GetBinCenter( refsyaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, refsyaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refsyaHisto.GetBinContent(1) ); // BG
    refsyaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFaligny = fgp0->GetParameter(1);
  }

  // finer alignment:

  if( REFaligniteration > 0 && fabs( newREFalignx - REFalignx ) < 0.1 ) {

    cout << endl << refdxHisto.GetTitle()
	 << " bin " << refdxHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refdxHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refdxHisto.GetBinCenter( refdxHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 8*refdxHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refdxHisto.GetBinContent(1) ); // BG
    refdxHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFalignx = REFalignx + fgp0->GetParameter(1);

  }

  if( REFaligniteration > 0 && fabs( newREFaligny - REFaligny ) < 0.1 ) {

    cout << endl << refdyHisto.GetTitle()
	 << " bin " << refdyHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refdyHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refdyHisto.GetBinCenter( refdyHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 5*refdyHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refdyHisto.GetBinContent(1) ); // BG
    refdyHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFaligny = REFaligny + fgp0->GetParameter(1);

    // dyvsx -> rot

    refdyvsx.Fit( "pol1", "q", "", -midx[iREF]+0.2, midx[iREF]-0.2 );
    TF1 * fdyvsx = refdyvsx.GetFunction( "pol1" );
    cout << endl << refdyvsx.GetTitle()
	 << ": extra rot " << fdyvsx->GetParameter(1) << endl;
    REFrot += fdyvsx->GetParameter(1);

    // dyvsty -> dz:

    refdyvsty.Fit( "pol1", "q", "", -0.002, 0.002 );
    TF1 * fdyvsty = refdyvsty.GetFunction( "pol1" );
    cout << endl << refdyvsty.GetTitle()
	 << ": z shift " << fdyvsty->GetParameter(1)
	 << " mm (add to REFz)"
	 << endl;
    REFz += fdyvsty->GetParameter(1);
  }

  // write new REF alignment:

  ofstream REFalignFile( outputDirectory+"/"+REFalignFileName.str() );

  REFalignFile << "# REF alignment for run " << run << endl;
  ++REFaligniteration;
  REFalignFile << "iteration " << REFaligniteration << endl;
  REFalignFile << "alignx " << newREFalignx << endl;
  REFalignFile << "aligny " << newREFaligny << endl;
  REFalignFile << "rot " << REFrot << endl;
  REFalignFile << "dz " << REFz - zz[5] << endl;

  REFalignFile.close();

  cout << endl << "wrote REF alignment iteration " << REFaligniteration
       << " to " << REFalignFileName.str()
       << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT alignment:

  double newDUTalignx = DUTalignx;
  double newDUTaligny = DUTaligny;

  if( cmsdxaHisto.GetMaximum() > cmssxaHisto.GetMaximum() ) {
    cout << endl << cmsdxaHisto.GetTitle()
	 << " bin " << cmsdxaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -10, 10 );
    double xpk = cmsdxaHisto.GetBinCenter( cmsdxaHisto.GetMaximumBin() );
    fgp0->SetParameter( 0, cmsdxaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, xpk );
    fgp0->SetParameter( 2, cmsdxaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, cmsdxaHisto.GetBinContent(1) ); // BG
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
    fgp0->SetParameter( 3, cmssxaHisto.GetBinContent(1) ); // BG
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
    fgp0->SetParameter( 3, cmsdyaHisto.GetBinContent(1) ); // BG
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
    fgp0->SetParameter( 3, cmssyaHisto.GetBinContent(1) ); // BG
    cmssyaHisto.Fit( "fgp0", "q", "", xpk-1, xpk+1 );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newDUTaligny = fgp0->GetParameter(1);
  }

  // finer alignment:

  if( DUTaligniteration > 0 && fabs( newDUTalignx - DUTalignx ) < 0.1 ) {

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

  }

  if( DUTaligniteration > 0 && fabs( newDUTaligny - DUTaligny ) < 0.1 ) {

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

    // dyvsx -> rot

    cmsdyvsx.Fit( "pol1", "q", "", -midx[iDUT]+0.2, midx[iDUT]-0.2 );
    TF1 * fdyvsx = cmsdyvsx.GetFunction( "pol1" );
    cout << endl << cmsdyvsx.GetTitle()
	 << ": extra rot " << fdyvsx->GetParameter(1) << endl;
    DUTrot += upsigny*fdyvsx->GetParameter(1);

    // dyvsy -> tilt:

    cmsdyvsy.Fit( "pol1", "q", "", -midy[iDUT]+0.2, midy[iDUT]-0.2 );
    TF1 * fdyvsy = cmsdyvsy.GetFunction( "pol1" );
    cout << endl << cmsdyvsy.GetTitle()
	 << ": slope " << fdyvsy->GetParameter(1)
	 << ", extra tilt " << fdyvsy->GetParameter(1)/wt/max(sa,0.01)
	 << " deg"
	 << endl;
    DUTtilt += fdyvsy->GetParameter(1)/wt/max(sa,0.01); // [deg] min 0.6 deg

    // dyvsty -> dz:

    cmsdyvsty.Fit( "pol1", "q", "", -0.002, 0.002 );
    TF1 * fdyvsty = cmsdyvsty.GetFunction( "pol1" );
    cout << endl << cmsdyvsty.GetTitle()
	 << ": z shift " << fdyvsty->GetParameter(1)
	 << " mm"
	 << endl;
    DUTz -= upsigny*fdyvsty->GetParameter(1);
  }

  // write new DUT alignment:

  ofstream DUTalignFile( outputDirectory+"/"+DUTalignFileName.str() );

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
       << " to " << DUTalignFileName.str()
       << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  cout << endl << histoFile->GetName() << endl;

  if(useMODasREF)  cout << "Using MOD for efficiency " << endl;
  else cout << "Using REF for efficiency " << endl;
  cout << endl;

  return 0;
}

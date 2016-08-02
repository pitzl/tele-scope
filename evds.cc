
// Daniel Pitzl, DESY Jul 2016
// event display 4 module planes, B-field, sagitta A C D

// evds -l 19 2321 # 2.0 GeV, 1.4 T, sagitta spacing

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TApplication.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TStyle.h> // gStyle
#include <TFile.h>
#include <TH2D.h>
#include <TF1.h>

#include <cstdlib> // atoi, drand48
#include <iostream> // cout
#include <iomanip> // setw
#include <string> // strings
#include <sstream> // stringstream
#include <fstream> // files
#include <vector>
#include <unistd.h> // usleep

//------------------------------------------------------------------------------
class MyMainFrame:public TGMainFrame
{
private:
  TGMainFrame * fMain;
  TRootEmbeddedCanvas *fEcanvas;
public:
  MyMainFrame( const TGWindow * p, UInt_t w, UInt_t h );
  //virtual ~ MyMainFrame(  );
  ~MyMainFrame(  );
  TCanvas *GetCanvas(  );
};

//------------------------------------------------------------------------------
using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int adc;
  double cal;
};

struct cluster {
  vector <pixel> vpix;
  int size;
  int sumA;
  double charge;
  double col,row;
  bool bigx, bigy;
  double x5;
  double y5;
  double z3;
};

//------------------------------------------------------------------------------
// globals:

pixel pb[66560]; // vector of pixels with hit, 16*4160 = 66560
int fNHit; // used in getClus

//------------------------------------------------------------------------------
MyMainFrame::MyMainFrame( const TGWindow * p, UInt_t w, UInt_t h )
  :TGMainFrame( p, w, h )
{
  cout << "MyMainFrame..." << endl;
  // Create a main frame:
  fMain = new TGMainFrame( p, w, h );

  fMain->SetWMPosition( 99, 0 ); // no effect

  // Create canvas widget:
  fEcanvas = new TRootEmbeddedCanvas( "Ecanvas", fMain, w, h );

  fMain->AddFrame( fEcanvas,
                   new TGLayoutHints( kLHintsExpandX | kLHintsExpandY, 1, 1, 1, 1 ) );

  // Set a name to the main frame:
  fMain->SetWindowName( "4-plane event display" );

  // Map all subwindows of main frame:
  fMain->MapSubwindows(  );

  // Initialize the layout algorithm:
  fMain->Resize( fMain->GetDefaultSize(  ) );

  // Map main frame:
  fMain->MapWindow(  );
}

MyMainFrame::~MyMainFrame(  )
{
  // Clean up used widgets: frames, buttons, layouthints
  fMain->Cleanup(  );
  cout << "MyMainFrame: Cleanup" << endl;
  delete fMain;
  //delete fEcanvas; // crash
}

TCanvas *MyMainFrame::GetCanvas(  )
{
  return ( fEcanvas->GetCanvas(  ) );
}

// ----------------------------------------------------------------------
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

    // added all I could. determine position and append it to the list o f clusters:

    c.sumA = 0;
    c.charge = 0;
    c.size = 0;
    c.col = 0;
    c.row = 0;
    double sumQ = 0;
    c.bigx = 0;
    c.bigy = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      c.sumA += p->adc; // Aout
      double Qpix = p->cal; // calibrated [Vcal]
      if( Qpix < 0 ) Qpix = 1; // DP 1.7.2012
      c.charge += Qpix;
      sumQ += Qpix;
      c.col += (*p).col*Qpix;
      c.row += (*p).row*Qpix;
      if( p->col ==  0 ) c.bigx = 1;
      if( p->col == 51 ) c.bigx = 1;
      if( p->row == 79 ) c.bigy = 1;
    }

    c.size = c.vpix.size();

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( ! c.charge == 0 ) {
      c.col /= sumQ;
      c.row /= sumQ;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with zero charge" << endl;
    }

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

  FileReader reader = FileReader( runnum.c_str(), "data/run$6R$X" );

  // further arguments:

  int lev = 99;

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] );

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  TFile* histoFile = new TFile( "evdB.root", "RECREATE" );

  cout << "ROOT application..." << endl;

  TApplication theApp( "comet", &argc, argv );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // set ROOT styles:

  gStyle->SetTextFont( 62 ); // 62 = Helvetica bold
  gStyle->SetTextAlign( 11 );

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.01, "y" );
  gStyle->SetTickLength( -0.01, "z" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.013, "y" );
  gStyle->SetLabelOffset( 0.022, "z" );

  gStyle->SetTitleOffset( 1.6, "x" );
  gStyle->SetTitleOffset( 1.6, "y" );
  gStyle->SetTitleOffset( 1.7, "z" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );
  gStyle->SetLabelFont( 62, "z" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );
  gStyle->SetTitleFont( 62, "z" );

  gStyle->SetTitleBorderSize( 0 ); // no frame around global title
  gStyle->SetTitleAlign( 13 ); // 13 = left top align
  gStyle->SetTitleX( 0.12 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetLineWidth( 1 ); // frames
  gStyle->SetHistLineColor( 4 ); // 4=blau
  gStyle->SetHistLineWidth( 3 );
  gStyle->SetHistFillColor( 5 ); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle( 1001 ); // 1001 = solid

  gStyle->SetFrameLineWidth( 2 );

  // statistics box:

  gStyle->SetOptStat( 10 );
  gStyle->SetStatFormat( "8.6g" ); // more digits, default is 6.4g
  gStyle->SetStatFont( 42 ); // 42 = Helvetica normal
  //  gStyle->SetStatFont(62); // 62 = Helvetica bold
  gStyle->SetStatBorderSize( 1 ); // no 'shadow'

  gStyle->SetStatX( 0.80 );
  gStyle->SetStatY( 0.95 );

  gStyle->SetPalette( 55 ); // rainbow colors

  gStyle->SetHistMinimumZero(  ); // no zero suppression

  gStyle->SetOptDate( 0 );

  cout << "open ROOT window..." << endl;
  MyMainFrame *myMF = new
    MyMainFrame( gClient->GetRoot(  ), 1400, 700 ); // 4 planes

  cout << "open Canvas..." << endl;
  TCanvas * c1 = myMF->GetCanvas(  );

  c1->SetBottomMargin( 0.11 );
  c1->SetLeftMargin( 0.06 );
  c1->SetRightMargin( 0.13 );

  gPad->Update(  ); // required

  //------------------------------------------------------------------------------
  // geometry:

  double pi = 4*atan(1);
  double wt = 180/pi;

  double zpos[4] = { -48, -16, 16, 48 }; // 32 mm spacing

  if( run >= 2308 && run < 2365 ) { // Sagitta spacing
    zpos[0] = -16 - 103;
    zpos[1] = -16;
    zpos[2] =  16;
    zpos[3] =  16 + 103;
  }

  double tricuty = 0.1; // [mm] aligned

  double tilt = 19.3; // always

  double turn = 27.7;
  if( run <= 2154 )
    turn = 0;
  if( run >= 2368 ) // showers
    turn = 0;

  double costilt = cos(tilt/wt);
  double sintilt = sin(tilt/wt);
  double costurn = cos(turn/wt);
  double sinturn = sin(turn/wt);

  double phi = 0;
  if( run >= 2187 ) phi =  -6;
  if( run >= 2192 ) phi =  -5;
  if( run >= 2199 ) phi =  -2;
  if( run >= 2203 ) phi = -15;
  if( run >= 2205 ) phi =  -2;
  if( run >= 2209 ) phi =  -4;
  if( run >= 2213 ) phi =  -3;
  if( run >= 2217 ) phi =  -2;
  if( run >= 2221 ) phi =  -1;
  if( run >= 2225 ) phi =   0;
  if( run >= 2229 ) phi =   1;
  if( run >= 2235 ) phi =   2;
  if( run >= 2239 ) phi =   3;
  if( run >= 2243 ) phi =   4;
  if( run >= 2250 ) phi = -15;
  if( run >= 2259 ) phi =  -6;
  if( run >= 2266 ) phi =  -2;
  if( run >= 2308 ) phi =  -4;
  if( run >= 2314 ) phi =  -6;
  if( run >= 2317 ) phi =   0;
  if( run >= 2321 ) phi =  -6;
  if( run >= 2326 ) phi =  -3;
  if( run >= 2331 ) phi =  -2;
  if( run >= 2336 ) phi = -10;
  if( run >= 2337 ) phi = -12;
  if( run >= 2338 ) phi = -15;
  if( run >= 2341 ) phi =  -2;
  if( run >= 2347 ) phi =  -2.2;
  if( run >= 2368 ) phi =  -2;

  double cosphi = cos(phi/wt);
  double sinphi = sin(phi/wt);

  // B field:

  double T = 1.35; // [T]

  //------------------------------------------------------------------------------
  // alignments:

  const int A = 0;
  const int B = 1;
  const int C = 2;
  const int D = 3;

  int aligniteration = 0;
  double alignx[4];
  double aligny[4];
  double fx[4];
  double fy[4];
  double tx[4];
  double p2x[4];
  double p3x[4];
  double ty[4];
  double p2y[4];
  double p3y[4];

  for( int ipl = 0; ipl < 4; ++ipl ) {
    alignx[ipl] = 0;
    aligny[ipl] = 0;
    fx[ipl] = 0;
    fy[ipl] = 0;
    tx[ipl] = 0;
    p2x[ipl] = 0;
    p3x[ipl] = 0;
    ty[ipl] = 0;
    p2y[ipl] = 0;
    p3y[ipl] = 0;
  }

  ostringstream alignFileName; // output string stream

  alignFileName << "align3D_" << run << ".dat";

  ifstream ialignFile( alignFileName.str() );

  cout << endl;
  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    cout << "no " << alignFileName.str() << ", will bootstrap" << endl;
    cout << endl;
  }
  else {

    cout << "read alignment from " << alignFileName.str() << endl;

    string Hash( "#" );
    string Iteration( "iteration" );
    string Plane( "plane" );
    string Alignx( "alignx" );
    string Aligny( "aligny" );
    string Rx( "fx" );
    string Ry( "fy" );
    string Tx( "tx" );
    string P2x( "p2x" );
    string P3x( "p3x" );
    string Ty( "ty" );
    string P2y( "p2y" );
    string P3y( "p3y" );

    int ipl = 0;

    while( ! ialignFile.eof() ) {

      string line;
      getline( ialignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == Hash ) // comments start with #
	continue;

      if( tag == Iteration ) 
	tokenizer >> aligniteration;

      if( tag == Plane )
	tokenizer >> ipl;

      if( ipl < 0 || ipl >= 4 ) {
	//cout << "wrong plane number " << ipl << endl;
	continue;
      }

      double val;
      tokenizer >> val;
      if(      tag == Alignx )
	alignx[ipl] = val;
      else if( tag == Aligny )
	aligny[ipl] = val;
      else if( tag == Rx )
	fx[ipl] = val;
      else if( tag == Ry )
	fy[ipl] = val;
      else if( tag == Tx )
	tx[ipl] = val;
      else if( tag == P2x )
	p2x[ipl] = val;
      else if( tag == P3x )
	p3x[ipl] = val;
      else if( tag == Ty )
	ty[ipl] = val;
      else if( tag == P2y )
	p2y[ipl] = val;
      else if( tag == P3y )
	p3y[ipl] = val;

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  ialignFile.close();

  //------------------------------------------------------------------------------
  //event loop:

  size_t nev = 0;
  int kev = 0;

  do {
    // Get next event:
    DetectorEvent evt = reader.GetDetectorEvent();

    if( evt.IsBORE() )
      eudaq::PluginManager::Initialize(evt);

    bool ldbg = 0;

    //if( nev%1000 == 0 ) ldbg = 1;

    nev++;

    if( ldbg ) cout << "Processing event " << nev << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    vector <cluster> cl[4];

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldbg ) cout << "PLANE " << plane.ID() << ": ";

      int mod = 0; // QUAD
      if( plane.ID() == 6 ) mod = 1; // TRP
      if( plane.ID() == 7 ) mod = 2; // DUT
      if( plane.ID() == 8 ) mod = 3; // REF

      int npx = 0;

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldbg ) 
	  std::cout << plane.GetX(ipix)
		    << " " << plane.GetY(ipix)
		    << " " << plane.GetPixel(ipix) << " ";

	int xm = plane.GetX(ipix); // global column 0..415
	int ym = plane.GetY(ipix); // global row 0..159
	int adc = plane.GetPixel(ipix); // ADC 0..255

	int roc = xm / 52; // 0..7

	// leave space for big pixels:

	int x = 1 + xm + 2*roc; // 1..52 per ROC
	int y = ym;
	if( ym > 79 ) y += 2;

	double cal = adc;

	// fill pixel block for clustering:
	pb[npx].col = x;
	pb[npx].row = y;
	pb[npx].adc = adc;
	pb[npx].cal = cal;
	++npx;

	// double big pixels:
	// 0+1
	// 2..51
	// 52+53

	if( xm%52 == 0 ) {
	  pb[npx].col = x-1; // double
	  pb[npx].row = y;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  ++npx;
	}
	else if( xm%52 == 51 ) {
	  pb[npx].col = x+1; // double
	  pb[npx].row = y;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  ++npx;
	}

	if( ym == 79 ) {
	  pb[npx].col = x; // double
	  pb[npx].row = 80;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  ++npx;
	}
	else if( ym == 80 ) {
	  pb[npx].col = x; // double
	  pb[npx].row = 81;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  ++npx;
	}

      } // pix
      
      if( ldbg ) cout << endl;

      // clustering:

      fNHit = npx; // for cluster search

      cl[mod] = getClus();

      if( ldbg ) cout << "clusters mod " << mod << " " << cl[mod].size() << endl;

    } // planes = mod

    if( cl[0].size() < 1 ) continue; // skip event
    if( cl[1].size() < 1 ) continue;
    if( cl[2].size() < 1 ) continue;
    if( cl[3].size() < 1 ) continue;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // book event histos:

    ++kev;

    cout << "run " << run
	 << ", trigger " << nev
	 << ", display " << kev
	 << endl;

    TH2D xview( Form( "ev%ix", kev ),
		Form( "x display %i trigger %i;z [mm];x [mm];y [mm]",
		      kev, (int) nev ),
		100, zpos[0]-21, zpos[3]+21, 120, -60, 60  );
		//		100, zpos[0]-21, zpos[3]+21, 120, -60, 60  );

    xview.GetYaxis()->SetTitleOffset(0.7);
    xview.GetZaxis()->SetTitleOffset(1.2);
    xview.SetMinimum( 0);
    xview.SetMaximum(20);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // local to global:

    for( int mod = 0; mod < 4; ++mod )

      for( vector<cluster>::iterator c = cl[mod].begin(); c != cl[mod].end(); ++c ) {

	// pixel (0,0) is top left
	// right handed coordinate system:
	// x to the right
	// y is up
	// z at you (along beam)

	// passive change of coordinate system (hits remain in space):

	double xm = c->col*0.15 - 32.325;
	double ym =-c->row*0.10 +  8.050; // invert
	double zm = 0;

	// tilt around local x:

	double x1 = xm;
	double y1 = costilt*ym + sintilt*zm;
	double z1 =-sintilt*ym + costilt*zm;

	// turn around y:

	double x2 = costurn*x1 + sinturn*z1;
	double y2 = y1;
	double z2 =-sinturn*x1 + costurn*z1;

	// spread along z:

	z2 += zpos[mod];

	// rotate plate around y :

	double x3 = cosphi*x2 + sinphi*z2;
	double y3 = y2;
	double z3 =-sinphi*x2 + cosphi*z2;

	// align:

	double x4 = x3 - alignx[mod];
	double y4 = y3 - aligny[mod];
	double x5 = x4 - y4*fx[mod] - tx[mod]*x4 - p2x[mod]*x4*x4 - p3x[mod]*x4*x4*x4;
	double y5 = y4 + x4*fy[mod] - ty[mod]*y4 - p2y[mod]*x4*x4 - p3y[mod]*x4*x4*x4;

	c->x5 = x5;
	c->y5 = y5;
	c->z3 = z3;

	xview.Fill( z3, x5, y5+10 );

      } // clusters

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplet ACB in y:

    int nADCy = 0;

    vector <TF1> tracks;

    for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {

      double xA = cA->x5;
      double yA = cA->y5;
      double zA = cA->z3;

      for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {

	double xD = cD->x5;
	double yD = cD->y5;
	double zD = cD->z3;

	// A-D vector:

	double dz2 = zD - zA;

	double txAD = (xD - xA) / dz2; // angle
	double tyAD = (yD - yA) / dz2; // angle

	// needed: intersect track with tilted plane

	for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	  double xC = cC->x5;
	  double yC = cC->y5;
	  double zC = cC->z3;

	  double xintC = xA + txAD * (zC-zA); // interpolate
	  double yintC = yA + tyAD * (zC-zA);

	  // tri ADC:

	  double dx3 = xC - xintC;
	  double dy3 = yC - yintC;

	  if( abs( dy3 ) > tricuty ) continue;

	  ++nADCy;
	  /*
	  cout << "  A " << distance( cl[A].begin(), cA )
	       << ", D " << distance( cl[D].begin(), cD )
	       << ", C " << distance( cl[C].begin(), cC )
	       << ": dy " << dy3
	       << endl;
	  */

	  // parabola through 3 points:

	  TF1 f1( "f1", "[1]+[2]*(x-[0])+[3]*(x-[0])*(x-[0])", zA-5, zD+5 );
	  f1.SetParameter( 0, zC );
	  f1.SetParameter( 1, xC );
	  double det = (zA-zC)*(zD-zC)*(zD-zC) - (zD-zC)*(zA-zC)*(zA-zC);
	  double b = ( (zD-zC)*(zD-zC)*(xA-xC) - (zA-zC)*(zA-zC)*(xD-xC) ) / det; // slope
	  double c = (-(zD-zC)        *(xA-xC) + (zA-zC)        *(xD-xC) ) / det; // curvature
	  f1.SetParameter( 2, b );
	  f1.SetParameter( 3, c );
	  f1.SetLineColor(1);
	  f1.SetLineWidth(1);
	  double R = 0.0005/c; // track radius [m]
	  double p = 0.3*R*T; // [GeV]
	  cout << "  ADC track radius " << R << " m"
	       << ", p " << p
	       << endl;

	  if( fabs(p) > 0.02 )
	    tracks.push_back( f1 );

	} // cl C

      } // cl D

    } // cl A

    cout << "  nADCy " << nADCy
	 << endl;

    if( nADCy == 0 ) continue;

    if( nADCy > 5 ) continue;

    xview.Draw( "colz" );
    //xview.Write();

    for( size_t ii = 0; ii < tracks.size(); ++ii )
      tracks[ii].Draw( "same" );

    xview.Draw( "colsame" ); // once more. on top

    c1->Update();

    //string input;
    //cout << "hit enter";
    //getline( cin, input, '\n' );

    usleep( 1900*1000 ); // sleep some micro-seconds

  } while( reader.NextEvent() && kev < lev );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  histoFile->Write();
  histoFile->Close();
  cout << "write " << histoFile->GetName() << endl;
  delete histoFile;

  delete myMF;

  return 0;
}

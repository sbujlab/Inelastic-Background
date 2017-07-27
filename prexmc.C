#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "math.h"
#include "stdlib.h"
#include "TLine.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TLegend.h"
#include "TRandom.h"
#include "TTree.h"
#include <stdlib.h>
#include <fstream>

//#include "msdist.h"
#include "msdist.cxx"

#include "exstate.h"
#include "leroseres.h"

//class msdist;

TF1 *fspence;

msdist *mydist;

double probextbrem( double E, double t, double Delta, double dmax = -1e9 );
double GenerateRad( double E, double t, double lmin, double lmax );
double raddeltaTsai( double E, double th, double Delta, double Z, double M );
double RadProfile(double E, double t, double loss);
double IRadProfile(double E, double t, double loss1, double loss2);
double raddeltaMY( double E, double th, double Delta );
double dspence( double *x, double *){ return -1.0*log(1.0 - x[0])/x[0]; }
double Spence(double x);

int _nE, _nth;
double _crs[70][150], _crs_de[70][150], _crs_dt[70][150], _crs_dedt[70][150];
double _A[70][150], _A_de[70][150], _A_dt[70][150], _A_dedt[70][150];
double _dA[70][150], _dA_de[70][150], _dA_dt[70][150], _dA_dedt[70][150];
double _thmin, _thmax;
double _Emin, _Emax;
void LoadData( char [], char [] );

double crsMott(double Z, double M, double E, double th );

double crsval( double E, double th );
double   Aval( double E, double th );
double  dAval( double E, double th );

double bcrsval( double E, double th );
double   bAval( double E, double th );

double getX0(double Z, double A);

int prexmc(const char inputfile[] = NULL){
    int i, j;
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    double E, Ef, th, thf, ph, phf, Q2, Q2r, Q2obs;
    double Er, Erf, Erfobs;
    double th0, ph0, th1, ph1, crs, crsr, A, rate, counts, Am;
    double Ar, Arm, q3v, FF;
    double radrate, radcounts;
    double dAdrr, dAmdrr;
    double zt;

    double pi = 3.1415927;
//    double HRSres = 1e-3; // 0.1% hardware momentum resolution
    double HRSres = 0; // Perfect resolution

    double Ebeam, Pe, current, duration;
    int nevt;

    char filename[255], xsfile[255], sxsfile[255];

    double t, Z, An, M, rho, omega;

    FILE *finp;

    if( inputfile == NULL ){
	printf("Opening prexmc.inp\n");
	finp = fopen("prexmc.inp", "r");
    } else {
	printf("Opening %s\n", inputfile);
	finp = fopen(inputfile, "r");
    }
    if( !finp ){ 
	printf("Could not open inputfile\n");
	return 1;
    }

    char dummy[255];

    fscanf(finp, "%s%lf%s", dummy, &Ebeam, dummy );
    fscanf(finp, "%s%lf%s", dummy, &Pe, dummy );
    fscanf(finp, "%s%lf%s", dummy, &current, dummy );
    fscanf(finp, "%s%lf%s", dummy, &duration, dummy );
    // Convert to s
    duration *= 24.0*3600.0;
    fscanf(finp, "%s%d%s\n", dummy, &nevt, dummy );
    fscanf(finp, "%s%lf%s", dummy, &t, dummy );
    fscanf(finp, "%s%lf", dummy, &Z);
    fscanf(finp, "%s%lf", dummy, &An);
    fscanf(finp, "%s%lf%s", dummy, &M, dummy );
    //  Convert to GeV
    M *= 0.9315;
    fscanf(finp, "%s%lf%s", dummy, &rho, dummy );
    fscanf(finp, "%s%s", dummy, xsfile );
    fscanf(finp, "%s%s", dummy, sxsfile );
    fscanf(finp, "%s%s", dummy, filename );
   
    printf("----------------------------------\n");
    printf("%10d events\n", nevt );
    printf("output\t%s\n\n", filename);

    printf("Running with parameters:\n");
    printf("Ebeam\t%5.3f GeV\n", Ebeam);
    printf("Pe   \t%4.2f\n", Pe);
    printf("I    \t%4.1f uA\n", current);
    printf("time \t%4.1f days\n", duration/(24.0*3600.0));
    printf("t    \t%5.3f radlen\n", t);
    printf("Z    \t%3.0f\n", Z);
    printf("A    \t%3.0f\n", An);
    printf("M    \t%6.2f GeV\n", M);
    printf("density\t%6.3f g/cm3\n", rho);
    printf("xsfile\t%s\n", xsfile);
    printf("sxsfile\t%s\n", sxsfile);


    printf("----------------------------------\n");

    LoadCa48ExState();
    LoadLeroseRes();

    /*
    double Ebeam   = 2.2;  // GeV
    double Pe      = 0.85; // Beam polarization
    double current = 50.0; // Beam current, uA
    double duration= 30.0*24.0*3600.0; // s
    int    nevt    = 2e6;

//    char filename[255] = "ca48fom.root";
    char filename[255] = "dummy.root";

    char xsfile[255] = "ca48_fsu.dat";

    double t = 0.05; // Radiation lengths of target
    double Z = 20.0; // Calcium

    double An= 48.0; // Calcium-48
    double M = 48.95*0.9315; // GeV, Calcium-48
    double rho = 1.822; // g/cm^3  Calcium-48
    
	
    //double An= 40.0; // Calcium-40
    //double M = 39.96*0.9315; // GeV, Calcium-40
    //double rho = 1.518; // g/cm^3  Calcium-48
    //char xsfile[255] = "ca40_fsu.dat";
    */

    double X0 = getX0( Z, An ); // g/cm2

    double thick = t*X0; // g/cm2

    printf("Target length = %f cm\n", thick/rho );


    double lumin = // Luminosity
	   current*1.0e-6/1.602e-19  // Current to e-/s
	*  thick*6.022e23/An;        // thickness to nucleons/cm^2


    printf("luminosity = %e Hz/cm2\n", lumin );

    /////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////

    // Assume that the spectrometer is set for Ebeam
    // acceptance is 0.6% below the main peak from PREX
//    double Delta = 0.006*Ebeam; // Energy acceptance of spectrometer
    double Delta = 0.003; // Energy acceptance of the detector

    double dmin = 1e-6; // Minimum radiated photon energy
    double Emin = 1.5; // GeV; Minimum energy in our table
//    double Emin = 0.0; // GeV; Minimum energy in our table
    double dmax;

    if( Emin > Ebeam ){ printf("Decrease the minimum radiative energy\n"); exit(1);}
    
    // true theta generation range
    double cthmax = cos(  1.1*pi/180.0 );
    double cthmin = cos( 10.0*pi/180.0 );
    /*
    double cthmax = cos( 80.0*pi/180.0 );
    double cthmin = cos( 100.0*pi/180.0 );
    */

    double V = (cthmax-cthmin)*2.0*pi;
    
    //  Cross section is in units of mb/str, differential in Omega
    //  Asymmetry is pure
    LoadData(xsfile, sxsfile);
/*
    ofstream ofs;
    ofs.open("output.txt");
    for(int i =0; i<70; i++){
	for(int j=0; j<150;j++){
	ofs << _crs[i][j] << " " << _A[i][j] << "\n";
	}
    }
    ofs.close();
*/
    // Initialize Spence Function
    fspence = new TF1("dspence", dspence, 0.0, 1e9, 0);

    mydist = new msdist(2.200, thick, An, Z);

//    printf("%f\n", exp(raddeltaTsai(2.2, 5.0*pi/180.0, Delta, Z, M)));
//    printf("%f\n", probextbrem(Ebeam, 0.1, Delta));

    gRandom->SetSeed(time(0));

    double rscale1, rscale2, smear;
    double radsupp;

    TTree *T = new TTree("T", "PV Calcium Studies");
    T->Branch("th",     &th,     "th/D");
    T->Branch("ph",     &ph,     "ph/D");
    T->Branch("thf",    &thf,    "thf/D");
    T->Branch("phf",    &phf,    "phf/D");
    T->Branch("crs",    &crs,    "crs/D");
    T->Branch("A",      &A,      "A/D");
    T->Branch("dAdrr",  &dAdrr,   "dAdrr/D");
    T->Branch("Am",     &Am,     "Am/D");
    T->Branch("dAmdrr",  &dAmdrr,   "dAmdrr/D");
    T->Branch("rate",   &rate,   "rate/D");
    T->Branch("counts", &counts, "counts/D");
    T->Branch("radsupp", &radsupp, "radsupp/D");

    T->Branch("Ar",      &Ar,      "Ar/D");
    T->Branch("Arm",     &Arm,     "Arm/D");

    T->Branch("Q2",    &Q2r,    "Q2/D");
    T->Branch("Q2obs",    &Q2obs,    "Q2obs/D");
    T->Branch("q3m",    &q3v,    "q3m/D");
    T->Branch("FF",    &FF,    "FF/D");
    T->Branch("Erf",    &Erf,    "Erf/D");
    T->Branch("Erfobs",    &Erfobs,    "Erfobs/D");
    T->Branch("crsr",    &crsr,    "crsr/D");
    T->Branch("radrate",   &radrate,   "radrate/D");
    T->Branch("radcounts", &radcounts, "radcounts/D");

    T->Branch("state", &j, "state/I");

    double x, y, z, xp, yp, zp;

    double eff_t, dv, Eex;
    smear = 0.0;

    for( i = 0; i < nevt; i++ ){
	if( (i%1000) == 0 ){ printf("Event %10d/%10d\n", i, nevt); }
	// Loop over discrete final states
	for( j = 0; j < nState(Z,An); j++ ){
	    Eex = ExEnergy(Z,An,j);

	    E  = Ebeam;
	    th = acos(gRandom->Uniform()*(cthmax-cthmin)+cthmin);

	    Ef = (E*M - Eex*(M+Eex))/(M + E*(1.0-cos(th)));
	    Q2 =  2.0*E*Ef*(1.0-cos(th));

	    ph = gRandom->Uniform()*2.0*pi;
	    zt = gRandom->Uniform();

	    // Base cross section in total detector acceptance
	    if( j == 0 ){
		crs = crsval(E, th)*exp(raddeltaTsai(E, th, Delta, Z, M))*(1.0-probextbrem(E,t,Delta))*V;
		radsupp = exp(raddeltaTsai(E, th, Delta, Z, M))*(1.0-probextbrem(E,t,Delta));
//		printf("%f %f %f\n", exp(raddeltaTsai(E, 4.3*3.14159/180.0, 0.006*2.2, Z, M))*(1.0-probextbrem(E,t,0.006*2.2)), exp(raddeltaTsai(E, 4.3*3.14159/180.0, 0.003, Z, M))*(1.0-probextbrem(E,t,0.003)) );

		A =   Aval(E, th);
		Am =  A*Pe;
		dAdrr =  dAval(E, th);
		dAmdrr =  dAdrr*Pe;
	    } else {
		crs    = 0.0;
		A      = 0.0;
		Am     = 0.0;
		dAdrr  = 0.0;
		dAmdrr = 0.0;
		radsupp = 0.0;
	    }

	    // Preradiation, effective radiator

	    // Vertex and vacuum corrections (going to assume loss is not so big)
	    // This should be subtracted from main contribution
	    dv = (-1.0/(137.0*3.14159))*(28.0/9.0 - 13.0/6.0*log(Q2/pow(0.000511,2.0)));

	    // We're being super pedantic and including recoil, etc.  It's all basically
	    // the same if we use the naive case
	    eff_t = -0.75*(raddeltaTsai(E, th, dmin, Z, M)-dv)/log(Ebeam/dmin)/2.0;
// Because we aren't sampling over the full range of tail, we
	    // will rescale the cross section
	    rscale1 = 1.0;

	    dmax = E-Emin;
	    if( dmax < dmin ){ printf("dmin > dmax\n"); continue; }

	    if( probextbrem(E,t*zt + eff_t, dmin) > gRandom->Uniform() ){
		// Radiative event!
		omega = GenerateRad(E, t*zt + eff_t, dmin, dmax );
		Er  = E - omega;
		rscale1 = 1.0-IRadProfile( E, t*zt + eff_t, dmax, E )/IRadProfile( E, t*zt + eff_t, dmin, E );
	    } else {
		Er = E;
	    }

	    Erf = (Er*M - Eex*(M+Eex))/(M + Er*(1.0-cos(th)));
	    Q2r =  2.0*Er*Erf*(1.0-cos(th));
	    q3v = sqrt(Ef*Ef - 2.0*Ef*Erf*cos(th) + Erf*Erf);


	    ///////////////////////////////////////  
	    //post radiation
	    rscale2 = 1.0;

	    dmax = Er-Emin;
	    if( probextbrem(Er,t*(1.0-zt) + eff_t, dmin) > gRandom->Uniform() ){
		// Radiative event!
		omega = GenerateRad(Er, t*(1.0-zt) + eff_t, dmin, dmax );
		Erf  = Erf - omega;
		rscale2 = 1.0-IRadProfile( Er, t*(1.0-zt) + eff_t, dmax, Er )/IRadProfile( Er, t*(1.0-zt) + eff_t, dmin, Er );
	    }	    

	    rate   = crs*lumin/nevt;
	    counts = rate*duration;

	    FF = -1e9;
	    if( j == 0 ){
		crsr = crsval(Er, th)*exp(dv)*V;
		FF =  crsval(Er, th)/crsMott(Z, M, Er, th);
	    } else {
		FF = exFF(Z,An,j,Q2r,th,q3v*q3v);
		crsr = crsMott(Z, M, Er, th)*exp(dv)*FF*V;
	    }
	    crsr *= rscale1*rscale2;
	    Ar =   Aval(Er, th);
	    Arm =  Ar*Pe;
	    radrate   = crsr*lumin/nevt;
	    radcounts = radrate*duration;


	    th0 = mydist->GenerateMS(E, thick*zt, An, Z);
	    ph0 = gRandom->Uniform()*2.0*pi;

	    th1 = mydist->GenerateMS(E, thick*(1.0-zt), An, Z);
	    ph1 = gRandom->Uniform()*2.0*pi;

	    //  Do rotations to get final angle after MS

	    // It goes  MS, inverse rotate to Z, rotate to hard, undo ms0 inverse,
	    // invert to Z, rotate to ms1, undo last invert

	    //  Initial vector after hard scattering 
	    x = sin(th)*cos(ph);
	    y = sin(th)*sin(ph);
	    z = cos(th);

	    //  Inverse of MS before hard
	    xp =  x*cos(th0) + z*sin(th0);
	    yp =  y;
	    zp = -x*sin(th0) + z*cos(th0);
	    x = xp; y = yp; z = zp;

	    xp =  x*cos(ph0) - y*sin(ph0);
	    yp =  x*sin(ph0) + y*cos(ph0);
	    zp =  z;
	    x = xp; y = yp; z = zp;

	    // th and ph after target scattering
	    thf = acos(z);
	    phf = atan2(y,x)+pi/2.0;

	    //  Coordinate after last leg of MS
	    x = sin(th1)*cos(ph1);
	    y = sin(th1)*sin(ph1);
	    z = cos(th1);

	    // Rotations for initial vector
	    xp =  x*cos(thf) + z*sin(thf);
	    yp =  y;
	    zp = -x*sin(thf) + z*cos(thf);
	    x = xp; y = yp; z = zp;

	    xp =  x*cos(phf) - y*sin(phf);
	    yp =  x*sin(phf) + y*cos(phf);
	    zp =  z;
	    x = xp; y = yp; z = zp;

	    // th and ph after everything
	    thf = acos(z);
	    phf = atan2(y,x)+pi/2.0;

	    //	printf("%f %f %f %f %e %e\n",th*180.0/pi, ph*180.0/pi, thf, phf, rate, A );

#ifdef __LEROSERES_H
	    smear = getlerose_ressmear()*Ebeam;
#else
	    smear = gRandom->Gaus(0.0, Ebeam*HRSres);
#endif
	    Erfobs = Erf + smear;
	    Q2obs = 2.0*E*Erfobs*(1.0-cos(thf));

	    T->Fill();
	}
    }

    TFile *fout = new TFile(filename, "RECREATE");
    T->Write("T", TObject::kOverwrite);
    fout->Close();

    return 0;
}

///////////////////////////////////////////////////////////

double getX0( double Z, double A ){
    double num = 716.4*A;
    double den = Z*(Z+1.0)*log(287.0/sqrt(Z));

    return num/den;
}

double RadProfile(double E, double t, double loss){
    double bt = 4.0*t/3.0;
    return 1./loss*(1.-loss/E+0.75*pow(loss/E,2))*pow(loss/E,bt);
}

double IRadProfile(double E, double t, double loss1, double loss2){
    double bt = 4.0*t/3.0;

    double I1 = pow( loss1/E, bt )/bt - pow( loss1/E, bt+1.0 )/(bt+1.0) + 0.75*pow(loss1/E, bt+2.0)/(bt+2.0);
    double I2 = pow( loss2/E, bt )/bt - pow( loss2/E, bt+1.0 )/(bt+1.0) + 0.75*pow(loss2/E, bt+2.0)/(bt+2.0);

    return I2-I1;
}

double GenerateRad( double E, double t, double lmin, double lmax ){
    double x = gRandom->Uniform();

    double ival = x*IRadProfile( E, t, lmin, lmax );

    const double conv = 1e-6; // Convergence criterion
    double delta = 1e9;
    int nmax = 500;
    int niter = 0;

    double x0 = lmin;
    double x1, ffp;

    // Newton's method to solve for loss value

    while(fabs(delta) > conv && niter < nmax ){
	ffp = (IRadProfile( E, t, x0, lmax ) - ival)/RadProfile(E, t, x0);
	x1 = x0 + ffp;
	delta = IRadProfile( E, t, x1, lmax ) - ival;
	x0 = x1;
	niter++;
    }

    if( niter >= nmax ){
	printf("Too many steps %d\n", niter);
	exit(1);
    }

    if( x0 < lmin || lmax < x1 ){
	printf("generated photon outside of given loss range\n");
	printf("Not %e < %e, %e < %e\n", lmin, x0, x1, lmax);
//	exit(1);
    }

    return x0;
}

double probextbrem( double E, double t, double Delta, double dmax ){
    double bt = 4.0*t/3.0;
    double num, den, prob;

    if( dmax < 0 ){
	num = 
	    1.0 - pow(Delta/E,bt) 
	    -   (bt/(bt+1.0))*( 1.0 - pow(Delta/E,bt+1.0))
	    +   0.75*(bt/(bt+2.0))*( 1.0 - pow(Delta/E,bt+2.0));

	den = TMath::Gamma(1.0 + bt);
	prob = num/den;
    } else {
	prob = probextbrem( E, t, Delta ) - probextbrem( E, t, dmax );
    }

    return prob;
}

double beta(double x);

double raddeltaTsai( double E, double th, double Delta, double Z, double M ){
    // Maddening delta from Tsai
    // Mo Tsai (1969), II.6
    //
    // Numerically checked this against their tables
    // accurate to < 0.1%

    // Assume this for ca48 elastic

    double m = 0.000511; // m_e in GeV

    /*
       double M = 0.938; // GeV
       double Z = 1.0;
       */

    double alpha = 1.0/137.0;
    double pi    = 3.1415927;

    double E1 = E;
    double E3 = E1*M/(M+E1*(1.0-cos(th)));
    double E4 = E1+M-E3;
    double p4 = sqrt( E4*E4 - M*M );

    double b4 = p4/E4;

    double eta = E1/E3;
    double Q2  = 2.0*E1*E3*(1.0-cos(th));
    double m2  = m*m;

    double line[7];

    //  7 lines, we omit -alpha/pi

    line[0] = (28.0/9.0) - (13.0/6.0)*log(Q2/m2) 
	+  ( log(Q2/m2) - 1.0 + 2.0*Z*log(eta))*( 2.0*log(E1/Delta) - 3.0*log(eta))
	-  Spence( (E3-E1)/E3 ) - Z*Z*log(E4/M);

    line[1] = Z*Z*log(M/(eta*Delta))*(log((1.0+b4)/(1.0-b4))/b4 - 2.0) 
	+ (Z*Z/b4)*( 
		0.5*log((1.0+b4)/(1.0-b4))*log((E4+M)/(2.0*M)) 
		- Spence( -1.0*sqrt( (E4-M)/(E4+M) )*sqrt((1.0+b4)/(1.0-b4))) 
		);

    line[2] = Z*(
	    Spence( -1.0*(M-E3)/(E1) ) - Spence( M*(M-E3)/(2.0*E3*E4 - M*E1) )
	    + Spence( 2.0*E3*(M-E3)/(2.0*E3*E4 - M*E1)) 
	    + log(fabs( (2.0*E3*E4 - M*E1)/(E1*(M-2.0*E3))) )*log(M/(2.0*E3))
	    );

    line[3] = -Z*(
	    Spence( -1.0*(E4-E3)/(E3) ) - Spence( M*(E4-E3)/(2.0*E1*E4 - M*E3) )
	    + Spence( 2.0*E1*(E4-E3)/(2.0*E1*E4 - M*E3)) 
	    + log(fabs( (2.0*E1*E4 - M*E3)/(E3*(M-2.0*E1))) )*log(M/(2.0*E1))
	    );

    line[4] = -Z*(
	    Spence( -1.0*(M-E1)/(E1) ) - Spence( (M-E1)/E1 )
	    + Spence( 2.0*(M-E1)/ M) 
	    + log(fabs( M/(2.0*E1-M) ))*log(M/(2.0*E1))
	    );

    line[5] = Z*(
	    Spence( -1.0*(M-E3)/(E3) ) - Spence( (M-E3)/E3 )
	    + Spence( 2.0*(M-E3)/ M) 
	    + log(fabs( M/(2.0*E3-M) ))*log(M/(2.0*E3))
	    );


    line[6] = 

	-1.0*Spence( (E1-E3)/E1 ) 
	+ (Z*Z/b4)*(
		Spence( sqrt( (E4-M)/(E4+M) )*sqrt( (1.0-b4)/(1.0+b4)) )
		-Spence( sqrt( (E4-M)/(E4+M) ) )
		+Spence( -sqrt( (E4-M)/(E4+M) ) )
		)
	;

    double sum = 0.0; int i;

    for( i = 0; i < 7; i++ ){
	//	printf("line %d  %f\n", i+1, line[i]);
	sum += line[i];
    }

    return -1.0*sum*alpha/pi; 
}

double raddeltaMY( double E, double th, double Delta ){
    E = th = Delta;
    return 0.0;
}


////////////////////////////////////////////////////////////////////
// Auxiliary functions

double Spence( double x ){
    if( x < -1.0 ){
	return -0.5*pow(log(fabs(x)),2.0) - 3.1415927*3.1415927/6.0 - fspence->Integral(0, 1.0/x);
    }
    if( x > 1.0 ){
	return -0.5*pow(log(fabs(x)),2.0) + 3.1415927*3.1415927/3.0 - fspence->Integral(0, 1.0/x);
    }

    return fspence->Integral(0.0,x);
}

double beta(double x){
    if( 1.0 - x > 0.0 ){
	return log(x)*log(x);
    }

    return 0;
}


////////////////////////////////////////////////////////////////////
// Cross section and asymmetry IO

void LoadData( char filename[],  char sfilename[]){
    FILE *f = fopen(filename, "r");
    FILE *fs = fopen(sfilename, "r");

    if(!f){
        fprintf(stderr, "Could not open %s\n", filename);
        exit(1);
    }
    if(!fs){
        fprintf(stderr, "Could not open %s\n", sfilename);
        exit(1);
    }


    double minE, maxE, minth, maxth;
    double dummy, sA;

    minE = 1e9;  maxE = 0.0;
    minth = 1e9; maxth = 0.0;

    int NE, Nth, dint;
    double E, th;

    fscanf(f, "%d%d", &NE, &Nth);
    fscanf(fs, "%d%d", &dint, &dint);

    int i,j;


    for( i = 0; i < NE; i++ ){
	fscanf(f,"%lf", &E );
	fscanf(fs,"%lf", &dummy );
	if( E < minE ){ minE = E; }
	if( E > maxE ){ maxE = E; }
	for( j = 0; j < Nth; j++ ){
	    fscanf(f,"%lf%lf%lf", &th, &_crs[i][j], &_A[i][j] );
	    if( th < minth ){ minth = th; }
	    if( th > maxth ){ maxth = th; }
	    fscanf(fs,"%lf%lf%lf", &dummy, &dummy, &sA );

	    // Streched is a 1% increase in r_n
	    // dA is dA/(dr/r)
	    _dA[i][j] = (sA-_A[i][j])*100.0;
	}
    }

    fclose(f);
    fclose(fs);

    /*
       double dE = (maxE - minE)/(NE-1);
       double dt = (maxth - minth)/(Nth-1);
       */
    // Keep these dimensionless
    double dE = 1.0;
    double dt = 1.0;

    for( i = 0; i < NE; i++ ){
	for( j = 0; j < Nth; j++ ){
	    if( i == 0 || j == 0 || i == NE-1 || j == Nth-1 ){
		_crs_de[i][j] = 0.0;
		_crs_dt[i][j] = 0.0;
		_crs_dedt[i][j] = 0.0;
		_A_de[i][j] = 0.0;
		_A_dt[i][j] = 0.0;
		_A_dedt[i][j] = 0.0;
		_dA_de[i][j] = 0.0;
		_dA_dt[i][j] = 0.0;
		_dA_dedt[i][j] = 0.0;

	    } else {
		_crs_de[i][j] = (_crs[i+1][j] - _crs[i-1][j])/(2.0*dE);
		_crs_dt[i][j] = (_crs[i][j+1] - _crs[i][j-1])/(2.0*dt);
		_crs_dedt[i][j] = (_crs[i+1][j+1] - _crs[i-1][j+1]
			- _crs[i+1][j-1] + _crs[i-1][j-1])/(4.0*dE*dt);

		_A_de[i][j] = (_A[i+1][j] - _A[i-1][j])/(2.0*dE);
		_A_dt[i][j] = (_A[i][j+1] - _A[i][j-1])/(2.0*dt);
		_A_dedt[i][j] = (_A[i+1][j+1] - _A[i-1][j+1]
			- _A[i+1][j-1] + _A[i-1][j-1])/(4.0*dE*dt);

		_dA_de[i][j] = (_dA[i+1][j] - _dA[i-1][j])/(2.0*dE);
		_dA_dt[i][j] = (_dA[i][j+1] - _dA[i][j-1])/(2.0*dt);
		_dA_dedt[i][j] = (_dA[i+1][j+1] - _dA[i-1][j+1]
			- _dA[i+1][j-1] + _dA[i-1][j-1])/(4.0*dE*dt);
	    }
	}
    }

    _nE = NE;
    _nth = Nth;

    _Emin = minE;
    _Emax = maxE;
    _thmin = minth;
    _thmax = maxth;

    return;
}

double interp( double E, double th, double grid[70][150], double grid_de[70][150], double grid_dt[70][150], double grid_dedt[70][150] ){
    double escale, tscale;
    int eidx, tidx;
    double ex, tx;

    th *= 180.0/3.1415927;
    E  *= 1000.0;

    escale = (E - _Emin)/(_Emax-_Emin);
    tscale = (th - _thmin)/(_thmax-_thmin);

    if( escale < 0.0 || escale >= 1.0 || tscale < 0.0 || tscale >= 1.0 ){ /*printf("OUT OF RANGE: E = %f (%f), th = %f (%f)\n", E, escale, th, tscale);*/ return 0.0; }

    eidx = (int) floor(escale*(_nE-1));
    tidx = (int) floor(tscale*(_nth-1));

    ex = escale*(_nE-1) - floor(escale*(_nE-1));
    tx = tscale*(_nth-1) - floor(tscale*(_nth-1));


    // Simple bilinear interpolation
    if( eidx == 0 || eidx == _nE-2 || tidx == 0 || tidx == _nth-2 ){
	double interp = grid[eidx][tidx]*(1.0-ex)*(1.0-tx)
	    + grid[eidx+1][tidx]*ex*(1.0-tx)
	    + grid[eidx][tidx+1]*(1.0-ex)*tx
	    + grid[eidx+1][tidx+1]*ex*tx;

	return interp;
    }
    // Simple bicubic spline interpolation

    double alpha[4][4], fvec[16];

    fvec[0]  = grid[eidx][tidx];
    fvec[1]  = grid[eidx+1][tidx];
    fvec[2]  = grid[eidx][tidx+1];
    fvec[3]  = grid[eidx+1][tidx+1];

    fvec[4]  = grid_de[eidx][tidx];
    fvec[5]  = grid_de[eidx+1][tidx];
    fvec[6]  = grid_de[eidx][tidx+1];
    fvec[7]  = grid_de[eidx+1][tidx+1];

    fvec[8]  = grid_dt[eidx][tidx];
    fvec[9]  = grid_dt[eidx+1][tidx];
    fvec[10] = grid_dt[eidx][tidx+1];
    fvec[11] = grid_dt[eidx+1][tidx+1];

    fvec[12] = grid_dedt[eidx][tidx];
    fvec[13] = grid_dedt[eidx+1][tidx];
    fvec[14] = grid_dedt[eidx][tidx+1];
    fvec[15] = grid_dedt[eidx+1][tidx+1];

    double splinesol[16][16] = {
	{1.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0},
	{-3.0, 3.0, 0.0, 0.0,  -2.0,-1.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0},
	{2.0,-2.0, 0.0, 0.0,  1.0, 1.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0},

	{0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  1.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0, -3.0, 3.0, 0.0, 0.0,  -2.0,-1.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,  0.0, 0.0, 0.0, 0.0,  2.0,-2.0, 0.0, 0.0,  1.0, 1.0, 0.0, 0.0},

	{-3.0, 0.0, 3.0, 0.0,  0.0, 0.0, 0.0, 0.0,  -2.0, 0.0, -1.0, 0.0,  0.0, 0.0, 0.0, 0.0},
	{ 0.0, 0.0, 0.0, 0.0, -3.0, 0.0, 3.0, 0.0,  0.0, 0.0, 0.0, 0.0, -2.0, 0.0, -1.0, 0.0 },
	{ 9.0,-9.0,-9.0, 9.0,  6.0, 3.0,-6.0,-3.0,  6.0,-6.0, 3.0,-3.0,  4.0, 2.0, 2.0, 1.0 },
	{-6.0, 6.0, 6.0,-6.0, -3.0,-3.0, 3.0, 3.0, -4.0, 4.0,-2.0, 2.0, -2.0,-2.0,-1.0,-1.0 },

	{2.0, 0.0,-2.0, 0.0,  0.0, 0.0, 0.0, 0.0,  1.0, 0.0, 1.0, 0.0,  0.0, 0.0, 0.0, 0.0},
	{0.0, 0.0, 0.0, 0.0,  2.0, 0.0,-2.0, 0.0,  0.0, 0.0, 0.0, 0.0,  1.0, 0.0, 1.0, 0.0},
	{-6.0, 6.0, 6.0,-6.0, -4.0,-2.0, 4.0, 2.0, -3.0, 3.0,-3.0, 3.0, -2.0,-1.0,-2.0,-1.0 },
	{ 4.0,-4.0,-4.0, 4.0,  2.0, 2.0,-2.0,-2.0,  2.0,-2.0, 2.0,-2.0,  1.0, 1.0, 1.0, 1.0 },
    };

    int i,j;
    for( i = 0; i < 16; i++ ){
	alpha[i%4][i/4] = 0.0;
	for( j = 0; j < 16; j++ ){
	    alpha[i%4][i/4] += splinesol[i][j]*fvec[j];
	}
    }

    double val = 0.0;

    for( i = 0; i < 4; i++ ){
	for( j = 0; j < 4; j++ ){
	    val += alpha[i][j]*pow(ex,i)*pow(tx,j);
	}
    }

    return val;
}


double crsMott( double Z, double M, double E, double th ){
    double Ef = E*M/(M+E*(1.0-cos(th)));
    double Q2 = 2.0*E*Ef*(1.0-cos(th));
    double cth2 = cos(th/2.0)*cos(th/2.0);

    double num = 4.0*Z*Z*(0.197*0.197)*Ef*Ef*Ef*cth2;
    double den = Q2*Q2*137.0*137.0*E;

    return num*1e-26/den; // Is in fm2, return in cm2
}

double crsval( double E, double th ){
    double val = interp(E,th,_crs, _crs_de, _crs_dt, _crs_dedt);
    //if( val < 0 ){ printf("NEGATIVE CROSS SECTION\n"); return 0.0; }

    // is in mb, return in units of cm^2
    return val*1.0e-27;
}

double Aval( double E, double th ){
    return interp(E,th,_A, _A_de, _A_dt, _A_dedt);
}

double dAval( double E, double th ){
    return interp(E,th,_dA, _dA_de, _dA_dt, _dA_dedt);
}

double biinterp( double E, double th, double grid[70][150], double [70][150], double [70][150], double [70][150] ){
    double escale, tscale;
    int eidx, tidx;
    double ex, tx;

    th *= 180.0/3.1415927;
    E  *= 1000.0;

    escale = (E - _Emin)/(_Emax-_Emin);
    tscale = (th - _thmin)/(_thmax-_thmin);

    if( escale < 0.0 || escale >= 1.0 || tscale < 0.0 || tscale >= 1.0 ){ /*printf("OUT OF RANGE\n");*/ return 0.0; }

    eidx = (int) floor(escale*(_nE-1));
    tidx = (int) floor(tscale*(_nth-1));

    ex = escale*(_nE-1) - floor(escale*(_nE-1));
    tx = tscale*(_nth-1) - floor(tscale*(_nth-1));

    // Simple bilinear interpolation
    double interp = grid[eidx][tidx]*(1.0-ex)*(1.0-tx)
	+ grid[eidx+1][tidx]*ex*(1.0-tx)
	+ grid[eidx][tidx+1]*(1.0-ex)*tx
	+ grid[eidx+1][tidx+1]*ex*tx;

    return interp;
}

double bcrsval( double E, double th ){
    // is in mb, return in cm^2
    return biinterp(E,th,_crs, _crs_de, _crs_dt, _crs_dedt)*1.0e-27;
}

double bAval( double E, double th ){
    return biinterp(E,th,_A, _A_de, _A_dt, _A_dedt);
}

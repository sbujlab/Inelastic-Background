#include "msdist.h"

msdist::msdist() {
    InitInternal();

    return;
}

msdist::msdist( double p, int nmat, double t[], double A[], double Z[] ){
    /*
       p    - electron momentum, [GeV]
       nmat - number of materials
       t    - Thickness [g/cm2]
       A    - Mass number
       Z    - Atomic number
       */
    InitInternal();
    msdist();
    Init( p, nmat, t, A, Z );
    return;
}

msdist::msdist( double p, double t, double A, double Z ){
    /*
       p    - electron momentum, [GeV]
       nmat - number of materials
       t    - Thickness [g/cm2]
       A    - Mass number
       Z    - Atomic number
       */
    InitInternal();
    Init( p, t, A, Z );

    return;
}

void msdist::InitInternal(){
    fInit = false;
    fErf2sig = TMath::Erf(2.0/sqrt(2.0));
    srand48(time(0));
    fNmat = 0;
}

void msdist::Init( double p, int nmat, double t[], double A[], double Z[] ){
    /* 
       Load materials and calculate necessary
       variables to generate distributions

       p    - electron momentum, [GeV]
       nmat - number of materials
       t    - Thickness [g/cm2]
       A    - Mass number
       Z    - Atomic number
    */

    int i;

    if( nmat >= MAT_MAX ){
	fprintf(stderr, "%s %s line %d: Too many materials.  Limited by MAT_MAX (%d)\n", 
		__FILE__, __FUNCTION__, __LINE__, MAT_MAX );
	return;
    }

    fNmat = nmat;

    fp = p;

    for( i = 0; i < nmat; i++ ){
	ft[i] = t[i];
	fA[i] = A[i];
	fZ[i] = Z[i];
    }

    double radsum = 0.0;
    double X0;

    for( i = 0; i < nmat; i++ ){
	X0  = 716.4*A[i]/(Z[i]*(Z[i]+1.0)*log(287.0/sqrt(Z[i])));
	radsum += t[i]/X0;
    }

    // First work out characteristic gaussian spread.
    // this is the PDG number, which is different from
    // the Moliere f0 width.  I think this number
    // accounts for the higher order terms in the sum,
    // so it's what we should use.
    double thpdg  = 13.6e-3*sqrt(radsum)*(1.0 + 0.038*log(radsum))/p;
    fthpdg = thpdg;

    // First calculate b

    
    double expb_num, expb_den;
    double bsum = 0.0;

    for( i = 0; i < fNmat; i++ ){
	expb_num = 6680.0*ft[i]*(fZ[i]+1.0)*pow(fZ[i],1.0/3.0);
	expb_den = fA[i]*(1.0+3.34*pow(fZ[i]/137.0,2.0));

	bsum += expb_num/expb_den;
    }

    double b = log( bsum );

    // Need to solve
    // B - log(B) = b
    // Use Newtons's method

    fB = solvelogeq(b);

    /////////////////////////////////////////

    // Change of variables

    double chi2, chi2_num, chi2_den;

    chi2 = 0.0;

    for( i = 0; i < fNmat; i++ ){
	// Use units of p of GeV
	// hbar.c = 0.197*1e-13 GeV.cm
	chi2_num = 4.0*3.14159*pow(0.197*1e-13, 2.0)*ft[i]*fZ[i]*(fZ[i]+1.0)*6.022e23;
	chi2_den = fp*fp*fA[i]*137.0*137.0;
	chi2  += chi2_num/chi2_den;
    }

    fchi2 = chi2;


    // This is the number from Moliere
    double th  = sqrt(fchi2*fB/2.0);
    th = th;

    fth = thpdg;

    double v0  = CalcMSDistPlane( 0.0    );
    double v2  = CalcMSDistPlane( 2.0*fth);
    double v10 = CalcMSDistPlane(10.0*fth);

    // Area under 2sigma gaussian
    double Agaus = v0*fErf2sig*sqrt(2.0*3.14159)*fth;

    // exponential parameters
    double l = log( v2/v10 )/(8.0*fth);
    fl = l;
    double C = v2*exp(l*2.0*fth);
    fC = C;

    // Area under tail envelope
    double Dt =  exp(-l*2.0*fth) - exp(-l*10.0*fth);
    fDt    = Dt;
    double Atail = C*Dt/l;

    ftailprob = Atail/(Atail+Agaus);

    fInit = true;

    return;
}

void   msdist::Init( double p, double t, double A, double Z ){
    /*
       p    - electron momentum, [GeV]
       nmat - number of materials
       t    - Thickness [g/cm2]
       A    - Mass number
       Z    - Atomic number
       */

    double tt[] = {t};
    double tA[] = {A};
    double tZ[] = {Z};

    Init( p, 1, tt, tA, tZ );
    return;
}

double msdist::J0(double x) {
    // Returns J0 for any real x
    // Stolen from ROOT in TMath.cxx
    
    double ax,z;
    double xx,y,result,result1,result2;
    const double p1  = 57568490574.0, p2  = -13362590354.0, p3 = 651619640.7;
    const double p4  = -11214424.18,  p5  = 77392.33017,    p6 = -184.9052456;
    const double p7  = 57568490411.0, p8  = 1029532985.0,   p9 = 9494680.718;
    const double p10 = 59272.64853,   p11 = 267.8532712;

    const double q1  = 0.785398164;
    const double q2  = -0.1098628627e-2,  q3  = 0.2734510407e-4;
    const double q4  = -0.2073370639e-5,  q5  = 0.2093887211e-6;
    const double q6  = -0.1562499995e-1,  q7  = 0.1430488765e-3;
    const double q8  = -0.6911147651e-5,  q9  = 0.7621095161e-6;
    const double q10 =  0.934935152e-7,   q11 = 0.636619772;

    if ((ax=fabs(x)) < 8) {
	y=x*x;
	result1 = p1 + y*(p2 + y*(p3 + y*(p4  + y*(p5  + y*p6))));
	result2 = p7 + y*(p8 + y*(p9 + y*(p10 + y*(p11 + y))));
	result  = result1/result2;
    } else {
	z  = 8/ax;
	y  = z*z;
	xx = ax-q1;
	result1 = 1  + y*(q2 + y*(q3 + y*(q4 + y*q5)));
	result2 = q6 + y*(q7 + y*(q8 + y*(q9 - y*q10)));
	result  = sqrt(q11/ax)*(cos(xx)*result1-z*sin(xx)*result2);
    }
    return result;
}

double msdist::solvelogeq(double b){
    // Newton's method to solve B - log(B) = b
    double err = 1e-4;

    double thisB = b;
    double lastB = 1e9;

    int n = 0;

    double f, df;
    
    // Fix at 100 iterations
    while( n < 100 && fabs(thisB-lastB)>err ){
	f  = thisB - log(thisB) - b;
	df = 1.0 - 1.0/thisB;

	lastB  = thisB;
	thisB -= f/df;

	n++;
    }

    return thisB;
}


double msdist::fn_integrand( double u, double th, int n ){
    // Check for bad values of logarithm
    if( log(u) != log(u) || !(log(u) > -1e9) ){ 
	return 0;
    }

    return u*J0(u*th)*exp(-0.25*u*u)*pow(0.25*u*u*log(0.25*u*u),n);
}

double msdist::intsimpson_fn( double th, int n ){
    if( n >= 5 ) {fprintf(stderr, "%s %s: %d:  Warning, integrating over integrand terms that are of too large n\n", 
	    __FILE__, __FUNCTION__, __LINE__ ); }

    /* Simpson's method of integration.
       We will choose the integration step
       to be dynamically generated based on th.
     */

    //  Zeros for J0 are spaced apart by at least  ~2.4
    //  so we will integrate at most in steps of 2.4/2

    double bess_step = 2.4/2.0/th;

    //  We want to integrate over the gaussian term at 
    //  most in 0.2 unit steps

    double gauss_step = 1.0;

    // Take the minimum
    double step = (bess_step < gauss_step? bess_step : gauss_step );

    // Integrate over to 8
    //  This should be good enough for n<3
    int  maxstep = (int) 8.0/step;
    int    nstep = 0;

    double stepsum;
    double sum = 0.0;

    for( nstep = 0; nstep < maxstep; nstep++ ){
	stepsum = (step/6.0)*(
	  	      fn_integrand(step*nstep, th, n)
	         +4.0*fn_integrand(step*(nstep+0.5), th, n)
	         +    fn_integrand(step*(nstep+1.0), th, n)
	       );

	sum += stepsum;
    }

    double fact = 1.0;
    int i;
    for( i = 1; i <= n; i++ ){
	fact *= i;
    }

    return sum/fact;
}

double msdist::CalcMSDistPlane( double theta, double p, double t, double A, double Z ){
    Init( p, t, A, Z );
    return CalcMSDistPlane(theta);
}

double msdist::CalcMSDistPlane( double theta, double p, int nmat, double t[], double A[], double Z[] ){
    Init( p, nmat, t, A, Z );
    return CalcMSDistPlane(theta);
}

double msdist::CalcMSDistPlane( double theta){
    /*
       p    - electron momentum, [GeV]
       nmat - number of materials
       t    - Thickness [g/cm2]
       A    - Mass number
       Z    - Atomic number
       */

    double th = fabs(theta)/sqrt(fchi2*fB);

    // Separately, we do three integrals from Moillere
    // Let's use Simpson's rule
    double f0 = 2.0*exp(-th*th);
    double f1 = intsimpson_fn( th, 1 );
    double f2 = intsimpson_fn( th, 2 );
//    double f3 = intsimpson_fn( th, 3 );
    double f3 = 0.0;

    return f0 + f1/fB + f2/pow(fB,2.0) + f3/pow(fB,3.0);
}

double msdist::CalcMSDist( double theta, double p, double t, double A, double Z ){
    Init( p, t, A, Z );
    return CalcMSDist(theta);
}

double msdist::CalcMSDist( double theta, double p, int nmat, double t[], double A[], double Z[] ){
    Init( p, nmat, t, A, Z );
    return CalcMSDist(theta);
}

double msdist::CalcMSDist( double theta){
    return CalcMSDistPlane(theta)*sin(fabs(theta));
}

double msdist::GenerateMSPlane(){
    /*
       Generate an event for a single momentum in a
       "plane projected" distribution (i.e. what you
       would measure if you took the angle in the a 
       single direction).


       We will assume perfect Gaussian up to 2
       sigma of that Gaussian.  This part is very fast

       For the tail events, we will use the sample/check
       method.  The envelope we sample is a decaying
       exponential from 2 - 10 sigma.  
     
       All the relevant parameters for this were set
       in Init()

     */
    if( !fInit ){
	return 0.0;
    }

    double trialv;

    // start rolling dice

    if( drand48() > ftailprob ){
	// Gaussian
	// Make sure we don't take more than two sigma here
	do {
	    trialv = sin(2.0*3.14159*drand48())*sqrt(-2.0*log(drand48()));
	}
	while( fabs(trialv) > 2.0 );

	return trialv*fth;
    } else {
	//  Now we have our long tail
	//  here we use sample and reject sampling from an exponential
	//  envelope
	//  We'll work with just positive for now and set the sign at
	//  the end
	//
	//  This has an efficiency of ~0.5, which is probably 
	//  pretty good since this are only 5% of the distribution
	do {
	    trialv = -log(exp(-fl*2.0*fth) - fDt*drand48()) /fl;
	} 
	while ( drand48() > CalcMSDistPlane( trialv )*
		exp(fl*trialv)/fC);


	// Choose side
	if( drand48() < 0.5 ){
	    trialv *= -1.0;
	}

	return trialv;
    }

    return -1e9;
}

double msdist::GenerateMSPlane( double p, int nmat, double t[], double A[], double Z[] ){
    /*
       p    - electron momentum, [GeV]
       nmat - number of materials
       t    - Thickness [g/cm2]
       A    - Mass number
       Z    - Atomic number
       */
    Init(p, nmat, t, A, Z);
    return GenerateMSPlane();
}

double msdist::GenerateMSPlane( double p, double t, double A, double Z ){
    /*
       p    - electron momentum, [GeV]
       nmat - number of materials
       t    - Thickness [g/cm2]
       A    - Mass number
       Z    - Atomic number
       */
    Init(p, t, A, Z);
    return GenerateMSPlane();
}


double msdist::GenerateMS(){
    // This returns the polar coordinate
    // theta for a single event

    double x = GenerateMSPlane();
    double y = GenerateMSPlane();
    
    return sqrt(x*x+y*y);
}

double msdist::GenerateMS( double p, int nmat, double t[], double A[], double Z[] ){
    /*
       p    - electron momentum, [GeV]
       nmat - number of materials
       t    - Thickness [g/cm2]
       A    - Mass number
       Z    - Atomic number
       */
    Init(p, nmat, t, A, Z);
    return GenerateMS();
}

double msdist::GenerateMS( double p, double t, double A, double Z ){
    /*
       p    - electron momentum, [GeV]
       nmat - number of materials
       t    - Thickness [g/cm2]
       A    - Mass number
       Z    - Atomic number
       */
    Init(p, t, A, Z);
    return GenerateMS();
}

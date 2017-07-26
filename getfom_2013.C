double acc[5000], ang[5000];
double thmin, thmax;
int Nacc;


double GetAcc(double th );
double GetBasicAcc(double th, double tc );

void getfom_2013(){
    thmin = 1e9; thmax = 0.0;
    int nread, i;
    Nacc = -1;
    double maxacc = 0.0;
    double phiacc = 2.9/9.093497;
   // double phiacc = 1.0;

    double cenang = 4.0;

    for( Nacc = 0; Nacc < 500; Nacc ++ ){

	ang[Nacc] = 1.0 + 9.0*((double) Nacc)/500.0;

	if( ang[Nacc] > 3.4 && ang[Nacc] < 4.6 ){
	    acc[Nacc] = 1.0;
	} else {
	    acc[Nacc] = 0.0;
	}

	if( ang[Nacc] < thmin ){ thmin = ang[Nacc]; }
	if( ang[Nacc] > thmax ){ thmax = ang[Nacc]; }

	if( acc[Nacc] > maxacc ){ maxacc = acc[Nacc]; }
    }


    TChain *T = new TChain("T");
    T->Add("ca48_2013_2200.root");

    double thf, phf, counts, rate, Am, Ef;
    int state;

    T->SetBranchAddress("thf", &thf);
    T->SetBranchAddress("phf", &phf);
    T->SetBranchAddress("counts", &counts);
    T->SetBranchAddress("rate", &rate);
    T->SetBranchAddress("A", &Am);

    int i;

    double ratesum = 0.0;
    double accsum = 0.0;
    double evsum = 0.0;
    double Nsum = 0.0;
    double Asum = 0.0;
    double thisacc;

    double totacc = (cos(1.0*3.14159/180)-cos(10.0*3.14159/180))*2.0*3.14159;

    for( i = 0; i < T->GetEntries(); i++ ){
	T->GetEntry(i);

	if( counts <= 0.0 ) continue;

	evsum += 1.0;

	if( thf != thf ) continue;

	thisacc = GetAcc(thf*180.0/3.1415927)*phiacc;

	if( thisacc > 0.0 && counts > 0.0){
	    accsum += thisacc;
	    Nsum += counts*thisacc;
	    ratesum += rate*thisacc;
	    Asum += Am*counts*thisacc*1e6;
//	    printf("%f %f %f\n", Nsum, ratesum, Asum );
	}
    }

    printf("rate = %e MHz\n", ratesum*1e-6 );
    printf("Am   = %f ppm\n", Asum/Nsum );
    printf("acc  = %f msr (of %f msr)\n", accsum*totacc*1e3/evsum, totacc*1e3);

    printf("accsum/evsum = %f/%f = %f\n", accsum, evsum, accsum/evsum  );

}

double GetAcc(double th ){
    double thscale = (th-thmin)/(thmax-thmin);

    if( thscale < 0 || thscale >= 1.0 ) return 0.0;

    int thidx = floor(thscale*(Nacc-1));
    double x = thscale*(Nacc-1) - floor(thscale*(Nacc-1));

    return acc[thidx]*(1.0-x) + acc[thidx+1]*x;
}

double GetBasicAcc( double th, double thc ){
    if( fabs(th-thc)<0.7 ){ return 1.0; }

    return 0.0;
}

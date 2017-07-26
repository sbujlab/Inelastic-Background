double acc[5000], ang[5000];
double thmin, thmax;
int Nacc;

double GetAcc(double th );
double GetBasicAcc(double th, double tc );

void getfom(){
    FILE *f = fopen("basicacc_fromdata.dat", "r");

    thmin = 1e9; thmax = 0.0;
    int nread, i;
    Nacc = -1;
    double maxacc = 0.0;

    double cenang = 5.4;
//    double cenang = 5.0;
    double phiscale = tan(55.0*3.14159/180.0)*tan(5.4*3.14159/180);
    // phi acceptance is ~55 degrees centered on 5.4 degrees
    
    double phiacc = atan(phiscale/tan(cenang*3.14159/180.0))/(2.0*3.14159);

    do {
	Nacc++;
	nread = fscanf(f, "%lf%lf", &ang[Nacc], &acc[Nacc] );
	// Trim off tail
	if( ang[Nacc] > 7.0 ){ acc[Nacc] = 0.0; }
	// Center on tail, offset is about 5.4 degrees
	ang[Nacc] -= 5.4;

	if( ang[Nacc] < thmin ){ thmin = ang[Nacc]; }
	if( ang[Nacc] > thmax ){ thmax = ang[Nacc]; }

	if( acc[Nacc] > maxacc ){ maxacc = acc[Nacc]; }
    } while( nread == 2 && !feof(f) );
    fclose(f);

    for( i = 0; i < Nacc; i++ ){
//	acc[i] /= maxacc;
    }

    TChain *T = new TChain("T");
    T->Add("dummy.root");

    double thf, phf, counts, rate, Am;

    T->SetBranchAddress("thf", &thf);
    T->SetBranchAddress("phf", &phf);
    T->SetBranchAddress("counts", &counts);
    T->SetBranchAddress("rate", &rate);
    T->SetBranchAddress("Am", &Am);

    int i;

    double ratesum = 0.0;
    double Nsum = 0.0;
    double Asum = 0.0;
    double thisacc;


    for( i = 0; i < T->GetEntries(); i++ ){
	T->GetEntry(i);

	thisacc = GetAcc(thf*180.0/3.1415927 - cenang)*phiacc;
//	thisacc = GetAcc(thf*180.0/3.1415927 - cenang);
//	thisacc = GetBasicAcc(thf*180.0/3.1415927, cenang);

	if( thisacc > 0.0 && counts > 0.0 ){
	    Nsum += counts*thisacc;
	    ratesum += rate*thisacc;
	    Asum += Am*counts*thisacc*1e6;
//	    printf("%f %f %f\n", Nsum, ratesum, Asum );
	}
    }

    printf("rate = %e\n", ratesum );
    printf("Am   = %f\n", Asum/Nsum );


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

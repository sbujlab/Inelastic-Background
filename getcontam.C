int Nacc;
double ACCOFF = 5.0;
double thmin = 1e9;
double thmax = -1e9;

double plow = 2.185;
double pcut = 2.1965;

double acc[500], ang[500];

//#define m_to_GeV 0.1276
#define m_to_GeV 0.1276

double GetAcc(double th ){
    double thscale = (th-thmin)/(thmax-thmin);

    if( thscale < 0 || thscale >= 1.0 ) return 0.0;

    int thidx = floor(thscale*(Nacc-1));

    double x = thscale*(Nacc-1) - floor(thscale*(Nacc-1));
//    printf("thscale %f (idx %d), %f %f, %f\n", thscale, thidx, acc[thidx], acc[thidx+1], x);

    return acc[thidx]*(1.0-x) + acc[thidx+1]*x;
}



void getcontam(){
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    TChain *T = new TChain("T");
    T->Add("ca48_2200_newres.root");

    TH1F *hel = new TH1F("hel", "hel", 150, pcut, 2.21);
    TH1F *hin = new TH1F("hin", "hin", 150, pcut, 2.21);

    TH1F *helf = new TH1F("helf", "helf", 150, (plow-2.2)/m_to_GeV, (2.21-2.2)/m_to_GeV);
    TH1F *hinf = new TH1F("hinf", "hinf", 150, (plow-2.2)/m_to_GeV, (2.21-2.2)/m_to_GeV);

    TH1F *hel2 = new TH1F("hel2", "hel2", 150, 2.17, 2.201);
    TH1F *hin2 = new TH1F("hin2", "hin2", 150, 2.17, 2.201);

    double Erfobs, thf, Erf, radrate, Q2;
    int state;

    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("Erfobs", true);
    T->SetBranchStatus("Erf", true);
    T->SetBranchStatus("thf", true);
    T->SetBranchStatus("state", true);
    T->SetBranchStatus("radrate", true);
    T->SetBranchStatus("Q2", true);

    T->SetBranchAddress("Erfobs", &Erfobs);
    T->SetBranchAddress("Erf", &Erf);
    T->SetBranchAddress("thf", &thf);
    T->SetBranchAddress("state", &state);
    T->SetBranchAddress("radrate", &radrate);
    T->SetBranchAddress("Q2", &Q2);

    int nread = 2;

    FILE *f = fopen("basicacc_fromdata.dat", "r");
    do {
	Nacc++;
	nread = fscanf(f, "%lf%lf", &ang[Nacc], &acc[Nacc] );
	if( nread != 2 ) continue;
	// Trim off tail

	if( ang[Nacc] > 7.0 ){ acc[Nacc] = 0.0; }
	// Center on tail, offset is about ACCOFF degrees
	ang[Nacc] -= ACCOFF;

	if( ang[Nacc] < thmin ){ thmin = ang[Nacc]; }
	if( ang[Nacc] > thmax ){ thmax = ang[Nacc]; }
    } while( nread == 2 && !feof(f) );
    fclose(f);

    int evt = 0;

    double thisacc;
    double Q2sum = 0.0;
    double Nsum = 0.0;
    for( evt = 0; evt < T->GetEntries(); evt++ ){
	T->GetEntry(evt);

	if( !(thf>0.0) ) continue;

	thisacc =  GetAcc(thf*180/3.14159-4.0)*radrate;

	if( !(radrate < 0) && !(radrate >= 0 )){
	    printf("Radrate nan...\n");
	    continue;
	}

	if( !(thisacc < 0) && !(thisacc>= 0 )){
	    printf("acc nan... %f\n", thf*180/3.14159);
	    continue;
	}

	if( Erfobs > plow ){
	    if( state==0 ){
//		printf("Filling %f (%f)\n", Erf, GetAcc(thf*180/3.14159-4.0)*radrate);
		if( Erfobs > pcut ){
		    hel->Fill(Erfobs,thisacc);
		    hel2->Fill(Erf, thisacc);
		}
		helf->Fill((Erfobs-2.2)/m_to_GeV,thisacc);
		Q2sum += thisacc*Q2;
		Nsum += thisacc;
	    } else {
		if( Erfobs > pcut ){
		    hin->Fill(Erfobs, thisacc);
		}
		hinf->Fill((Erfobs-2.2)/m_to_GeV,thisacc);
		hin2->Fill(Erf, thisacc);
	    }
	}
    }

    printf("Avg Q2 = %f\n", Q2sum/Nsum);

    helf->SetTitle("Elastic/Inelastic {}^{48}Ca Focal Plane x Distribution, 2.2 GeV, Anticipated HRS Resolution");
    helf->GetXaxis()->SetTitle("x [m]");
    helf->GetXaxis()->CenterTitle();

    hin->SetLineColor(kRed);
    hinf->SetLineColor(kRed);

    TCanvas *c = new TCanvas();
    c->SetLogy();
    helf->Draw();
    hinf->Draw("same");

    double cont = hin->Integral()/(hin->Integral()+hel->Integral());
    printf("%f\n", cont);

    TLegend *leg = new TLegend(0.195, 0.17, 0.45, 0.33);
    leg->AddEntry(helf, "Elastic", "l");
    leg->AddEntry(hinf, "Inelastic", "l");
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->Draw();

    TPaveLabel *lab = new TPaveLabel(0.60, 0.77, 0.88, 0.86, Form("Cont. = %3.1f%%", cont*100.0),  "NDC");
    lab->SetFillStyle(0);
    lab->SetBorderSize(0);
    lab->Draw();

    TLine *l = new TLine( (pcut-2.2)/m_to_GeV, helf->GetMinimum(), (pcut-2.2)/m_to_GeV, helf->GetMaximum());
    l->SetLineWidth(2);
    l->SetLineColor(kGreen);
    l->Draw();

    /*
    c->Print("20120430/contspec_full.pdf");
    c->Print("20120430/contspec_full.png");
    */

    //////////////////////////////////////////////////////////////////
    TCanvas *c2 = new TCanvas();
    c2->SetLogy();

    hel2->SetTitle("Elastic/Inelastic {}^{48}Ca Spectrum, 2.2 GeV, Accepted Events");
    hel2->GetXaxis()->SetTitle("p [GeV]");
    hel2->GetXaxis()->CenterTitle();

    hin2->SetLineColor(kRed);

    hel2->Draw();
    hin2->Draw("same");

    TLegend *leg2 = new TLegend(0.195, 0.67, 0.45, 0.83);
    leg2->AddEntry(helf, "Elastic", "l");
    leg2->AddEntry(hinf, "Inelastic", "l");
    leg2->SetFillStyle(0);
    leg2->SetBorderSize(0);
    leg2->Draw();

    /*
    c2->Print("20120430/accespec.pdf");
    c2->Print("20120430/accespec.png");
    */
}

#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TMultiGraph.h"

double acc[5000], ang[5000];
double thmin, thmax;
int Nacc;

double GetAcc(double th );
double GetBasicAcc(double th, double tc );

#define NENE 58
//#define NENE 4
#define minang 2.0
#define maxang 8.0

#define ACCOFF 5.0

//#define polerr 0.010
#define syserr 0.018

#define NRADIUS 3.436 //fm

#define CURRENT 100.0


void energyscan(){
    gROOT->SetStyle("Plain");
    gStyle->SetMarkerStyle(7);

    gStyle->SetFillStyle(0);

    FILE *f = fopen("basicacc_fromdata.dat", "r");

    thmin = 1e9; thmax = 0.0;
    int nread, i;
    Nacc = -1;
    double maxacc = 0.0;

    double phiscale = tan(55.0*3.14159/180.0)*tan(ACCOFF*3.14159/180);
    // phi acceptance is ~55 degrees centered on ACCOFF degrees
    

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

	if( acc[Nacc] > maxacc ){ maxacc = acc[Nacc]; }
    } while( nread == 2 && !feof(f) );


    fclose(f);

    for( i = 0; i < Nacc; i++ ){
	acc[i] /= maxacc;
    }

    /////////////////////////////////////////////////////////////

    double thf, phf, counts, rate, Am, dAmdrr;
    double phiacc, thisang, thisene;

    int eneidx, nucl;
    double ratesum = 0.0;
    double Nsum = 0.0;
    double Asum = 0.0;
    double dAAdrrsum = 0.0;
    double thisacc;

    TChain *T;

    double enearr[2][NENE], ratearr[2][NENE], asymarr[2][NENE], errarr[2][NENE];
    double dAAdrrarr[2][NENE], drrarr[2][NENE];
    double drarr[2][NENE];
    double drParr[2][NENE];
    double drrParr[2][NENE];

    TMultiGraph *mgrate = new TMultiGraph();
    TMultiGraph *mgerr  = new TMultiGraph();
    TMultiGraph *mgasym = new TMultiGraph();
    TMultiGraph *mgdAAdrr = new TMultiGraph();
    TMultiGraph *mgdrr = new TMultiGraph();
     TMultiGraph *mgdr = new TMultiGraph();



    TLegend *lrate = new TLegend(0.6, 0.7, 0.86, 0.86);
    TLegend *lasym = new TLegend(0.2, 0.7, 0.46, 0.86);
    TLegend *lerr  = new TLegend(0.2, 0.7, 0.46, 0.86);
    TLegend *ldAAdrr= new TLegend(0.2, 0.7, 0.46, 0.86);
    TLegend *ldrr= new TLegend(0.2, 0.7, 0.56, 0.86);
    TLegend *ldr= new TLegend(0.14, 0.14, 0.48, 0.36);


    for( nucl = 0; nucl < 1; nucl++ ){
	for( eneidx = 0; eneidx < NENE; eneidx++ ){
	    thisang = 4.0;
	    thisene = 0.6 + 0.05*eneidx;


	    T = new TChain("T");

	    if( nucl == 0 ){
		T->Add(Form("output_old/ca48_%.0f.root", thisene*1000.0));
	    } else {
		T->Add(Form("output_old/ca40_%.0f.root", thisene*1000.0));
	    }

	    T->SetBranchAddress("thf", &thf);
	    T->SetBranchAddress("phf", &phf);
	    T->SetBranchAddress("counts", &counts);
	    T->SetBranchAddress("rate", &rate);
	    T->SetBranchAddress("Am", &Am);
	    T->SetBranchAddress("dAmdrr", &dAmdrr);


	    phiacc = atan(phiscale/tan(thisang*3.14159/180.0))/(2.0*3.14159);
//	    phiacc *= pow(1.5/2.1,2.0);
	    phiacc *= pow(2.4/2.85,2.0);
	    phiacc *= 0.85;  // Kludge for output_old to account for reduced momentum
	    		    // acceptance from more relatistic inelastic cut


	    ratesum = 0.0;
	    Nsum = 0.0;
	    Asum = 0.0;
	    dAAdrrsum = 0.0;

	    for( i = 0; i < T->GetEntries(); i++ ){
		T->GetEntry(i);

		counts *= 2.0; // 2 HRS

		counts *= CURRENT/70.0;
		rate   *= CURRENT/70.0;

		thisacc = GetAcc(thf*180.0/3.1415927 - thisang)*phiacc;

		if( thisacc > 0.0 && counts > 0.0 ){
		    Nsum += counts*thisacc;
		    ratesum += rate*thisacc;
		    Asum += Am*counts*thisacc;
		    dAAdrrsum += dAmdrr*counts*thisacc/Am;
		    //	    printf("%f %f %f\n", Nsum, ratesum, Asum );
		}

	    }

	    printf("E = %e\n", thisene);
	    printf("rate = %e, counts = %e\n", ratesum, Nsum );
	    printf("Am   = %f\n", 1e6*Asum/Nsum );
	    printf("dA/A   = %f\n", sqrt(Nsum)/Asum);
	    printf("dA/A / dr/r   = %f\n", dAAdrrsum/Nsum);

	    enearr[nucl][eneidx] = thisene;
	    asymarr[nucl][eneidx] = Asum*1e6/Nsum;
	    ratearr[nucl][eneidx] = ratesum;
	    errarr[nucl][eneidx] = 100.0*sqrt(Nsum)/Asum;
	    dAAdrrarr[nucl][eneidx] =  -dAAdrrsum*0.01/Nsum;

	    drrarr[nucl][eneidx] = 1e-2*errarr[nucl][eneidx]/fabs(dAAdrrarr[nucl][eneidx]);

	    drrParr[nucl][eneidx] =  sqrt(
		    pow(1e-2*errarr[nucl][eneidx],2.0)+
		    pow(syserr,2.0)
		    )
		/fabs(dAAdrrarr[nucl][eneidx]);

	    drarr[nucl][eneidx] = drrarr[nucl][eneidx]*0.01*NRADIUS;
	    drParr[nucl][eneidx] = drrParr[nucl][eneidx]*0.01*NRADIUS;

	    printf("dr/r   = %f%%\n", drrarr[nucl][eneidx]);
	    printf("dr   = %f fm\n", drarr[nucl][eneidx]);


	}

	TGraph *grate = new TGraph(NENE, enearr[nucl], ratearr[nucl]);
	TGraph *gasym = new TGraph(NENE, enearr[nucl], asymarr[nucl]);
	TGraph *gerr  = new TGraph(NENE, enearr[nucl], errarr[nucl]);
	TGraph *gdAAdrr  = new TGraph(NENE, enearr[nucl], dAAdrrarr[nucl]);
	TGraph *gdrr  = new TGraph(NENE, enearr[nucl], drrarr[nucl]);
	TGraph *gdrrP  = new TGraph(NENE, enearr[nucl], drrParr[nucl]);

	 TGraph *gdr  = new TGraph(NENE, enearr[nucl], drarr[nucl]);
	 TGraph *gdrP = new TGraph(NENE, enearr[nucl], drParr[nucl]);


	if( nucl == 0 ){
	    grate->SetTitle("Ca48");
	    gasym->SetTitle("Ca48");
	    gerr->SetTitle("Ca48");
	    gdAAdrr->SetTitle("Ca48");
	    gdrr->SetTitle("Ca48");
	    gdrrP->SetTitle(Form("Ca48 - w/ #deltasys = %3.1f%%", syserr*100));

	    gdr->SetTitle("Ca48");
	    gdrP->SetTitle(Form("Ca48 - w/ #deltasys = %3.1f%%", syserr*100));


	    grate->SetMarkerColor(kRed);
	    gasym->SetMarkerColor(kRed);
	    gerr->SetMarkerColor(kRed);
	    gdAAdrr->SetMarkerColor(kRed);
	    gdrr->SetMarkerColor(kRed);
	    gdrrP->SetMarkerColor(kRed);
	    gdr->SetMarkerColor(kRed);
	    gdrP->SetMarkerColor(kRed);

	    grate->SetLineColor(kRed);
	    gasym->SetLineColor(kRed);
	    gerr->SetLineColor(kRed);
	    gdAAdrr->SetLineColor(kRed);
	    gdrr->SetLineColor(kRed);
	    gdrrP->SetLineColor(kRed);
	    gdr->SetLineColor(kRed);
	    gdrP->SetLineColor(kRed);

	}

	if( nucl == 1 ){
	    grate->SetTitle("Ca40");
	    gasym->SetTitle("Ca40");
	    gerr->SetTitle("Ca40");
	    gdAAdrr->SetTitle("Ca40");
	    gdrr->SetTitle("Ca40");

	    grate->SetMarkerColor(kBlue);
	    gasym->SetMarkerColor(kBlue);
	    gerr->SetMarkerColor(kBlue);
	    gdrr->SetMarkerColor(kBlue);

	    grate->SetLineColor(kBlue);
	    gasym->SetLineColor(kBlue);
	    gerr->SetLineColor(kBlue);
	    gdAAdrr->SetLineColor(kBlue);
	    gdrr->SetLineColor(kBlue);

	}

	grate->SetLineWidth(2);
	gasym->SetLineWidth(2);
	gerr->SetLineWidth(2);
	gdAAdrr->SetLineWidth(2);
	gdrr->SetLineWidth(2);
	gdrr->SetLineStyle(7);
	gdrrP->SetLineWidth(2);
	gdr->SetLineWidth(2);
	gdr->SetLineStyle(7);
	gdrP->SetLineWidth(2);

	lrate->AddEntry(grate, grate->GetTitle(), "L");
	lasym->AddEntry(gasym, gasym->GetTitle(), "L");
	lerr->AddEntry(gerr, gerr->GetTitle(), "L");
	ldAAdrr->AddEntry(gdAAdrr, gdAAdrr->GetTitle(), "L");
	ldrr->AddEntry(gdrr, gdrr->GetTitle(), "L");
	ldrr->AddEntry(gdrrP, gdrrP->GetTitle(), "L");

	ldr->AddEntry(gdr, gdr->GetTitle(), "L");
	ldr->AddEntry(gdrP, gdrP->GetTitle(), "L");

	mgrate->Add(grate, "C");
	mgasym->Add(gasym, "C");
	mgerr->Add(gerr, "C");
	mgdAAdrr->Add(gdAAdrr, "C");
	mgdrr->Add(gdrr, "C");
	mgdrr->Add(gdrrP, "C");
	mgdr->Add(gdr, "C");
	mgdr->Add(gdrP, "C");

    }


    TCanvas *crate = new TCanvas();
    crate->SetLogy();
    crate->SetGridx();
    crate->SetGridy();
    mgrate->Draw("A");

    mgrate->SetTitle("Rate vs. Energy, {}^{48}Ca, #theta = 4#circ, 100#muA, 1 HRS");
    mgrate->GetXaxis()->SetTitle("E_{beam} [GeV]");
    mgrate->GetXaxis()->CenterTitle();
    mgrate->GetYaxis()->SetTitle("Rate [Hz]");
    mgrate->GetYaxis()->CenterTitle();

    mgrate->Draw("A");

    lrate->SetBorderSize(0);
    //lrate->Draw();
    crate->Print("20120530/escan_rate.png");

    ////////////////////////////////////////////////////////////////

    TCanvas *casym = new TCanvas();
    casym->SetGridx();
    casym->SetGridy();

    mgasym->Draw("A");
    mgasym->SetTitle("Measured asymmetry vs. Energy, {}^{48}Ca, #theta = 4#circ");
    mgasym->GetXaxis()->SetTitle("E_{beam} [GeV]");
    mgasym->GetXaxis()->CenterTitle();
    mgasym->GetYaxis()->SetTitle("A_{m} [ppm]");
    mgasym->GetYaxis()->CenterTitle();

    mgasym->Draw("A");

    lasym->SetBorderSize(0);
    //lasym->Draw();
    casym->Print("20120530/escan_asym.png");

    ////////////////////////////////////////////////////////////////

    TCanvas *cerr  = new TCanvas();
    cerr->SetGridx();
    cerr->SetGridy();
    mgerr->Draw("A");

    mgerr->SetTitle("Relative statistical uncertainty vs. Energy, {}^{48}Ca, #theta = 4#circ, 100#muA, 2 HRS");
    mgerr->GetXaxis()->SetTitle("E_{beam} [GeV]");
    mgerr->GetXaxis()->CenterTitle();
    mgerr->GetYaxis()->SetTitle("#deltaA/A [%]");
    mgerr->GetYaxis()->CenterTitle();

    mgerr->Draw("A");

    lerr->SetBorderSize(0);
    //lerr->Draw();
    cerr->Print("20120530/escan_err.png");

    ////////////////////////////////////////////////////////////////

    TCanvas *cdAAdrr  = new TCanvas();
    cdAAdrr->SetGridx();
    cdAAdrr->SetGridy();
    mgdAAdrr->Draw("A");

    mgdAAdrr->SetTitle("dA/A for {}^{48}Ca 1% change in R vs. Energy, #theta = 4#circ");
    mgdAAdrr->GetXaxis()->SetTitle("E_{beam} [GeV]");
    mgdAAdrr->GetXaxis()->CenterTitle();
    mgdAAdrr->GetYaxis()->SetTitle("#deltaA/A");
    mgdAAdrr->GetYaxis()->CenterTitle();

    mgdAAdrr->Draw("A");

    ldAAdrr->SetBorderSize(0);
    //ldAAdrr->Draw();
    cdAAdrr->Print("20120530/escan_dAAdrr.png");

    ////////////////////////////////////////////////////////////////

    TCanvas *cdrr  = new TCanvas();
    cdrr->SetGridx();
    cdrr->SetGridy();
    mgdrr->Draw("A");

    mgdrr->SetTitle("dR/R vs. Energy, {}^{48}Ca, #theta = 4#circ, 100#muA, 2 HRS");
    mgdrr->GetXaxis()->SetTitle("E_{beam} [GeV]");
    mgdrr->GetXaxis()->CenterTitle();
    mgdrr->GetYaxis()->SetTitle("#deltaR/R [%]");
    mgdrr->GetYaxis()->CenterTitle();
    mgdrr->SetMinimum(0.0);
    mgdrr->SetMaximum(5.0);

    mgdrr->Draw("A");

    ldrr->SetBorderSize(0);
    ldrr->Draw();

    cdrr->Print("20120530/escan_drr.png");
    //
    ////////////////////////////////////////////////////////////////

    TCanvas *cdr  = new TCanvas();
    cdr->SetGridx();
    cdr->SetGridy();
    mgdr->Draw("A");

    mgdr->SetTitle("dR vs. Energy, {}^{48}Ca, #theta = 4#circ, 100#muA, 2 HRS");
    mgdr->GetXaxis()->SetTitle("E_{beam} [GeV]");
    mgdr->GetXaxis()->CenterTitle();
    mgdr->GetYaxis()->SetTitle("#deltaR [fm]");
    mgdr->GetYaxis()->CenterTitle();
    mgdr->SetMinimum(0.0);
    mgdr->SetMaximum(0.1);

    mgdr->Draw("A");

    ldr->SetBorderSize(0);
    ldr->Draw();

    cdr->Print("20120530/escan_dr.png");

    return;
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

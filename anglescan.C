#include "TChain.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TROOT.h"
#include "TH1.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TLegend.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double acc[5000], ang[5000];
double thmin, thmax;
int Nacc;

double GetAcc(double th );
double GetBasicAcc(double th, double tc );


//#define USEBASICACC
#define accwidth 2.00

#define NANG 60
//#define NANG 4
#define minang 2.0
#define maxang 8.0

#define ACCOFF 5.0

#define CURRENT 100.0
#define radlen 0.05
//#define radlen 0.10

//#define polerr 0.005
//#define polerr 0.010
#define syserr 0.018

#define GENMIN 1.0
#define GENMAX 10.0

#define NRADIUS 3.436 //fm

#define JULIETTE_PLOTS 

void anglescan(){
    gROOT->SetStyle("Plain");
    gStyle->SetMarkerStyle(7);
    gStyle->SetFillStyle(0);

#ifdef JULIETTE_PLOTS
    gStyle->SetTitleXSize(0.06);
    gStyle->SetTitleXOffset(0.75);
    gStyle->SetTitleYSize(0.06);
    gStyle->SetTitleYOffset(0.75);
    gStyle->SetLabelSize(0.05, "X");
    gStyle->SetLabelOffset(0.00, "X");
    gStyle->SetLabelSize(0.05, "Y");
    gStyle->SetLabelOffset(0.00, "Y");

    gStyle->SetTitleW(0.95);
    gStyle->SetTitleH(0.10);
    gStyle->SetTitleFontSize(0.12);
#endif//JULIETTE_PLOTS

    FILE *f = fopen("basicacc_fromdata.dat", "r");

    thmin = 1e9; thmax = 0.0;
    int nread, i;
    Nacc = -1;
    double maxacc = 0.0;

    double phiscale = tan(55.0*3.14159/180.0)*tan(ACCOFF*3.14159/180);
    // phi acceptance is ~55 degrees centered on ACCOFF degrees
    
    char accstr[255];
#ifdef USEBASICACC
    sprintf(accstr, "Flat acceptance, #Delta#theta = #pm%5.3f deg", accwidth/2.0);
#else
    sprintf(accstr, "Scaled PREX acceptance");
#endif

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

    double thf, phf, counts, rate, Am, dAmdrr, th;
    double phiacc, thisang;

    int angidx, nucl;
    double ratesum = 0.0;
    double Nsum = 0.0;
    double Asum = 0.0;
    double dAAdrrsum = 0.0;
    double sangsum = 0.0;
    double thisacc;

    double sgen = (acos(GENMIN*3.14159/180.0) - acos(GENMAX*3.14159/180.0))*2.0*3.14159;

    printf("Generated solid angle is %f msr\n", sgen*1e3);

    TChain *T;

    double angarr[NANG], ratearr[NANG], asymarr[NANG], errarr[NANG];
    double drrParr[NANG];
    double dAAdrrarr[NANG], drrarr[NANG], sangarr[NANG];
    double drarr[NANG], drParr[NANG];

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
    TLegend *ldrr= new TLegend(0.15, 0.63, 0.62, 0.83);
    TLegend *ldr= new TLegend(0.412, 0.16, 0.882, 0.33);


    for( nucl = 0; nucl < 1; nucl++ ){

	T = new TChain("T");

	if( nucl == 0 ){
//	    T->Add("output/ca48_2200.root");
	    if( radlen < 0.075 ){
		T->Add("ca48_2200_newdelta.root");
	    } else {
		T->Add("ca48_10_2200.root");
	    }
	} else {
	    T->Add("output/ca40_2200.root");
	}

	T->SetBranchAddress("thf", &thf);
	T->SetBranchAddress("th", &th);
	T->SetBranchAddress("phf", &phf);
	T->SetBranchAddress("counts", &counts);
	T->SetBranchAddress("rate", &rate);
	T->SetBranchAddress("Am", &Am);
	T->SetBranchAddress("dAmdrr", &dAmdrr);

	for( angidx = 0; angidx < NANG; angidx++ ){
	    thisang = (maxang-minang)*((double) angidx)/((double) NANG) +minang;

	    phiacc = atan(phiscale/tan(thisang*3.14159/180.0))/(2.0*3.14159);
	//    phiacc = atan(phiscale/tan(ACCOFF*3.14159/180.0))/(2.0*3.14159);

//	    phiacc *= 1.5; // Increase for potential acceptance 

	//  From moving the target back
//	    phiacc *= pow(1.5/2.1,2.0);
//	    // 2.4 is from PREX survey
//	    // +45 cm is from Juliette on what John used for the septum
	    phiacc *= pow(2.4/2.85,2.0);

	    ratesum = 0.0;
	    Nsum = 0.0;
	    Asum = 0.0;
	    dAAdrrsum = 0.0;
	    sangsum = 0.0;

	    for( i = 0; i < T->GetEntries(); i++ ){
		T->GetEntry(i);

		if( th < GENMIN*3.14159/180.0-1e-4 || th > GENMAX*3.14159/180.0+1e-4 ){
		    printf("WARNING:  Generated value may be out of range.\nSolid angle calculations could be wrong!!!\n !  %f < %f < %f\n\n", GENMIN, th*180.0/3.14159, GENMAX);
		}
		/*
		counts *= 0.51/0.61;
		rate   *= 0.51/0.61;
		*/

		counts *= 2.0;

		// 100 uA
		/*
		counts *= CURRENT/70.0;
		rate   *= CURRENT/70.0;
		*/
		counts *= CURRENT/100.0;
		rate   *= CURRENT/100.0;

	
#ifdef USEBASICACC
		thisacc = GetBasicAcc(thf*180.0/3.1415927, thisang)*phiacc;
#else
		thisacc = GetAcc(thf*180.0/3.1415927 - thisang)*phiacc;
#endif

		if( thisacc > 0.0 && counts > 0.0 ){
		    Nsum += counts*thisacc;
		    ratesum += rate*thisacc;
		    Asum += Am*counts*thisacc;
		    dAAdrrsum += dAmdrr*counts*thisacc/Am;
		    sangsum += 1.0*thisacc;
		    //	    printf("%f %f %f\n", Nsum, ratesum, Asum );
		}
	    }

	    printf("-----------------------\n");
	    printf("th = %e\n", thisang);
	    printf("rate = %e, counts = %e\n", ratesum, Nsum );
	    printf("Am   = %f\n", 1e6*Asum/Nsum );
	    printf("dA/A   = %f\n", sqrt(Nsum)/Asum);
	    printf("dA/A / dr/r   = %f\n", dAAdrrsum*0.01/Nsum);
	    printf("solidang = %f msr\n", sangsum*sgen*1e3/T->GetEntries());

	    angarr[angidx] = thisang;
	    asymarr[angidx] = Asum*1e6/Nsum;
	    ratearr[angidx] = ratesum;
	    errarr[angidx] = 100.0*sqrt(Nsum)/Asum;
	    sangarr[angidx] = sangsum*4.0*3.14159*1e-3/T->GetEntries();
	    dAAdrrarr[angidx] =  -dAAdrrsum*0.01/Nsum;
	    drrarr[angidx] = 1e-2*errarr[angidx]/fabs(dAAdrrarr[angidx]);
	    drrParr[angidx] = sqrt(
			pow(1e-2*errarr[angidx],2.0)+
			pow(syserr,2.0)
		    )
		/fabs(dAAdrrarr[angidx]);

	    drarr[angidx] = drrarr[angidx]*0.01*NRADIUS;
	    drParr[angidx] = drrParr[angidx]*0.01*NRADIUS;

	    printf("dr/r   = %f%%\n", drrarr[angidx]);
	    printf("dr   = %f fm\n", drarr[angidx]);
	}

	TGraph *grate = new TGraph(NANG, angarr, ratearr);
	TGraph *gasym = new TGraph(NANG, angarr, asymarr);
	TGraph *gerr  = new TGraph(NANG, angarr, errarr);
	TGraph *gdAAdrr  = new TGraph(NANG, angarr, dAAdrrarr);
	TGraph *gdrr  = new TGraph(NANG, angarr, drrarr);
	TGraph *gdrrP = new TGraph(NANG, angarr, drrParr);

	TGraph *gdr  = new TGraph(NANG, angarr, drarr);
	TGraph *gdrP = new TGraph(NANG, angarr, drParr);

	/*
	mgrate->Add(grate, "P");
	mgasym->Add(gasym, "P");
	mgerr->Add(gerr, "P");
	*/
	mgrate->Add(grate, "C");
	mgasym->Add(gasym, "C");
	mgerr->Add(gerr, "C");
	mgdAAdrr->Add(gdAAdrr, "C");
	mgdrr->Add(gdrr, "C");
	mgdrr->Add(gdrrP, "C");

	mgdr->Add(gdr, "C");
	mgdr->Add(gdrP, "C");

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
	    gdrrP->SetTitle(Form("Ca40 - w/ #deltasys = %3.1f%%", syserr*100));
	    gdr->SetTitle("Ca40");
	    gdrP->SetTitle(Form("Ca40 - w/ #deltasys = %3.1f%%", syserr*100));

	    grate->SetMarkerColor(kBlue);
	    gasym->SetMarkerColor(kBlue);
	    gerr->SetMarkerColor(kBlue);
	    gdAAdrr->SetMarkerColor(kBlue);
	    gdrr->SetMarkerColor(kBlue);
	    gdrrP->SetMarkerColor(kBlue);
	    gdr->SetMarkerColor(kBlue);
	    gdrP->SetMarkerColor(kBlue);

	    grate->SetLineColor(kBlue);
	    gasym->SetLineColor(kBlue);
	    gerr->SetLineColor(kBlue);
	    gdAAdrr->SetLineColor(kBlue);
	    gdrr->SetLineColor(kBlue);
	    gdrrP->SetLineColor(kBlue);
	    gdr->SetLineColor(kBlue);
	    gdrP->SetLineColor(kBlue);
	}


	grate->SetLineWidth(2);
	gasym->SetLineWidth(2);
	gerr->SetLineWidth(2);
	gdAAdrr->SetLineWidth(2);
	gdrr->SetLineWidth(2);
	gdrrP->SetLineWidth(2);
	gdrr->SetLineStyle(7);
	gdr->SetLineWidth(2);
	gdrP->SetLineWidth(2);
	gdr->SetLineStyle(7);

	lrate->AddEntry(grate, grate->GetTitle(), "L");
	lasym->AddEntry(gasym, gasym->GetTitle(), "L");
	lerr->AddEntry(gerr, gerr->GetTitle(), "L");
	ldAAdrr->AddEntry(gdAAdrr, gdAAdrr->GetTitle(), "L");
	ldrr->AddEntry(gdrr, gdrr->GetTitle(), "L");

	ldrr->AddEntry(gdrrP, gdrrP->GetTitle(), "L");
	ldr->AddEntry(gdr, gdr->GetTitle(), "L");

	ldr->AddEntry(gdrP, gdrP->GetTitle(), "L");


    }


    TCanvas *cmain = new TCanvas();
    cmain->Divide(2,2);

    TCanvas *crate = new TCanvas();

    cmain->cd(1);
    crate->SetLogy();
    crate->SetGridx();
    crate->SetGridy();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    mgrate->Draw("A");

    mgrate->SetTitle(Form("Rate vs. Central Angle, %2.0f%% rad, E = 2.2 GeV, %3.0f#muA, 1 HRS, %s", radlen*100.0, CURRENT, accstr));
    mgrate->GetXaxis()->SetTitle("#theta [deg]");
    mgrate->GetXaxis()->CenterTitle();
    mgrate->GetYaxis()->SetTitle("Rate [Hz]");
    mgrate->GetYaxis()->CenterTitle();

    mgrate->Draw("A");

    lrate->SetBorderSize(0);
//    lrate->Draw();
    crate->cd();
    mgrate->Draw("A");
    crate->Print("20120919/angscan_rate.pdf");

    ////////////////////////////////////////////////////////////////

    TCanvas *casym = new TCanvas();
    cmain->cd(2);
    casym->SetGridx();
    casym->SetGridy();
    gPad->SetGridx();
    gPad->SetGridy();

    mgasym->Draw("A");
    mgasym->SetTitle(Form("Measured asymmetry vs. Central Angle, E = 2.2 GeV, %s", accstr));
    mgasym->GetXaxis()->SetTitle("#theta [deg]");
    mgasym->GetXaxis()->CenterTitle();
    mgasym->GetYaxis()->SetTitle("A_{m} [ppm]");
    mgasym->GetYaxis()->CenterTitle();
    mgasym->SetMinimum(0.0);

    mgasym->Draw("A");

    lasym->SetBorderSize(0);
//    lasym->Draw();
    casym->cd();
    mgasym->Draw("A");
    casym->Print("20120919/angscan_asym.pdf");

    
    ////////////////////////////////////////////////////////////////

    cmain->cd(3);
    TCanvas *cerr  = new TCanvas();
    cerr->SetGridx();
    cerr->SetGridy();
    gPad->SetGridx();
    gPad->SetGridy();
    mgerr->Draw("A");

    mgerr->SetTitle(Form("Relative statistical uncertainty vs. Central Angle, %2.0f%% rad, E = 2.2 GeV, %3.0f#muA, 2 HRS, %s", radlen*100, CURRENT, accstr));
    mgerr->GetXaxis()->SetTitle("#theta [deg]");
    mgerr->GetXaxis()->CenterTitle();
    mgerr->GetYaxis()->SetTitle("#deltaA/A [%]");
    mgerr->GetYaxis()->CenterTitle();
    mgerr->SetMinimum(0.0);

    mgerr->Draw("A");

    lerr->SetBorderSize(0);
 //   lerr->Draw();
    cerr->cd();
    mgerr->Draw("A");
    cerr->Print("20120919/angscan_err.pdf");
    
    ////////////////////////////////////////////////////////////////

    cmain->cd(4);
    TCanvas *cdAAdrr  = new TCanvas();
    cdAAdrr->SetGridx();
    cdAAdrr->SetGridy();
    gPad->SetGridx();
    gPad->SetGridy();
    mgdAAdrr->Draw("A");

    mgdAAdrr->SetTitle(Form("dA/A for 1%% change in R vs. Central Angle, E = 2.2 GeV, %s", accstr));
    mgdAAdrr->GetXaxis()->SetTitle("#theta [deg]");
    mgdAAdrr->GetXaxis()->CenterTitle();
    mgdAAdrr->GetYaxis()->SetTitle("#deltaA/A");
    mgdAAdrr->GetYaxis()->CenterTitle();

    mgdAAdrr->Draw("A");
    ldAAdrr->SetBorderSize(0);
  //  ldAAdrr->Draw();
    cdAAdrr->cd();
    mgdAAdrr->Draw("A");
    cdAAdrr->Print("20120919/angscan_dAAdrr.pdf");

    
    if( radlen < 0.075 ){
	cmain->Print("20120919/angscan.pdf");
    } else {
	cmain->Print("20120919/angscan_10rad.pdf");
    }
    ////////////////////////////////////////////////////////////////

    TCanvas *cdrr  = new TCanvas();
    cdrr->SetGridx();
    cdrr->SetGridy();
    mgdrr->Draw("A");

    mgdrr->SetTitle(Form("dR/R vs. Central Angle, %2.0f%% rad, E = 2.2 GeV, %3.0f#muA, 2 HRS, %s", radlen*100, CURRENT, accstr));
    mgdrr->GetXaxis()->SetTitle("#theta [deg]");
    mgdrr->GetXaxis()->CenterTitle();
    mgdrr->GetYaxis()->SetTitle("#deltaR/R [%]");
    mgdrr->GetYaxis()->CenterTitle();
    mgdrr->SetMinimum(0.0);
    mgdrr->SetMaximum(5.0);

    mgdrr->Draw("A");

    ldrr->SetBorderSize(0);
    ldrr->Draw();

//    cdrr->Print("20120919/angscan_drr_nogeo.png");
    if( radlen < 0.075 ){
	cdrr->Print("20120919/angscan_drr.pdf");
    } else {
	cdrr->Print("20120919/angscan_drr_10rad.png");
    }

    ////////////////////////////////////////////////////////////////////////
    TCanvas *cdr  = new TCanvas();
    cdr->SetGridx();
    cdr->SetGridy();
    mgdr->Draw("A");

    mgdr->SetTitle(Form("dR vs. Central Angle, %2.0f%% rad, E = 2.2 GeV, %3.0f#muA, 2 HRS, %s", radlen*100, CURRENT, accstr));
    mgdr->GetXaxis()->SetTitle("#theta [deg]");
    mgdr->GetXaxis()->CenterTitle();
    mgdr->GetYaxis()->SetTitle("#deltaR [fm]");
    mgdr->GetYaxis()->CenterTitle();
    mgdr->SetMinimum(0.0);
    mgdr->SetMaximum(0.1);

    mgdr->Draw("A");

    ldr->SetBorderSize(0);
    ldr->Draw();

//    cdrr->Print("20120919/angscan_drr_nogeo.png");
    if( radlen < 0.075 ){
	cdr->Print("20120919/angscan_dr.pdf");
    } else {
	cdr->Print("20120919/angscan_dr_10rad.png");
    }

    TCanvas *cpanel = new TCanvas();
    cpanel->Divide(2,2);
    cpanel->cd(1);
    mgrate->Draw("A");
    cpanel->cd(2);
    mgasym->Draw("A");
    cpanel->cd(3);
    mgdAAdrr->Draw("A");
    cpanel->cd(4);
    mgdrr->Draw("A");


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
    if( fabs(th-thc)<accwidth/2.0 ){ return 1.0; }
    return 0.0;
}

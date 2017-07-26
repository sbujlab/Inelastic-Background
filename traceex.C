void traceex(){
    gROOT->SetStyle("Plain");
    gStyle->SetOptStat(0);

    TChain *T = new TChain("T");
    T->Add("ca48_2200.root");
    TCanvas *c =  new TCanvas();
    c->Divide(2,2);

    c->cd(1);
    gPad->SetLogy();
    T->Draw("Erf >> h1(100, 2.19, 2.201)", "radcounts*(Erf>(2.2-0.01) && abs(th*180/3.14159-3)<0.5)");
    h1->SetMinimum(1e11);
    h1->SetTitle("{}^{48}Ca 5% rad, E = 2.2 GeV, #theta = 3#pm0.5#circ");
    h1->GetXaxis()->SetTitle("E_{f} [GeV]");
    h1->GetXaxis()->CenterTitle();

    T->SetLineColor(kRed);
    T->Draw("Erf", "radcounts*(Erf>(2.2-0.010) && state>0 && abs(th*180/3.14159-3.0)<0.5 )", "same");

    c->cd(2);
    gPad->SetLogy();
    T->SetLineColor(kBlack);
    T->Draw("Erf >> h2(100, 2.19, 2.201)", "radcounts*(Erf>(2.2-0.01) && abs(th*180/3.14159-4)<0.5)");
    h2->SetMinimum(1e11);
    h2->SetTitle("{}^{48}Ca 5% rad, E = 2.2 GeV, #theta = 4#pm0.5#circ");
    h2->GetXaxis()->SetTitle("E_{f} [GeV}");
    h2->GetXaxis()->CenterTitle();
    T->SetLineColor(kRed);
    T->Draw("Erf", "radcounts*(Erf>(2.2-0.010) && state>0 && abs(th*180/3.14159-4.0)<0.5 )", "same");

    c->cd(3);
    gPad->SetLogy();
    T->SetLineColor(kBlack);
    T->Draw("Erf >> h3(100, 2.19, 2.201)", "radcounts*(Erf>(2.2-0.01) && abs(th*180/3.14159-5)<0.5)");
    h3->SetMinimum(1e11);
    h3->SetTitle("{}^{48}Ca 5% rad, E = 2.2 GeV, #theta = 5#pm0.5#circ");
    h3->GetXaxis()->SetTitle("E_{f} [GeV]");
    h3->GetXaxis()->CenterTitle();
    T->SetLineColor(kRed);
    T->Draw("Erf", "radcounts*(Erf>(2.2-0.010) && state>0 && abs(th*180/3.14159-5.0)<0.5 )", "same");

    ////////////////////////////////////////
    c->cd(4);
    gPad->SetLogy();
    T->SetLineColor(kBlack);
    T->Draw("Erfobs >> h2s(100, 2.185, 2.21)", "radcounts*(Erf>(2.2-0.01) && abs(th*180/3.14159-4)<0.5)");
    h2s->SetMinimum(1e11);
    h2s->SetTitle("{}^{48}Ca 5% rad, E = 2.2 GeV, #theta = 4#pm0.5#circ, #deltaE = 10^{-3}");
    h2s->GetXaxis()->SetTitle("E_{f} [GeV]");
    h2s->GetXaxis()->CenterTitle();
    T->SetLineColor(kRed);
    T->Draw("Erfobs", "radcounts*(Erf>(2.2-0.010) && state>0 && abs(th*180/3.14159-4.0)<0.5 )", "same");

}

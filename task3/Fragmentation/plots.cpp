void plots(){
   double BF_bupsik = 1.010E-3;
    double BF_bupsik_err = 0.028 / 1.010;
    double BF_bspsiphi = 1.08E-3 * 0.492;
    double BF_bspsiphi_err = sqrt(pow(0.08 / 1.08, 2) + pow(0.5 / 49.2, 2));
      double Factor = 0.93518/0.489;
      std::cout<<"factor" << Factor<<std::endl;
  
 double labels[] = {12, 13, 14, 15, 16, 18, 20, 23, 26, 29, 34, 45, 70};
    int n_labels = sizeof(labels) / sizeof(labels[0]);


TFile *fInput = TFile::Open("output.root");
TH1D *htemp =(TH1D*) fInput->Get("fragm_ratio");
const Double_t ptBin[13] = {12,13,14,15,16,18,20,23,26,29,34,45,70};

TH1D *hFragmFrac = new TH1D("hFragmFrac", "hFragmFrac", 12, ptBin);
TH1D *hFragmFrac_cms = new TH1D("hFragmFrac_cms", "hFragmFrac_cms", 12, ptBin);
TLatex text;

double PaperMeas[12] = {0.1314,0.1196,0.1165,0.1154,0.1135,0.1106,0.1105,0.1110,0.1091,0.1095,0.1088,0.1117};
double PaperMeasStat[12] = {2.1,1.6,1.3,1.2,0.8,0.8,0.7,0.8,0.9,0.9,0.9,1.3};
double PaperMeasSys[12] = {3.1,2.7,2.4,2.6,2.6,2.8,2.9,2.6,3.2,2.3,2.8,2.6};
for(int i=0; i<13; i++){
    hFragmFrac->SetBinContent(i+1, htemp->GetBinContent(i+1));
    hFragmFrac->SetBinError(i+1, 0.079*htemp->GetBinContent(i+1));
    std::cout<<hFragmFrac->GetBinError(i+1)<<std::endl;
    hFragmFrac_cms->SetBinContent(i+1, PaperMeas[i]);
   hFragmFrac_cms->SetBinError(i+1, sqrt(TMath::Power(PaperMeasStat[i]*PaperMeas[i]*0.01,2) + TMath::Power(PaperMeasSys[i]*PaperMeas[i]*0.01,2)));
}
hFragmFrac_cms->Scale(Factor);

TCanvas *c_fsfu = new TCanvas("c_fsfu", "c_fsfu", 1500,2000);
c_fsfu->Divide(1,2);
TPad *padMain = (TCanvas*) c_fsfu->cd(1);
  padMain->Draw();
  padMain->SetTopMargin(0.075);
  padMain->SetBottomMargin(0.);
  padMain->SetLeftMargin(0.15);
  padMain->SetRightMargin(0.03);
  hFragmFrac->GetYaxis()->SetRangeUser(0.,0.40);
//  padMain->SetLogy();
  padMain->SetLogx(0);  
 // padMain->SetTicks(1, 1);

  padMain->cd();   

   hFragmFrac->SetFillColor(kWhite); // Fill color for the error boxes
  hFragmFrac->SetFillStyle(1001); // Solid fill
    hFragmFrac->SetLineColor(kBlack);
    hFragmFrac->Draw("ep");

hFragmFrac->SetMarkerStyle(kFullCircle);
hFragmFrac->SetTitle("");
hFragmFrac->GetYaxis()->SetTitle("f_{s}/f_{u}");
hFragmFrac->GetYaxis()->SetTitleSize(0.07);
hFragmFrac->GetYaxis()->SetTitleOffset(0.8);
hFragmFrac->SetMarkerSize(1.5);

  TLatex textEta;
#ifdef 
 textEta.DrawLatexNDC(0.8,0.8, "|#eta|<2.5");
 TLatex textDecay;
 textDecay.DrawLatexNDC(0.2,0.8, "B_{s}->J/#Psi K");
 TLatex textDecay2;
 textDecay2.DrawLatexNDC(0.2,0.7,"B^{+}->J/#Psi #Phi");
hFragmFrac->SetStats(0);

TPad *padRef = (TCanvas*) c_fsfu->cd(2);
  padRef->SetTopMargin(0.);
  padRef->SetBottomMargin(0.3);
  padRef->SetLeftMargin(0.15);
  padRef->SetRightMargin(0.03);
  padRef->SetLogy(0);
  padRef->SetLogx(0);  
//  padRef->SetTicks(1, 1);
  padRef->Draw();
  padRef->cd();
      for (int i = 1; i <= n_labels; ++i) {
        double bin_center = hFragmFrac->GetBinCenter(i);
        TText *label = new TText(bin_center, -0.1, Form("%g", labels[i - 1]));
        label->SetTextAlign(22); // Center align
        label->SetTextSize(0.03); // Adjust text size as needed
        label->Draw();
    }

hFragmFrac_cms->Draw("ep ");
hFragmFrac_cms->SetMarkerStyle(kFullCircle);
hFragmFrac_cms->SetLineColor(kBlack);
hFragmFrac_cms->SetMarkerColor(kBlack);
hFragmFrac_cms->SetMarkerSize(1.5);
hFragmFrac_cms->GetXaxis()->SetTitle("p_{T}");
hFragmFrac_cms->SetStats(0);
TLatex textY;
textY.DrawLatexNDC(0.8,0.8, "|y|<2.5");
hFragmFrac_cms->GetYaxis()->SetRangeUser(0.15,0.30);
hFragmFrac_cms->SetTitle("");
c_fsfu->SaveAs("FragmFractio.pdf");




}
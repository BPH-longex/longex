 #define ptfixed
//#define etafixed

void plots2(){

   double BF_bupsik = 1.010E-3;
    double BF_bupsik_err = 0.028 / 1.010;
    double BF_bspsiphi = 1.08E-3 * 0.492;
    double BF_bspsiphi_err = sqrt(pow(0.08 / 1.08, 2) + pow(0.5 / 49.2, 2));
      double Factor = (BF_bupsik/BF_bspsiphi);
      std::cout<<"factor" << Factor<<std::endl;
  
 double labels[] = {10,12,14,16,18,20,24,30,40,60,120};
    int n_labels = sizeof(labels) / sizeof(labels[0]);

TFile *fInput;
TH1D *hFragmFrac;
#ifdef ptfixed
fInput = TFile::Open("output2_eta.root");
#elif defined etafixed
fInput = TFile::Open("output2_eta.root");
#endif


TH1D *htemp =(TH1D*) fInput->Get("fragm_ratio");

const Double_t ptBin[11] = {10,12,14,16,18,20,24,30,40,60};
const Double_t etaBin[7] = {0.,0.2,0.4,0.8,1.2,2.5,4};
#ifdef ptfixed
hFragmFrac = new TH1D("hFragmFrac", "hFragmFrac", 5, etaBin);
#elif defined etafixed
hFragmFrac = new TH1D("hFragmFrac", "hFragmFrac", 9, ptBin);
#endif
TH1D *hFragmFrac_cms = new TH1D("hFragmFrac_cms", "hFragmFrac_cms", 9, ptBin);
TLatex text;

double PaperMeas[12] = {0.1314,0.1196,0.1165,0.1154,0.1135,0.1106,0.1105,0.1110,0.1091,0.1095,0.1088,0.1117};
double PaperMeasStat[12] = {2.1,1.6,1.3,1.2,0.8,0.8,0.7,0.8,0.9,0.9,0.9,1.3};
double PaperMeasSys[12] = {3.1,2.7,2.4,2.6,2.6,2.8,2.9,2.6,3.2,2.3,2.8,2.6};
#ifdef ptfixed
for(int i=0; i<9; i++){
#else
for(int i=0; i<13; i++){
#endif
  if(!std::isnan(htemp->GetBinContent(i+1))){
    hFragmFrac->SetBinContent(i+1, htemp->GetBinContent(i+1));
    hFragmFrac->SetBinError(i+1, htemp->GetBinError(i+1));
    }
    std::cout<<hFragmFrac->GetBinContent(i+1)<<std::endl;
 //   hFragmFrac_cms->SetBinContent(i+1, PaperMeas[i]);
 //  hFragmFrac_cms->SetBinError(i+1, sqrt(TMath::Power(PaperMeasStat[i]*PaperMeas[i]*0.01,2) + TMath::Power(PaperMeasSys[i]*PaperMeas[i]*0.01,2)));
}
//hFragmFrac_cms->Scale(Factor);

TCanvas *c_fsfu = new TCanvas("c_fsfu", "c_fsfu", 2000,1500);
TPad *padMain = (TCanvas*) c_fsfu->cd();
  padMain->Draw();
  padMain->SetTopMargin(0.075);
  padMain->SetLeftMargin(0.15);
  padMain->SetRightMargin(0.03);
  hFragmFrac->GetYaxis()->SetRangeUser(0.,0.40);
//  padMain->SetLogy();
  padMain->SetLogx(0);  
 // padMain->SetTicks(1, 1);


  padMain->cd();   

hFragmFrac->SetFillColor(68); // Fill color for the error boxes
hFragmFrac->SetFillStyle(1111); // Solid fill
hFragmFrac->SetLineColor(kBlack);
hFragmFrac->SetMarkerColor(kBlue);

hFragmFrac->Draw("E2 ][ p");
hFragmFrac->SetMarkerSize(2);
 TLatex textCMS;
 textCMS.SetTextSize(0.04);
 textCMS.DrawLatexNDC(0.15,0.94,"CMS Data(Private Work), #sqrt{s}=13 TeV");
hFragmFrac->SetMarkerStyle(kFullCircle);
hFragmFrac->SetTitle("");
hFragmFrac->GetYaxis()->SetTitle("R_{s}");
hFragmFrac->GetYaxis()->SetTitleSize(0.08);
hFragmFrac->GetYaxis()->SetTitleOffset(0.8);
#ifdef ptfixed
hFragmFrac->GetXaxis()->SetTitle("#eta");
#else
hFragmFrac->GetXaxis()->SetTitle("p_{T}");
#endif
hFragmFrac->GetXaxis()->SetTitleOffset(0.55);
hFragmFrac->GetXaxis()->SetTitleSize(0.05);
hFragmFrac->SetMarkerSize(2);

TLatex textResult;
textResult.SetTextSize(0.03);
#ifdef etafixed
textResult.DrawLatexNDC(0.5,0.6, "Our Result: R_{s} = 0.130 #pm 0.008 ");
#else
textResult.DrawLatexNDC(0.5,0.6, "Our Result: R_{s} = 0.14 #pm 0.01 ");
#endif

TF1 *fit = new TF1("fit", "[0] + [1]*x", 0, 300);
fit->SetParameter(0, 0.13);
fit->SetParameter(1, 0);
//hFragmFrac->Fit("fit", "", "0", 0,1.2);


  TLatex textEta;
#ifdef ptfixed
 textEta.SetTextSize(0.04);
 textEta.DrawLatexNDC(0.6,0.8, "10 GeV <|p_{T}|<120 GeV");
#else
 textEta.DrawLatexNDC(0.8,0.8, "|#eta|<2.5");
#endif
 TLatex textDecay;
 textDecay.DrawLatexNDC(0.2,0.85, "B_{s}->J/#Psi #Phi");
 TLatex textDecay2;
 textDecay2.DrawLatexNDC(0.2,0.75,"B^{+}->J/#Psi K^{+}");
hFragmFrac->SetStats(0);
#ifdef etafixed
c_fsfu->SaveAs("FragmFractio2_pt.pdf");
#else
c_fsfu->SaveAs("FragmFractio2_eta.pdf");
#endif

/*

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
*/




}
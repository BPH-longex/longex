void plots(){

TFile *fInput = TFile::Open("output.root");
TH1D *hFragmFrac =(TH1D*) fInput->Get("fragm_ratio");
TCanvas *c_fsfu = new TCanvas("c_fsfu", "c_fsfu", 2000,1500);
hFragmFrac->Draw("ep");
c_fsfu->SaveAs("FragmFractio.pdf");




}
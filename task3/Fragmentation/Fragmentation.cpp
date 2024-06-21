#include "common.h"
#include "def_variables.h"
#include "FitJpsiK.h"
#include "FitJpsiPhi.h"
#include "ComputeFragmFrac.h"

void Fragmentation(){

for(int i=0; i<(number_BinPt-1); i++){
    FitJpsiPhi(ptBin[i], ptBin[i+1], 0, 2.5, i,false);
    FitJpsiK(ptBin[i], ptBin[i+1], 0, 2.5, i, false);
}

/*
for(int i=0; i<(number_BinEta-1); i++){
    FitJpsiPhi(20,120, etaBin[i], etaBin[i+1], i, true);
   FitJpsiK(20,120, etaBin[i], etaBin[i+1], i, true);
}
*/




for(int i=0; i<number_BinPt; i++){
    std::cout<<"*********************************************"<<std::endl;
    std::cout<<" FIXED ETA: " << std::endl;
    std::cout<<"******* Yield for B-->Jpsi Phi************" << std::endl;
    std::cout<<" Bin: "<<i+1 << "   Yield: " << hYield_phi_fixedEta->GetBinContent(i+1) <<std::endl;
    std::cout<<"         Effi: "<<hEffi_phi_fixedEta->GetBinContent(i+1);
    std::cout<<"******* Yield for B-->Jpsi K ************" << std::endl;
    std::cout<<" Bin: "<<i+1 << "   Yield: " << hYield_K_fixedEta->GetBinContent(i+1) <<std::endl;
    std::cout<<"         Effi: "<<hEffi_K_fixedEta->GetBinContent(i+1);
}

for(int ii=0; ii<number_BinEta; ii++){
    std::cout<<"*********************************************"<<std::endl;
    std::cout<<" FIXED PT: " << std::endl;
    std::cout<<"******* Yield for B-->Jpsi Phi************" << std::endl;
    std::cout<<" Bin: "<<ii+1 << "   Yield: " << hYield_phi_fixedPt->GetBinContent(ii+1) <<std::endl;
    std::cout<<"         Effi: "<<hEffi_phi_fixedPt->GetBinContent(ii+1);
    std::cout<<"******* Yield for B-->Jpsi K ************" << std::endl;
    std::cout<<" Bin: "<<ii+1 << "   Yield: " << hYield_K_fixedPt->GetBinContent(ii+1) <<std::endl;
    std::cout<<"         Effi: "<<hEffi_K_fixedPt->GetBinContent(ii+1);
}


// B->Jpsi phi //
TCanvas *c_fixedEta_Phi = new TCanvas();
hYield_phi_fixedEta->Draw("ep");
hYield_phi_fixedEta->GetYaxis()->SetTitle("N(B_{s}->J/#Psi#Phi)");
hYield_phi_fixedEta->GetXaxis()->SetTitle("p_T");
hYield_phi_fixedEta->GetXaxis()->SetNdivisions(10);
for(int i=0; i<number_BinPt; i++){
    hYield_K_fixedEta->GetXaxis()->SetBinLabel(i+1,(std::to_string(ptBin[i])).c_str());
}
c_fixedEta_Phi->SetLogy(1);
c_fixedEta_Phi->SaveAs("Yield_phi_fixedEta.png");

TCanvas *c_fixedPt_Phi = new TCanvas();
hYield_phi_fixedPt->Draw("ep");
hYield_phi_fixedPt->GetYaxis()->SetTitle("N(B_{s}->J/#Psi#Phi)");
hYield_phi_fixedPt->GetXaxis()->SetTitle("#eta");
hYield_phi_fixedPt->GetXaxis()->SetNdivisions(8);
for(int i=0; i<number_BinEta; i++){
    hYield_phi_fixedPt->GetXaxis()->SetBinLabel(i+1,(std::to_string(etaBin[i])).c_str());
}
c_fixedPt_Phi->SetLogy(1);
c_fixedPt_Phi->SaveAs("Yield_phi_fixedPt.png");

TCanvas *c_fixedPt_Effi_Phi = new TCanvas();
hEffi_phi_fixedPt->Draw("ep");
c_fixedPt_Effi_Phi->SaveAs("Effi_fixedPt_Phi.png");
TCanvas *c_fixedEta_Effi_Phi = new TCanvas();
hEffi_phi_fixedEta->Draw("ep");
c_fixedEta_Effi_Phi->SaveAs("Effi_fixedEta_Phi.png");


// ********************* //
// B->Jpsi K //
TCanvas *c_fixedEta_K = new TCanvas();
hYield_K_fixedEta->Draw("ep");
hYield_K_fixedEta->GetYaxis()->SetTitle("N(B_{s}->J/#Psi K)");
hYield_K_fixedEta->GetXaxis()->SetTitle("p_T");
hYield_phi_fixedEta->GetXaxis()->SetNdivisions(10);
c_fixedEta_K->SetLogy(1);
for(int i=0; i<number_BinPt; i++){
    hYield_K_fixedEta->GetXaxis()->SetBinLabel(i+1,(std::to_string(ptBin[i])).c_str());
}
c_fixedEta_K->SaveAs("Yield_K_fixedEta.png");

TCanvas *c_fixedPt_K = new TCanvas();
hYield_K_fixedPt->Draw("ep");
hYield_K_fixedPt->GetYaxis()->SetTitle("N(B_{s}->J/#Psi K)");
hYield_K_fixedPt->GetXaxis()->SetTitle("#eta");
hYield_phi_fixedPt->GetXaxis()->SetNdivisions(8);
for(int i=0; i<number_BinEta; i++){
    hYield_K_fixedPt->GetXaxis()->SetBinLabel(i+1,(std::to_string(etaBin[i])).c_str());
}
c_fixedPt_K->SetLogy(1);
c_fixedPt_K->SaveAs("Yield_K_fixedPt.png");

TCanvas *c_fixedPt_Effi_K = new TCanvas();
hEffi_K_fixedPt->Draw("ep");
c_fixedPt_Effi_K->SaveAs("Effi_fixedPt_K.png");
TCanvas *c_fixedEta_Effi_K = new TCanvas();
hEffi_K_fixedEta->Draw("ep");
c_fixedEta_Effi_K->SaveAs("Effi_fixedEta_K.png");

ComputeFragmFrac(hEffi_K_fixedEta , hEffi_phi_fixedEta, hYield_K_fixedEta, hYield_phi_fixedEta, false); // effiK, effiPhi, YieldK, YieldPhi, 
//ComputeFragmFrac(hEffi_phi_fixedEta, hYield_phi_fixedEta, false);
}
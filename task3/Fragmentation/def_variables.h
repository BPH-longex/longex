const int number_BinPt=11;
const int number_BinEta = 7;

TH1D *hEffi_K_fixedEta = new TH1D("hEffi_K","hEffi_K", number_BinPt, 0, number_BinPt);
TH1D *hEffi_phi_fixedEta = new TH1D("hEffi_K","hEffi_K", number_BinPt, 0, number_BinPt);
TH1D *hYield_K_fixedEta =new TH1D("hYield_K_fixedEta","hYield_K_fixedEta", number_BinPt, 0, number_BinPt);
TH1D *hYield_phi_fixedEta = new TH1D("hYield_phi_fixedEta","hYield_phi_fixedEta", number_BinPt, 0, number_BinPt);

TH1D *hEffi_K_fixedPt = new TH1D("hEffi_K","hEffi_K", number_BinEta, 0, number_BinEta);
TH1D *hEffi_phi_fixedPt = new TH1D("hEffi_K","hEffi_K", number_BinEta, 0, number_BinEta);
TH1D *hYield_K_fixedPt =new TH1D("hYield_K_fixedPt","hYield_K_fixedPt", number_BinEta, 0, number_BinEta);
TH1D *hYield_phi_fixedPt = new TH1D("hYield_phi_fixedPt","hYield_phi_fixedPt", number_BinEta, 0, number_BinEta);

TH1D *hFrag_s_fixedEta = new TH1D("hFrag_s_fixedEta", "hFrag_s_fixedEta", number_BinPt, 0, number_BinPt);
TH1D *hFrag_s_fixedPt = new TH1D("hFrag_s_fixedPt", "hFrag_s_fixedPt", number_BinEta, 0, number_BinEta);

TH1D *hFrag_u_fixedEta = new TH1D("hFrag_u_fixedEta", "hFrag_u_fixedEta", number_BinPt, 0, number_BinPt);
TH1D *hFrag_u_fixedPt = new TH1D("hFrag_u_fixedPt", "hFrag_u_fixedPt", number_BinEta, 0, number_BinEta);

TH1F *fragm_ratio;
TH1F *fs;
TH1F *fu;


double ptBin[number_BinPt] = {10,12,14,16,18,20,24,30,40,60,120};
double etaBin[number_BinEta] = {0.,0.2,0.4,0.8,1.2,2.5,4};

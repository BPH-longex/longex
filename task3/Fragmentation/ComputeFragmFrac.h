void ComputeFragmFrac(TH1D *hEffiK, Th1D *hEffiPhi, TH1D *hYieldPsi, TH1D *hYieldPhi, Bool_t FixedPt){
    using namespace RooFit;	
    uint cate = 2;
    double pt_min = 20.0;
    double pt_max = 30.0;
    double eta_min = 0.;
    double eta_max = 2.5;
    
    //****** branching fraction ********//
    double BF_bupsik = 1.010E-3;
    double BF_bupsik_err = 0.028/1.010;
    double BF_bspsiphi = 1.08E-3 * 0.492;
    double BF_bspsiphi_err = sqrt(pow(0.08/1.08,2) + pow(0.5/49.2,2));

    TH1D *fs = new TH1D("fs", "fs", number_BinPt, 0, number_BinPt);
    TH1D *fu = new TH1D("fu", "fu", number_etaBin, 0, number_etaBin);

    /*
     From task3_1.cc
     Category: 2
     pt range: 20, 30
     |eta| range: 0, 2.5
     Selection efficiency: 0.32762 +- 0.0014842
     Observed yield: 106740 +- 782.594     

     From task3_2.cc
     Category: 2
     pt range: 20, 30
     |eta| range: 0, 2.5
     Selection efficiency: 0.39079 +- 0.00154296
     Observed yield: 7593.83 +- 207.071     */

    double ret_bupsik[4], ret_bspsiphi[4];
    
    ret_bupsik[0] = 0.32762;
    ret_bupsik[1] = 0.0014842/ret_bupsik[0];
    ret_bupsik[2] = 106740;
    ret_bupsik[3] = 782.594/ret_bupsik[2];
    
    ret_bspsiphi[0] = 0.39079;
    ret_bspsiphi[1] = 0.00154296/ret_bspsiphi[0];
    ret_bspsiphi[2] = 7593.83;
    ret_bspsiphi[3] = 207.071/ret_bspsiphi[2];
    
    double baseeff_bupsik = 0.001267;
    double baseeff_bspsiphi = 0.0005476;

    if(!fixedPt){
     for(int ptBin=0; ptBin<number_etaBin; ptBin++){
        fs->SetBinContent(ptBin+1,hYieldPhi->GetBinContent(ptBin+1)/(hEffiPhi->GetBinContent(ptBin+1)/BF_bspsiphi; 
        fu->SetBinContent(ptBin+1, hYieldK->GetBinContent(ptBin+1)/(hEffiK->GetBinContent(ptBin+1)/BF_bspsiphi; 
     }
    }
    else if(fixedPt){
        for(int etaBin=0; etaBin<number_EtaPt; ptBin++){
        fs->SetBinContent(etaBin+1, hYieldPhi->GetBinContent(etaBin+1)/(hEffiPhi->GetBinContent(etaBin+1)/BF_bspsiphi/baseeff_bspsiphi); 
        fu->SetBinContent(etaBin+1, YieldK->GetBinContent(etaBin+1)/(hEffiK->GetBinContent(etaBin+1)/BF_bspsiphi/baseeff_bupsik; 
     }
    }
    TH1D *fragm_ratio;
    fragm_ratio =(TH1D*) fs->Clone("fragm_ratio");
    fragm_ratio->Divide(fu);
    for(int Bin=0; Bin<fragm_ratio->GeNbinsX()+1; Bin++){
    fragm_ratio->SetBinError(Bin+1, fragm_ratio->GetBinContent(Bin+1)*sqrt(pow(BF_bupsik_err,2)+pow(BF_bspsiphi_err,2));
    }
 
    
    
}


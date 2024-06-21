void ComputeFragmFrac(TH1D *hEffiK, TH1D *hEffiPhi, TH1D *hYieldK, TH1D *hYieldPhi, Bool_t fixedPt){

    using namespace RooFit;
    uint cate = 2;
    double pt_min = 20.0;
    double pt_max = 30.0;
    double eta_min = 0.;
    double eta_max = 2.5;

    double baseeff_bupsik = 0.001267;
    double baseeff_bspsiphi = 0.0005476;

    //****** branching fraction ********//
    double BF_bupsik = 1.010E-3;
    double BF_bupsik_err = 0.028 / 1.010;
    double BF_bspsiphi = 1.08E-3 * 0.492;
    double BF_bspsiphi_err = sqrt(pow(0.08 / 1.08, 2) + pow(0.5 / 49.2, 2));
    // ********************************* //

    if(!fixedPt){
    fragm_ratio = new TH1F("fragm_ratio", "fragm_ratio", number_BinPt, 0, number_BinPt);
    fs = new TH1F("fs", "fs", number_BinPt, 0, number_BinPt);
    fu = new TH1F("fu", "fu", number_BinPt, 0, number_BinPt);
    }



    if(!fixedPt) {
        for (int ptBin = 0; ptBin < number_BinPt; ptBin++) {
            double effiPhi = hEffiPhi->GetBinContent(ptBin + 1);
            double yieldPhi = hYieldPhi->GetBinContent(ptBin + 1);
            double effiK = hEffiK->GetBinContent(ptBin + 1);
            double yieldK = hYieldK->GetBinContent(ptBin + 1);

            std::cout<<yieldK<<"    " << effiK << "   " << yieldPhi << "  " << effiPhi << std::endl;


            if (!std::isnan(effiPhi) && !std::isnan(yieldPhi)) {
                fs->SetBinContent(ptBin + 1, yieldPhi/effiPhi/BF_bspsiphi/baseeff_bspsiphi);
            }
            if (!std::isnan(effiK) && !std::isnan(yieldK)) {
                fu->SetBinContent(ptBin + 1, yieldK/effiK/BF_bupsik/baseeff_bupsik);
            }
        }
    } else if (fixedPt) {
        for (int etaBin = 0; etaBin < number_BinEta; etaBin++) {
            double effiPhi = hEffiPhi->GetBinContent(etaBin + 1);
            double yieldPhi = hYieldPhi->GetBinContent(etaBin + 1);
            double effiK = hEffiK->GetBinContent(etaBin + 1);
            double yieldK = hYieldK->GetBinContent(etaBin + 1);

            std::cout<<yieldK<<"    " << effiK << "   " << yieldPhi << "  " << effiPhi << std::endl;

            if (!std::isnan(effiPhi) && !std::isnan(yieldPhi)) {
                fs->SetBinContent(etaBin + 1, yieldPhi / effiPhi / BF_bspsiphi / baseeff_bspsiphi);
            }
            if (!std::isnan(effiK) && !std::isnan(yieldK)) {
                fu->SetBinContent(etaBin + 1, yieldK / effiK / BF_bupsik / baseeff_bupsik);
            }
        }
    }


    for (int Bin = 0; Bin < (number_BinPt-1); Bin++) {
        fragm_ratio->SetBinContent(Bin+1, fs->GetBinContent(Bin+1)/fu->GetBinContent(Bin+1));
        fragm_ratio->SetBinError(Bin + 1, fragm_ratio->GetBinContent(Bin + 1) * sqrt(pow(BF_bupsik_err, 2) + pow(BF_bspsiphi_err, 2)));
        std::cout<<fragm_ratio->GetBinContent(Bin+1)<<std::endl;
    }
    TCanvas *c_fragmRatioPt = new TCanvas();
    fragm_ratio->Draw("ep");
    fragm_ratio->SetMarkerStyle(kFullCircle);
    fragm_ratio->SetTitle("");
    fragm_ratio->SetStats(0);
    //fragm_ratio->GetXaxis()->SetRangeUser(0,10);
    fragm_ratio->GetYaxis()->SetTitle("#frac{f_{s}}{f_{u}}(p_{T})");
    fragm_ratio->GetXaxis()->SetTitle("p_T");
    TLatex textPhi;
    textPhi.DrawLatexNDC(0.65, 0.7, "f_{s}=f(b->B_{s}) ");
    TLatex textK;
    textK.DrawLatexNDC(0.65,0.8, "f_{u}=f(b->B^{+})");
    c_fragmRatioPt->SaveAs("PtDependenceFragment.png");

}

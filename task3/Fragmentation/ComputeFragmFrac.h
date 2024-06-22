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

    double Factor = BF_bupsik/(BF_bspsiphi*0.489);

    if(!fixedPt){
    fragm_ratio = new TH1F("fragm_ratio", "fragm_ratio", number_BinPt, 0, number_BinPt);
    fs = new TH1F("fs", "fs", number_BinPt, 0, number_BinPt);
    fu = new TH1F("fu", "fu", number_BinPt, 0, number_BinPt);
    }
    else if(fixedPt){
    fragm_ratio = new TH1F("fragm_ratio", "fragm_ratio", number_BinEta, 0, number_BinEta);
    fs = new TH1F("fs", "fs", number_BinEta, 0, number_BinEta);
    fu = new TH1F("fu", "fu", number_BinEta, 0, number_BinEta);
    }



    if(!fixedPt) {
        for (int ptBin = 0; ptBin < number_BinPt; ptBin++) {
            double effiPhi = hEffiPhi->GetBinContent(ptBin + 1);
            double yieldPhi = hYieldPhi->GetBinContent(ptBin + 1);
            double effiK = hEffiK->GetBinContent(ptBin + 1);
            double yieldK = hYieldK->GetBinContent(ptBin + 1);

            std::cout<<yieldK<<"    " << effiK << "   " << yieldPhi << "  " << effiPhi << std::endl;


            if (!std::isnan(effiPhi) && !std::isnan(yieldPhi)) {
                fs->SetBinContent(ptBin + 1, yieldPhi/effiPhi /baseeff_bspsiphi);
                fs->SetBinError(ptBin+1, fs->GetBinContent(ptBin+1)*sqrt(pow(hEffiPhi->GetBinError(ptBin+1)/hEffiPhi->GetBinContent(ptBin+1),2) + pow(hYieldPhi->GetBinError(ptBin+1)/hYieldPhi->GetBinContent(ptBin+1),2)));
            }
            if (!std::isnan(effiK) && !std::isnan(yieldK)) {
                fu->SetBinContent(ptBin + 1, yieldK/effiK / baseeff_bupsik);
                fu->SetBinError(ptBin+1, fu->GetBinContent(ptBin+1)*sqrt(pow(hEffiK->GetBinError(ptBin+1)/hEffiK->GetBinContent(ptBin+1),2) + pow(hYieldK->GetBinError(ptBin+1)/hYieldK->GetBinContent(ptBin+1),2)));

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
                fs->SetBinContent(etaBin + 1, yieldPhi / effiPhi / baseeff_bspsiphi );
                fs->SetBinError(etaBin+1, fs->GetBinContent(etaBin+1)*sqrt(pow(hEffiPhi->GetBinError(etaBin+1)/hEffiPhi->GetBinContent(etaBin+1),2) + pow(hYieldPhi->GetBinError(etaBin+1)/hYieldPhi->GetBinContent(etaBin+1),2)));

            }
            if (!std::isnan(effiK) && !std::isnan(yieldK)) {
                fu->SetBinContent(etaBin + 1, yieldK / effiK / baseeff_bupsik);
                fs->SetBinError(etaBin+1, fs->GetBinContent(etaBin+1)*sqrt(pow(hEffiK->GetBinError(etaBin+1)/hEffiK->GetBinContent(etaBin+1),2) + pow(hYieldK->GetBinError(etaBin+1)/hYieldK->GetBinContent(etaBin+1),2)));

            }
        std::cout<<fs->GetBinContent(etaBin+1)<<std::endl;
        }
    }

    if(!fixedPt){
    for (int Bin = 0; Bin < (number_BinPt-1); Bin++) {
        fragm_ratio->SetBinContent(Bin+1, fs->GetBinContent(Bin+1)/fu->GetBinContent(Bin+1));
        fragm_ratio->SetBinError(Bin + 1, fragm_ratio->GetBinContent(Bin + 1) * sqrt(pow(BF_bupsik_err, 2) + pow(BF_bspsiphi_err, 2)));
        std::cout<<fragm_ratio->GetBinContent(Bin+1)<<std::endl;
    }
    }
    else if(fixedPt){
        for (int Bin = 0; Bin < (number_BinEta-1); Bin++) {
        fragm_ratio->SetBinContent(Bin+1, fs->GetBinContent(Bin+1)/fu->GetBinContent(Bin+1));
        fragm_ratio->SetBinError(Bin + 1, fragm_ratio->GetBinContent(Bin + 1) * sqrt(pow(BF_bupsik_err, 2) + pow(BF_bspsiphi_err, 2)));
        std::cout<<fragm_ratio->GetBinContent(Bin+1)<<std::endl;
    }
    }
    

  //  char const *range[number_BinPt-1] = {"10", "11", "12","13","14","15","16","18","20","23","26","29","34","45","70","80","100","120"};
   // fragm_ratio->Scale(1./Factor);
    TCanvas *c_fragmRatioPt = new TCanvas();
    fragm_ratio->Draw("ep BOX");
    fragm_ratio->SetMarkerStyle(kFullCircle);
    fragm_ratio->SetTitle("");
    gStyle->SetErrorX(0);
    TLatex textCMS;
    textCMS.DrawLatexNDC(0.15, 0.92, "CMS Data (Private Work)");
    //fragm_ratio->GetXaxis()->SetNdivisions(11);
    fragm_ratio->SetStats(0);
    //fragm_ratio->GetXaxis()->SetRangeUser(0,11);
    auto *x = fragm_ratio->GetXaxis();
   // for(int i=1; i<number_BinPt; i++){
   //     x->ChangeLabel(i+1, -1,-1,-1,-1,-1, range[i-1]);
   // }
    fragm_ratio->GetYaxis()->SetRangeUser(0, 0.40);
    fragm_ratio->GetYaxis()->SetTitle("f_{s}/f_{u}");
    fragm_ratio->GetXaxis()->SetTitle("p_{T}");
    TLatex textPhi;
    textPhi.DrawLatexNDC(0.65, 0.7, "#font[52]{f_{s}=f(b->B_{s})} ");
    TLatex textK;
    textK.DrawLatexNDC(0.65,0.8, "#font[52]{f_{u}=f(b->B^{+})}");
    //fragm_ratio->Fit("pol1");
    c_fragmRatioPt->SaveAs("PtDependenceFragment.png");

  


}

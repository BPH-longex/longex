{
    using namespace RooFit;

    RooRealVar m("m","",5.0,5.8);
    RooRealVar wgt("wgt","",1.,0.,1000.);
    RooDataSet *rds_data = new RooDataSet("rds_data","",RooArgSet(m,wgt),"wgt");

    TFile *fin = new TFile("/eos/user/c/cmsdas/2024/long-ex-bph/bupsikData.root");
    TTree *tin = (TTree*)fin->Get("bupsikData");

    unsigned int cate_t;
    float wgt_t,m_t;
    tin->SetBranchAddress("cate",&cate_t);
    tin->SetBranchAddress("wgt",&wgt_t);
    tin->SetBranchAddress("m",&m_t);
    
    for(int evt=0; evt<tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (cate_t!=0) continue;
        if (m_t<5.0 || m_t>=5.8) continue;
        m.setVal(m_t);
        wgt.setVal(wgt_t);
        rds_data->add(RooArgSet(m,wgt),wgt_t);
    }
    delete fin;

    // Parameters dump from from task 2a 
    // sig_frac = 0.926022; 
    // sig_mean1 = 5.27971;    
    // sig_mean2 = 5.26646;    
    // sig_sigma1 = 0.0156007;
    // sig_sigma2 = 0.0627436;
    
    
    RooRealVar sig_shift("sig_shift","",0.,-0.02,0.02);
    RooRealVar sig_scale("sig_scale","",1.,0.8,1.2);

    RooRealVar sigmc_mean1("sigmc_mean1","",5.27971);
    RooRealVar sigmc_mean2("sigmc_mean2","",5.26646);
    RooRealVar sigmc_sigma1("sigmc_sigma1","",0.0156007);
    RooRealVar sigmc_sigma2("sigmc_sigma2","",0.0627436);

    //RooRealVar sig_mean1("sig_mean1","", 5.28, 5.2, 5.4);
    //RooRealVar sig_mean2("sig_mean2","", 5.28, 5.2, 5.4);
    //RooRealVar sig_sigma1("sig_sigma1","", 0.030, 0.005, 0.060);
    //RooRealVar sig_sigma2("sig_sigma2","", 0.080, 0.040, 0.200);
    RooAddition sig_mean1("sig_mean1", "", RooArgList(sigmc_mean1, sig_shift));
    RooAddition sig_mean2("sig_mean2", "", RooArgList(sigmc_mean2, sig_shift));
    RooProduct sig_sigma1("sig_sigma1", "", RooArgList(sigmc_sigma1, sig_scale));
    RooProduct sig_sigma2("sig_sigma2", "", RooArgList(sigmc_sigma2, sig_scale));

    RooRealVar sig_frac("sig_frac","",0.9,0.5,1.0);
    RooGaussian sig_g1("sig_g1", "", m, sig_mean1, sig_sigma1);
    RooGaussian sig_g2("sig_g2", "", m, sig_mean2, sig_sigma2);
    RooAddPdf pdf_sig("pdf_sig", "", RooArgList(sig_g1,sig_g2), RooArgList(sig_frac));

    RooRealVar comb_coeff("comb_coeff","",-1.2,-10.,10.);
    RooExponential pdf_comb("pdf_comb","",m,comb_coeff);

    RooRealVar jpsix_scale("jpsix_scale","",0.02,0.001,0.08);
    RooRealVar jpsix_shift("jpsix_shift","",5.13,5.12,5.16);
    RooGenericPdf pdf_jpsix("pdf_jpsix","","TMath::Erfc((@0-@1)/@2)",RooArgList(m,jpsix_shift,jpsix_scale));

    RooRealVar n_sig("n_sig","",100000,0.,1E8);
    RooRealVar n_comb("n_comb","",80000,0.,1E6);
    RooRealVar n_jpsix("n_jpsix","",20000,0.,1E5);
    RooAddPdf model("model","",RooArgList(pdf_sig,pdf_comb,pdf_jpsix),RooArgList(n_sig,n_comb,n_jpsix));

    model.fitTo(*rds_data, Extended(true), SumW2Error(true));

    RooPlot *frame = m.frame();
    rds_data->plotOn(frame, Name("rds_data"));
    //pdf_sig.plotOn(frame, Name("sig_g1"), Components(RooArgList(sig_g1)), LineColor(kRed));
    model.plotOn(frame, Name("pdf_comb"), Components(RooArgList(pdf_comb)), LineColor(kGreen));
    model.plotOn(frame, Name("pdf_sig"), Components(RooArgList(pdf_sig)), LineColor(kRed));
    model.plotOn(frame, Name("pdf_jpsix"), Components(RooArgList(pdf_jpsix)), LineColor(kBlue));
    model.plotOn(frame);

    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);

    canvas->Divide(1,2);

    canvas->cd(1);
    frame->Draw();

    canvas->cd(2);
    RooHist* residHist = frame->pullHist();
    RooPlot *frame_resid = m.frame();
    frame_resid->addPlotable(residHist, "P");
  
    frame->GetYaxis()->SetTitle("Entries / 0.01 GeV");
    frame->GetXaxis()->SetTitle("M(#mu#muK) [GeV]");

    frame_resid->Draw(); 

    canvas->SaveAs("out.png");

}

{
    using namespace RooFit;

    unsigned int cate = 0;

    RooRealVar m("m","",5.0,5.8);
    RooRealVar wgt("wgt","",1.,0.,1000.);

    // --- DATA --- //

    RooDataSet *rds_data = new RooDataSet("rds_data", "", RooArgSet(m,wgt), "wgt");

    TFile *fin = new TFile("/eos/user/c/cmsdas/2024/long-ex-bph/bupsikData.root");
    TTree *tin = (TTree*)fin->Get("bupsikData");

    unsigned int cate_t;
    float wgt_t,m_t;
    tin->SetBranchAddress("cate",&cate_t);
    tin->SetBranchAddress("wgt",&wgt_t);
    tin->SetBranchAddress("m",&m_t);
    
    for(int evt=0; evt<tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (cate_t!=cate) continue;
        if (m_t<5.0 || m_t>=5.8) continue;
        m.setVal(m_t);
        wgt.setVal(wgt_t);
        rds_data->add(RooArgSet(m,wgt),wgt_t);
    }
    delete fin;


    // --- MC --- //

    RooDataSet *rds_mc = new RooDataSet("rds_mc", "", RooArgSet(m));

    fin = new TFile("/eos/user/c/cmsdas/2024/long-ex-bph/bupsikMc.root");
    tin = (TTree*)fin->Get("bupsikMc");

    tin->SetBranchAddress("cate",&cate_t);
    tin->SetBranchAddress("m",&m_t);
    
    for(int evt=0; evt<tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (cate_t!=cate) continue;
        if (m_t<5.0 || m_t>=5.8) continue;
        m.setVal(m_t);
        rds_mc->add(RooArgSet(m));
    }
    delete fin;

    // --- param from MC --- //  

	RooRealVar sigmc_mean1("sigmc_mean1","",5.28,5.2,5.4);
	RooRealVar sigmc_mean2("sigmc_mean2","",5.28,5.2,5.4);
	RooRealVar sigmc_sigma1("sigmc_sigma1","",0.030,0.005,0.060);
	RooRealVar sigmc_sigma2("sigmc_sigma2","",0.080,0.040,0.200);
	RooRealVar sig_frac("sig_frac","",0.9,0.5,1.0);
	RooGaussian sigmc_g1("sigmc_g1","",m,sigmc_mean1,sigmc_sigma1);
	RooGaussian sigmc_g2("sigmc_g2","",m,sigmc_mean2,sigmc_sigma2);
	RooAddPdf pdf_sigmc("pdf_sig","",RooArgList(sigmc_g1,sigmc_g2),RooArgList(sig_frac));

	pdf_sigmc.fitTo(*rds_mc);

    // --- draw MC and fit --- //
    RooPlot *frame_mc = m.frame();
    rds_mc->plotOn(frame_mc, Name("rds_mc"));
    pdf_sigmc.plotOn(frame_mc, Name("sigmc_g1"), Components(RooArgList(sigmc_g1)), LineColor(kRed));
    pdf_sigmc.plotOn(frame_mc, Name("sigmc_g2"), Components(RooArgList(sigmc_g2)), LineColor(kGreen));
    pdf_sigmc.plotOn(frame_mc);

    TCanvas* canvas_mc = new TCanvas("canvas_mc", "", 600, 600);

    canvas_mc->Divide(1,2);

    canvas_mc->cd(1);
    frame_mc->Draw();

    // Drwaing the pull below
    canvas_mc->cd(2);
    RooHist* residHist_mc = frame_mc->pullHist();
    RooPlot *frame_mc_resid = m.frame();
    frame_mc_resid->addPlotable(residHist_mc, "P");
  
    frame_mc->GetYaxis()->SetTitle("Entries / 0.01 GeV");
    frame_mc->GetXaxis()->SetTitle("M(#mu#muK) [GeV]");

    frame_mc_resid->Draw(); 

    canvas_mc->SaveAs("mcout.pdf");
 

    // --- fit the data -- //

    RooRealVar sig_shift("sig_shift","",0.,-0.02,0.02);
    RooRealVar sig_scale("sig_scale","",1.,0.8,1.2);

    RooAddition sig_mean1("sig_mean1", "", RooArgList(sigmc_mean1, sig_shift));
    RooAddition sig_mean2("sig_mean2", "", RooArgList(sigmc_mean2, sig_shift));
    RooProduct sig_sigma1("sig_sigma1", "", RooArgList(sigmc_sigma1, sig_scale));
    RooProduct sig_sigma2("sig_sigma2", "", RooArgList(sigmc_sigma2, sig_scale));

    RooGaussian sig_g1("sig_g1", "", m, sig_mean1, sig_sigma1);
    RooGaussian sig_g2("sig_g2", "", m, sig_mean2, sig_sigma2);
    RooAddPdf pdf_sig("pdf_sig", "", RooArgList(sig_g1,sig_g2), RooArgList(sig_frac));

    RooRealVar comb_coeff("comb_coeff", "", -1.2, -10., 10.);
    RooExponential pdf_comb("pdf_comb", "", m, comb_coeff);

    RooRealVar jpsix_scale("jpsix_scale", "", 0.02, 0.001, 0.08);
    RooRealVar jpsix_shift("jpsix_shift", "", 5.13, 5.12, 5.16);
    RooGenericPdf pdf_jpsix("pdf_jpsix", "", "TMath::Erfc((@0-@1)/@2)", RooArgList(m,jpsix_shift, jpsix_scale));

    RooRealVar n_sig("n_sig","", 100000, 0., 1E8);
    RooRealVar n_comb("n_comb", "", 80000, 0., 1E6);
    RooRealVar n_jpsix("n_jpsix", "", 20000, 0., 1E5);
    RooAddPdf model("model", "", RooArgList(pdf_sig, pdf_comb, pdf_jpsix), RooArgList(n_sig, n_comb, n_jpsix));

    model.fitTo(*rds_data, Extended(true), SumW2Error(true));

    RooWorkspace wspace("wspace`:", "wspace");
    //wspace.import(*rds_data, model, comb_coeff, pdf_comb, sig_g1, sig_g2, pdf_sig,n_jpsix, n_comb, n_sig);
    wspace.import(*rds_data);
    wspace.import(model);
    wspace.import(comb_coeff);
    wspace.import(pdf_comb);
    wspace.import(sig_g1);
    wspace.import(sig_g2);
    wspace.import(pdf_sig);
    wspace.import(n_jpsix);
    wspace.import(n_comb);
    wspace.import(n_sig);

    TFile *fout = new TFile("jpsik.root","RECREATE");
    wspace.Write();
    fout->Close();

    // --- draw data and fit --- //
    RooPlot *frame_data = m.frame();
    rds_data->plotOn(frame_data, Name("rds_data"));
    model.plotOn(frame_data, Name("pdf_sig"), Components(RooArgList(pdf_sig)), LineColor(kRed));
    model.plotOn(frame_data, Name("pdf_comb"), Components(RooArgList(pdf_comb)), LineColor(kGreen));
    model.plotOn(frame_data, Name("pdf_jpsix"), Components(RooArgList(pdf_jpsix)), LineColor(kMagenta));
    model.plotOn(frame_data);

    TCanvas* canvas_data = new TCanvas("canvas_data", "", 600, 600);

    canvas_data->Divide(1,2);

    canvas_data->cd(1);
    frame_data->Draw();

    // Drawing the pull below
    canvas_data->cd(2);
    RooHist* residHist_data = frame_data->pullHist();
    RooPlot *frame_data_resid = m.frame();
    frame_data_resid->addPlotable(residHist_data, "P");
  
    frame_data->GetYaxis()->SetTitle("Entries / 0.01 GeV");
    frame_data->GetXaxis()->SetTitle("M(#mu#muK) [GeV]");

    frame_data_resid->Draw(); 

    canvas_data->SaveAs("dataout.pdf");

}

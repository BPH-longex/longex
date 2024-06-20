{
    using namespace RooFit;

	unsigned int cate_t = 0;
	RooRealVar m("m","",5.0,5.8);
	RooRealVar cate("cate","",-1,20);

	TFile *fin = new TFile("/eos/user/c/cmsdas/2024/long-ex-bph/bupsikMc.root");
	TTree *tin = (TTree*)fin->Get("bupsikMc");

	RooDataSet* rds_mc = new RooDataSet("rds_mc", "rds_mc", tin, RooArgSet(cate, m));

	rds_mc = (RooDataSet*)rds_mc->reduce(m, Form("m >= 5.0 && m <= 5.8 && cate==%d", cate_t));

	delete fin;

	RooRealVar sig_mean1("sig_mean1","",5.28,5.2,5.4);
	RooRealVar sig_mean2("sig_mean2","",5.28,5.2,5.4);
	RooRealVar sig_sigma1("sig_sigma1","",0.030,0.005,0.060);
	RooRealVar sig_sigma2("sig_sigma2","",0.080,0.040,0.200);
	RooRealVar sig_frac("sig_frac","",0.9,0.5,1.0);
	RooGaussian sig_g1("sig_g1","",m,sig_mean1,sig_sigma1);
	RooGaussian sig_g2("sig_g2","",m,sig_mean2,sig_sigma2);
	RooAddPdf pdf_sig("pdf_sig","",RooArgList(sig_g1,sig_g2),RooArgList(sig_frac));

	pdf_sig.fitTo(*rds_mc);

    RooPlot *frame = m.frame();
    rds_mc->plotOn(frame, Name("rds_mc"));
    pdf_sig.plotOn(frame, Name("sig_g1"), Components(RooArgList(sig_g1)), LineColor(kRed));
    pdf_sig.plotOn(frame, Name("sig_g2"), Components(RooArgList(sig_g2)), LineColor(kGreen));
    pdf_sig.plotOn(frame);

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

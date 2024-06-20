#include "common.h"

using namespace RooFit;

void fill_dataset(RooWorkspace* wspace, unsigned int cate = 0){
    double bdt_min = 0.8;

    vector<TString> decay;
    vector<double> yield, yield_err;

    TString input_path = "/eos/user/c/cmsdas/2024/long-ex-bph/";

    // Bs -> K mu nu
    decay.push_back("bskmunuMcBg");
    yield.push_back(effyield[cate][N_bskmunuMcBg]);
    yield_err.push_back(effyield[cate][dN_bskmunuMcBg]);

    // B0 -> pi mu mu
    decay.push_back("bdpimunuMcBg");
    yield.push_back(effyield[cate][N_bdpimunuMcBg]);
    yield_err.push_back(effyield[cate][dN_bdpimunuMcBg]);

    // B+ -> pi mu mu
    decay.push_back("bupimumuMcBg");
    yield.push_back(effyield[cate][N_bupimumuMcBg]);
    yield_err.push_back(effyield[cate][dN_bupimumuMcBg]);

    RooRealVar m("m","",4.9,5.9);
    RooRealVar wgt("wgt","",1.,0.,1000.);
    RooDataSet *rds = new RooDataSet("rds","",RooArgSet(m,wgt),"wgt");

    double sum_weight = 0.;
    double sum_weight_err = 0.;
    rds->weightError(RooAbsData::SumW2);

    double weight_max = 0.; // get the largest weight value across samples
    for(int proc=0; proc<(int)decay.size(); proc++) {
        TFile *fin = new TFile(input_path + decay[proc]+".root");
        TTree *tin = (TTree*)fin->Get(decay[proc]);
        double weight = yield[proc]/(double)tin->GetEntries(Form("cate==%d",cate));
        if (weight>weight_max) weight_max = weight;
        delete fin;
    }

    for(int proc=0; proc<(int)decay.size(); proc++) {
        TFile *fin = new TFile(input_path + decay[proc]+".root");
        TTree *tin = (TTree*)fin->Get(decay[proc]);

        unsigned int cate_t;
        float m_t,bdt_t;
        tin->SetBranchAddress("cate",&cate_t);
        tin->SetBranchAddress("m",&m_t);
        tin->SetBranchAddress("bdt",&bdt_t);

        double weight = yield[proc]/(double)tin->GetEntries(Form("cate==%d",cate));
	    double weight_err = yield_err[proc]/(double)tin->GetEntries(Form("cate==%d",cate))/weight_max;

	    std::cout<<"weight "<<weight<<" weight_err"<<weight_err<<std::endl;

        for(int evt=0; evt<tin->GetEntries(); evt++) {
            tin->GetEntry(evt);
            if (cate_t!=cate) continue;
            if (bdt_t<=bdt_min) continue;
            m.setVal(m_t);
	        wgt.setVal(weight/weight_max); // rescale the event with largest weight to be 1
            rds->add(RooArgSet(m,wgt), weight/weight_max);

            sum_weight += weight;
            sum_weight_err += weight_err; // systematics; linear sum
        }

        delete fin;
    }

    cout << "Category: " << cate << endl;
    cout << "BDT min: " << bdt_min << endl;
    cout << "Sum of weights: " << sum_weight << " +- " << sum_weight_err << endl;

    wspace->import(m);
    wspace->import(*rds);
}

void fit(RooWorkspace* wspace, unsigned int cate = 0){

    RooDataSet *rds = (RooDataSet*)wspace->data("rds");
    RooRealVar m = *(RooRealVar*)wspace->var("m");

    RooNDKeysPdf pdf("pdf", "", m, *rds, "ma", 2);

    // pdf.fitTo(*rds, SumW2Error(true));

    wspace->import(pdf);
}

void draw(RooWorkspace* wspace, unsigned int cate = 0){

    // reload workspace
    RooDataSet *rds = (RooDataSet*)wspace->data("rds");
    RooRealVar m = *(RooRealVar*)wspace->var("m");
    RooNDKeysPdf pdf = *(RooNDKeysPdf*)wspace->pdf("pdf");

    RooPlot *frame = m.frame(Title(" "), Bins(70));
    rds->plotOn(frame, Name("t_rds"), MarkerSize(0.8));
    pdf.plotOn(frame, Name("t_pdf"), LineWidth(3));

    TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
    canvas->Divide(1,2);
    canvas->cd(1);
    canvas->SetMargin(0.15,0.06,0.13,0.07);

    frame->GetYaxis()->SetTitleOffset(1.50);
    frame->GetYaxis()->SetTitle("Entries / 0.01 GeV");
    frame->GetXaxis()->SetTitleOffset(1.15);
    frame->GetXaxis()->SetLabelOffset(0.01);
    frame->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
    frame->GetXaxis()->SetTitleSize(0.043);
    frame->GetYaxis()->SetTitleSize(0.043);
    frame->Draw();

    TLegend* leg = new TLegend(0.58,0.77,0.93,0.92);
    leg->SetFillStyle(0);
    leg->SetLineWidth(0);
    leg->SetHeader(Form("Category %d",cate));
    leg->AddEntry(frame->findObject("t_rds"),"Simluation","EP");
    leg->AddEntry(frame->findObject("t_pdf"),"PDF","L");
    leg->Draw();
    canvas->cd(2);

    // Create a RooHist object to plot the ratio between fit and data
    RooHist* pullHist = frame->residHist();
    RooPlot* frame_pull = m.frame();
    frame_pull->addPlotable(pullHist,"P");
    pullHist->SetMarkerStyle(20);
    pullHist->SetMarkerSize(0.8);
    pullHist->SetLineWidth(1);
    pullHist->SetLineColor(kBlack);
    pullHist->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
    pullHist->GetYaxis()->SetTitle("Pull");
    pullHist->GetYaxis()->SetTitleOffset(0.5);
    pullHist->GetYaxis()->SetTitleSize(0.043);
    pullHist->GetXaxis()->SetTitleSize(0.043);
    pullHist->GetXaxis()->SetLabelSize(0.035);
    pullHist->GetYaxis()->SetLabelSize(0.035);
    frame_pull->Draw();

    canvas->Print("task_4_2.pdf");
    canvas->Print("task_4_2.png");
}

void task_4_2(unsigned int cate = 0)
{
    gROOT->SetBatch();

    RooWorkspace *wspace = new RooWorkspace("wspace","wspace");

    fill_dataset(wspace, cate);
    fit(wspace, cate);
    draw(wspace, cate);

    // save RooFit workspace
    TFile *fout = new TFile("wspace_4_2.root","RECREATE");
    wspace->Write();
    fout->Close();

    return;

}
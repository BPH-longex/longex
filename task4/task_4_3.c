#include "common.h"

using namespace RooFit;

namespace signal_pdf {

    void fill_dataset(RooWorkspace* wspace, uint cate = 0, double bdt_min = 0.8){

        RooRealVar m("m","",4.9,5.9);
        RooDataSet *rds = new RooDataSet("rds","",RooArgSet(m));

        TFile *fin = new TFile("/eos/user/c/cmsdas/2024/long-ex-bph/bsmmMc.root");
        TTree *tin = (TTree*)fin->Get("bsmmMc");

        unsigned int cate_t;
        float m_t,bdt_t;
        tin->SetBranchAddress("cate",&cate_t);
        tin->SetBranchAddress("m",&m_t);
        tin->SetBranchAddress("bdt",&bdt_t);

        int evtcount[2] = {0,0};
        for(int evt=0; evt<tin->GetEntries(); evt++) {
            tin->GetEntry(evt);
            if (cate_t!=cate) continue;
            evtcount[0]++;
            if (bdt_t<=bdt_min) continue;
            evtcount[1]++;
            m.setVal(m_t);
            rds->add(RooArgSet(m));
        }

        double eff = (double)evtcount[1]/(double)evtcount[0] * effyield[cate][Eff_bsmmMc];
        double eff_error = effyield[cate][dEff_bsmmMc]/effyield[cate][Eff_bsmmMc]*eff; // ignore MC statistical error

        cout << "Category: " << cate << endl;
        cout << "BDT min: " << bdt_min << endl;
        cout << "Selection efficiency: " << eff << " +- " << eff_error << endl;

        wspace->import(m);
        wspace->import(*rds, Rename(Form("rds_signalMc_%d", cate)));
        wspace->import(RooRealVar(Form("Eff_bs_%d",cate),"",eff));

        delete fin;
    }

    void fit(RooWorkspace* wspace, uint cate = 0, double bdt_min = 0.8){

        RooRealVar m = *(RooRealVar*)wspace->var("m");
        RooDataSet *rds = (RooDataSet*)wspace->data(Form("rds_signalMc_%d", cate));

        RooRealVar bs_mean1("bs_mean1","",5.37,5.2,5.5);
        RooRealVar bs_mean2("bs_mean2","",5.37,5.2,5.5);
        RooRealVar bs_sigma1("bs_sigma1","",0.030,0.005,0.060);
        RooRealVar bs_sigma2("bs_sigma2","",0.080,0.040,0.200);
        RooRealVar bs_cbalpha("bs_cbalpha","",1.,0.,4.);
        RooRealVar bs_cbn("bs_cbn","",1.,0.,4.);
        RooRealVar bs_frac("bs_frac","",0.7,0.,1.);
        RooGaussian bs_gaus("bs_gaus", "", m, bs_mean1, bs_sigma1);
        RooCBShape bs_cbline("bs_cbline", "", m, bs_mean2, bs_sigma2, bs_cbalpha, bs_cbn);
        RooAddPdf pdf("pdf", "", RooArgList(bs_gaus, bs_cbline), RooArgList(bs_frac));

        // pdf.fitTo(*rds);

        wspace->import(pdf, RenameAllNodes(Form("signalMc_%d", cate)), RenameAllVariablesExcept(Form("signalMc_%d", cate), "m"));
    }

    void draw(RooWorkspace* wspace, uint cate = 0, double bdt_min = 0){

        RooDataSet *rds = (RooDataSet*)wspace->data(Form("rds_signalMc_%d", cate));
        RooRealVar m = *(RooRealVar*)wspace->var("m");
        RooAddPdf pdf = *(RooAddPdf*)wspace->pdf(Form("pdf_signalMc_%d", cate));

        RooPlot *frame = m.frame(Title(" "),Bins(100));
        rds->plotOn(frame, Name("t_rds"));
        pdf.plotOn(frame, Name("t_pdf"), LineWidth(3));

        TCanvas* canvas = new TCanvas("canvas", "", 600, 600);
        canvas->SetMargin(0.15,0.06,0.13,0.07);
        canvas->Divide(1,2);
        canvas->cd(1);
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
        RooHist *hpull = frame->residHist();
        RooPlot *frame2 = m.frame(Title(" "),Bins(100),Range(5.2,5.5));
        frame2->addPlotable(hpull,"P");
        frame2->Draw();

        string pdf_name="task_4_3_"+ to_string(cate) +".pdf";
        string png_name="task_4_3_"+ to_string(cate) +".png";
        canvas->Print(pdf_name.c_str());
        canvas->Print(png_name.c_str());
    }

}

void task_4_3(unsigned int cate=0, double bdt_min=0.8) {

    RooWorkspace* wspace = new RooWorkspace("wspace","wspace");

    gROOT->SetBatch(kTRUE);

    for(int i = 0; i < 2; i++){
        signal_pdf::fill_dataset(wspace, i, bdt_min);
        signal_pdf::fit(wspace, i, bdt_min);
        signal_pdf::draw(wspace, i, bdt_min);
    }


    // save workspace, including dataset, 
    string wp_name = "wspace_signal_" + to_string(cate) + ".root";
    TFile *fout = new TFile(wp_name.c_str(), "RECREATE");
    wspace->Write();
    fout->Close();
}


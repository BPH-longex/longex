#include "task_4_1.C"
#include "task_4_2.C"
#include "task_4_3.c"
#include "task_4_4.C"

#define N_CATE 8

const double bdt_cut[] = {
   0.993425,
   0.991260,
   0.995055,
   0.993425,
   0.997204,
   0.997898,
   0.997204,
   0.996281
};

using namespace RooFit;

void build_pdf_peak(RooWorkspace* wspace, uint cate = 0, double bdt_min = 0.8){
    peaking_bkg_pdf::fill_dataset(wspace, cate, bdt_min);
    peaking_bkg_pdf::fit(wspace, cate, bdt_min);
}

void build_pdf_semi(RooWorkspace* wspace, uint cate = 0, double bdt_min = 0.8){
    semilep_bkg_pdf::fill_dataset(wspace, cate, bdt_min);
    semilep_bkg_pdf::fit(wspace, cate, bdt_min);
}

void build_pdf_signal(RooWorkspace* wspace, uint cate = 0, double bdt_min = 0.8){
    signal_pdf::fill_dataset(wspace, cate, bdt_min);
    signal_pdf::fit(wspace, cate, bdt_min);
}

void build_pdf_norm(RooWorkspace* wspace, uint cate = 0){
    norm_pdf::fill_dataset_data(wspace, cate);
    norm_pdf::fill_dataset_mc(wspace, cate);
    norm_pdf::fit(wspace, cate);
}

void build_pdf_comb(RooWorkspace *wspace, int cate, double bdt_min)
{
    RooRealVar m("m","",4.9,5.9);

    RooRealVar comb_B1(Form("comb_B1_%d",cate), "", 0.5, 0. , 1);
    RooFormulaVar comb_B2(Form("comb_B2_%d",cate), "", "1.-@0", RooArgList(comb_B1));
    RooBernstein pdf(Form("pdf_comb_%d",cate), "", m, RooArgList(comb_B1, comb_B2));
    
    TFile *fin = new TFile("/eos/user/c/cmsdas/2024/long-ex-bph/bmmData-blind.root");
    TTree *tin = (TTree*)fin->Get("bmmData");
    
    double n_comb_guess = (double)tin->GetEntries(Form("cate==%d&&bdt>%g&&m>5.45", cate, bdt_min));
    n_comb_guess *= 1.0/0.45; // scale to full mass region
    delete fin;
 
    wspace->import(pdf);
    wspace->import(RooRealVar(Form("n_comb_%d",cate),"",n_comb_guess,0.,n_comb_guess*10.));
}

//void task_4_5(unsigned int cate=0, double bdt_min = 0.8){
void task_4_5(){

    RooWorkspace* wspace = new RooWorkspace("wspace","wspace");

    // set batch mode
    gROOT->SetBatch(kTRUE);
    
    for(int cate=0; cate<N_CATE; cate++) {
        build_pdf_peak(wspace, cate, bdt_cut[cate]);
        build_pdf_semi(wspace, cate, bdt_cut[cate]);
        build_pdf_signal(wspace, cate, bdt_cut[cate]);
        build_pdf_norm(wspace, cate);
        build_pdf_comb(wspace, cate, bdt_cut[cate]);
    }

    // save workspace
    wspace->writeToFile("task_4_5.root");
}

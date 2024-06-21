#include "task_4_1.C"
#include "task_4_2.C"
#include "task_4_3.c"
#include "task_4_4.C"

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

void build_pdf_norm(RooWorkspace* wspace, uint cate = 0, double bdt_min = 0.8){
    norm_pdf::fill_dataset_data(wspace, cate, bdt_min);
    norm_pdf::fill_dataset_mc(wspace, cate, bdt_min);
    norm_pdf::fit(wspace, cate, bdt_min);
}

void task_4_5(unsigned int cate=0, double bdt_min = 0.8){

    RooWorkspace* wspace = new RooWorkspace("wspace","wspace");

    // set batch mode
    gROOT->SetBatch(kTRUE);

    // build_pdf_peak(wspace, cate, bdt_min);
    build_pdf_semi(wspace, cate, bdt_min);
    build_pdf_signal(wspace, cate, bdt_min);
    build_pdf_norm(wspace, cate, bdt_min);

    // save workspace
    wspace->writeToFile("task_4_5.root");

}
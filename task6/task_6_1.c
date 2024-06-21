#include "common.h"

void task_6_1(string filename, float cut=0, bool profile=false) {

    std::vector<float> bdt_cuts = {0.993425,
                                   0.991260,
                                   0.995055,
                                   0.993425,
                                   0.997204,
                                   0.997898,
                                   0.997204,
                                   0.996281};

    TFile *fin_wspace = new TFile(filename.c_str());
    RooWorkspace *wspace = (RooWorkspace *)fin_wspace->Get("wspace");
    RooArgSet vars = wspace->allVars();
    for (auto var : vars){
        auto name = TString(var->GetName());
        // && !name.Contains("comb_B1")
        if (name!="m" && !name.Contains("n_comb") && !name.Contains("comb_B1")){
            wspace->var(name)->setConstant(true);
            wspace->var(name)->Print("v");
        }
    }

    RooRealVar *m = wspace->var("m");
    RooArgSet nset(*m);
    RooCategory cate("cate", "");
    for (int idx = 0; idx < N_Categories; idx++)
        cate.defineType(Form("c%d", idx), idx);

    RooSimultaneous *model = new RooSimultaneous("model", "", cate);

    // Main POI: Bs->mumu branching fraction
    RooRealVar *BF_bs = new RooRealVar("BF_bs", "", 3.57E-9, 0., 2E-8);

    // BF(B+ -> J/psi K+) = (1.010 +- 0.028) E-3 (PDG)
    // BF(J/psi -> mu+mu-) = (5.961 +- 0.033) E-2 (PDG)

    //!NEW ones!!!!
    //Br(jpsiK+) (pdg2023) = (1.020+/-0.019)10^{-3}
    //Br(jpsi->mumu)= (5.961+/-0.033)10^{-2}
    RooRealVar *BF_bu_mean = new RooRealVar("BF_bu_mean", "", 1.020E-3 * 5.961E-2);
    RooRealVar *BF_bu = new RooRealVar("BF_bu", "", 1.020E-3*5.961E-2);//, 0.5E-8, 2*(1.020E-3*5.961E-2));
    // fs/fu = 0.252 +- 0.012 (PDG) +- 0.015 (energy/pt dependence)
    RooRealVar *fs_over_fu = new RooRealVar("fs_over_fu", "", 0.252);

    //!Nuisance
    RooRealVar* Bf_bu_sigma= new RooRealVar("Bf_bu_sigma", "", 1.216044E-06);

    RooGaussian *constr_BF_bu = new RooGaussian("constr_BF_bu", "",  *BF_bu, *BF_bu_mean, *Bf_bu_sigma);


    RooFormulaVar *N_bs[N_Categories];
    RooAddPdf *pdf_sum[N_Categories];
    for (int idx = 0; idx < N_Categories; idx++) {

        RooRealVar *Eff_bs = wspace->var(Form("Eff_bs_%d", idx));
        RooRealVar *Eff_bu = wspace->var(Form("Eff_norm_%d", idx));
        RooRealVar *N_bu = wspace->var(Form("n_sig_norm_%d", idx));
        Eff_bu->setVal(effyield[idx][Eff_bupsikMc]);
        N_bu->setVal(effyield[idx][N_bupsikData]);

        N_bs[idx] = new RooFormulaVar(Form("N_bs_%d", idx), "", "@0*@1*@2*@3/@4/@5",
                                      RooArgList(*BF_bs, *N_bu, *fs_over_fu, *Eff_bs, *Eff_bu, *BF_bu));
        RooRealVar *N_peak = wspace->var(Form("N_peak_%d", idx));
        RooRealVar *N_semi = wspace->var(Form("N_semilep_%d", idx));
        RooRealVar *N_comb = wspace->var(Form("n_comb_%d", idx));

        RooArgList pdf_list;
        pdf_list.add(*wspace->pdf(Form("pdf_signalMc_%d", idx)));
        pdf_list.add(*wspace->pdf(Form("pdf_peakingBkg_%d", idx)));
        pdf_list.add(*wspace->pdf(Form("pdf_semilepBkg_%d", idx)));
        pdf_list.add(*wspace->pdf(Form("pdf_comb_%d", idx)));

        RooArgList N_list;
        N_list.add(*N_bs[idx]);
        N_list.add(*N_peak);
        N_list.add(*N_semi);
        N_list.add(*N_comb);

        pdf_sum[idx] = new RooAddPdf(Form("pdf_sum_%d", idx), "", pdf_list, N_list);
        model->addPdf(*pdf_sum[idx], Form("c%d", idx));
    }

    RooDataSet *rds_data = new RooDataSet("rds_data", "", RooArgSet(*m, cate));

    //TFile *fin_data = new TFile("/eos/user/c/cmsdas/2024/long-ex-bph/bmmSoup10.root");

    //TTree *tin = (TTree *)fin_data->Get("bmmSoup10_100");

    TFile *fin_data = new TFile("/eos/user/c/cmsdas/2024/long-ex-bph/bmmData.root");
    TTree *tin = (TTree *)fin_data->Get("bmmData");

    uint cate_t;
    float m_t, bdt_t;
    tin->SetBranchAddress("cate", &cate_t);
    tin->SetBranchAddress("m", &m_t);
    tin->SetBranchAddress("bdt", &bdt_t);
    double bdt_min;
    for (int evt = 0; evt < tin->GetEntries(); evt++) {
        tin->GetEntry(evt);
        if (cut==0){
            bdt_min = bdt_cuts[cate_t];
        }else{
            bdt_min = cut;
        }
        if (bdt_t <= bdt_min)
            continue;
        cate.setIndex(cate_t);
        m->setVal(m_t);
        rds_data->add(RooArgSet(*m, cate));
    }

    RooFitResult* result=model->fitTo(*rds_data,Extended(true),Minos(RooArgSet(*BF_bs)),Save(true)
    //,ExternalConstraints(RooArgSet(*constr_BF_bu))
    );





    TCanvas *canvas = new TCanvas("canvas", "", 1200, 600);
    canvas->Divide(4, 2);

    RooPlot *frame[N_Categories];
    for (int idx = 0; idx < N_Categories; idx++) {
        TVirtualPad *pad = canvas->cd(idx + 1);
        pad->SetMargin(0.15, 0.06, 0.13, 0.07);

        frame[idx] = m->frame(Title(" "), Bins(25));
        rds_data->plotOn(frame[idx], Cut(Form("cate==%d", idx)), MarkerSize(0.8));

        double norm = pdf_sum[idx]->expectedEvents(&nset);
        pdf_sum[idx]->plotOn(frame[idx], Normalization(norm, RooAbsReal::NumEvent), LineWidth(3));
        pdf_sum[idx]->plotOn(frame[idx], Normalization(norm, RooAbsReal::NumEvent), Components(RooArgSet(*wspace->pdf(Form("pdf_signalMc_%d", idx)))), DrawOption("F"), FillColor(kRed), FillStyle(3365));
        pdf_sum[idx]->plotOn(frame[idx], Normalization(norm, RooAbsReal::NumEvent), Components(RooArgSet(*wspace->pdf(Form("pdf_peakingBkg_%d", idx)))), DrawOption("F"), FillColor(kViolet - 4), FillStyle(3344));
        pdf_sum[idx]->plotOn(frame[idx], Normalization(norm, RooAbsReal::NumEvent), Components(RooArgSet(*wspace->pdf(Form("pdf_semilepBkg_%d", idx)))), DrawOption("L"), LineColor(kGreen - 3), LineStyle(2));
        frame[idx]->GetYaxis()->SetTitleOffset(1.50);
        frame[idx]->GetYaxis()->SetTitle("Entries / 0.04 GeV");
        frame[idx]->GetXaxis()->SetTitleOffset(1.15);
        frame[idx]->GetXaxis()->SetLabelOffset(0.01);
        frame[idx]->GetXaxis()->SetTitle("M(#mu#mu) [GeV]");
        frame[idx]->GetXaxis()->SetTitleSize(0.043);
        frame[idx]->GetYaxis()->SetTitleSize(0.043);
        frame[idx]->Draw();
    }
    result->Print("v");
    canvas->Print("task_6_1.pdf");
    canvas->Print("task_6_1.png");


    if (profile){
        std::vector<float> centers;
        std::vector<float> NLL;
        RooAbsReal* nll = model->createNLL(*rds_data,Extended(true));
        TH1D *frameNLL = new TH1D("frameNLL","",35, 1.4E-9, 4.6E-9);
        for(int bin=1; bin<=frameNLL->GetNbinsX(); bin++) {
            BF_bs->setVal(frameNLL->GetBinCenter(bin));
            BF_bs->setConstant(true);
            model->fitTo(*rds_data,Extended(true));
            frameNLL->SetBinContent(bin,(nll->getVal()-result->minNll())*2.);
            centers.push_back(frameNLL->GetBinCenter(bin));
            NLL.push_back((nll->getVal()-result->minNll())*2.);
        }

        for (int i=0; i<centers.size(); i++){
            cout << "BF_bs: " << centers[i] << " NLL: " << NLL[i] << endl;
        }

        TCanvas* canvasNLL = new TCanvas("canvasNLL", "", 600, 600);
        frameNLL->Draw("c");
        canvasNLL->Print("task_6_1_nll.pdf");
    }

    RooWorkspace *wspace_new = (RooWorkspace *)fin_wspace->Get("wspace");
    wspace_new->import(*model);
    wspace_new->writeToFile("outfile.root");

}

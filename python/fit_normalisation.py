import ROOT
from ROOT import RooFit, RooRealVar, RooGaussian, RooAddPdf, RooExponential, RooGenericPdf, RooArgList, RooAddPdf, RooDataSet, TFile, TTree

import argparse
import mplhep as hep
from pdfs import pdf_mc_sig, pdf_norm_data
from read_data import read_data, read_mc
def parse_arguments():
    parser = argparse.ArgumentParser(description="PyROOT fitting script with user-defined arguments.")
    
    parser.add_argument('--data_in', type=str, required=True, help="Input data file (ROOT file).")
    parser.add_argument('--tree_name', type=str, required=True, help="Name of the tree in the ROOT file.")
    #parser.add_argument('--cuts', type=str, required=True, help="Cuts to apply to the data (as a string).")
    parser.add_argument('--mc', action='store_true', help="Flag to indicate if the input data is Monte Carlo simulation.")
    parser.add_argument('--output_path', type=str, required=True, help="Output path to save the results.")
    
    return parser.parse_args()


def fit(data_in, tree_name, mc, output_path):
    # Define the variables
    m = RooRealVar("m", "m", 5.0, 5.8)

    if not mc: 
       wgt = RooRealVar("wgt", "wgt", 1.0, 0.0, 1000.0)
       variables = [m, wgt]
       model = pdf_norm_data(m)
       rds = read_data(data_in, tree_name, variables)

    if mc:
       cate = RooRealVar("cate","",-1,20)
       variables = [m, cate]
       model = pdf_norm_sig(m)
       rds = read_mc(data_in, tree_name, variables)
     
    # Fit the model to the data
    fit_result = model.fitTo(rds, RooFit.Extended(True), RooFit.SumW2Error(True))

    # Plot the results
    frame = m.frame(RooFit.Title("Fit result"))
    rds.plotOn(frame)
    model.plotOn(frame)
    if mc:
        model.plotOn(frame, RooFit.Components("sig_g1"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(kRed))
        model.plotOn(frame, RooFit.Components("sig_g2"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(kGreen))
        model.plotOn(frame, RooFit.Components("pdf_sig"), RooFit.LineColor(ROOT.kBlue))


    if not mc:
        model.plotOn(frame, RooFit.Components("pdf_comb"), RooFit.LineStyle(ROOT.kDashed))
        model.plotOn(frame, RooFit.Components("pdf_jpsix"), RooFit.LineStyle(ROOT.kDotted))
        model.plotOn(frame, RooFit.Components("pdf_sig"), RooFit.LineColor(ROOT.kRed))

    # Pull plot
    pulls = frame.pullHist()
    pull_frame = m.frame(ROOT.RooFit.Title("Pull Plot"))
    pull_frame.addPlotable(pulls, "P")

    # Save the fit result and pull plot using mplhep CMS style
    mplhep.style.use("CMS")

    # Fit result plot
    c = ROOT.TCanvas("c", "c", 800, 600)
    frame.Draw()
    c.SaveAs(f"{output_path}/fit_result_root.png")

    # Convert ROOT frame to matplotlib
    fig, ax = plt.subplots()
    mplhep.histplot(
        np.array([frame.getObject(i).getY() for i in range(frame.numItems()) if "h_data" in frame.getObject(i).GetName()]),
        bins=np.array([frame.getObject(i).getX() for i in range(frame.numItems()) if "h_data" in frame.getObject(i).GetName()]),
        histtype="step",
        label="Data",
        ax=ax
    )
    ax.errorbar(
        np.array([frame.getObject(i).getX() for i in range(frame.numItems()) if "h_data" in frame.getObject(i).GetName()]),
        np.array([frame.getObject(i).getY() for i in range(frame.numItems()) if "h_data" in frame.getObject(i).GetName()]),
        yerr=np.array([frame.getObject(i).getErrorYhigh() for i in range(frame.numItems()) if "h_data" in frame.getObject(i).GetName()]),
        fmt='o', color='black'
    )

    # Apply the CMS style
    mplhep.cms.label("Preliminary", ax=ax)
    plt.legend()
    plt.savefig(f"{output_path}/fit_result.png")

    # Pull plot
    fig, ax = plt.subplots()
    mplhep.histplot(
        np.array([pull_frame.getObject(i).getY() for i in range(pull_frame.numItems())]),
        bins=np.array([pull_frame.getObject(i).getX() for i in range(pull_frame.numItems())]),
        histtype="step",
        label="Pull",
        ax=ax
    )

    # Apply the CMS style
    mplhep.cms.label("Preliminary", ax=ax)
    plt.legend()
    plt.savefig(f"{output_path}/pull_plot.png")



if __name__ == "__main__":
    args = parse_arguments()
    
    fit(**vars(args))

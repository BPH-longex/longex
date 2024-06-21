import ROOT
from ROOT import RooFit, RooRealVar, RooGaussian, RooAddPdf, RooExponential, RooGenericPdf, RooArgList, RooAddPdf, RooDataSet, TFile, TTree, RooWorkspace
import argparse
import mplhep as hep
from pdfs import pdf_mc_sig, pdf_norm_data, pdf_mc_sig_me
from read_data import read_data, read_mc
import matplotlib.pyplot as plt
import numpy as np
import json
from ROOT import RooFit, RooGaussian, RooRealVar, RooArgList, RooAddPdf, RooExponential, RooGenericPdf
import yaml

with open('configuration.yaml', 'r') as file:
    categories = yaml.safe_load(file)["categories"]

def extract_data_points(frame):
    x_data, y_data, y_errors_l, y_errors_h = [], [], [], []
    for i in range(int(frame.numItems())):
        if frame.getObject(i).IsA().InheritsFrom("RooHist"):
            hist = frame.getObject(i)
            for j in range(hist.GetN()):
                x, y = hist.GetX()[j], hist.GetY()[j]
                ey_high = hist.GetErrorYhigh(j)
                ey_low = hist.GetErrorYlow(j)
                x_data.append(x)
                y_data.append(y)
                y_errors_l.append(ey_low)
                y_errors_h.append(ey_high)
    return np.array(x_data), np.array(y_data), np.array(y_errors_l), np.array(y_errors_h)


def extract_pdf_curve(frame, name):
    x_curve, y_curve = [], []
    for i in range(int(frame.numItems())):
        if frame.getObject(i).IsA().InheritsFrom("RooCurve") and frame.getObject(i).GetName() == name:
            print(frame.getObject(i).GetName())
            curve = frame.getObject(i)
            for j in range(curve.GetN()):
                x, y = curve.GetX()[j], curve.GetY()[j]
                x_curve.append(x)
                y_curve.append(y)
    return np.array(x_curve), np.array(y_curve)

def calculate_pulls(y_data, y_fit, y_errors_l, y_errors_h):
    y_errors = [ (x+y)/2 for x, y in zip( y_errors_l, y_errors_h)]
    pulls = []#(y_data - y_fit) / y_errors
    pulls_err_l, pulls_err_h = [], []
    for i in range(len(y_data)):
        resid = (y_data[i] - y_fit[i])
        norm = y_errors_l[i] if resid > 0  else y_errors_h[i]
        pull= resid/norm
        pulls.append(pull)
        pulls_err_l.append(y_errors_l[i]/norm)
        pulls_err_h.append(y_errors_h[i]/norm)
    
    return pulls, pulls_err_l, pulls_err_h

def parse_arguments():
    parser = argparse.ArgumentParser(description="PyROOT fitting script with user-defined arguments.")
    
    parser.add_argument('--data_in', type=str, required=True, help="Input data file (ROOT file).")
    parser.add_argument('--tree_name', type=str, required=True, help="Name of the tree in the ROOT file.")
    parser.add_argument('--parameters', type=str, default=None, help="Parameters To fix")
    parser.add_argument('--cat', type=str, required=True, help="Category.")
    parser.add_argument('--mc', action='store_true', help="Flag to indicate if the input data is Monte Carlo simulation.")
    parser.add_argument('--output_path', type=str, required=True, help="Output path to save the results.")
    
    return parser.parse_args()


def fit(data_in, tree_name, parameters, cat, mc, output_path):
    ws = RooWorkspace(f"ws_{cat}")
    # Define the variables
    m = RooRealVar("m", "m", 5.0, 5.8)
    me = RooRealVar("me", "me", 0., 1.)
    cate = RooRealVar("cate","cate",-1,20)
    variables = [m, me, cate]
 
    if not mc: 
       wgt = RooRealVar("wgt", "wgt", 1.0, 0.0, 1000.0)
       variables.append(wgt)
       read_data(data_in, tree_name, variables, ws, categories[cat]["val"])
       pdf_norm_data(m, ws, shape="dcb")
    if mc:
       read_mc(data_in, tree_name, variables, ws, categories[cat]["val"]) 
       pdf_mc_sig(m, ws, shape = "dcb")
   
   
    ws.Print("v")
    model = ws.pdf("model")
    rds = ws.data("rds_data")
    # Fit the model to the data
    if parameters != None:
        print("fixing parameters")
        with open(parameters, 'r') as file:
            parameters = json.load(file)
        for p in parameters.keys():
            if p in ["alpha1", "alpha2", "n1", "n2", "sig_mean", "sig_sigma"]:
                ws.var(p).setVal(parameters[p]["value"])
                ws.var(p).setConstant(True)

    if not mc:
        fit_result = model.fitTo(rds, RooFit.Extended(True), RooFit.SumW2Error(True), RooFit.Save(True))
    if mc: 
        fit_result = model.fitTo(rds, RooFit.Save(True), ConditionalObservables={me})

    fit_result.Print("v")
    fit_result.correlationMatrix().Print("v")

    # Plot the results
    frame = m.frame(RooFit.Title("Fit result"))
    rds.plotOn(frame)
    model.plotOn(frame)
    if mc:
        model.plotOn(frame, RooFit.Components("sig_g1"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kRed))
        model.plotOn(frame, RooFit.Components("sig_g2"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen))
        model.plotOn(frame, RooFit.Components("model"), RooFit.LineColor(ROOT.kBlue))


    if not mc:
        model.plotOn(frame, RooFit.Components("pdf_comb"), RooFit.LineStyle(ROOT.kDashed))
        model.plotOn(frame, RooFit.Components("pdf_jpsix"), RooFit.LineStyle(ROOT.kDotted))
        model.plotOn(frame, RooFit.Components("pdf_sig"), RooFit.LineColor(ROOT.kRed))
        model.plotOn(frame, RooFit.Components("model"), RooFit.LineColor(ROOT.kBlue))
    # Pull plot
    pulls = frame.pullHist()
    pull_frame = m.frame()
    pull_frame.addPlotable(pulls, "P")

    # Save the fit result and pull plot using hep CMS style
    hep.style.use("CMS")

    # Fit result plot
    c = ROOT.TCanvas("c", "c", 800, 600)
    frame.Draw()
    c.SaveAs(f"{output_path}/fit_result_root.png")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), gridspec_kw={'height_ratios': [3, 1]})
    hep.cms.label('$B^{+} \\rightarrow J/\psi K^{+}$ '+categories[cat]["eta"], year = categories[cat]["year"], data=mc, ax = ax1, fontsize=12.)



    x_data, y_data, y_errors_l, y_errors_h  = extract_data_points(frame) 
    x_pdf, y_pdf = extract_pdf_curve(frame, "model_Norm[m]_Comp[model]")
    ax1.set_xlim(5.0, 5.8)

    ax1.errorbar(x_data, y_data, yerr=[y_errors_l, y_errors_h], fmt='o', label='Data points', zorder=1)


    if not mc:

       x1_pdf, y1_pdf = extract_pdf_curve(frame, "model_Norm[m]_Comp[model]")
       x2_pdf, y2_pdf = extract_pdf_curve(frame, "model_Norm[m]_Comp[model]")
 
    if mc:

       x1_pdf, y1_pdf = extract_pdf_curve(frame, "model_Norm[m]_Comp[sig_g1]")
       x2_pdf, y2_pdf = extract_pdf_curve(frame, "model_Norm[m]_Comp[sig_g2]")
       ax1.plot(x1_pdf, y1_pdf, color = 'red', linestyle='-', linewidth = 3, label = "Gauss 1", zorder=2)
       ax1.plot(x2_pdf, y2_pdf, color = 'green', linestyle='-', linewidth = 3, label = "Gauss 2", zorder=2)

    ax1.plot(x_pdf, y_pdf, color = 'black', linewidth = 3, label = "", zorder=3)
    params = fit_result.floatParsFinal()
    fit_params = {param.GetName(): (param.getVal(), param.getError()) for param in params}

    # Prepare the annotation text
    param_text = '\n'.join([f'{name}: {val:.3f} ± {err:.3f}' for name, (val, err) in fit_params.items()])
    ax1.annotate(param_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12,
             horizontalalignment='left', verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))

    ax1.set_ylabel("N")
    # Apply the CMS style
    y_errors = (y_errors_l + y_errors_h) / 2
    y_fit = np.interp(x_data, x_pdf, y_pdf)
    pulls, pulls_err_l, pulls_err_h = calculate_pulls(y_data, y_fit, y_errors_l, y_errors_h)
    pulls_err = [ (x+y)/2 for x, y in zip(pulls_err_l,pulls_err_h)]

    ax2.errorbar(x_data, pulls, yerr=[pulls_err_l, pulls_err_h], fmt='o', label='Pulls', capsize=2)
    ax2.set_xlim(5.0, 5.8)
    ax2.axhline(0, color='red', linestyle='--')
    ax2.set_xlabel('m [GeV]')
    ax2.set_ylabel('Pull')
    ax2.legend()

    plt.tight_layout()
    plt.legend()
    #plt.show()
    plt.savefig(f"{output_path}/fit_result.png")
    plt.clf()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), gridspec_kw={'height_ratios': [3, 1]})
    hep.cms.label('$B^{+} \\rightarrow J/\psi K^{+}$ '+categories[cat]["eta"], year = categories[cat]["year"], data=mc, ax = ax1, fontsize=12.)




    ax1.errorbar(x_data, y_data, yerr=[y_errors_l, y_errors_h], fmt='o', label='Data points', zorder=1)


    if not mc:

       x1_pdf, y1_pdf = extract_pdf_curve(frame, "model_Norm[m]_Comp[model]")
       x2_pdf, y2_pdf = extract_pdf_curve(frame, "model_Norm[m]_Comp[model]")
 
    if mc:

       x1_pdf, y1_pdf = extract_pdf_curve(frame, "model_Norm[m]_Comp[sig_g1]")
       x2_pdf, y2_pdf = extract_pdf_curve(frame, "model_Norm[m]_Comp[sig_g2]")
       ax1.plot(x1_pdf, y1_pdf, color = 'red', linestyle='-', linewidth = 3, label = "Gauss 1", zorder=2)
       ax1.plot(x2_pdf, y2_pdf, color = 'green', linestyle='-', linewidth = 3, label = "Gauss 2", zorder=2)
    ax1.set_xlim(5.0, 5.8)
    ax1.set_yscale('log')
    ax1.plot(x_pdf, y_pdf, color = 'black', linewidth = 3, label = "", zorder=3)
    params = fit_result.floatParsFinal()
    fit_params = {param.GetName(): (param.getVal(), param.getError()) for param in params}

    # Prepare the annotation text
    param_text = '\n'.join([f'{name}: {val:.3f} ± {err:.3f}' for name, (val, err) in fit_params.items()])
    ax1.annotate(param_text, xy=(0.05, 0.95), xycoords='axes fraction', fontsize=12,
             horizontalalignment='left', verticalalignment='top', bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))

    ax1.set_ylabel("N")
    # Apply the CMS style
    y_errors = (y_errors_l + y_errors_h) / 2
    y_fit = np.interp(x_data, x_pdf, y_pdf)
    pulls, pulls_err_l, pulls_err_h = calculate_pulls(y_data, y_fit, y_errors_l, y_errors_h)
    pulls_err = [ (x+y)/2 for x, y in zip(pulls_err_l,pulls_err_h)]
    ax2.set_xlim(5.0, 5.8)
    ax2.errorbar(x_data, pulls, yerr=[pulls_err_l, pulls_err_h], fmt='o', label='Pulls', capsize=2)
    ax2.axhline(0, color='red', linestyle='--')
    ax2.set_xlabel('m [GeV]')
    ax2.set_ylabel('Pull')
    ax2.legend()

    plt.tight_layout()
    plt.legend()
    #plt.show()
    plt.savefig(f"{output_path}/fit_result_log.png")
 

    fit_results = {}
    for i in range(params.getSize()):
        param = params.at(i)
        fit_results[param.GetName()] = {
            "value": param.getVal(),
            "error": param.getError()
        }

    with open(output_path+"/parameters.json", "w") as json_file:
        json.dump(fit_results, json_file, indent=4)

if __name__ == "__main__":
    args = parse_arguments()
    
    fit(**vars(args))

from ROOT import TFile, TTree
from ROOT import RooFit, RooDataSet, RooArgSet, RooRealVar, RooArgList   
import ROOT


def read_data(data_in, tree_name, variables, ws, cat):
    # Extract the necessary variables from the list
    m = None
    wgt = None
    cate = None
    for var in variables:
        if var.GetName() == "m":
            m = var
        elif var.GetName() == "wgt":
            wgt = var
        elif var.GetName() == "cate":
            cate = var
    if m is None or wgt is None or cate is None:
        raise ValueError("The list of RooRealVar must include 'm' and 'wgt' and 'cate' variables.")
 


    # Load the data from the ROOT file
    fin = TFile(data_in)
    tin = fin.Get(tree_name)

    rds_data = RooDataSet("rds_data", "", RooArgSet(m, wgt, cate), RooFit.Import(tin), RooFit.WeightVar(wgt))
    cate_t = 0
    rds_data.reduce(f"m >= 5.0 && m <= 5.8 && cate=={cat}")

    fin.Close()
    getattr(ws, 'import')(rds_data)


def read_mc(data_in, tree_name, variables, ws, cat):
    m = None
    cate = None
    me  = None
    for var in variables:
        if var.GetName() == "m":
            m = var
        elif var.GetName() == "cate":
            cate = var
        elif var.GetName() == "me":
            me = var
    if m is None or cate is None or me is None:
        raise ValueError("The list of RooRealVar must include 'm' and 'cate' and 'me' variables.")
    
 
    cate_t = 0
   
    fin = TFile(data_in)
    tin = fin.Get(tree_name)

    rds_mc   = RooDataSet("rds_data", "rds_data", RooArgSet(cate, m, me), RooFit.Import(tin))

    rds_mc = rds_mc.reduce(f"m >= 5.0 && m <= 5.8 && cate=={cat}")
    fin.Close()
    getattr(ws, 'import')(rds_mc)


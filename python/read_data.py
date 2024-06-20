from ROOT import TFile, TTree
from ROOT import RooFit, RooDataSet, RooArgSet, RooRealVar, RooArgList   
import ROOT


def read_data(data_in, tree_name, variables):
    # Extract the necessary variables from the list
    m = None
    wgt = None
    for var in variables:
        if var.GetName() == "m":
            m = var
        elif var.GetName() == "wgt":
            wgt = var

    if m is None or wgt is None:
        raise ValueError("The list of RooRealVar must include 'm' and 'wgt' variables.")
 


    # Load the data from the ROOT file
    fin = TFile(data_in)
    tin = fin.Get(tree_name)

    cate_t = ROOT.std.vector('unsigned int')()
    wgt_t = ROOT.std.vector('float')()
    m_t = ROOT.std.vector('float')()

    tin.SetBranchAddress("cate", cate_t)
    tin.SetBranchAddress("wgt", wgt_t)
    tin.SetBranchAddress("m", m_t)

    rds_data = RooDataSet("rds_data", "", RooArgSet(m, wgt), RooFit.WeightVar(wgt), "cate==0 && m > 5.0 && <= 5.8")

    fin.Close()
    return rds_data


def read_mc(data_in, tree_name, variables):
    m = None
    cate = None
    for var in variables:
        if var.GetName() == "m":
            m = var
        elif var.GetName() == "cate":
            cate = var

    if m is None or cate is None:
        raise ValueError("The list of RooRealVar must include 'm' and 'cate' variables.")
    
 
    cate_t = 0
   
    fin = TFile(data_in)
    tin = fin.Get(tree_name)

    rds_mc   = RooDataSet("rds_mc", "rds_mc", tin,  RooArgSet(cate, m))

    rds_mc = rds_mc.reduce(m, Form("m >= 5.0 && m <= 5.8 && cate==%d", cate_t))
    fin.Close()
    return rds_mc
 

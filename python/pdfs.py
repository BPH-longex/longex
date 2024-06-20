
from ROOT import RooFit, RooGaussian, RooRealVar, RooArgList, RooAddPdf, RooExponential, RooGenericPdf

def pdf_mc_sig(m, shape = "gauss"):
    sig_mean1 = RooRealVar("sig_mean1", "", 5.28, 5.2, 5.4)
    sig_mean2 = RooRealVar("sig_mean2", "", 5.28, 5.2, 5.4)
    sig_sigma1 = RooRealVar("sig_sigma1", "", 0.030, 0.005, 0.060)
    sig_sigma2 = RooRealVar("sig_sigma2", "", 0.080, 0.040, 0.200)
    sig_frac = RooRealVar("sig_frac", "", 0.9, 0.5, 1.0)

    # Define the signal Gaussians
    sig_g1 = RooGaussian("sig_g1", "", m, sig_mean1, sig_sigma1)
    sig_g2 = RooGaussian("sig_g2", "", m, sig_mean2, sig_sigma2)

    # Define the combined signal PDF
    pdf_sig = RooAddPdf("pdf_sig", "", RooArgList(sig_g1, sig_g2), RooArgList(sig_frac))
    return pdf_sig




def pdf_norm_data(m, shape="gauss_err"):

    sig_mean1 = RooRealVar("sig_mean1", "", 5.27971)
    sig_mean2 = RooRealVar("sig_mean2", "", 5.26646)
    sig_sigma1 = RooRealVar("sig_sigma1", "", 0.0156007)
    sig_sigma2 = RooRealVar("sig_sigma2", "", 0.0627436)
    sig_frac = RooRealVar("sig_frac", "", 0.926022)

    sig_g1 = RooGaussian("sig_g1", "", m, sig_mean1, sig_sigma1)
    sig_g2 = RooGaussian("sig_g2", "", m, sig_mean2, sig_sigma2)

    pdf_sig = RooAddPdf("pdf_sig", "", RooArgList(sig_g1, sig_g2), RooArgList(sig_frac))

    # Combinatorial background
    comb_coeff = RooRealVar("comb_coeff", "", -1.2, -10., 10.)
    pdf_comb = RooExponential("pdf_comb", "", m, comb_coeff)

    # J/psi X background
    jpsix_scale = RooRealVar("jpsix_scale", "", 0.02, 0.001, 0.08)
    jpsix_shift = RooRealVar("jpsix_shift", "", 5.13, 5.12, 5.16)
    pdf_jpsix = RooGenericPdf("pdf_jpsix", "", "TMath::Erfc((@0-@1)/@2)", RooArgList(m, jpsix_shift, jpsix_scale))

    # Yields
    n_sig = RooRealVar("n_sig", "", 100000, 0., 1E8)
    n_comb = RooRealVar("n_comb", "", 80000, 0., 1E6)
    n_jpsix = RooRealVar("n_jpsix", "", 20000, 0., 1E5)

    # Combined model
    model = RooAddPdf("model", "", RooArgList(pdf_sig, pdf_comb, pdf_jpsix), RooArgList(n_sig, n_comb, n_jpsix))
    return model


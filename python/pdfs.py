
from ROOT import RooFit, RooGaussian, RooRealVar, RooArgList, RooAddPdf, RooExponential, RooGenericPdf, RooCrystalBall, RooFormulaVar, RooProdPdf
from ROOT import RooWorkspace
def pdf_mc_sig(m, ws, shape = "gauss"):
    if shape == "gauss":
        sig_mean1 = RooRealVar("sig_mean1", "", 5.28, 5.2, 5.4)
        sig_mean2 = RooRealVar("sig_mean2", "", 5.28, 5.2, 5.4)
        sig_sigma1 = RooRealVar("sig_sigma1", "", 0.030, 0.005, 0.060)
        sig_sigma2 = RooRealVar("sig_sigma2", "", 0.080, 0.040, 0.200)
        sig_frac = RooRealVar("sig_frac", "", 0.9, 0.5, 1.0)

        # Define the sigmal Gaussians
        sig_g1 = RooGaussian("sig_g1", "", m, sig_mean1, sig_sigma1)
        sig_g2 = RooGaussian("sig_g2", "", m, sig_mean2, sig_sigma2)

        # Define the combined sigmal PDF
        model = RooAddPdf("model", "", RooArgList(sig_g1, sig_g2), RooArgList(sig_frac))
    if shape == "dcb":
        mean = RooRealVar("sig_mean", "", 5.28, 5.2, 5.4)
        sigma = RooRealVar("sig_sigma", "width", 0.030, 0.005, 0.060)
        alpha1 = RooRealVar("alpha1", "alpha1", 1.5, 0.1, 10)
        n1 = RooRealVar("n1", "n1", 2, 0.1, 10)
        alpha2 = RooRealVar("alpha2", "alpha2", 1.5, 0.1, 10)
        n2 = RooRealVar("n2", "n2", 2, 0.1, 10)
        model = RooCrystalBall("model", "", m, mean, sigma, alpha1, n1, alpha2, n2) 


    if shape == "dcb_gauss":
        mean = RooRealVar("sig_mean", "", 5.28, 5.2, 5.4)
        sigma = RooRealVar("sig_sigma", "width", 0.030, 0.005, 0.060)
        alpha1 = RooRealVar("alpha1", "alpha1", 1.5, 0.1, 10)
        n1 = RooRealVar("n1", "n1", 2, 0.1, 10)
        alpha2 = RooRealVar("alpha2", "alpha2", 1.5, 0.1, 10)
        n2 = RooRealVar("n2", "n2", 2, 0.1, 10)
        model = RooCrystalBall("model", "", m, mean, sigma, alpha1, n1, alpha2, n2) 
        sig_mean1 = RooRealVar("sig_mean1", "", 5.28, 5.2, 5.4)
        sig_mean2 = RooRealVar("sig_mean2", "", 5.28, 5.2, 5.4)
        sig_sigma1 = RooRealVar("sig_sigma1", "", 0.030, 0.005, 0.060)
        sig_sigma2 = RooRealVar("sig_sigma2", "", 0.080, 0.040, 0.200)
        sig_frac = RooRealVar("sig_frac", "", 0.9, 0.5, 1.0)

        # Define the sigmal Gaussians
        sig_g1 = RooGaussian("sig_g1", "", m, sig_mean1, sig_sigma1)
        sig_g2 = RooGaussian("sig_g2", "", m, sig_mean2, sig_sigma2)

        # Define the combined sigmal PDF
        model = RooAddPdf("model", "", RooArgList(sig_g1, sig_g2), RooArgList(sig_frac))
 

    getattr(ws, 'import')(model)
    

def pdf_mc_sig_me(m, me, ws, shape = "gauss"):
    if shape == "gauss":
        sig_mean1 = RooRealVar("sig_mean1", "", 5.28, 5.2, 5.4)
        sig_mean2 = RooRealVar("sig_mean2", "", 5.28, 5.2, 5.4)
        sig_sigma1 = RooRealVar("sig_sigma1", "", 0.030, 0.005, 0.060)
        sig_sigma2 = RooRealVar("sig_sigma2", "", 0.080, 0.040, 0.200)
        sig_frac = RooRealVar("sig_frac", "", 0.9, 0.5, 1.0)

        # Define the sigmal Gaussians
        sig_g1 = RooGaussian("sig_g1", "", m, sig_mean1, sig_sigma1)
        sig_g2 = RooGaussian("sig_g2", "", m, sig_mean2, sig_sigma2)

        # Define the combined sigmal PDF
        model = RooAddPdf("model", "", RooArgList(sig_g1, sig_g2), RooArgList(sig_frac))
    if shape == "dcb":
        mean = RooRealVar("sig_mean", "", 5.28, 5.2, 5.4)
        scale = RooRealVar("scale", "", 0.01, 0., 1.)
        bias = RooRealVar("bias", "", 0., -0.1, 0.1)
        sigma = RooFormulaVar("sig_sigma", "[0]*[1]+[2]", RooArgList(scale, me, bias))
        alpha1 = RooRealVar("alpha1", "alpha1", 1.5, 0.1, 10)
        n1 = RooRealVar("n1", "n1", 2, 0.1, 10)
        alpha2 = RooRealVar("alpha2", "alpha2", 1.5, 0.1, 10)
        n2 = RooRealVar("n2", "n2", 2, 0.1, 10)
        dcb = RooCrystalBall("dcb", "", m, mean, sigma, alpha1, n1, alpha2, n2)
        model = RooProdPdf("model","",dcb);

    getattr(ws, 'import')(model)
 


def pdf_norm_data(m, ws, shape="gauss_err"):

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
    getattr(ws, 'import')(model)
   

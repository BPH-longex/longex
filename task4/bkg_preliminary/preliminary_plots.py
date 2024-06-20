import ROOT
import matplotlib as mpl
import matplotlib.pyplot as plt
import mplhep as hep
import uproot

# set CMS style
hep.style.use("CMS")
plt.style.use(hep.style.CMS)

# Styling options
mpl.rcParams['patch.linewidth'] = 2 #for step hist plots

# set batch mode
ROOT.gROOT.SetBatch(ROOT.kTRUE)

input_path = "/eos/user/c/cmsdas/2024/long-ex-bph/"

# background categories
bkg_categories = {
    "peaking" : {
        "decays" : ["bdmmMc", "bstohhMcBg", "bdtohhMcBg"],
        "decays_latex" : [r"$B^0 \rightarrow \mu^+ \mu^-$", r"$B_s \rightarrow h h\ (h = K, \pi)$", r"$B^0 \rightarrow h h \ (h = K, \pi)$"]
    }, 

    "semileptonic" : {
        "decays" : ["bskmunuMcBg", "bdpimunuMcBg", "bupimumuMcBg"],
        "decays_latex" : [r"$B_s \rightarrow K^+ \mu^- \nu_\mu$", r"$B^0 \rightarrow \pi^+ \mu^- \nu_\mu$", r"$B^+ \rightarrow \pi^+ \mu^+ \mu^-$"]
    }
}

for category in bkg_categories:

    decays = bkg_categories[category]["decays"]
    decays_latex = bkg_categories[category]["decays_latex"]

    fig, axs = plt.subplots(1, 3, figsize=(30, 10))

    for ax, proc, proc_latex in zip(axs, decays, decays_latex):
        fin = uproot.open(input_path + proc + ".root")
        tin = fin[proc]
        data = tin["m"].array(library="np")
        ax.hist(tin["m"].array(library="np"), bins=50, histtype='step', label=f"Entries = {data.size}")
        ax.set_title(proc_latex, pad = 40)
        ax.set_xlabel("m [GeV]")
        ax.set_ylabel("Entries")
        hep.cms.label(loc=0, ax=ax, label = "Private Work", rlabel = "", fontsize = 20)
        ax.legend()

    plt.tight_layout()
    
    outname = f"{category}_preliminary_plots"
    plt.savefig(f"{outname}.png")
    plt.savefig(f"{outname}.pdf")
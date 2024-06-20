# %%
import ROOT
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
import hist
from cycler import cycler


filename = "task2/jpsik.root"
workspace_name = "wspace"
dataset_name = "rds_data"

#! THE TOTAL MUST BE THE FIRST
pdf_names = ["model","pdf_comb", "sig_g1", "sig_g2", "pdf_sig", "pdf_jpsix"]
labels = ["Total", "Combinatorial", "Gauss1", "Gauss2", "GaussSum", "ErrorFunc"]
hist_range = (5, 5.8)


#filename="task2/norm_wspace.root"
#pdf_names=[""]
#labels=[""]

#!MUST BE EQUAL TO THE ONE SETTED IN THE FITTING SCRIPT
bins = 100
x_sampling=400

acab_palette = (
    "#1f77b4",
    "#ff7f0e",
    "#2ca02c",
    "#d62728",
    "#9467bd",
    "#8c564b",
    "#e377c2",
    "#7f7f7f",
    "#bcbd22",
    "#17becf",
)


hep.styles.cms.CMS["patch.linewidth"] = 2
hep.styles.cms.CMS["lines.linewidth"] = 3
hep.styles.cms.CMS["axes.prop_cycle"] = cycler("color", acab_palette)
hep.styles.cms.CMS["legend.frameon"] = True
# hep.styles.cms.CMS["figure.autolayout"] = True
hep.style.CMS["axes.grid"] = True
hep.set_style("CMS")
# Initialize the ROOT framework
ROOT.gROOT.SetBatch(True)


# %%
def pdf(x, frame, data, model):
    data.plotOn(frame)

    components = ROOT.RooFit.Components(ROOT.RooArgSet(model))
    model = workspace.pdf(pdf_names[0])
    model.plotOn(frame, components)

    # Extract data points and fit curve from the frame
    fit_curve = frame.getObject(1)
    y = np.array([fit_curve.Eval(xi) for xi in x])
    return y


# Open the ROOT file containing the workspace
file = ROOT.TFile.Open(filename)

# Load the workspace
workspace = file.Get(workspace_name)

# Close the ROOT file
file.Close()

# Get dataset and model from workspace
data = workspace.data(dataset_name)

# %%
# Create a frame for plotting


np_data = data.to_numpy()

x = np.linspace(*hist_range, x_sampling)
h = hist.Hist(hist.axis.Regular(bins, *hist_range), storage=hist.storage.Weight())
h.fill(np_data["m"], weight=np_data["wgt"])

fig, ax = plt.subplots(2, 1, gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
plt.subplots_adjust(hspace=0)
h.plot(histtype="errorbar", color="black", ax=ax[0], linewidth=2, label="Data")

#! Fight the garbage collector with his own weapons
observables = []
frames = []
models = []
pdfs = []
for pdf_name, label in zip(pdf_names, labels):
    observable = workspace.var("m")
    frame = observable.frame()
    model = workspace.pdf(pdf_name)
    observables.append(observable)
    frames.append(frame)
    models.append(model)
    pdfs.append(pdf(x, frame, data, model))
    ax[0].plot(x, pdfs[-1], label=label)

ax[0].set_xlim(*hist_range)
residuals = h.values() - pdf(h.axes[0].centers, frames[0], data, models[0])
ax[1].errorbar(h.axes[0].centers, residuals, fmt="o", yerr=np.sqrt(h.variances()))
ax[1].axhline(0, color="red", linestyle="--")
ax[1].set_xlim(*hist_range)

hep.cms.text("Private Work", ax=ax[0])
hep.cms.lumitext("13 TeV", ax=ax[0])

ax[0].set_ylabel("Events")
ax[1].set_xlabel("$m_{\mu \mu}$ [GeV]")

ax[0].legend()

fig.savefig(filename.split(".root")[0] + ".pdf")

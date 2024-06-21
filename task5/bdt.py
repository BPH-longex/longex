# %%
from coffea.nanoevents import NanoEventsFactory, NanoAODSchema
from collections import OrderedDict
import matplotlib.pyplot as plt
import numpy as np
import mplhep as hep
from cycler import cycler
import pandas as pd
import awkward as ak
import hist

path = "/eos/user/c/cmsdas/2024/long-ex-bph/"

NanoAODSchema.error_missing_event_ids = False


signal =NanoEventsFactory.from_root(path + "bsmmMc.root", schemaclass=NanoAODSchema).events().compute()

bkg=NanoEventsFactory.from_root(path + "bmmData-blind.root", schemaclass=NanoAODSchema).events().compute()


df_eff=pd.read_csv("effyield.csv")
#%%
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


hep.styles.cms.CMS["patch.linewidth"] = 3
hep.styles.cms.CMS["lines.linewidth"] = 4
hep.styles.cms.CMS["axes.prop_cycle"] = cycler("color", acab_palette)
hep.styles.cms.CMS["legend.frameon"] = True
# hep.styles.cms.CMS["figure.autolayout"] = True
hep.style.CMS["axes.grid"] = True
hep.style.CMS["axes.axisbelow"]=True
hep.style.use("CMS")

bins=35
func = lambda x: np.arctanh(x)


def significance(ns,nb):
    return np.sqrt(2*((ns+nb)*np.log(1+ns/nb)-ns))

best_cut = []
for cat in range(8):
    fig, ax = plt.subplots(2, 1, gridspec_kw={"height_ratios": [3, 1]}, sharex=True)
    plt.subplots_adjust(hspace=0)

    sig_weight=(df_eff[df_eff["Category"]==cat][" N_bsmmMc"].iloc[0])/len(signal[signal.cate==cat])*np.ones(len(signal[signal.cate==cat]))

    bkg_weight=np.ones(len(bkg[bkg.cate==cat]))*0.25/0.45

    ax[0].hist(func(signal[signal["cate"]==cat]["bdt"]),bins=bins,range=(0,5),weights=sig_weight,histtype="step",color="red", label="Signal",hatch="/")

    ax[0].hist(func(bkg[bkg["cate"]==cat]["bdt"]),bins=bins,range=(0,5),weights=bkg_weight,histtype="stepfilled",color="dodgerblue",label="Background")

    sig_h=hist.Hist(hist.axis.Regular(bins, 0,5, name="atanh(BDT)"))
    bkg_h=hist.Hist(hist.axis.Regular(bins, 0,5, name="atanh(BDT)"))

    sig_h.fill(func(signal[signal["cate"]==cat]["bdt"]),weight=sig_weight)

    bkg_h.fill(func(bkg[bkg["cate"]==cat]["bdt"]),weight=bkg_weight)
    cuts=sig_h.axes[0].edges[:-1]
    n_sig_arr=np.array([sig_h.integrate(0,hist.loc(cut),hist.loc(5)) for cut in cuts])
    n_bkg_arr=np.array([bkg_h.integrate(0,hist.loc(cut),hist.loc(5)) for cut in cuts])


    ax[0].set_yscale("log")
    ax[0].legend()
    ax[0].set_ylabel("Events")
    hep.cms.text("Private Work",ax=ax[0])
    hep.cms.lumitext(f"Category {cat} (13 TeV)",ax=ax[0])

    z=significance(n_sig_arr,n_bkg_arr)



    ax[1].stairs(z,sig_h.axes[0].edges,linewidth=2,color="red")
    z=np.nan_to_num(z,0)
    z[z>1e100]=0
    maxidx=np.argmax(z)
    maxcut=cuts[maxidx]
    best_cut.append(maxcut)
    ax[1].axvline(maxcut,linestyle="--",color="royalblue",label=f"BDT={np.tanh(maxcut):.4f}")

    ax[1].set_ylabel("Z")
    ax[1].legend()
    ax[1].set_xlabel("atanh(BDT)")
    ax[1].set_xlim(-0.2,5)
    ax[0].set_xlim(-0.2,5)
    ax[0].set_ylim(2e-2,1e4)
    fig.savefig(f"fig/bdt_cat{cat}.pdf")


# %%
for idx,cut in enumerate(best_cut):
    print(f"Cat:{idx} Score:{np.tanh(cut):4f}")
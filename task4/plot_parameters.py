import json
import yaml
import matplotlib.pyplot as plt
import mplhep as hep
import os
import pandas as pd
import numpy as np

with open('../python/configuration.yaml', 'r') as file:
    categories = yaml.safe_load(file)["categories"]


category_code = {
        '0': "2016BFC",
        '1': "2016BFF",
        '2': "2016GHC",
        '3': "2016GHF",
        '4': "2017C",
        '5': "2017F",
        '6': "2018C",
        '7': "2018F"
        }

def read_json_files(directory, variable):
    data = {"category": [], "value":[], "error":[]}
    with open(os.path.join(directory, "parameters.json"), 'r') as file:
        content = json.load(file)
    for c in content:
        if variable not in c: continue
        cat_inx = c[-1]
        category = category_code[cat_inx]
        data["category"].append(str(categories[category]["year"]) + "\n" + categories[category]["eta"])
        data["value"].append(content[c]["value"])
        data["error"].append(content[c]["error"])
                
    return pd.DataFrame(data)


def plot_parameters(data, parameter, mc, output_path):
    
    titles = data["category"]
    values = data["value"]
    errors = data["error"]
    plt.figure()
    x_pos = np.arange(len(titles))  # The label locations
    hep.style.use("CMS")


    # Create the bar plot
    fig, ax = plt.subplots()
    plt.errorbar(x_pos, np.array(values), yerr=np.array(errors), fmt='o', linestyle='', markersize=5)
    if "norm" not in parameter:
        hep.cms.label('$B^{+} \\rightarrow \mu^{-}/\psi \mu^{+}$ ', data=True, ax = ax)
    else:
        hep.cms.label('$B^{+} \\rightarrow J/\psi K^{+}$ ', data=True, ax = ax)
 
    # Add labels and title
    ax.set_xlabel('Category')
    ax.set_ylabel(parameter)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(titles, fontsize=10)

    # Show the plot
    #plt.show()
    plt.savefig(output_path+'/'+parameter+'_per_category.png')



output_path = "./results/"
for v in ["N_peak", "N_semilep", "N_comb", "bs_mean1_signalMC", "bs_sigma1_signal_Mc", "Eff_bs", "bs_mean2_signalMc", "bs_frac_signalMc", "bs_sigma2_signalMc", "Eff_norm"]:
    dataset = read_json_files("./", v)
    plot_parameters(dataset, v, False, output_path)
 

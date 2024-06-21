import os
import json
import pandas as pd
import matplotlib.pyplot as plt
import argparse 
import re
import yaml
import mplhep as hep 
import numpy as np
with open('configuration.yaml', 'r') as file:
    categories = yaml.safe_load(file)["categories"]

def parse_arguments():
    parser = argparse.ArgumentParser(description="PyROOT fitting script with user-defined arguments.")
    
    parser.add_argument('--parameter', type=str, default=None, help="Parameter to plot")
    parser.add_argument('--results', type=str, required=True, help="Results")
    parser.add_argument('--mc', action='store_true', help="Flag to indicate if the input data is Monte Carlo simulation.")
    parser.add_argument('--output_path', type=str, required=True, help="Output path to save the results.")
    
    return parser.parse_args()

def list_subdirectories(directory_path):
    all_entries = os.listdir(directory_path)
        
    # Filter out only the directories
    subdirectories = [directory_path+'/'+entry for entry in all_entries if os.path.isdir(os.path.join(directory_path, entry))]
        
    return subdirectories

def extract_category_name(directory_path):
    # Define the regex pattern for the category name (e.g., "2016HG")
    pattern = r'\b\d{4}[A-Z]{1,}\b'

    # Search for the pattern in the directory path
    match = re.search(pattern, directory_path)

    if match:
        return match.group(0)
    else:
        return None

# Step 1: Reading JSON files
def read_json_files(directory, variable):
    data = {"category": [], "value":[], "error":[]}
    directories = list_subdirectories(directory)
    for d in directories:
        category = extract_category_name(d)
        with open(os.path.join(d, "parameters.json"), 'r') as file:
            print(os.path.join(d, "parameters.json"))
            content = json.load(file)
            data["category"].append(str(categories[category]["year"]) + "\n" + categories[category]["eta"])
            data["value"].append(content[variable]["value"])
            data["error"].append(content[variable]["error"])
                
    return pd.DataFrame(data)

# Step 2: Plotting the data
def plot_parameters(data, parameter, mc, output_path):
    
    titles = data["category"]
    values = data["value"]
    errors = data["error"]
    plt.figure()
    x_pos = np.arange(len(titles))  # The label locations
    hep.style.use("CMS")


    # Create the bar plot
    fig, ax = plt.subplots()
    ax.bar(x_pos, values, yerr=errors, align='center', alpha=0.7, ecolor='black', capsize=10)
    hep.cms.label('$B^{+} \\rightarrow J/\psi K^{+}$ ', data=mc, ax = ax, fontsize=12.)


    # Add labels and title
    ax.set_xlabel('Category')
    ax.set_ylabel(parameter)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(titles)

    # Show the plot
    plt.show()
    plt.savefig(output_path+'/'+parameter+'_per_category.png')

if __name__ == "__main__":
    args = parse_arguments()
    parameter = args.parameter
    results_path = args.results
    mc = args.mc
    output_path = args.output_path
    data = read_json_files(results_path, parameter)
    plot_parameters(data, parameter, mc, output_path)

# Example usage
directory = './results/plots/'
main(directory)

import ROOT
import json

# Function to extract parameter values from RooWorkspace
def extract_parameters_from_workspace(workspace_file, workspace_name):
    # Open the ROOT file
    file = ROOT.TFile.Open(workspace_file)
    if not file or file.IsZombie():
        raise FileNotFoundError(f"File {workspace_file} could not be opened")

    # Retrieve the workspace
    workspace = file.Get(workspace_name)
    if not workspace:
        raise ValueError(f"Workspace {workspace_name} not found in file {workspace_file}")

    # Get all parameters from the workspace
    params = workspace.allVars()
    param_dict = {}

    # Loop over all parameters and store their values
    param_iter = params.createIterator()
    param = param_iter.Next()
    while param:
        param_dict[param.GetName()] = {"value" : param.getVal(), "error": param.getError()}
        param = param_iter.Next()

    # Close the file
    file.Close()

    return param_dict

# Function to save parameters to JSON file
def save_parameters_to_json(param_dict, json_file):
    with open(json_file, 'w') as f:
        json.dump(param_dict, f, indent=4)

# Main function
def main(workspace_file, workspace_name, json_file):
    param_dict = extract_parameters_from_workspace(workspace_file, workspace_name)
    save_parameters_to_json(param_dict, json_file)
    print(f"Parameters saved to {json_file}")

# Example usage
workspace_file = 'task_4_5.root'
workspace_name = 'wspace'
json_file = 'parameters.json'

main(workspace_file, workspace_name, json_file)

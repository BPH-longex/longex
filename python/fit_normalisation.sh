#!/bin/bash

# Define the categories in an array
categories=("2016BFC" "2016BFF" "2016GHC" "2016GHF" "2017C" "2017F" "2018C" "2018F") 

for cat in "${categories[@]}"
do
  # Create directory for the current category
  mkdir -p ./results/mc/$cat/
  mkdir -p ./results/data/$cat/

  # Run the first command
  python3 fit_normalisation.py --data_in /eos/user/c/cmsdas/2024/long-ex-bph/bupsikMc.root --mc --cat $cat --tree_name bupsikMc --output_path ./results/mc/$cat/

  # Run the second command
  python3 fit_normalisation.py --data_in /eos/user/c/cmsdas/2024/long-ex-bph/bupsikData.root --tree_name bupsikData --parameters ./results/mc/$cat/parameters.json --cat $cat --output_path ./results/data/$cat/

  python3 plot_fit_results.py --parameter n1 --results ./results/mc/ --output_path ./results/mc/

done



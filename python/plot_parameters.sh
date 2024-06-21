#!/bin/bash
parameters=("n1" "n2" "alpha1" "alpha2" "sig_mean" "sig_sigma")
for par in "${parameters[@]}"
do
   python3 plot_fit_results.py --parameter $par --results ./results/mc/ --mc --output_path ./results/mc/
done

#!/bin/bash
parameters=("scale" "bias" "jpsix_scale" "jpsix_shift" "n_sig" "n_comb" "n_jpsix"  )
for par in "${parameters[@]}"
do
   python3 plot_fit_results.py --parameter $par --results ./results/data/ --output_path ./results/data/
done



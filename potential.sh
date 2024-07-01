#!/bin/bash

rm ./potential_loop/data.npz
rm -r ./gfortcache

python_script="potentialV8.py"

## Type of problemns: Type = 1 (Equidist); Type = 2 (qbinary); Type = 3 (ternary);##
typeproblem=1


## Types of initialization
# typeinit = 1 Just one initialization and one optimization (For additional initializations and optimizations, add 2, 3, etc. to ninit_list) 
# Warning: if typeinit = 2 or 3, then make ninit_list=(1).
# typeinit = 2 N perturbations of the ground truth (The value of N must be provided) 
# typeinit = 3 N random initializations (The value of N must be provided) 
typeinit=3

N=100

nsites_list=(9 10 12)
ninit_list=(1)
nsources_list=(1)
nmesh_list=(128)
noise_coeff_list=(0.02)


# Loop through the arguments and run the Python script
for nsites in "${nsites_list[@]}"; do
for ninit in "${ninit_list[@]}"; do
for nsources in "${nsources_list[@]}"; do
for nmesh in "${nmesh_list[@]}"; do
for noise_coeff in "${noise_coeff_list[@]}"; do
    echo "Running $python_script with argument: $arg1"
    python3 potentialV8.py $nsites $ninit $nsources $nmesh $noise_coeff $typeproblem $typeinit $N
done
done
done
done
done

python3 createlatex.py

cd ./potential_loop

pdflatex results.tex



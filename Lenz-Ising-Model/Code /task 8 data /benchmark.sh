#!/bin/bash
set -e

for lattice_size in 40 60 80 100 
do 
    echo "Running with $lattice_size times $lattice_size grid  "

    OMP_NUM_THREADS=10 ./ising_8.exe $lattice_size 500000 "results_zoomed_in${lattice_size}.csv"
done

echo "All gridsizes finished!"
#!/bin/bash
set -e

for threads in 1 2 4 6 8 10 
do 
    echo "Running with $threads threads "

    OMP_NUM_THREADS=$threads ./ising_omp.exe 20 100000 "results_${threads}.csv"
done

echo "All benchmarks finished!"
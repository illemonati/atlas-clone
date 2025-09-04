#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate 5 samples in Hardy–Weinberg Equilibrium and write vcf file
$atlas simulate --vcf --type HW --sampleSize 3 --logFile simulate.out

# Calculate inbreeding coefficient
# Use only 1 burnin and 200 iterations for speed (keep the default values for real data)
$atlas inbreeding --vcf ATLAS_simulations.vcf.gz --numBurnin 1 --iterations 200 --logFile inbreeding.out

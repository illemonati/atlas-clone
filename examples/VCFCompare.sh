#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate 2 samples in Hardy–Weinberg Equilibrium and write vcf file
$atlas simulate --vcf --type HW --sampleSize 2 --depth 2,5,8 --chrLength 2321 --ploidy 2 --logFile simulate.out

# Compare the samples
$atlas VCFCompare --vcf ATLAS_simulations.vcf.gz --samples "ind_0,ind_1" --logFile VCFCompare.out

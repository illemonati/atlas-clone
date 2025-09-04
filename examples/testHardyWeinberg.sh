#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate 5 samples in Hardy–Weinberg Equilibrium and write vcf file
$atlas simulate --vcf --type HW --sampleSize 5 --logFile simulate.out

# Test for Hardy Weinberg equilibrium
$atlas testHardyWeinberg --vcf ATLAS_simulations.vcf.gz --logFile testHardyWeinberg.out

#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate 5 samples in Hardy–Weinberg Equilibrium and write vcf file
$atlas simulate --vcf --type HW --sampleSize 5 --logFile simulate.out

# Create allele-count file
$atlas alleleCounts --vcf ATLAS_simulations.vcf.gz --logFile alleleCounts.out

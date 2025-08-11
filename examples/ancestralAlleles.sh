#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate 5 samples in Hardy–Weinberg Equilibrium and write vcf file
$atlas simulate --vcf --type HW --sampleSize 5 --logFile simulate.out

# Create allele-count file with alleles
$atlas alleleCounts --vcf ATLAS_simulations.vcf.gz --outFormat withAlleles --logFile alleleCounts.out

# Create fasta-file with ancetral alleles using the calculated allele counts
$atlas ancestralAlleles --alleleCounts ATLAS_simulations_alleleCounts.txt.gz --logFile ancestralAlleles.out

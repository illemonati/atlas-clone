#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate 5 samples in Hardy–Weinberg Equilibrium and write vcf file
$atlas simulate --vcf --type HW --sampleSize 5 --logFile simulate.out

# Create bed file
echo "chr1 0 3333" > bed.bed
echo "chr1 5000 8888" >> bed.bed

# Create polymorphic Window
$atlas polymorphicWindows --vcf ATLAS_simulations.vcf.gz --regions bed.bed --logFile polymorphicWindows.out

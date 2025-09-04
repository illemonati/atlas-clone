#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate 5 BAM files in Hardy–Weinberg equilibrium
$atlas simulate --type HW --sampleSize 5 --logFile simulate.out

# Create GLF files
for f in *.bam; do
	$atlas GLF --bam $f
done
samples=$(ls -1 *.glf.gz | paste -s -d ',' -)

# Create site allele frequencies file
$atlas saf --glf $samples --fasta ATLAS_simulations.fasta --logFile saf.out

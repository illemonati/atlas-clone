#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate 3 BAM Files in Hardy–Weinberg Equilibrium
$atlas simulate --type HW --sampleSize 3 --logFile simulate.out

# Create GLF files
for f in *.bam; do
	$atlas GLF --bam $f
done
samples=$(ls -1 *.glf.gz | paste -s -d ',' -)

# Calcuate genetic distance between the samples
$atlas geneticDist --glf $samples --logFile geneticDist.out

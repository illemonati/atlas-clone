#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File with specific HKY85 genotype distribution
$atlas simulate --type "HKY85" --mu 0.71 --thetaG 0.0002 --thetaR 0.003 --logFile simulate.out 

# Estimate HKY85 genotype distribution
$atlas HKY85 --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --logFile HKY85.out

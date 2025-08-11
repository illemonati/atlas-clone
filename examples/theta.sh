#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File with specific theta
$atlas simulate  --theta 0.00123 --logFile simulate.out 

# Estimate Felsenstein genotype distribution
$atlas theta --bam ATLAS_simulations.bam --logFile theta.out

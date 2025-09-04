#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# exponential PMD pattern 
pmd="CT5:0.2*exp(-0.3*p)+0.07;GA3:0.2*exp(-0.25*p)+0.06"
# recal pattern, linear quality transformation and position dependency
recal="quality:polynomial[0.9];position:polynomial[0.02]"

# Simulate a BAM File with PMD and recal
$atlas simulate --pmd $pmd --recal $recal --logFile simulate.out 

# Estimate errors
$atlas estimateErrors --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --logFile estimateErrors.out

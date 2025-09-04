#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File with specific theta
$atlas simulate --chrLength "10000,20000" --theta "0.001,0.002" --logFile simulate.out 

# Create region files
echo "chr1 1 10000" > region1.bed
echo "chr2 1 10000" > region2.bed

# 
$atlas thetaRatio --bam ATLAS_simulations.bam --region1 region1.bed --region2 region2.bed --logFile theta.out

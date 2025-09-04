#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate 5 BAM Files in Hardy–Weinberg Equilibrium
$atlas simulate --type HW --sampleSize 5 --logFile simulate.out

# Create pileup file of all bam files and name it multiBam
$atlas pileup --bam "$(ls *.bam)" --logFile pileup.out --out multiBam

# Create bed file from pileup output with 5 <= depth < 30
$atlas pileupToBed --pileup multiBam_pileup.txt.gz --depth "[10,30)" --logFile pileupToBed.out

#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File
$atlas simulate --logFile simulate.out

# Create pileup file
$atlas pileup --bam ATLAS_simulations.bam --logFile pileup.out

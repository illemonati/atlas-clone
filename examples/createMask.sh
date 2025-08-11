#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File
$atlas simulate --logFile simulate.out

# Create mask for depth within [5, 10]
$atlas createMask --type depth --minDepth 5 --maxDepth 10 --bam ATLAS_simulations.bam --logFile createMask.out

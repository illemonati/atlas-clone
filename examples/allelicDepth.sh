#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File
$atlas simulate --logFile simulate.out

# Create allelicDepth histogram, limit depth to 100
$atlas allelicDepth --readUpToDepth 100 --bam ATLAS_simulations.bam --logFile allelicDepth.out

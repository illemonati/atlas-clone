#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File
$atlas simulate --logFile simulate.out

# Create read statistics
$atlas BAMDiagnostics --bam ATLAS_simulations.bam --logFile BAMDiagnostics.out

#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File
$atlas simulate --logFile simulate.out

# Generating a PSMC Input file
$atlas PSMC --bam ATLAS_simulations.bam --logFile PSMC.out

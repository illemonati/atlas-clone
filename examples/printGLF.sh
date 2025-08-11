#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File
$atlas simulate --logFile simulate.out

# Create GLF File
$atlas GLF --bam ATLAS_simulations.bam --logFile GLF.out

# Create txt-glf-file
$atlas printGLF --glf ATLAS_simulations.glf.gz --logFile printGLF.out

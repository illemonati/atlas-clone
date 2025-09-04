#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File with depth 15
$atlas simulate --logFile simulate.out --depth 15

# downsample to depth 3.14
$atlas downsample --bam ATLAS_simulations.bam --depth 3.14 --logFile downsample.out 

# simplify name of downsampled bam-file
mv ATLAS_simulations_downsampled_*.bam downsampled.bam
mv ATLAS_simulations_downsampled_*.bam.bai downsampled.bam.bai

# Calculate depth of downsampled BAM File
$atlas BAMDiagnostics --bam downsampled.bam --logFile BAMDiagnostics.out

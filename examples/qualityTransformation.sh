#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File with actual quality = 0.9*BaseQuality
$atlas simulate --logFile simulate.out --recal "quality:polynomial[0.9]"

# Create Quality Transformation
$atlas qualityTransformation --bam ATLAS_simulations.bam --logFile qualityTransformation.out

#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File
$atlas simulate --logFile simulate.out

# Assess softclipped bases
$atlas assessSoftClipping --bam ATLAS_simulations.bam --logFile assessSoftClipping.out

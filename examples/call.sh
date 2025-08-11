#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File
$atlas simulate --logFile simulate.out

# Call sites
$atlas call --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --logFile call.out

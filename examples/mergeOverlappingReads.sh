#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File with paired end reads
$atlas simulate --logFile simulate.out --seqType "paired"

# Merge reads
$atlas mergeOverlappingReads --bam ATLAS_simulations.bam --logFile mergeOverlappingReads.out

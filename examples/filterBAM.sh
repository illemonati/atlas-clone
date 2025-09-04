#! /bin/bash

# Set atlas path
atlas=$(dirname "$0")/../build/atlas

# Simulate a BAM File with 3 diploid chromosomes and 1 haploid chromosome
$atlas simulate --logFile simulate.out --ploidy "2{3},1"

# Keep only the 3 diploid chromosomes and MQ>50
$atlas filterBAM --bam ATLAS_simulations.bam --chr "chr1,chr2,chr3" --filterMQ "[50,256]" --logFile filterBAM.out

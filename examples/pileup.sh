#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --logFile simulate.out

# test bam-index with .bai instead of .bam.bai suffix
mv ATLAS_simulations.bam.bai ATLAS_simulations.bai

$atlas --task pileup --bam ATLAS_simulations.bam --fields depth,alleles --window 100000 --printAll --logFile pileup.out

#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

# test bam-index with .bai instead of .bam.bai suffix
mv ATLAS_simulations.bam.bai ATLAS_simulations.bai

$atlas --task pileup --bam ATLAS_simulations.bam --fields depth,alleles,mates,strands --window 100000 --histograms depth,allelicDepth,contexts,qualities --readUpToDepth 100 --printAll --fixedSeed 0 --out pileup --logFile pileup.out

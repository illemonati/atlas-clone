#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --logFile simulate.out
$atlas --task pileup --bam ATLAS_simulations.bam --fields depth,alleles --window 100000 --printAll --logFile pileup.out

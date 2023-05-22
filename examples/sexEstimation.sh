#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task sexEstimation --wholeGenome --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --fixedSeed 0 --out sexEstimation --logFile sexEstimation.out

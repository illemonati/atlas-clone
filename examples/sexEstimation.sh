#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --logFile simulate.out
$atlas --task sexEstimation --wholeGenome --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --fixedSeed 0 --out sexEstimation --logFile sexEstimation.out

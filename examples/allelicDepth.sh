#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --logFile simulate.out
$atlas --task allelicDepth --readUpToDepth 100 --bam ATLAS_simulations.bam --fixedSeed 0 --out allelicDepth --logFile allelicDepth.out

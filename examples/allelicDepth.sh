#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task allelicDepth --readUpToDepth 100 --bam ATLAS_simulations.bam --fixedSeed 0 --out allelicDepth --logFile allelicDepth.out

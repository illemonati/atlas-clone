#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task allelicDepth --readUpToDepth 97 \
	   --bam ATLAS_simulations.bam --window 4567 \
	   --fixedSeed 0 --out allelicDepth --logFile allelicDepth.out

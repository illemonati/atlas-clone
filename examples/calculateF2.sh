#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 11

$atlas --task calculateF2 --vcf ATLAS_simulations.vcf.gz \
	   --filterDepth [2,] --maxMissing 0.1 --minMAF 0.01 \
	   --fixedSeed 0 --out calculateF2 --logFile calculateF2.out


#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --chrLength 111111

$atlas --task genotypeDistribution --window 15345 --HKY85 \
	   --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta \
	   --out genotypeDistribution --logFile genotypeDistribution.out 

$atlas --task theta  --window 15345 --bam ATLAS_simulations.bam \
	   --out theta --logFile theta.out

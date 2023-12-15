#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 2 --chrLength 1111

$atlas --task alleleFreq --vcf ATLAS_simulations.vcf.gz \
	   --fixedSeed 0 --out alleleFreq --logFile alleleFreq.out

$atlas --task alleleFreq --likelihoods --vcf ATLAS_simulations.vcf.gz \
	   --fixedSeed 0 --out alleleFreqLKs --logFile alleleFreqLKs.out

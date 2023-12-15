#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 47 --chrLength 1212

$atlas inbreeding --numBurnin 1 --iterations 97 \
	   --vcf ATLAS_simulations.vcf.gz \
	   --fixedSeed 1 --out inbreeding --logFile inbreeding.out

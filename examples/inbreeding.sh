#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 47 --chrLength 1212 --fixedSeed 129

out="inbreeding"
$atlas inbreeding --numBurnin 1 --iterations 97 --vcf simulate.vcf.gz \
	   --fixedSeed 1 --out $out --logFile $out.out 2> $out.eout

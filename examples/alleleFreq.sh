#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 2 --chrLength 1111 --fixedSeed 19

out="alleleFreq"
$atlas --task alleleFreq --vcf simulate.vcf.gz \
	   --fixedSeed 18 --out $out --logFile $out.out 2> $out.eout

out="alleleFreqLKs"
$atlas --task alleleFreq --likelihoods --vcf simulate.vcf.gz \
	   --fixedSeed 17 --out $out --logFile $out.out 2> $out.eout

#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 11 --fixedSeed 35

out="calculateF"
$atlas --task calculateF2 --vcf simulate.vcf.gz \
	   --filterDepth [2,] --maxMissing 0.1 --minMAF 0.01 \
	   --fixedSeed 36 --out $out --logFile $out.out 2> $out.eout


#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 11 \
  --chrLength 13217 --ploidy 2 --fixedSeed 242

out="testHardyWeinberg"
$atlas --task testHardyWeinberg --vcf simulate.vcf.gz \
	   --fixedSeed 245 --out $out --logFile $out.out 2> $out.eout

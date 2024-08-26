#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 2 --fixedSeed 11

out="alleleCounts"
$atlas --task alleleCounts --vcf simulate.vcf.gz \
	   --fixedSeed 12 --out $out --logFile $out.out 2> $out.eout

out="alleleCountsSAF"
$atlas --task alleleCounts --dosaf --vcf simulate.vcf.gz \
	   --fixedSeed 13 --out $out --logFile $out.out 2> $out.eout

out="transform"
$atlas --task alleleCounts --transform alleleCounts_alleleCounts.txt.gz \
	   --fixedSeed 14 --out $out --logFile $out.out 2> $out.eout

#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --fixedSeed 28

out="alleleCounts"
$atlas --task alleleCounts --vcf simulate.vcf.gz --outFormat withAlleles \
	   --fixedSeed 27 --out $out --logFile $out.out 2> $out.eout

out="ancestralAlleles"
$atlas --task ancestralAlleles --alleleCounts alleleCounts_alleleCounts.txt.gz \
	   --fastaIndex simulate.fasta.fai \
	   --minorCountMaximum 2 --totalCountMinimum 4 \
	   --fixedSeed 26 --out $out --logFile $out.out 2> $out.eout

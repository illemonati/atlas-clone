#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --writeTrueGenotypes --fixedSeed 162

rm simulate.bam.bai # will automatically recreate it

# prepare alleles file
zcat simulate_trueGenotypes.vcf.gz | awk 'BEGIN{print "Chr", "Pos", "Allele1"}$1!~/#/{print $1, $2, $4}' > alleles.txt

out="mutationLoad"
$atlas --task mutationLoad --bam simulate.bam \
	   --alleles alleles.txt --window 4567 \
	   --fixedSeed 168 --out $out --logFile $out.out 2> $out.eout

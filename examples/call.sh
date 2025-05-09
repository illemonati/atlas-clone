#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 44 --writeTrueGenotypes

# prepare alleles file
gunzip -c simulate_trueGenotypes.vcf.gz | awk 'BEGIN{print "Chr", "Pos", "Ref", "Alt"}$1!~/#/{if ($5 == ".") {print $1, $2, $4} else {print $1, $2, $4, $5}}' > alleles.txt

out="callRef"
$atlas --task call --window 4567 \
	   --bam simulate.bam --fasta simulate.fasta \
	   --fixedSeed 42 --out $out --logFile $out.out 2> $out.eout

out="callAll"
$atlas --task call --window 4567 --alleles alleles.txt \
	   --bam simulate.bam \
	   --fixedSeed 42 --out $out --logFile $out.out 2> $out.eout

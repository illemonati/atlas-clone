#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 44 --writeTrueGenotypes

# prepare alleles file
#special cases:
echo "Chr Pos Ref Alt" > alleles.txt
echo "chr1 40 N" >> alleles.txt
echo "chr1 41 N N" >> alleles.txt
echo "chr1 42 N T" >> alleles.txt
echo "chr1 43 T N" >> alleles.txt
gunzip -c simulate_trueGenotypes.vcf.gz | tail -n +50 | awk '$1!~/#/{if ($5 == ".") {print $1, $2, $4} else {print $1, $2, $4, $5}}' >> alleles.txt

out="callRef"
$atlas --task call --window 4567 \
	   --bam simulate.bam --fasta simulate.fasta \
	   --fixedSeed 42 --out $out --logFile $out.out 2> $out.eout

out="callAll"
$atlas --task call --alleles alleles.txt --allowRefN \
	   --bam simulate.bam \
	   --fixedSeed 42 --out $out --logFile $out.out 2> $out.eout

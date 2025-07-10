#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 2 --chrLength 1111 --fixedSeed 19

out="alleleFreq"
$atlas --task alleleFreq --vcf simulate.vcf.gz \
	   --fixedSeed 18 --out $out --logFile $out.out 2> $out.eout

out="alleleFreqLKs"
$atlas --task alleleFreq --likelihoods --vcf simulate.vcf.gz \
	   --fixedSeed 17 --out $out --logFile $out.out 2> $out.eout

echo "sample population" > samples.txt
echo "ind_0 pop_0" >> samples.txt
echo "ind_1 pop_1" >> samples.txt

out="alleleFreqComp"
$atlas --task alleleFreq --compare --chr chr1 \
	   --vcf simulate.vcf.gz --samples samples.txt \
	   --fixedSeed 16 --out $out --logFile $out.out 2> $out.eout

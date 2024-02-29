#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --chrLength 132109,143210 --ploidy 2,1 --depth 8,11

echo "chr1 0 5" > diplo.bed
echo "chr2 5 10" > haplo.bed
for i in {1..13000}; do
	echo "chr1 ${i}0 ${i}5" >> diplo.bed
	echo "chr2 ${i}5 $((i+1))0" >> haplo.bed
done

name=windows
$atlas --task genotypeDistribution --window 15345 \
	   --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta \
	   --fixedSeed 0 --out $name --logFile $name.out

name=diplo
$atlas --task genotypeDistribution --window 15345 --regions diplo.bed --ploidy 2 \
	   --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta \
	   --fixedSeed 0 --out $name --logFile $name.out

name=haplo
$atlas --task genotypeDistribution --window 15345 --regions haplo.bed --ploidy 1 \
	   --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta \
	   --fixedSeed 0 --out $name --logFile $name.out

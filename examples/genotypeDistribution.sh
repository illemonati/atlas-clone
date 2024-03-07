#! /bin/bash

. $(dirname $0)/find_atlas

recal="intercept[0.1];quality:polynomial[0.9,0.01];position:polynomial[0.01]"
pmd="CT5:0.1*exp(-0.2*p)+0.001;GA3:0.1*exp(-0.2*p)+0.001"
pmd=""

. $(dirname $0)/simulate --chrLength 5021091,50321 --ploidy 2,1 --depth 11,8 --recal $recal --pmd $pmd --baseQuality "categorical(12,22,27,32,37,41)"

echo "chr1 0 5" > diplo.bed
echo "chr2 5 10" > haplo.bed
for i in {1..500001}; do
	echo "chr1 ${i}0 ${i}5" >> diplo.bed
	echo "chr2 ${i}5 $((i+1))0" >> haplo.bed
done

probs="0.5"
name=windows
$atlas --task genotypeDistribution --window 15345 --prob "1.0,$probs" \
	   --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta \
	   --recal $recal --pmd $pmd \
	   --fixedSeed 0 --out $name --logFile $name.out

probs="0.5,0.2,0.1,0.05,0.02,0.01"
name=diploRaw
$atlas --task genotypeDistribution --genomeWide --prob "$probs" \
	   --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta \
	   --regions diplo.bed --ploidy 2 \
	   --fixedSeed 0 --out $name --logFile $name.out

name=diploEE
$atlas --task genotypeDistribution  --genomeWide --prob "$probs" \
	   --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta \
	   --regions diplo.bed --ploidy 2 --recal $recal --pmd $pmd \
	   --fixedSeed 0 --out $name --logFile $name.out

name=thetaRaw
$atlas --task theta --genomeWide --prob "$probs" \
	   --bam ATLAS_simulations.bam --regions diplo.bed  \
	   --fixedSeed 0 --out $name --logFile $name.out

name=thetaEE
$atlas --task theta --genomeWide --prob "$probs" \
	   --bam ATLAS_simulations.bam --regions diplo.bed --recal $recal --pmd $pmd \
	   --fixedSeed 0 --out $name --logFile $name.out

# no downsampling
name=haplo
$atlas --task genotypeDistribution --genomeWide \
	   --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta \
	   --regions haplo.bed --ploidy 1 --recal $recal --pmd $pmd \
	   --fixedSeed 0 --out $name --logFile $name.out

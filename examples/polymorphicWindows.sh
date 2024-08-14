#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate_vcf --sampleSize 11

echo "chr1 1000 10000" > bed.bed
echo "chr2 2 3" >> bed.bed
echo "chr3 0 10" >> bed.bed
echo "chr3 100 110" >> bed.bed
echo "chr3 200 210" >> bed.bed
echo "chr3 300 310" >> bed.bed
echo "chr3 400 410" >> bed.bed
echo "chr3 500 510" >> bed.bed

name="polymorphicWindows"
$atlas --task polymorphicWindows --vcf ATLAS_simulations.vcf.gz --regions bed.bed \
	   --fixedSeed 0 --out $name.out --logFile $name.out


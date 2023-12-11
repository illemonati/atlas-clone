#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type HW --sampleSize 19

echo "chr1	0	4567" >> window.txt
echo "chr1	4567	9134" >> window.txt
echo "chr1	9134	11111" >> window.txt
echo "chr2	0	4567" >> window.txt
echo "chr2	4567	5432" >> window.txt
echo "chr3	0	4567" >> window.txt
echo "chr3	4567	9134" >> window.txt
echo "chr3	9134	12345" >> window.txt

$atlas --task pileup --bam ATLAS_simulations_ind1.bam --window window.txt  --printAll \
	   --histograms depth,allelicDepth,contexts,qualities --readUpToDepth 97 \
	   --fixedSeed 0 --out singleBam --logFile singleBam.out

bams=$(ls *.bam)
$atlas --task pileup --fields depth --bam "$bams"  \
	   --fixedSeed 0 --out multiBam --logFile multiBam.out

#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --type HW --sampleSize 19 --fixedSeed 177

echo "chr1	0	4567" >> window.txt
echo "chr1	4567	9134" >> window.txt
echo "chr1	9134	11111" >> window.txt
echo "chr2	0	4567" >> window.txt
echo "chr2	4567	5432" >> window.txt
echo "chr3	0	4567" >> window.txt
echo "chr3	4567	9134" >> window.txt
echo "chr3	9134	12345" >> window.txt

echo "chr1 0 1" > bed.bed
echo "chr1 2 3" >> bed.bed
echo "chr1 4 5" >> bed.bed
echo "chr1 6 7" >> bed.bed
echo "chr1 8 9" >> bed.bed
echo "chr1 10 111" >> bed.bed
echo "chr1 200 333" >> bed.bed
echo "chr1 400 555" >> bed.bed
echo "chr1 600 777" >> bed.bed
echo "chr1 800 999" >> bed.bed
echo "chr1 1000 11111" >> bed.bed

out="default"
$atlas --task pileup \
	   --bam simulate_ind1.bam --fasta simulate.fasta \
	   --fixedSeed 171 --out $out --logFile $out.out 2> $out.eout

out="printAll"
$atlas --task pileup --printAll \
	   --bam simulate_ind2.bam --fasta simulate.fasta \
	   --window window.txt  --readUpToDepth 97 \
	   --histograms depth,allelicDepth,contexts,qualities \
	   --fixedSeed 173 --out $out --logFile $out.out 2> $out.eout

out="regions"
$atlas --task pileup --printAll \
	   --bam simulate_ind3.bam --fasta simulate.fasta \
	   --regions bed.bed --histograms depth \
	   --fixedSeed 175 --out $out --logFile $out.out 2> $out.eout

out="multiBam"
bams=$(ls *.bam)
$atlas --task pileup --fields "depth,bases,sampleBases" --bam "$bams"  \
	   --fixedSeed 179 --out $out --logFile $out.out 2> $out.eout

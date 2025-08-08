#! /bin/bash

. $(dirname $0)/find_atlas

. $(dirname $0)/simulate  --fixedSeed 176 --chrLength 1234{7} --depth 3 --ploidy 2

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

echo "chr1 0 10000" > region.bed
echo "chr2 100 1000" >> region.bed
echo "chr3 0 20" >> region.bed
echo "chr3 100 120" >> region.bed
echo "chr3 500 876" >> region.bed
echo "chr3 1000 2000" >> region.bed

out="default"
$atlas --task pileup \
	   --bam simulate.bam --fasta simulate.fasta \
	   --fixedSeed 171 --out $out --logFile $out.out 2> $out.eout

out="filter"
$atlas --task filterBAM --bam simulate.bam --chr "chr1,chr7" --downsampleReads 0.3 \
        --fixedSeed 172 --out $out --logFile $out.out 2> $out.eout

out="filtered"
$atlas --task pileup \
	   --bam filter_filtered.bam --fasta simulate.fasta \
	   --fixedSeed 173 --out $out --logFile $out.out 2> $out.eout

out="bases"
$atlas --task pileup --shuffleSites --fields "bases" --downsampleSites 0.5 \
	   --bam simulate.bam --fasta simulate.fasta \
	   --fixedSeed 174 --out $out --logFile $out.out 2> $out.eout

out="printAll"
$atlas --task pileup --printAll --regions region.bed --mask bed.bed \
	   --bam simulate.bam --fasta simulate.fasta \
	   --window window.txt  --readUpToDepth 97 \
	   --histograms depth,allelicDepth,contexts,qualities,transitions,prevBases \
	   --fixedSeed 175 --out $out --logFile $out.out 2> $out.eout

out="regions"
$atlas --task pileup --printAll \
	   --bam simulate.bam --fasta simulate.fasta \
	   --regions bed.bed --histograms depth \
	   --fixedSeed 176 --out $out --logFile $out.out 2> $out.eout

. $(dirname $0)/simulate --type HW --sampleSize 19 --fixedSeed 177
out="multiBam"
bams=$(ls *.bam)
$atlas --task pileup --fields "depth,bases,sampleBases" --bam "$bams" \
	   --fixedSeed 178 --out $out --logFile $out.out 2> $out.eout

out="p2b"
$atlas --task pileupToBed --pileup multiBam_pileup.txt.gz --depth "[7,30)" \
	   --fixedSeed 179 --out $out --logFile $out.out 2> $out.eout

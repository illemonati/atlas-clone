#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --chrLength 15111 --depth 8 --ploidy 2 --numReadGroups 1

rm ATLAS_simulations.bam.bai # will automatically recreate it

echo "chr1 0 1000" > bed.bed
echo "chr1 2000 3000" >> bed.bed
echo "chr1 4000 5000" >> bed.bed
echo "chr1 6000 7000" >> bed.bed
echo "chr1 8000 9000" >> bed.bed
echo "chr1 10000 11000" >> bed.bed
echo "chr1 12000 13000" >> bed.bed

# theta
$atlas --task theta --bam ATLAS_simulations.bam --printAll \
	   --extraVerbose --window 4567 --regions bed.bed \
	   --fixedSeed 0 --out standard --logFile standard.out

# theta with downsample
$atlas --task theta --bam ATLAS_simulations.bam --prob 1,0.5 \
	   --extraVerbose --window 4567 \
	   --fixedSeed 0 --out downsample --logFile downsample.out 

# theta genomewide
$atlas --task theta --genomeWide --bootstraps 3  --window 4567 \
	   --bam ATLAS_simulations.bam --extraVerbose  --regions bed.bed \
	   --fixedSeed 0 --out genomewide --logFile genomewide.out

# theta genomewide with downsample
$atlas --task theta  --prob 1,0.5 --genomeWide --bootstraps 3  --window 4567 \
	   --bam ATLAS_simulations.bam --extraVerbose \
	   --fixedSeed 0 --out genomewide_downsample --logFile genomewide_downsample.out

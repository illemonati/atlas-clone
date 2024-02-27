#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --chrLength 15111 --depth 8 --ploidy 2 --numReadGroups 1

rm ATLAS_simulations.bam.bai # will automatically recreate it

echo "chr1 0 1" > bed.bed
for i in {2..15110..2}; do
	echo "chr1 $i $((i+1))" >> bed.bed
done

# theta
name=standard
$atlas --task theta --bam ATLAS_simulations.bam --printAll \
	   --extraVerbose --window 4567  \
	   --fixedSeed 0 --out $name --logFile $name.out

# theta with downsample
name=downsample
$atlas --task theta --bam ATLAS_simulations.bam --prob 1,0.5 \
	   --extraVerbose --window 4567 \
	   --fixedSeed 0 --out $name --logFile $name.out

# theta genomewide
name=genomewide
$atlas --task theta --genomeWide --bootstraps 3  --window 4567 \
	   --bam ATLAS_simulations.bam --extraVerbose  --regions bed.bed \
	   --fixedSeed 0 --out $name --logFile $name.out

# theta genomewide with downsample
name=genomewide_downsample
$atlas --task theta --genomeWide --bootstraps 3  --window 4567 --prob 1,0.5 \
	   --bam ATLAS_simulations.bam --extraVerbose  --regions bed.bed \
	   --fixedSeed 0 --out $name --logFile $name.out

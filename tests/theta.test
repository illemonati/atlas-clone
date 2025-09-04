#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --numReadGroups 1 --fixedSeed 266 \
  --chrLength 15111 --depth 8 --ploidy 2 

rm simulate.bam.bai # will automatically recreate it

echo "chr1 0 1" > bed.bed
for i in {2..15110..2}; do
	echo "chr1 $i $((i+1))" >> bed.bed
done

# theta
out=standard
$atlas --task theta --bam simulate.bam --printAll \
	   --extraVerbose --window 4567  --filterDepth 4,16 \
	   --fixedSeed 263 --out $out --logFile $out.out 2> $out.eout

# theta with downsample
out=downsample
$atlas --task theta --bam simulate.bam --depth "20,5" \
	   --extraVerbose --window 4567 \
	   --fixedSeed 260 --out $out --logFile $out.out 2> $out.eout

# theta genomewide
out=genomewide
$atlas --task theta --genomeWide --bootstraps 3 --window 4567 \
	   --bam simulate.bam --extraVerbose  --regions bed.bed \
	   --fixedSeed 269 --out $out --logFile $out.out 2> $out.eout

# theta genomewide with downsample
out=genomewide_downsample
$atlas --task theta --genomeWide --bootstraps 3  --window 4567 --prob 1,0.5 \
	   --bam simulate.bam --extraVerbose  --regions bed.bed \
	   --fixedSeed 266 --out $out --logFile $out.out 2> $out.eout

#! /bin/bash

. $(dirname $0)/find_atlas

echo "chr1 1 10000" > region1.bed
echo "chr2 1 10000" > region2.bed

. $(dirname $0)/simulate --fixedSeed 250

out="thetaRatio"
$atlas --task thetaRatio --bam simulate.bam \
	   --region1 region1.bed --region2 region2.bed \
	   --fixedSeed 255 --out $out --logFile $out.out 2> $out.eout

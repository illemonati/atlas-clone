#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 66

out="downsample_reads"
$atlas --task downsample --bam simulate.bam --prob 0.5,0.3,0.1 \
	   --fixedSeed 67 --out $out --logFile $out.out 2> $out.eout

out="downsample_separate"
$atlas --task downsample --separateReads \
	   --bam simulate.bam --prob 0.5,0.3,0.1 \
	   --fixedSeed 68 --out $out --logFile $out.out 2> $out.eout

out="downsample_bases"
$atlas --task downsample --downsampleBases \
	   --bam simulate.bam --prob 0.5,0.3,0.1 \
	   --fixedSeed 69 --out $out --logFile $out.out 2> $out.eout

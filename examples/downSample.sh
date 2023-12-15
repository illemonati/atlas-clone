#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task downsample --bam ATLAS_simulations.bam --prob 0.5,0.3,0.1 \
	   --fixedSeed 0 --out downsample_reads --logFile downsample_reads.out

$atlas --task downsample --separateReads --bam ATLAS_simulations.bam --prob 0.5,0.3,0.1 \
	   --fixedSeed 0 --out downsample_separate --logFile downsample_separate.out

$atlas --task downsample --downsampleBases --bam ATLAS_simulations.bam --prob 0.5,0.3,0.1 \
	   --fixedSeed 0 --out downsample_bases --logFile downsample_bases.out

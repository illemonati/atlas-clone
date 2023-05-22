#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --out simulate 
$atlas --task qualityTransformation --bam simulate.bam --fixedSeed 0 --out qualityTransformation --logFile qualityTransformation.out

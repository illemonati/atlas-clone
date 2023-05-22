#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate 
$atlas --task qualityTransformation --bam ATLAS_simulations.bam --fixedSeed 0 --out qualityTransformation --logFile qualityTransformation.out

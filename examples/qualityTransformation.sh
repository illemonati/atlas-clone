#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --out simulate --logFile simulate.out 
$atlas --task qualityTransformation --bam simulate.bam --fixedSeed 0 --out qualityTransformation --logFile qualityTransformation.out

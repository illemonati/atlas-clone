#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --out simulate --logFile simulate.out 
$atlas --task assessSoftClipping --bam simulate.bam --fixedSeed 0 --out assessSoftClipping --logFile assessSoftClipping.out

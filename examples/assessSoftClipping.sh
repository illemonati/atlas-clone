#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task assessSoftClipping --bam simulate.bam --fixedSeed 0 --out assessSoftClipping --logFile assessSoftClipping.out

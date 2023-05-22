#! /bin/bash

. $(dirname $0)/find_atlas

echo "chr1 1 10000" > region1
echo "chr2 1 10000" > region2

. $(dirname $0)/simulate

$atlas --task thetaRatio --bam ATLAS_simulations.bam --region1 region1 --region2 region2 --fixedSeed 0 --out thetaRatio --logFile thetaRatio.out

#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --chrLength 100000 --fixedSeed 0 --logFile simulate.out
$atlas --task downsample --bam ATLAS_simulations.bam --prob 0.5,0.3,0.1 --fixedSeed 0 --logFile downsample.out

#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --logFile simulate.out

$atlas --task downsample --bam ATLAS_simulations.bam --prob 0.5,0.3,0.1 --fixedSeed 0 --out downsample --logFile downsample.out
$atlas --task downsample --separateReads --bam ATLAS_simulations.bam --prob 0.5,0.3,0.1 --fixedSeed 0 --out downsampleS --logFile downsampleS.out
$atlas --task downsample --writeN --bam ATLAS_simulations.bam --prob 0.5,0.3,0.1 --fixedSeed 0 --out downsampleN --logFile downsampleN.out

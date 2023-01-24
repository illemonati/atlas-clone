#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --logFile simulate.out
rm ATLAS_simulations.bam.bai
$atlas --task theta --bam ATLAS_simulations.bam --prob 1,0.5,0.2 --fixedSeed 0 --extraVerbose --logFile theta.out

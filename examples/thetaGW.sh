#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --logFile simulate.out
rm ATLAS_simulations.bam.bai
$atlas --task thetaGenomeWide --bootstraps 10 --bam ATLAS_simulations.bam --fixedSeed 0 --extraVerbose --logFile thetaGW.out

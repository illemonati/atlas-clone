#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --logFile simulate.out
$atlas --task PSMC --bam ATLAS_simulations.bam --fixedSeed 0 --logFile PSMC.out

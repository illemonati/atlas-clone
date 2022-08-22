#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --logFile simulate.out
$atlas --task theta --bam ATLAS_simulations.bam --printAll --fixedSeed 0 --extraVerbose --logFile theta.out

#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --logFile simulate.out
$atlas --task GLF --bam ATLAS_simulations.bam --printAll --logFile GLF.out

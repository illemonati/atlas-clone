#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task PSMC --bam ATLAS_simulations.bam --fixedSeed 0 --logFile PSMC.out

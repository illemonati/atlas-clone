#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task GLF --bam ATLAS_simulations.bam --printAll --fixedSeed 0 --out GLF --logFile GLF.out

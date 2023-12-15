#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task GLF --printAll --window 7777 \
	   --bam ATLAS_simulations.bam \
	   --fixedSeed 0 --out GLF --logFile GLF.out

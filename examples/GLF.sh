#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task GLF --printAll --window 7777 \
	   --bam ATLAS_simulations.bam \
	   --fixedSeed 0 --out GLF --logFile GLF.out

$atlas --task printGLF --glf GLF.glf.gz \
	   --fixedSeed 0 --logFile printGLF.out | tail -n +14 > GLF.txt

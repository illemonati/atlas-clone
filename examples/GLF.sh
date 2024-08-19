#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task GLF --printAll --window 7777 \
	   --bam ATLAS_simulations.bam \
	   --fixedSeed 0 --out GLF --logFile GLF.out

$atlas --task printGLF --glf GLF.glf.gz \
	   --fixedSeed 0 --logFile printGLF.out | tail -n +14 | head -n -4 > GLF.txt

mv GLF.glf.idx GLF.glf.idx.old

$atlas --task indexGLF --glf GLF.glf.gz \
	   --fixedSeed 0 --logFile indexGLF.out

$atlas --task indexGLF --glf GLF.glf.gz --check \
	   --fixedSeed 0 --logFile checkIndex.out

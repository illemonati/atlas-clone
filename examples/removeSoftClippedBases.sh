#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate
$atlas --task removeSoftClippedBases --bam ATLAS_simulations.bam --fixedSeed 0 --out removeSoftClippedBases --logFile removeSoftClippedBases.out 

for f in *.bam; do
	name=${f%.bam}
	$atlas --task assessSoftClipping --bam $f --fixedSeed 0 --out $name --logFile assess_$name.out
done

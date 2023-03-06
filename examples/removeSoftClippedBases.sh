#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --fixedSeed 0 --out simulate --logFile simulate.out 
$atlas --task removeSoftClippedBases --bam simulate.bam --fixedSeed 0 --out removeSoftClippedBases --logFile removeSoftClippedBases.out 

for f in *.bam; do
	name=${f%.bam}
	$atlas --task assessSoftClipping --bam $f --fixedSeed 0 --out $name --logFile assess_$name.out
done

#! /bin/bash

. $(dirname $0)/find_atlas

$atlas --task simulate --seqType single --fixedSeed 0 --out single --logFile single_simulate.out 
$atlas --task simulate --seqType paired --fixedSeed 0 --out paired --logFile paired_simulate.out 

for f in *.bam; do
	name=${f%.bam}
	$atlas --task readOverlap --bam $f --fixedSeed 0 --out $name --logFile $name.out
done

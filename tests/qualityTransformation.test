#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 199

out="qualityTransformation"
$atlas --task qualityTransformation --bam simulate.bam \
	   --fixedSeed 190 --out $out --logFile $out.out 2> $out.eout

#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 30

out="assessSoftClipping"
$atlas --task assessSoftClipping --bam simulate.bam \
	   --fixedSeed 31 --out $out --logFile $out.out 2> $out.eout

#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 277

out="trimSoftClips"
$atlas --task trimSoftClips --bam simulate.bam \
	   --fixedSeed 270 --out $out --logFile $out.out 2> $out.eout

for f in *.bam; do
	out=${f%.bam}
	$atlas --task assessSoftClipping --bam $f \
		   --fixedSeed 271 --out $out --logFile $out.out 2> $out.eout
done

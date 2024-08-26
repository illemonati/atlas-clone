#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --seqType single \
	--fixedSeed 200 --out single --logFile single_simulate.out 

. $(dirname $0)/simulate --seqType paired \
	--fixedSeed 201 --out paired --logFile paired_simulate.out 

for f in *.bam; do
	out=${f%.bam}
	$atlas --task readOverlap --bam $f \
		   --fixedSeed 202 --out $out --logFile $out.out 2> $out.eout
done

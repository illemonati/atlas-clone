#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 1

out="GLF"
$atlas --task GLF --printAll --window 7777 \
	   --bam simulate.bam \
	   --fixedSeed 2 --out $out --logFile $out.out 2> $out.eout

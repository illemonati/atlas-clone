#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 2

out="PSMC"
$atlas --task PSMC --bam simulate.bam --window 4567 \
	   --fixedSeed 1 --out $out --logFile $out.out 2> $out.eout

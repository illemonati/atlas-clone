#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 4321

echo "chr1 40 50" > mask.bed

out="GLF"
$atlas --task GLF --bam simulate.bam --haploid chr2 \
	   --fixedSeed 321 --out $out --logFile $out.out 2> $out.eout

out="printGLF"
$atlas --task printGLF --glf GLF.glf.gz --mask mask.bed \
	   --fixedSeed 3 --out $out --logFile $out.out 2> $out.eout

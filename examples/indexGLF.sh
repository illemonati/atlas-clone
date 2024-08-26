#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 1111

out="GLF"
$atlas --task GLF \
	   --bam simulate.bam \
	   --fixedSeed 1112 --out $out --logFile $out.out 2> $out.eout

mv GLF.glf.idx GLF.glf.idx.old

out="indexGLF"
$atlas --task indexGLF --glf GLF.glf.gz \
	   --fixedSeed 1113 --out $out --logFile $out.out 2> $out.eout

out="checkGLF"
$atlas --task indexGLF --glf GLF.glf.gz --check \
	   --fixedSeed 1114 --out $out --logFile $out.out 2> $out.eout

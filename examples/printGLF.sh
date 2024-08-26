#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 4321

out="GLF"
$atlas --task GLF --bam simulate.bam \
	   --fixedSeed 321 --out $out --logFile $out.out 2> $out.eout

out="printGLF"
$atlas --task printGLF --glf GLF.glf.gz \
	   --fixedSeed 3 --out $out --logFile $out.out 2> $out.eout \
	   | tail -n +14 | head -n -4 > GLF.txt

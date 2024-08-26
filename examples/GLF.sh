#! /bin/bash

# `--fixedSeed = N` is needed to have reproducable results in regression test

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 1

out="GLF"
$atlas --task GLF --printAll --window 7777 \
	   --bam simulate.bam \
	   --fixedSeed 2 --out $out --logFile $out.out 2> $out.eout

out="printGLF"
$atlas --task printGLF --glf GLF.glf.gz \
	   --fixedSeed 3 --out $out --logFile $out.out 2> $out.eout \
	   | tail -n +14 | head -n -4 > GLF.txt

mv GLF.glf.idx GLF.glf.idx.old

out="indexGLF"
$atlas --task indexGLF --glf GLF.glf.gz \
	   --fixedSeed 4 --out $out --logFile $out.out 2> $out.eout

out="checkGLF"
$atlas --task indexGLF --glf GLF.glf.gz --check \
	   --fixedSeed 5 --out $out --logFile $out.out 2> $out.eout

#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 50

out="depth"
$atlas --task createMask --type depth --bam simulate.bam --window 4567 \
	   --fixedSeed 51 --out $out --logFile $out.out 2> $out.eout

out="nonRef"
$atlas --task createMask --type nonRef --window 7654 \
	   --bam simulate.bam --fasta simulate.fasta \
	   --fixedSeed 52 --out $out --logFile $out.out 2> $out.eout

out="invariant"
$atlas --task createMask --type invariant --window 5674 --bam simulate.bam \
	   --fixedSeed 53 --out $out --logFile $out.out 2> $out.eout

out="variant"
$atlas --task createMask --type variant --window 6745 --bam simulate.bam \
	   --fixedSeed 54 --out $out --logFile $out.out 2> $out.eout

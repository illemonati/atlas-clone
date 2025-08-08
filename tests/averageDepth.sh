#! /bin/bash

. $(dirname $0)/find_atlas

out=simulate
$atlas --task simulate --depth 10,5 --chrLength 12345 --ploidy "2,1" \
	   --fixedSeed 0 --out $out --logFile $out.out 2> $out.eout

out=ploidy22
$atlas --task averageDepth --bam simulate.bam \
	   --fixedSeed 1 --out $out --logFile $out.out 2> $out.eout

out=ploidy21
$atlas --task averageDepth --bam simulate.bam --haploid "chr2" \
	   --fixedSeed 2 --out $out --logFile $out.out 2> $out.eout

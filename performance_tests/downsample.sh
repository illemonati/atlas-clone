#! /bin/bash

. $(dirname $0)/find_atlas

bname=$(basename $0)
name=${bname%.sh}
out=$name/$name

$atlas --task downsample \
	   --bam BAM/BAM.bam --prob 0.5,0.3,0.1 \
	   --fixedSeed 0 --out $out --logFile $out.out 2> $out.err > /dev/null

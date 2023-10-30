#! /bin/bash

. $(dirname $0)/find_atlas

bname=$(basename $0)
name=${bname%.sh}
out=$name/$name

$atlas --task filterBAM \
	   --bam BAM/BAM.bam \
	   --filterMQ 0,37 --filterReadLength 0,77 \
	   --filterFragmentLength 0,111 --filterSoftClips \
	   --fixedSeed 0 --out $out --logFile $out.out 2> $out.err > /dev/null

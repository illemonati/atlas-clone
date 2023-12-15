#! /bin/bash

. $(dirname $0)/find_atlas

bname=$(basename $0)
name=${bname%.sh}
out=$name/$name

echo "SimReadGroup1 SimReadGroup2" > mergeRG/rgs.txt

$atlas --task mergeRG \
	   --bam BAM/BAM.bam --readGroups mergeRG/rgs.txt \
	   --fixedSeed 0 --out $out --logFile $out.out 2> $out.err > /dev/null

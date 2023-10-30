#! /bin/bash

. $(dirname $0)/find_atlas

bname=$(basename $0)
name=${bname%.sh}
out=$name/$name

echo "readGroup seqType seqCycles" > splitMerge/rgs_paired.txt
echo "SimReadGroup1 paired 200" >> splitMerge/rgs_paired.txt

$atlas --task splitMerge --mergingMethod middle \
	   --bam BAM/BAM.bam --readGroupSettings splitMerge/rgs_paired.txt \
	   --fixedSeed 0 --out $out --logFile $out.out 2> $out.err > /dev/null 

#! /bin/bash

. $(dirname $0)/find_atlas

bname=$(basename $0)
name=${bname%.sh}
out=$name/$name

$atlas --task HKY85 --genomeWide \
	   --bam BAM/BAM.bam --fasta BAM/BAM.fasta --chr chr1 \
	   --fixedSeed 0 --out $out --logFile $out.out  2> $out.err > /dev/null

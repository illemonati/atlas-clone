#! /bin/bash

. $(dirname $0)/find_atlas

bname=$(basename $0)
name=${bname%.sh}
out=$name/$name

$atlas --task PMD --bam BAM/BAM.bam --fasta BAM/BAM.fasta --pmd "doubleStrand:Exponential:Empiric" --fixedSeed 0 --out $out --logFile $out.out 2> $out.err > /dev/null

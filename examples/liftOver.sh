#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --fixedSeed 44 --chrLength 1000{3} --depth 0 --out liftOverFasta

(
	echo "chr1 200 201 Info1"
	echo "chr1 400 401 Info2"
	echo "chr2 50 51 Info3"
	echo "chr2 950 951 Info4"
	echo "chr3 200 210 Info5"
) > liftOver.bed 
	

#out="call"
#$atlas --task call --window 4567 \
#	   --bam simulate.bam --fasta simulate.fasta \
#	   --fixedSeed 42 --out $out --logFile $out.out 2> $out.eout

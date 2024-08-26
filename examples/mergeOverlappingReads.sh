#! /bin/bash

. $(dirname $0)/find_atlas

# paired end
out="paired"
. $(dirname $0)/simulate --numReadGroups 10 --fixedSeed 140 \
  --seqType paired --out $out --logFile $out.out

echo "readGroup seqType seqCycles" > rgs_paired.txt
echo "SimReadGroup1 paired 100" >> rgs_paired.txt
echo "SimReadGroup3 paired 100" >> rgs_paired.txt
echo "SimReadGroup5 paired 100" >> rgs_paired.txt
echo "SimReadGroup7 paired 100" >> rgs_paired.txt
echo "SimReadGroup9 paired 100" >> rgs_paired.txt


for name in "middle" "firstMate" "secondMate" "highestQuality" "randomRead"; do
	out="$name"
	$atlas --task mergeOverlappingReads  --mergingMethod $name \
		   --bam paired.bam --readGroupSettings rgs_paired.txt\
		   --fixedSeed 141 --out $out --logFile $out.out 2> $out.eout

	out="${name}_2nd"
	$atlas --task mergeOverlappingReads  --mergingMethod $name \
		   --bam ${name}_merged.bam --readGroupSettings rgs_paired.txt\
		   --fixedSeed 142 --out $out --logFile $out.out 2> $out.eout

	if ! diff -q <(samtools view ${name}_merged.bam) <(samtools view ${name}_2nd_merged.bam) > /dev/null; then
		samtools view ${name}_merged.bam > ${name}_merged.sam
		samtools view ${name}_2nd_merged.bam > ${name}_2nd_merged.sam
		>&2 echo "${name}_merged.bam and  ${name}_2nd_merged.bam differ"
	fi
done

# single end
out="single"
. $(dirname $0)/simulate --numReadGroups 10 --fixedSeed 143 \
  --seqType single --out $out --logFile $out.out

echo "readGroup seqType seqCycles" > rgs_single.txt
echo "SimReadGroup1 single 100" >> rgs_single.txt
echo "SimReadGroup3 single 100" >> rgs_single.txt
echo "SimReadGroup5 single 100" >> rgs_single.txt
echo "SimReadGroup7 single 100" >> rgs_single.txt
echo "SimReadGroup9 single 100" >> rgs_single.txt

out="merge"
$atlas --task mergeOverlappingReads \
	   --bam single.bam --readGroupSettings rgs_single.txt \
	   --fixedSeed 144 --out $out --logFile $out.out 2> $out.eout

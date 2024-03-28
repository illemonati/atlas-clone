#! /bin/bash

. $(dirname $0)/find_atlas

# paired end
. $(dirname $0)/simulate --numReadGroups 10 --seqType paired --out paired --logFile paired.out
echo "readGroup seqType seqCycles" > rgs_paired.txt
echo "SimReadGroup1 paired 100" >> rgs_paired.txt
echo "SimReadGroup3 paired 100" >> rgs_paired.txt
echo "SimReadGroup5 paired 100" >> rgs_paired.txt
echo "SimReadGroup7 paired 100" >> rgs_paired.txt
echo "SimReadGroup9 paired 100" >> rgs_paired.txt


for name in "middle" "firstMate" "secondMate" "highestQuality" "randomRead"; do
	$atlas --task splitMerge  --mergingMethod $name \
		   --bam paired.bam --readGroupSettings rgs_paired.txt\
		   --fixedSeed 0 --out $name --logFile $name.out

	$atlas --task splitMerge  --mergingMethod $name \
		   --bam ${name}_splitMerged.bam --readGroupSettings rgs_paired.txt\
		   --fixedSeed 0 --out ${name}_2nd --logFile ${name}_2nd.out

	if ! diff -q <(samtools view ${name}_splitMerged.bam) <(samtools view ${name}_2nd_splitMerged.bam) > /dev/null; then
		samtools view ${name}_splitMerged.bam > ${name}_splitMerged.sam
		samtools view ${name}_2nd_splitMerged.bam > ${name}_2nd_splitMerged.sam
		>&2 echo "${name}_splitMerged.bam and  ${name}_2nd_splitMerged.bam differ"
	fi
done

# single end
. $(dirname $0)/simulate --numReadGroups 10 --seqType single --out single --logFile single.out

echo "readGroup seqType seqCycles" > rgs_single.txt
echo "SimReadGroup1 single 100" >> rgs_single.txt
echo "SimReadGroup3 single 100" >> rgs_single.txt
echo "SimReadGroup5 single 100" >> rgs_single.txt
echo "SimReadGroup7 single 100" >> rgs_single.txt
echo "SimReadGroup9 single 100" >> rgs_single.txt

$atlas --task splitMerge --bam single.bam --readGroupSettings rgs_single.txt \
	   --fixedSeed 0 --out split --logFile split.out

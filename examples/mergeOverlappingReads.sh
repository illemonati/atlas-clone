#! /bin/bash

. $(dirname $0)/find_atlas

echo '{"RG1":{"seqType": "paired"}, "RG2":{"seqType": "single"}, "RG3":{"seqType": "paired"}, "RG4":{"seqType": "single"}, "RG5":{"seqType": "paired"}}' > sim.json

# paired end
out="simulate"
$atlas --task simulate --RGInfo sim.json --chrLength 1111111 \
	   --fixedSeed 140 --out $out --logFile $out.out 2> $out.eout

echo "RG1 RG2" > rgs.txt

out="mergeRG"
$atlas --task mergeRG --bam simulate.bam --readGroups rgs.txt \
	   --fixedSeed 141 --out $out --logFile $out.out 2> $out.eout

bam=mergeRG_mergedRG.bam

out="simulate_bd"
$atlas --task BAMDiagnostics --bam $bam \
	   --fixedSeed 141 --logFile $out.out 2> $out.eout

for name in "middle" "firstMate" "secondMate" "highestQuality" "randomRead"; do
	out="$name"
	$atlas --task mergeOverlappingReads  --mergingMethod $name --bam $bam  \
		   --fixedSeed 142 --out $out --logFile $out.out 2> $out.eout

	out="${name}_bd"
	$atlas --task BAMDiagnostics --bam ${name}_merged.bam \
		   --fixedSeed 143 --logFile $out.out 2> $out.eout

	out="${name}_tt"
	$atlas --task transitionTable --bam ${name}_merged.bam --fasta simulate.fasta \
		   --fixedSeed 144 --out $out --logFile $out.out 2> $out.eout

	out="${name}_2nd"
	$atlas --task mergeOverlappingReads --mergingMethod $name --bam ${name}_merged.bam \
		   --fixedSeed 145 --out $out --logFile $out.out 2> $out.eout

	if ! diff -q <(samtools view ${name}_merged.bam) <(samtools view ${name}_2nd_merged.bam) > /dev/null; then
		samtools view ${name}_merged.bam > ${name}_merged.sam
		samtools view ${name}_2nd_merged.bam > ${name}_2nd_merged.sam
		>&2 echo "${name}_merged.bam and  ${name}_2nd_merged.bam differ"
	fi
done

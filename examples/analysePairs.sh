#! /bin/bash

. $(dirname $0)/find_atlas

echo '{"RG1":{"seqType": "paired"}}' > sim.json

# paired end
out="simulate"
$atlas --task simulate --RGInfo sim.json --chrLength 111111 \
	   --fixedSeed 1 --out $out --logFile $out.out 2> $out.eout

name=analysePairs
out=$name
$atlas --task analysePairs --bam simulate.bam \
	   --fixedSeed 2 --out $out --logFile $out.out 2> $out.eout

name=newMerge
out=$name
$atlas --task newMerge --bam simulate.bam \
	   --fixedSeed 3 --out $out --logFile $out.out 2> $out.eout

name=analyseNewMerge
out=$name
$atlas --task analysePairs --bam newMerge_merged.bam \
	   --fixedSeed 4 --out $out --logFile $out.out 2> $out.eout

name=oldMerge
out=$name
$atlas --task mergeOverlappingReads --bam simulate.bam \
	   --fixedSeed 3 --out $out --logFile $out.out 2> $out.eout

name=analyseOldMerge
out=$name
$atlas --task analysePairs --bam oldMerge_merged.bam \
	   --fixedSeed 4 --out $out --logFile $out.out 2> $out.eout

#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --duplicationRate 0.113 --fixedSeed 77

bam="simulate.bam"
out="original"
$atlas --task BAMDiagnostics --bam $bam --perChromosome \
	   --fixedSeed 71 --out $out --logFile $out.out 2> $out.eout

out="original_markDup"
$atlas --task BAMDiagnostics --bam $bam --markDuplicates --perChromosome \
	   --fixedSeed 73 --out $out --logFile $out.out 2> $out.eout

out="markDup"
#mark but keep
$atlas --task filterBAM --bam $bam --markDuplicates --keepDuplicates \
	   --fixedSeed 75 --out $out --logFile $out.out 2> $out.eout

out="markDupAgn"
$atlas --task filterBAM --bam $bam --markDuplicates --RGagnostic --keepDuplicates \
	   --fixedSeed 77 --out $out --logFile $out.out 2> $out.eout

bam="markDup_filtered.bam"
out="BD_markDup"
$atlas --task BAMDiagnostics --bam $bam --perChromosome \
	   --fixedSeed 79 --out $out --logFile $out.out 2> $out.eout

out="BD_markDup_keep"
$atlas --task BAMDiagnostics --bam $bam --keepDuplicates --perChromosome \
	   --fixedSeed 72 --out $out --logFile $out.out 2> $out.eout

out="BD_reset_markDup"
$atlas --task BAMDiagnostics --bam $bam --resetDuplicates --markDuplicates \
	   --fixedSeed 74 --out $out --logFile $out.out 2> $out.eout

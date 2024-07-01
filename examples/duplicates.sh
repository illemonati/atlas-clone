#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --duplicationRate 0.113

bam="ATLAS_simulations.bam"
name="BD_original"
$atlas --task BAMDiagnostics --bam $bam \
	   --out $name --perChromosome --logFile $name.out

name="BD_original_markDup"
$atlas --task BAMDiagnostics --bam $bam --markDuplicates \
	   --out $name --perChromosome --logFile $name.out

name="markDup"
#mark but keep
$atlas --task filterBAM --bam $bam \
	   --markDuplicates --keepDuplicates \
	   --out $name --logFile $name.out

name="markDupAgn"
$atlas --task filterBAM --bam $bam \
	   --markDuplicates --RGagnostic --keepDuplicates \
	   --out $name --logFile $name.out

bam="markDup_filtered.bam"
name="BD_markDup"
$atlas --task BAMDiagnostics --bam $bam \
	   --out $name --perChromosome --logFile $name.out

name="BD_markDup_keep"
$atlas --task BAMDiagnostics --bam $bam --keepDuplicates \
	   --out $name --perChromosome --logFile $name.out

name="BD_reset_markDup"
$atlas --task BAMDiagnostics --bam $bam \
	   --resetDuplicates --markDuplicates \
	   --out $name --perChromosome --logFile $name.out

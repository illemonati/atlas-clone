#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate

$atlas --task createMask --type depth --bam ATLAS_simulations.bam --fixedSeed 0 --out createMask_depth --logFile createMask_depth.out
$atlas --task createMask --type nonRef --bam ATLAS_simulations.bam --fasta ATLAS_simulations.fasta --fixedSeed 0 --out createMask_nonRef --logFile createMask_nonRef.out
$atlas --task createMask --type invariant --bam ATLAS_simulations.bam --fixedSeed 0 --out createMask_invariant --logFile createMask_invariant.out
$atlas --task createMask --type variant --bam ATLAS_simulations.bam --fixedSeed 0 --out createMask_variant --logFile createMask_variant.out

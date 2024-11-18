#! /bin/bash

. $(dirname $0)/find_atlas

recal="intercept[.0];quality:polynomial[0.95,0.01];position:polynomial[0.02]"
pmd="CT5:0.2*exp(-0.3*p)+0.01;GA3:0.2*exp(-0.3*p)+0.01"

delta=0.01
k="111"
L="200$k"

. $(dirname $0)/simulate --recal $recal --pmd $pmd  \
  --type "HKY85" --mu 0.55 --thetaG 0.00033 --thetaR 0.015 \
  --chrLength $L,$L --depth 20 --ploidy 2,1 --numReadGroups 3 --seqType paired \
  --baseQuality 'categorical(12:1237,22:845,27:1912,32:21069,37:34246,41:339557)' \
  --fragmentLength 'normal(40,10)[30,101]' \
  --mappingQuality 'normal(60,10)[30,60]' --fixedSeed 0

out="BAMDiagnostics"
$atlas --task BAMDiagnostics --bam simulate.bam --mergeInput \
	   --fixedSeed 1 --out $out --logFile $out.out 2> $out.eout

out="merge"
$atlas --task mergeOverlappingReads  \
	   --bam simulate.bam --readGroupSettings BAMDiagnostics_mergeInput.txt \
	   --fixedSeed 2 --out $out --logFile $out.out 2> $out.eout

depth="20,10,5,3,2,1"
out=HK85_unmerged_raw
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam simulate.bam --fasta simulate.fasta  --chr "chr1" \
	   --fixedSeed 3 --out $out --logFile $out.out 2> $out.eout

out=HK85_unmerged_truth
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam simulate.bam --fasta simulate.fasta  --chr "chr1" \
	   --recal $recal --pmd $pmd \
	   --fixedSeed 4 --out $out --logFile $out.out 2> $out.eout

out=HK85_merged_raw
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam merge_merged.bam --fasta simulate.fasta  --chr "chr1" \
	   --fixedSeed 3 --out $out --logFile $out.out 2> $out.eout

out=HK85_merged_truth
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam merge_merged.bam --fasta simulate.fasta  --chr "chr1" \
	   --recal $recal --pmd $pmd \
	   --fixedSeed 4 --out $out --logFile $out.out 2> $out.eout

out=theta_unmerged_truth
$atlas --task theta --genomeWide \
	   --bam simulate.bam --chr "chr1" \
	   --readUpToDepth 2  \
	   --recal $recal --pmd $pmd \
	   --fixedSeed 5 --out $out --logFile $out.out 2> $out.eout

out=theta_merged_truth
$atlas --task theta --genomeWide \
	   --bam merge_merged.bam --chr "chr1" \
	   --readUpToDepth 2 --sample upToDepth \
	   --recal $recal --pmd $pmd \
	   --fixedSeed 6 --out $out --logFile $out.out 2> $out.eout

out=ee_unmerged
$atlas --task estimateErrors --minDeltaLL $delta \
	   --bam simulate.bam --fasta simulate.fasta  --chr "chr1,chr2" \
	   --poolRecal "all" --poolPMD "all" \
	   --fixedSeed 7 --out $out --logFile $out.out 2> $out.eout

out=ee_merged
$atlas --task estimateErrors --minDeltaLL $delta \
	   --bam simulate.bam --fasta simulate.fasta  --chr "chr1,chr2" \
	   --poolRecal "all" --poolPMD "all" \
	   --fixedSeed 8 --out $out --logFile $out.out 2> $out.eout

out=HK85_unmerged_ee
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam simulate.bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo  ee_unmerged_RGInfo.json \
	   --fixedSeed 9 --out $out --logFile $out.out 2> $out.eout

out=HK85_merged_ee
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam merge_merged.bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo  ee_merged_RGInfo.json \
	   --fixedSeed 10 --out $out --logFile $out.out 2> $out.eout

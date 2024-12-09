#! /bin/bash

. $(dirname $0)/find_atlas

echo '{
  "RG1": {
  "readGroup": "RG1",
  "frequency": "1.0",
  "seqType": "paired",
  "seqCycles": "100",
  "fragmentLength": "normal(50,10)[30,101]",
  "baseQuality": "categorical(12:1237,22:845,27:1912,32:21069,37:34246,41:339557)",
  "mappingQuality": "normal(70,10)[30,60]",
  "softClipping": "poisson(0.5)[0,20]",
  "recal": {
    "Mate1": {
	  "intercept": 0.0,
	  "position": {"polynomial": [0.02]},
	  "quality": {"polynomial": [0.95,0.01]},
	  "rho": [
		[0.0, 0.275, 0.45, 0.275],
		[0.06, 0.0, 0.01, 0.93],
		[0.93, 0.01, 0.0, 0.06],
		[0.275, 0.45, 0.275, 0.0]]
      },
    "Mate2": {
	  "intercept": 2.0,
	  "position": {"polynomial": [0.02]},
	  "quality": {"polynomial": [0.95,0.01]},
	  "rho": [
		[0.0, 0.35, 0.25, 0.4],
		[0.16, 0.0, 0.01, 0.81],
		[0.81, 0.01, 0.0, 0.16],
		[0.4, 0.25, 0.35, 0.0]]
      }
    },
    "pmd": "CT5:0.27*exp(-0.3*p)+0.01"
  }
}' > workflow.json

delta=0.001

N=50
if [[ "$#" -eq 1 ]]; then
    N=$1
fi
k="111"
L="${N}$k"

out="simulate"
$atlas --task simulate --RGInfo "workflow.json" \
	   --type "HKY85" --mu 0.55 --thetaG 0.00033 --thetaR 0.015 \
	   --chrLength $L --depth 10 --ploidy 2 \
	   --fixedSeed 1 --out $out --logFile $out.out 2> $out.eout

out="BAMDiagnostics"
$atlas --task BAMDiagnostics --bam simulate.bam --mergeInput \
	   --fixedSeed 1 --out $out --logFile $out.out 2> $out.eout

out="merge"
$atlas --task mergeOverlappingReads  \
	   --bam simulate.bam --readGroupSettings BAMDiagnostics_mergeInput.txt \
	   --fixedSeed 2 --out $out --logFile $out.out 2> $out.eout

bam=merge_merged.bam

depth="10,5,3,2,1"
out=HK85upTo_raw
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --fixedSeed 3 --out $out --logFile $out.out 2> $out.eout

out=HK85upTo_truth
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo workflow.json \
	   --fixedSeed 4 --out $out --logFile $out.out 2> $out.eout

probs="0.5,0.2,0.1,0.05,0.02,0.01"
out=HK85reads_raw
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample reads \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --fixedSeed 5 --out $out --logFile $out.out 2> $out.eout

out=HK85reads_truth
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample reads \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo workflow.json \
	   --fixedSeed 6 --out $out --logFile $out.out 2> $out.eout

out=HK85sites_raw
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample sites \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --fixedSeed 6 --out $out --logFile $out.out 2> $out.eout

out=HK85sites_truth
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample sites \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo workflow.json \
	   --fixedSeed 6 --out $out --logFile $out.out 2> $out.eout

out=ee
$atlas --task estimateErrors --minDeltaLL $delta \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --poolRecal "all" --poolPMD "all" \
	   --fixedSeed 7 --out $out --logFile $out.out 2> $out.eout

ee=ee_RGInfo.json

out=HK85upTo_ee
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo $ee \
	   --fixedSeed 8 --out $out --logFile $out.out 2> $out.eout

out=HK85reads_ee
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample reads \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo $ee \
	   --fixedSeed 9 --out $out --logFile $out.out 2> $out.eout

out=HK85sites_ee
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample sites \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo $ee \
	   --fixedSeed 6 --out $out --logFile $out.out 2> $out.eout

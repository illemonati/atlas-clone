#! /bin/bash

. $(dirname $0)/find_atlas

echo '{
  "RG1": {
  "readGroup": "RG1",
  "seqType": "paired",
  "seqCycles": "100",
  "fragmentLength": "normal(50,50)[30,1025]",
  "baseQuality": "categorical(12:1,22:1,27:2,32:2,37:3,41:100)",
  "mappingQuality": "normal(70,10)[30,60]",
  "softClipping": "poisson(0.5)[0,20]",
  "recal": {
    "Mate1": {
	  "intercept": -1.0,
	  "position": {"polynomial": [0.02]},
	  "quality": {"polynomial": [0.95,0.01]},
      "fragmentLength": {"empiric": [[32, 1], [64, 0.5], [128, 0], [256, -0.5], [512, -1]]},
      "context": {"empiric": [[0, -1], [1, -1], [2, 1], [3, 1], [4, 0]]},
	  "rho": [
		[0.0, 0.9, 0.05, 0.05],
		[0.33, 0.0, 0.33, 0.33],
		[0.33, 0.33, 0.0, 0.33],
		[0.45, 0.45, 0.1, 0.0]]
      },
    "Mate2": {
	  "intercept": 0.0,
	  "position": {"polynomial": [0.04]},
	  "quality": {"polynomial": [0.95,0.01]},
      "fragmentLength": {"empiric": [[32, 1], [64, 0.5], [128, 0], [256, -0.5], [512, -1]]},
	  "rho": [
		[0.0, 0.1, 0.45, 0.45],
		[0.33, 0.0, 0.33, 0.33],
		[0.33, 0.33, 0.0, 0.33],
		[0.05, 0.05, 0.9, 0.0]]
      }
    },
    "pmd": "CT5:0.2*exp(-0.3*p)+0.05;GA5:0.05*exp(0*p)+0"
  },
  "RG2": {
  "readGroup": "RG2",
  "seqType": "single",
  "seqCycles": "150",
  "fragmentLength": "normal(200,200)[30,1025]",
  "baseQuality": "categorical(12:1,22:2,27:3,32:4,37:5,41:100)",
  "mappingQuality": "normal(60,10)[30,60]",
  "softClipping": "poisson(0.2)[0,5]",
  "recal": {
	  "intercept": -1.0,
	  "position": {"polynomial": [0.02]},
	  "quality": {"polynomial": [0.95,0.02]},
	  "fragmentLength": {"empiric": [[32, 1], [64, 0.5], [128, 0], [256, -0.5], [512, -1]]},
	  "rho": [
		[0.0, 0.2, 0.3, 0.5],
		[0.1, 0.0, 0.2, 0.7],
		[0.3, 0.4, 0.0, 0.4],
		[0.2, 0.2, 0.6, 0.0]]
	},
    "pmd": "CT5:0.2*exp(-0.3*p)+0.05;GA5:0.05*exp(0*p)+0"
  }
}' > workflow.json

delta=0.1
N=50

if [ $1 ]; then
    N=$1
	if [ $2 ]; then
		delta=$2
	fi
fi
k="111"
L="${N}$k"

out="simulate"
$atlas --task simulate --RGInfo "workflow.json" \
	   --type "HKY85" --mu 0.55 --thetaG 0.00033 --thetaR 0.015 \
	   --chrLength $L --depth 10 --ploidy 2 \
	   --fixedSeed 1 --out $out --logFile $out.out 2> $out.eout

out="pileup"
$atlas --task pileup --onlySummaries --histograms depth,transitions,strandMate \
	   --bam simulate.bam --fasta simulate.fasta \
	   --fixedSeed 2 --out $out --logFile $out.out 2> $out.eout

out="merge"
$atlas --task mergeOverlappingReads --bam simulate.bam  \
	   --fixedSeed 3 --out $out --logFile $out.out 2> $out.eout

bam=merge_merged.bam
depth="10,5,3,2,1"
probs="0.5,0.2,0.1,0.05,0.02,0.01"

out=HK85upTo_raw
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --fixedSeed 10 --out $out --logFile $out.out 2> $out.eout

out=HK85upTo_truth
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo workflow.json \
	   --fixedSeed 11 --out $out --logFile $out.out 2> $out.eout

out=HK85reads_raw
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample reads \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --fixedSeed 20 --out $out --logFile $out.out 2> $out.eout

out=HK85reads_truth
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample reads \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo workflow.json \
	   --fixedSeed 21 --out $out --logFile $out.out 2> $out.eout

out=HK85sites_raw
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample sites \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --fixedSeed 30 --out $out --logFile $out.out 2> $out.eout

out=HK85sites_truth
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample sites \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo workflow.json \
	   --fixedSeed 31 --out $out --logFile $out.out 2> $out.eout

out=ee
$atlas --task estimateErrors --minDeltaLL $delta --minData 1000 \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --fixedSeed 4 --out $out --logFile $out.out 2> $out.eout

ee=ee_RGInfo.json

out=HK85upTo_ee
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --depth $depth --sample upToDepth \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo $ee \
	   --fixedSeed 12 --out $out --logFile $out.out 2> $out.eout

out=HK85reads_ee
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample reads \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo $ee \
	   --fixedSeed 22 --out $out --logFile $out.out 2> $out.eout

out=HK85sites_ee
$atlas --task HKY85 --minDeltaLL $delta --genomeWide \
	   --prob $probs --sample sites \
	   --bam $bam --fasta simulate.fasta  --chr "chr1" \
	   --RGInfo $ee \
	   --fixedSeed 32 --out $out --logFile $out.out 2> $out.eout

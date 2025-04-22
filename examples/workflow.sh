#! /bin/bash

. $(dirname $0)/find_atlas

echo '{
  "paired": {
  "readGroup": "paired",
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
      "context": {"empiric": [[0, -1], [1, -1], [2, 1], [3, 1]]},
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
  },"single": {
  "readGroup": "single",
  "seqType": "single",
  "seqCycles": "100",
  "fragmentLength": "normal(100,70)[30,100]",
  "baseQuality": "categorical(12:1,22:1,27:2,32:2,37:3,41:100)",
  "mappingQuality": "normal(70,10)[30,60]",
  "softClipping": "poisson(0.5)[0,20]",
  "recal": {
	  "intercept": -1.0,
	  "position": {"polynomial": [0.02]},
	  "quality": {"polynomial": [0.95,0.01]},
      "fragmentLength": {"empiric": [[32, 1], [64, 0.5], [128, 0], [256, -0.5], [512, -1]]},
      "context": {"empiric": [[0, -1], [1, -1], [2, 1], [3, 1]]},
	  "rho": [
		[0.0, 0.9, 0.05, 0.05],
		[0.33, 0.0, 0.33, 0.33],
		[0.33, 0.33, 0.0, 0.33],
		[0.45, 0.45, 0.1, 0.0]]
    },
    "pmd": "CT5:0.2*exp(-0.3*p)+0.05;GA5:0.05*exp(0*p)+0"
  }
}' > workflow.json

delta=0.1
N=10

if [ $1 ]; then
    N=$1
	if [ $2 ]; then
		delta=$2
	fi
fi
k="111"
L="${N}$k"
probs="0.5,0.2,0.1,0.05"

out="simulate"
$atlas --task simulate --RGInfo "workflow.json" \
	--type "HKY85" --mu 0.55 --thetaG 0.00033 --thetaR 0.0075 \
	--chrLength $L --depth 10 --ploidy 2 --refBias 0.51 \
	--fixedSeed 1 --out $out --logFile $out.out 2> $out.eout
bam=$out.bam
fasta=$out.fasta

out="merge"
$atlas --task mergeOverlappingReads --mergingMethod middle --bam $bam \
		--fixedSeed 3 --out $out --logFile $out.out 2> $out.eout
bam=${out}_merged.bam 

out="pileup"
$atlas --task pileup --onlySummaries --histograms depth,transitions \
	--bam $bam --fasta $fasta \
	--fixedSeed 2 --out $out --logFile $out.out 2> $out.eout

out="bd"
$atlas --task BAMDiagnostics --bam $bam --writeMates \
	--fixedSeed 2 --out $out --logFile $out.out 2> $out.eout

out=HKY85_raw
$atlas --task HKY85 --minDeltaLL $delta --genomeWide 10 \
	--prob $probs --sample reads \
	--bam $bam --fasta $fasta  --chr "chr1" \
	--fixedSeed 20 --out $out --logFile $out.out 2> $out.eout

out=HKY85_truth
$atlas --task HKY85 --minDeltaLL $delta --genomeWide 10 \
	--prob $probs --sample reads \
	--bam $bam --fasta $fasta  --chr "chr1" \
	--RGInfo workflow.json \
	--fixedSeed 21 --out $out --logFile $out.out 2> $out.eout

out=ee
$atlas --task estimateErrors --minDeltaLL $delta \
	--bam $bam --fasta $fasta  --chr "chr1" \
	--fixedSeed 4 --out $out --logFile $out.out 2> $out.eout

ee=ee_RGInfo.json

out=HKY85_ee
$atlas --task HKY85 --minDeltaLL $delta --genomeWide 10 \
	--prob $probs --sample reads \
	--bam $bam --fasta $fasta  --chr "chr1" \
	--RGInfo $ee \
	--fixedSeed 22 --out $out --logFile $out.out 2> $out.eout

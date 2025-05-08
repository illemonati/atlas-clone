#! /bin/bash

k="321"
. $(dirname $0)/find_atlas

runSim() {
	out="${1}$2"
	$atlas --task simulate --numReadGroups 2 --type $2 --seqType paired --seqCycles 200 \
		   --chrLength 50$k{2},30$k,40$k,60$k --ploidy 2{3},1,2 --depth 10,0,8,0.5{2} \
		   --baseFreq 0.5,0.3,0.2,0 $3 \
		   --fragmentLength 'fixed(500)' --baseQuality 'binomial(95,0.01)[0,93]' \
		   --mappingQuality 'normal(60,10)[1,255]' --softClipping 'poisson(20)[0,50]' \
		   --pmd 'CT5:0.2*exp(-0.2*p)+0.02;GA3:0.3' --frequency 0.2 --fixedSeed 234 \
		   --out $out --logFile $out.out 2> $out.eout
}

for t in "one" "pair"; do
	runSim "" $t
done

for t in  "HKY85" "SFS" "HW" ; do
	runSim "single" $t "--sampleSize 1"
done

for t in "HKY85" "SFS" "HW" ; do
	runSim "multi" $t "--sampleSize 3"
done

out="fromFasta"
$atlas --task simulate --type HKY85 --fasta singleHKY85.fasta \
		   --out $out --logFile $out.out 2> $out.eout

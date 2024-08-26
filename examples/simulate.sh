#! /bin/bash

k="321"
. $(dirname $0)/find_atlas

. $(dirname $0)/simulate --type pair --seqType paired --seqCycles 200 \
	--chrLength 50$k{2},30$k,40$k,60$k --ploidy 2{3},1,2 --depth 10,8{2},5{2} \
	--phi 0.1{8},0.2 --baseFreq 0.5,0.3,0.2,0 --refDiv 0.5 \
	--fragmentLength 'fixed(500)' --baseQuality 'binomial(95,0.01)[0,93]' \
	--mappingQuality 'normal(60,10)[1,255]' --softClipping 'poisson(20)[0,50]' \
	--pmd 'CT5:0.2*exp(-0.2*p)+0.02;GA3:0.3' --frequency 0.2 --fixedSeed 234

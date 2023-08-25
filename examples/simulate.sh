#! /bin/bash

. $(dirname $0)/find_atlas
. $(dirname $0)/simulate --chrLength 50000{2},30000,40000,60000 --ploidy 2{3},1,2 --depth 10,8{2},5{2} --type pair --phi 0.1{8},0.2 --baseFreq 0.5,0.3,0.2,0 --refDiv 0.5 --seqType paired --seqCycles 200 --fragmentLength 'fixed(500)' --baseQuality 'binomial(95,0.01)[0,93]' --mappingQuality 'normal(60,10)[1,255]' --softClipping 'poisson(20)[0,50]' --pmd 'CT5:0.2*exp(-0.2*p)+0.02;GA3:0.3' --frequency 0.2

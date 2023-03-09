#! /bin/bash

. $(dirname $0)/find_atlas

mkdir BAMDiagnostics

file=BAMDiagnostics/times
if [ ! -f "$file" ]; then
        echo -e "BAMDiagnostics" > BAMDiagnostics/times
fi


timeFor10Runs=0
for i in {1..10}; do
start=`date +%s.%N`

$atlas --task BAMDiagnostics --diagnosticsPerChromosome --bam simulate/bamFile.bam --fixedSeed 0 --out BAMDiagnostics/out --logfile BAMDiagnostics/diagnostics.out

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo  "$timeFor10Runs" >> BAMDiagnostics/times
#add paste simulate/times BAMDiagnostics/times with header etc
#also count runs, put it in nice table like
#run    simulate        BAMDiagnostics
#1      1000.2          1200.6
#2      2000.7          1185.6

#! /bin/bash

. $(dirname $0)/find_atlas

mkdir GLF
file=GLF/times
if [ ! -f "$file" ]; then
	echo -e "GLF" > GLF/times
fi

timeFor10Runs=0
for i in {1..10}; do
start=`date +%s.%N`

$atlas --task GLF --bam simulate/bamFile.bam --printAll --logFile GLF/GLF.out --out GLF/GLF

end=`date +%s.%N`
runtime=$( echo "$end - $start" | bc -l )
timeFor10Runs=$(echo "$timeFor10Runs+$runtime" | bc -l)
done
echo -e "$timeFor10Runs" >> GLF/times

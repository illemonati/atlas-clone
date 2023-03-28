#!/bin/bash

cd $(dirname $0)

file=performanceTable

if [ ! -f $file ]; then
	./bamSimulate.sh
	cp simulate/BAMtimes performanceTable
	./vcfSimulate.sh
	paste performanceTable simulate/VCFtimes > tmp
	mv tmp performanceTable
	for i in *.sh; do
        	if [ $i != "bamSimulate.sh" ] && [ $i != "runPerformanceTests.sh" ] && [ $i != "vcfSimulate.sh" ]; then
         	       ./$i
         	        paste performanceTable ${i/.sh/}/times > tmp
         	        mv tmp performanceTable
        	fi
	done

	echo "total" > total	
	tail -n+2 performanceTable |  awk '{for(i=1;i<=NF;i++) t+=$i; printf "%.5f" t, printf "\n"; t=0}' >> total
	paste performanceTable total > tmp
	mv tmp performanceTable
	rm total

	lineNum=$(wc -l performanceTable | cut -d" " -f1)
	num=1
	echo "run" > runCount
	while [ $num -lt $lineNum ]; do
                echo "$num" >> runCount
                num=$(($num+1))
        done
	paste runCount performanceTable > tmp
	mv tmp performanceTable 
	rm runCount
else
	./bamSimulate.sh
	tail -1 simulate/BAMtimes > newLine
	./vcfSimulate.sh
	tail -1 simulate/VCFtimes > tmp
	paste newLine tmp > tmp2
	mv tmp2 newLine
	rm tmp
	for i in *.sh; do
        	if [ $i != "bamSimulate.sh" ] && [ $i != "runPerformanceTests.sh" ] && [ $i != "vcfSimulate.sh" ]; then
         	       ./$i
			tail -1 ${i/.sh/}/times > tmp
			paste newLine tmp > tmp2	
			mv tmp2 newLine
			rm tmp
        	fi
	done
	
	cat newLine |  awk '{for(i=1;i<=NF;i++) t+=$i; printf "%.5f", t; printf "\n"}' > tmp
	paste newLine tmp > tmp2
	mv tmp2 newLine
	rm tmp

	runNum=$(wc -l performanceTable | cut -d" " -f1)
	echo "$runNum" > runNum
	paste runNum newLine > tmp
	mv tmp newLine
	rm runNum

	cat performanceTable newLine > tmp
	mv tmp performanceTable
	rm newLine
fi

lineNum=$(wc -l $file | cut -d" " -f1)
if [ $lineNum -gt 11 ]; then
	tail -11 $file | head -10 > tmp
	awk '{for (o=2;o<=NF;o++){sum[o]+=$o;a[o,NR]=$o}}END{for (o=2;o<=NF;o++){for(i=1;i<=NR;i++)y+=(a[o,i]-(sum[o]/NR))^2;printf "%.3f", sum[o]/NR; printf "\t";printf "%.3f", sqrt(y/(NR-1));printf "\n"; y=0}}' tmp > averagesAndSD
	rm tmp
	tail -1 $file > newValues
	colNum=$(wc -l averagesAndSD | cut -d" " -f1)
	i=1
	while [ $i -le $colNum ]; do
		average=$(head -"$i" averagesAndSD | tail -1 | cut -f1)
		SD=$(head -"$i" averagesAndSD | tail -1 | cut -f2)
		column=$(($i+1))
		compValue=$(awk -v i=$column '{print $i}' newValues)
		comparison=$(echo "$compValue - $average" | bc -l | tr -d "-")
		if (( $(echo "$comparison > $SD" | bc -l) )); then
			taskName=$(head -1 $file | cut -f$column)
			echo "$taskName deviates from the average of the last 10 runs by over one standard deviation!"
		fi
		i=$(($i+1))
	done
	rm newValues averagesAndSD
fi


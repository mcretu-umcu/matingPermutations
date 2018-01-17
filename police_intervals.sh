#!/bin/bash





combFile=$1
vcfFile=$2
pedFile=$3
tempDir=$4
horse=$5
logs=$6
jobN=$7
rt=$8


Ncombs=$(wc -l $combFile | awk '{print $1}')

for (( i=1; i<=Ncombs; i++))
do
	if [ ! -f $tempDir/simmilarity.combination.interval.$i ]; then
		pythonArgs=$(awk -v lineNR=$i 'NR == lineNR {print}' $combFile)
		echo "python $horse $vcfFile $tempDir/simmilarity.combination.interval.$SGE_TASK_ID $pedFile $pythonArgs"
		echo "module load python/2.7.10; python $horse $vcfFile $tempDir/simmilarity.combination.interval.$i $pedFile $pythonArgs " | qsub -N $jobN -l h_rt=$rt,h_vmem=5G -cwd -e $logs/police.interval.${jobN}.e -o $logs/police.interval.${jobN}.o
	fi
done



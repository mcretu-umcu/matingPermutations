#!/bin/bash





combFile=$1
pedFile=$2
tempDir=$3
horse=$4
logs=$5
jobN=$6
rt=$7
pedFam=$8


Ncombs=$(wc -l $combFile | awk '{print $1}')



echo "%%% Police is on the roll... %%%"
for (( i=1; i<=Ncombs; i++))
do
	if [ ! -f $tempDir/simmilarity.combination.$i ]; then
		pythonArgs=$(awk -v lineNR=$i 'NR == lineNR {print}' $combFile)
		echo "python $horse $pedFile $tempDir/simmilarity.combination.$SGE_TASK_ID $pythonArgs"
		echo "started"
		echo "module load python/2.7.10; python $horse $pedFile $tempDir/simmilarity.combination.$i $pedFam $pythonArgs " | qsub -N $jobN -l h_rt=$rt,h_vmem=5G -cwd -e $logs/police.${jobN}.e -o $logs/police.${jobN}.o
	fi
done



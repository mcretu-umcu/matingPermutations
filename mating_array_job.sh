#!/bin/bash





combFile=$1
pedFile=$2
tempDir=$3
horse=$4
pedFam=$5

pythonArgs=$(awk -v lineNR=$SGE_TASK_ID 'NR == lineNR {print}' $combFile)


module load python/2.7.10
echo "python $horse $pedFile $tempDir/simmilarity.combination.$SGE_TASK_ID $pedFam $pythonArgs"
echo "started"
python $horse $pedFile $tempDir/simmilarity.combination.$SGE_TASK_ID $pedFam $pythonArgs






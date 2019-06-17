#!/bin/bash 
stamp="$5"
output="$6"
cd ..
docker run --name=lumpy -v $(pwd)/:/in/ -w /in/ erictdawson/lumpy-sv lumpyexpress -B "$1".bam -S "$1".splitters.bam -D "$1".discordants.bam -o "$output"/lumpy/lumpyoutput_"$stamp".vcf

./lumpy/lumpy_awk_filter.sh "$output"/lumpy/lumpyoutput_"$stamp".vcf 100 "$output"/lumpy/lumpyoutput_CUTOFF_"$stamp".vcf

#rm ./lumpy/output/lumpyoutput_"$stamp".vcf

START=$(docker inspect --format='{{.State.StartedAt}}' lumpy)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' lumpy)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $FTIME seconds 2>&1 | tee "$output"/lumpy/lumpy_runtime_"$stamp".txt

docker container rm lumpy

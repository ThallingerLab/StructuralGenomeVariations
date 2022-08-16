#!/bin/bash 
stamp="$5"
bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
output="$6"
cd ..
mkdir "$output"/softsv/softsv_output_"$stamp"
docker run --name=softsv -v $(pwd)/:/in/ -w /in/ chrishah/softsv:1.4.2 SoftSV -i "$bam".bam -o "$output"/softsv/softsv_output_"$stamp"

./softsv/softsv_awk_filter.sh "$output"/softsv/softsv_output_"$stamp" 100

START=$(docker inspect --format='{{.State.StartedAt}}' softsv)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' softsv)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

TIME1=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $TIME1 seconds 2>&1 | tee "$output"/softsv/softsv_runtime_"$stamp".txt

docker container rm softsv

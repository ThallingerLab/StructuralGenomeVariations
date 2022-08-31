#!/bin/bash 

threads=8

bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
outdir="$5"
threads="$6"
tools_dir="$7"

log_eval $PWD "docker run --name=softsv -v $(pwd):$(pwd) -w $outdir chrishah/softsv:1.4.2 SoftSV \
  -i $bam -o $outdir"

log_eval $PWD "docker run --rm --name=python3.8 -v $(pwd):$(pwd) -w $outdir python:3.8-buster python3 $tools_dir/softsv/SoftSVtoVCF.py \
  -i $outdir -r $fasta -c LT962478.1:2263458 LT962479.1:1827941"
# At some point i should replace the fixed contig entries
#./softsv/softsv_awk_filter.sh "$output"/softsv/softsv_output_"$stamp" 100

START=$(docker inspect --format='{{.State.StartedAt}}' softsv)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' softsv)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

TIME1=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $TIME1 seconds 2>&1 | tee "$outdir"/softsv_runtime.txt

docker container rm softsv

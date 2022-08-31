#!/bin/bash

threads=8

bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
outdir="$5"
threads="$6"

log_eval $PWD "docker run --name=wham -v $(pwd):$(pwd) -w $(pwd) gatksv/wham:8645aa whamg \
-a $bam \
2> $outdir/wham.err > $outdir/wham.vcf"

START=$(docker inspect --format='{{.State.StartedAt}}' wham)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' wham)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $FTIME seconds 2>&1 | tee "$outdir/wham_runtime.txt"

docker container rm wham

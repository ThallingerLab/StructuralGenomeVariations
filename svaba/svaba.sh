#!/bin/bash

threads=8

bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
outdir="$5"
threads="$6"
tools_dir="$7"
timing="$8"

log_eval $PWD "docker run --name=svaba -v $(pwd):$(pwd) -w $outdir iarcbioinfo/svaba-nf:v1.0 svaba run \
-t $bam \
-a $(basename $outdir) \
-G $fasta \
--germline"

START=$(docker inspect --format='{{.State.StartedAt}}' svaba)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' svaba)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $FTIME seconds 2>&1
echo -e "${outdir}\t${FTIME}" | tee -a "$timing"

docker container rm svaba

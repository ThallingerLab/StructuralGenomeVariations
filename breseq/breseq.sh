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


# Think about using bwa BAM files, running bowtie on all files will take forever, especially with so many settings ?!

log_eval $PWD "docker run --name=breseq -u 1001:1001 -v $(pwd)/:$(pwd) -w $outdir pvstodghill/breseq:0.35.7__2021-08-03 breseq \
  -o $outdir -r $fasta \
  --brief-html-output \
  --no-javascript \
  -j $threads \
  $fastq1 $fastq2"

START=$(docker inspect --format='{{.State.StartedAt}}' breseq)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' breseq)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $FTIME seconds
echo -e "${outdir}\t${FTIME}" | tee -a "$timing"

docker container rm breseq

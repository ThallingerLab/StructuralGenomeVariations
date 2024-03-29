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

log_eval $PWD "docker run --rm --name=python3.8 -v $(pwd):$(pwd) -v $tools_dir:$tools_dir -w $outdir python:3.8-buster python3 $tools_dir/breseq/breseqToVCF.py \
  -i $outdir/data/output.gd -r $fasta -c LT962478.1:2263458 LT962479.1:1827941"

START=$(docker inspect --format='{{.State.StartedAt}}' breseq)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' breseq)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $FTIME seconds
echo -e "${outdir}\t${FTIME}" | tee -a "$timing"

docker container rm breseq

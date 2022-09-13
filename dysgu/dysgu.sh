#!/bin/bash
bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
outdir="$5"
threads="$6"
tools_dir="$7"
timing="$8"

log_eval $PWD "docker run --name=dysgu -v $(pwd)/:$(pwd) -w $outdir zeunas/dysgu:1.3.12 dysgu run \
  --mode pe \
  -o $outdir/dysgu.vcf \
  --max-cov -1 \
  --diploid False \
  $fasta \ #Reference
  $outdir/tmp \ #Working directory
  $bam"

START=$(docker inspect --format='{{.State.StartedAt}}' dysgu)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' dysgu)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

docker container rm dysgu

echo final time: $FTIME seconds 2>&1
echo "${outdir}\t$FTIME}" | tee -a "$timing"

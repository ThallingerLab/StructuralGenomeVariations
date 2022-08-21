#!/bin/bash

threads=8

bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
outdir="$5"
threads="$6"

log_eval $PWD "docker run --name=gridss -v $(pwd):$(pwd) -w $outdir gridss/gridss:2.13.2 gridss \
 --reference $fasta \
 --output ${outdir}/svs.vcf \
 --assembly ${outdir}/assembly.bam \
 --threads $threads \
 $bam"

START=$(docker inspect --format='{{.State.StartedAt}}' gridss)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' gridss)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

docker container rm gridss

echo final time: $FTIME seconds 2>&1 | tee "$outdir"/gridss_runtime.txt

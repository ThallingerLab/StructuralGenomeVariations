#!/bin/bash

bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
outdir="$5"

cd ..
mkdir "$output"/gridss/gridss_out_"$stamp"
docker run --name=gridss -v $(pwd)/:/in/ -w /in/ gridss/gridss:2.13.2 OUTPUT="${outdir}/svs.vcf" INPUT="$bam" REFERENCE_SEQUENCE="$fasta" ASSEMBLY="${outdir}/assembly.bam"

START=$(docker inspect --format='{{.State.StartedAt}}' gridss)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' gridss)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $FTIME seconds 2>&1 | tee "$output"/gridss/gridss_runtime_"$stamp".txt

docker container rm gridss

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

IS=$(docker run --rm -v $(pwd):$(pwd) -w $(pwd) staphb/samtools:1.15 samtools stats $bam | grep "insert size average" | sed 's/[^0-9.]*//g')

config="${outdir}/pindel.config"

base=$(basename $bam)
base=${base/.bam/}

log_eval $PWD "echo $bam $IS $base > $config"

log_eval $PWD "docker run --name=pindel1 -v $(pwd):$(pwd) -w $outdir shuangbroad/pindel:v0.2.5b8 pindel \
-f $fasta \
-i $config \
-c ALL \
-T $threads \
-o ${base}"

START=$(docker inspect --format='{{.State.StartedAt}}' pindel1)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' pindel1)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME1=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

log_eval $PWD "docker run --name=pindel2 -v $(pwd):$(pwd) -w $outdir shuangbroad/pindel:v0.2.5b8 pindel2vcf \
-P ${base} \
-d 20220505 \
-r $fasta \
-R CBS7435 \
-co"

START=$(docker inspect --format='{{.State.StartedAt}}' pindel2)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' pindel2)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME2=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

FTIME=$(($FTIME1+$FTIME2))

echo final time: $FTIME seconds 2>&1
echo -e "${outdir}\t${FTIME}" | tee -a "$timing"

docker container rm pindel1 pindel2

#!/bin/bash
#Usage: bdmax.sh <bamfile> <output directory> <cutoff value> <calculate_lower>
#example : ./bdmax.sh lumpy_BSY_red.bam output_bdmax 60 3
#cutoff  = 60
#calculate_lower should be 3 for BSY10dKU70 and 2 for BSY11dKU70 (only the full version works with breakdancer)

threads=8

bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
outdir="$5"
threads="$6"
tools_dir="$7"
timing="$8"

log_eval $PWD "docker run --name=bdmax1 -v $(pwd):$(pwd) -w $outdir vrohnie/bdmax:v1.4.5 perl \
  /root/breakdancer-1.4.5/perl/bam2cfg.pl \
  $bam > $outdir/bdmaxconfig"

START=$(docker inspect --format='{{.State.StartedAt}}' bdmax1)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' bdmax1)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

TIME1=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo first time: $TIME1 seconds

docker container rm bdmax1

log_eval $PWD "docker run --name=bdmax2 -v $(pwd):$(pwd) -w $outdir vrohnie/bdmax:v1.4.5 breakdancer-max \
  $outdir/bdmaxconfig > $outdir/bdmax.tsv"

log_eval $PWD "docker run --rm --name=python3.8 -v $(pwd):$(pwd) -v $tools_dir:$tools_dir -w $outdir python:3.8-buster python3 $tools_dir/bdmax/bdmaxToVCF.py \
  -i $outdir/bdmax.tsv -r $fasta -c LT962478.1:2263458 LT962479.1:1827941"

START=$(docker inspect --format='{{.State.StartedAt}}' bdmax2)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' bdmax2)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

TIME2=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo second time: $TIME2 seconds
FTIME=$(($TIME1+$TIME2))

echo final time: $FTIME seconds
echo -e "${outdir}\t${FTIME}" | tee -a "$timing"

docker container rm bdmax2



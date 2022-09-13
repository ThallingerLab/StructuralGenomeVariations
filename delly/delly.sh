#!/bin/bash
bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
outdir="$5"
threads="$6"
tools_dir="$7"
timing="$8"

log_eval $PWD "docker run --name=delly -v $(pwd)/:$(pwd) -w $outdir dellytools/delly delly call \
  -g $fasta \
  $bam > $outdir/delly.vcf"

START=$(docker inspect --format='{{.State.StartedAt}}' delly)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' delly)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo first time: "$FTIME"

docker container rm delly

#docker run --name=bcftools -v $(pwd)/delly/output/:/in/ -w /in/ halllab/bcftools:v1.9 bcftools view -i'FILTER="PASS" && MAPQ>40 && (PE>5 || SR>5)' delly_output_"$stamp".vcf > ./delly/output/delly_output_BCFTOOLS_"$stamp".vcf.vcf

#START=$(docker inspect --format='{{.State.StartedAt}}' bcftools)
#STOP=$(docker inspect --format='{{.State.FinishedAt}}' bcftools)

#START_TIMESTAMP=$(date --date=$START +%s)
#STOP_TIMESTAMP=$(date --date=$STOP +%s)

#TIME2=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

#echo seocnd time: "$TIME2"

#docker container rm bcftools

#FTIME=$(($TIME1+$TIME2))

echo final time: $TIME1 seconds 2>&1
echo "${outdir}\t$FTIME}" | tee -a "$timing"

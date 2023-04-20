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

CLOVE_VCF="$outdir/clove.vcf"
GRIDSS_VCF="${outdir/clove/gridss}/svs.vcf"
DELLY_VCF="${outdir/clove/delly}/delly.vcf"

CLOVE="/app/custom_scripts/clovebiotech_v1.0.1.jar"

if [ -s "$GRIDSS_VCF" ] && [ -s "$DELLY_VCF" ]; then
  log_eval $PWD "docker run --name=clove -v $(pwd):$(pwd) -w $outdir vrohnie/micronap:v1.2.1 java -jar $CLOVE \
   -i $GRIDSS_VCF GRIDSS \
   -i $DELLY_VCF DELLY2 \
   -b $bam \
   -o ${CLOVE_VCF}.temp \
   -r 25000"

   log_eval $PWD "cat $CLOVE_VCF.temp | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > $CLOVE_VCF"
fi

START=$(docker inspect --format='{{.State.StartedAt}}' gridss)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' gridss)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)


FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

docker container rm clove

echo final time: $FTIME seconds 2>&1
echo -e "${outdir}\t${FTIME}" | tee -a "$timing"


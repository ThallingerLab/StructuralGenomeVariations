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
GRIDSS_VCF="${outdir/clovebiotech/gridss}/svs.vcf"
DELLY_VCF="${outdir/clovebiotech/delly}/delly.vcf"

CLOVE="/bin/clove-0.17-jar-with-dependencies.jar"

if [ -s "$GRIDSS_VCF" ] && [ -s "$DELLY_VCF" ]; then
  log_eval $PWD "docker run --name=clove -v $(pwd):$(pwd) -w $outdir vrohnie/clove:0.17 java -jar $CLOVE \
   -i $GRIDSS_VCF GRIDSS \
   -i $DELLY_VCF DELLY2 \
   -b $bam \
   -o ${CLOVE_VCF}.temp"

  START=$(docker inspect --format='{{.State.StartedAt}}' clove)
  STOP=$(docker inspect --format='{{.State.FinishedAt}}' clove)

  START_TIMESTAMP=$(date --date=$START +%s)
  STOP_TIMESTAMP=$(date --date=$STOP +%s)


  FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

  docker container rm clove

  cat $CLOVE_VCF.temp | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1V -k2,2n"}' > $CLOVE_VCF

  echo final time: $FTIME seconds 2>&1
  echo -e "${outdir}\t${FTIME}" | tee -a "$timing"

fi


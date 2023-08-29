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
GRIDSS_Filtered_VCF="${outdir/clove/gridss}/svs_filtered.vcf"
DELLY_VCF="${outdir/clove/delly_0.7.9}/delly_0.7.9.vcf"

#CLOVE="/bin/clove-0.15-jar-with-dependencies.jar"
CLOVE="/bin/clove.jar"

#if [ -s "$GRIDSS_VCF" ] && [ -s "$DELLY_VCF" ]; then
if [ -s "$GRIDSS_VCF" ]; then
  log_eval $PWD "egrep '^#|\[|\]' $GRIDSS_VCF > $GRIDSS_Filtered_VCF"

  # shellcheck disable=SC1065
  dp=$(basename $outdir)
  dp=${dp/clove_/}

  dp_std=$(echo "scale=2 ; $dp / 5" | bc)

  log_eval $PWD "docker run --user 0:0 --name=clove -v $(pwd):$(pwd) -w $outdir vrohnie/clove:lastes java -jar $CLOVE \
   -i $GRIDSS_Filtered_VCF GRIDSS \
   -b $bam \
   -o ${CLOVE_VCF}.temp \
   -c $dp $dp_std"
#   -i $DELLY_VCF DELLY2 \

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


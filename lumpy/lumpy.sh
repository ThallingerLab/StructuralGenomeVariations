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

#docker run --name=lumpy -v $(pwd)/:/in/ -w /in/ szarate/lumpy-sv:v0.3.0 lumpyexpress -B "$1".bam -S "$1".splitters.bam -D "$1".discordants.bam -o "$output"/lumpy/lumpyoutput_"$stamp".vcf

## It is advised to call lumpy via smoove
log_eval $PWD "docker run --name=lumpy -v $(pwd):$(pwd) -w $outdir brentp/smoove:v0.2.7 smoove call \
--fasta $fasta \
--name $(basename $outdir) \
--genotype \
--outdir $outdir \
$bam"

#./lumpy/lumpy_awk_filter.sh "$output"/lumpy/lumpyoutput_"$stamp".vcf 100 "$output"/lumpy/lumpyoutput_CUTOFF_"$stamp".vcf

#rm ./lumpy/output/lumpyoutput_"$stamp".vcf

START=$(docker inspect --format='{{.State.StartedAt}}' lumpy)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' lumpy)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $FTIME seconds 2>&1
echo -e "${outdir}\t${FTIME}" | tee -a "$timing"

docker container rm lumpy

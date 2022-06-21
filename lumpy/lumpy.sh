#!/bin/bash 
stamp="$5"
output="$6"
cd ..

docker pull szarate/lumpy-sv:v0.3.0

docker run --name=lumpy -v $(pwd)/:/in/ -w /in/ szarate/lumpy-sv:v0.3.0 lumpyexpress -B "$1".bam -S "$1".splitters.bam -D "$1".discordants.bam -o "$output"/lumpy/lumpyoutput_"$stamp".vcf

## It is advised to call lumpy via smoove
docker run -it -v $(pwd):/in -w /in  brentp/smoove:v0.2.7 smoove call --fasta CBS7435_3-4.fasta --name CHR-1 --genotype --outdir bam/HSXt_l150m550s165f100_EF bam/HSXt_l150m550s165f100_EF/CHR-1_sorted.bam

./lumpy/lumpy_awk_filter.sh "$output"/lumpy/lumpyoutput_"$stamp".vcf 100 "$output"/lumpy/lumpyoutput_CUTOFF_"$stamp".vcf

#rm ./lumpy/output/lumpyoutput_"$stamp".vcf

START=$(docker inspect --format='{{.State.StartedAt}}' lumpy)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' lumpy)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $FTIME seconds 2>&1 | tee "$output"/lumpy/lumpy_runtime_"$stamp".txt

docker container rm lumpy

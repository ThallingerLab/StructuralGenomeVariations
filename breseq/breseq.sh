#!/bin/bash
bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
stamp="$5"
output="$6"
cd ..

# Think about using bwa BAM files, running bowtie on all files will take forever, especially with so many settings ?!

docker run --name=breseq -v $(pwd)/:/in/ -w /in/ pvstodghill/breseq:0.35.7__2021-08-03 breseq -o "$output"/breseq/breseq_output_"$stamp" -r "$fasta".fasta "$fastq1".fastq "$fastq2".fastq

START=$(docker inspect --format='{{.State.StartedAt}}' breseq)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' breseq)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $FTIME seconds 2>&1 | tee "$output"/breseq_runtime_"$stamp".txt

docker container rm breseq

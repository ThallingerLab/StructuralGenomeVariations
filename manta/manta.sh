# Usage: manta.sh <bamfile> <reference fasta> <output directory> <number of available cores>
#

threads=8

bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
outdir="$5"
threads="$6"
tools_dir="$7"
timing="$8"

log_eval $PWD "docker run --name=manta1 -v $(pwd):$(pwd) -w $outdir kfdrc/manta:1.6.0 /manta-1.6.0.centos6_x86_64/bin/configManta.py \
  --bam=$bam \
  --referenceFasta=$fasta \
  --runDir=$outdir"

START=$(docker inspect --format='{{.State.StartedAt}}' manta1)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' manta1)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

TIME1=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo first time: $TIME1 seconds

docker container rm manta1


log_eval $PWD "docker run --name=manta2 -v $(pwd):$(pwd) -w $(pwd) kfdrc/manta:1.6.0 python $outdir/runWorkflow.py \
  -m local -j $threads"

START=$(docker inspect --format='{{.State.StartedAt}}' manta2)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' manta2)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

TIME2=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo second time: $TIME2 seconds

docker container rm manta2

FTIME=$(($TIME1+$TIME2))

echo final time: $FTIME seconds 2>&1
echo "${outdir}\t$FTIME}" | tee -a "$timing"

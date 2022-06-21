# Usage: manta.sh <bamfile> <reference fasta> <output directory> <number of available cores>
#
bam="$1"
fasta="$2"
fastq1="$3"
fastq2="$4"
cores=1
stamp="$5"
output="$6"

cd ..
docker run --name=manta1 -v $(pwd):/in/ -w /in/ kfdrc/manta:1.6.0  ../manta-1.6.0.centos6_x86_64/bin/configManta.py --bam="$bam".bam --referenceFasta="$fasta".fasta --runDir="$output"/manta/output_"$stamp"


START=$(docker inspect --format='{{.State.StartedAt}}' manta1)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' manta1)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

TIME1=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo first time: $TIME1 seconds

docker container rm manta1


docker run --name=manta2 -v $(pwd):/in/ -w /in/ kfdrc/manta python "$output"/manta/output_"$stamp"/runWorkflow.py -m local -j $cores

START=$(docker inspect --format='{{.State.StartedAt}}' manta2)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' manta2)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

TIME2=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo second time: $TIME2 seconds

docker container rm manta2

FTIME=$(($TIME1+$TIME2))

echo final time: $FTIME seconds 2>&1 | tee "$output"/manta/manta_runtime_"$stamp".txt

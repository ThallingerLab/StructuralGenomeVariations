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
docker run --name=metasv -v $(pwd):/in/ -w /in/ zhouanbo/metasv:1.0  ../manta-1.6.0.centos6_x86_64/bin/configManta.py --bam="$bam".bam --referenceFasta="$fasta".fasta --runDir="$output"/manta/output_"$stamp"

START=$(docker inspect --format='{{.State.StartedAt}}' metasv)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' metasv)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

FTIME=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo final time: $FTIME seconds 2>&1 | tee "$output"/metasv/metasv_runtime_"$stamp".txt

docker container rm metasv

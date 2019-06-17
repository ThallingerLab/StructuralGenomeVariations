#!/bin/bash
#Usage: bdmax.sh <bamfile> <output directory> <cutoff value> <calculate_lower>
#example : ./bdmax.sh lumpy_BSY_red.bam output_bdmax 60 3
#cutoff  = 60
#calculate_lower should be 3 for BSY10dKU70 and 2 for BSY11dKU70 (only the full version works with breakdancer)
bamfile="$1"
cutoff=41
stamp="$5"
output="$6"
conf="bdmaxconfig"

cd ..

if [[ $1 == *"BSY11"* ]]
then
	c=2
else
	c=3
fi

if [[ $1 == *"BSY11_red"* ]]
then
    conf="BSY11_red_config"
else
    docker run --name=bdmax1 -v $(pwd):/in/ -w /in/ molecular/breakdancer perl ../home/bio/breakdancer-1.4.5/perl/bam2cfg.pl -c "$c" "$bamfile".bam > ./bdmax/bdmaxconfig
fi


START=$(docker inspect --format='{{.State.StartedAt}}' bdmax1)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' bdmax1)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

TIME1=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo first time: $TIME1 seconds

docker container rm bdmax1

echo $conf

docker run --name=bdmax2 -v $(pwd):/in/ -w /in/ molecular/breakdancer breakdancer-max -y $cutoff ./bdmax/"$conf" > "$output"/bdmax/bdmax_"$stamp".vcf

START=$(docker inspect --format='{{.State.StartedAt}}' bdmax2)
STOP=$(docker inspect --format='{{.State.FinishedAt}}' bdmax2)

START_TIMESTAMP=$(date --date=$START +%s)
STOP_TIMESTAMP=$(date --date=$STOP +%s)

TIME2=$(($STOP_TIMESTAMP-$START_TIMESTAMP))

echo second time: $TIME2 secondssss
FTIME=$(($TIME1+$TIME2))

echo final time: $FTIME seconds 2>&1 | tee "$output"/bdmax/bdmax_runtime_"$stamp".txt

docker container rm bdmax2



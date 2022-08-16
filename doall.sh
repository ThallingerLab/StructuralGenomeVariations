#!/bin/bash
#tests all tools with all alignments we had deemed necessary 
#usage: ./doall.sh <output name>
#example: ./doall.sh newoutput
#						creates a directory called "newoutput" that contains the results of all four alignments for each tool
#if no outputn name is provided, a timestamp is used as a default

echo "calling tools with predeterminded parameters"

declare -a params=("BSY10_red"
		   "BSY11_red"
                   "BSY10_full"
                   "BSY11_full")

stamp="$(date +'%Y_%m_%d-%H_%M_%S')"
oname="$1"

if [ -z "$1" ]
then
    echo "No name supplied, setting time stamp as default"
    oname="$stamp"
fi

mkdir output_"$oname"
mkdir output_"$oname"/BSY11
mkdir output_"$oname"/BSY11/red
mkdir output_"$oname"/BSY11/full

mkdir output_"$oname"/BSY10
mkdir output_"$oname"/BSY10/red
mkdir output_"$oname"/BSY10/full

declare -a tools=("bdmax" "delly" "gridss" "lumpy" "manta" "softsv" "breseq" "metasv")
for tool in "${tools[@]}"
do
	mkdir output_"$oname"/BSY11/red/"$tool"
	mkdir output_"$oname"/BSY11/full/"$tool"
	mkdir output_"$oname"/BSY10/red/"$tool"
	mkdir output_"$oname"/BSY10/full/"$tool"
done

ref="reference/cbs_7435_full"

for p in "${params[@]}"
do
	echo "calling with $p"

	if [[ $p == *"BSY11"* ]]
	then
		if [[ $p == *"_red"* ]]; then
			output=output_"$oname"/BSY11/red
			bam="BSY11/BSY11_red"
			fastq1="BSY11/BSY11_1_trimmed_red"
			fastq2="BSY11/BSY11_2_trimmed_red"
		else
			output=output_"$oname"/BSY11/full
			bam="BSY11/BSY11_full"
			fastq1="BSY11/BSY11_1_trimmed"
			fastq2="BSY11/BSY11_2_trimmed"
		fi
	else
		if [[ $p == *"_red"* ]]; then
			output=output_"$oname"/BSY10/red
			bam="BSY10/BSY10_red"
			fastq1="BSY10/BSY10_1_trimmed_red"
			fastq2="BSY10/BSY10_2_trimmed_red"
		else
			output=output_"$oname"/BSY10/full
			bam="BSY10/BSY10_full"
			fastq1="BSY10/BSY10_1_trimmed"
			fastq2="BSY10/BSY10_2_trimmed"
		fi
	fi

    	bash ./doit.sh "$bam" "$ref" "$fastq1" "$fastq2" "$oname" "$output" 2>&1
done

echo "doall done, fingers crossed"

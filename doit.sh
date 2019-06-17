#!/bin/bash
#tests all tools usings a single alignment
#usage: ./doit.sh <bamfile name> <reference file name> <fastq1 file name> <fastq2 file name> <output name> <output folder>
# all file names are provided WITHOUT FILE EXTENSION
#example: ./doit.sh BSY11_trimmed_red cbs_7435_full 100549_1_trimmed 100549_2_trimmed test1 testoutput
echo "testing all tools..."

docker kill $(docker ps -q)
docker rm $(docker ps -a -q)

declare -a tools=("bdmax" "delly" "gridss" "lumpy" "manta" "softsv" "breseq")
#("bdmax" "breseq" "delly" "gridss" "lumpy" "manta" "metasv" "softsv")

dir="$(pwd)"

for tool in "${tools[@]}"
do
	stamp="$(date +'%Y_%d_%m-%H_%M_%S')"
	echo "╔══════════════════════════════════════════════════════════════╗"
	echo "║                        starting $tool                        ║"
	echo "╚══════════════════════════════════════════════════════════════╝"
    	cd "$dir"/"$tool"
        ./$tool.sh "$1" "$2" "$3" "$4" "$5" "$6" 2>&1 | tee ../"$6"/"$tool"/"$tool"_out_"$5".txt
	echo "#####################   DONE WITH $tool   ######################"
done

echo "That's all, folks!"

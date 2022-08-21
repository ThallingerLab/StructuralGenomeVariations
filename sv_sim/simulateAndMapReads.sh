#!/bin/bash

OPTSTRING="hr:s:i:o:t:c:"

usage()
{
	echo -e  "You did it wrong"
}

declare SWITCH
threads=8

# Examine individual options
while getopts "$OPTSTRING" SWITCH; do
	case $SWITCH in

		r) ref="$OPTARG"
		ref=$(readlink -e "$ref")
		echo "Reference = $ref"
		;;

		s) settings="$OPTARG"
		echo "Settings = $settings"
		;;

    i) base_dir="$OPTARG"
    base_dir=$(readlink -e "$base_dir")
		echo "Basedir = $base_dir"
		;;

    o) out_dir="$OPTARG"
    out_dir=$(readlink -e "$out_dir")
		echo "Outdir = $out_dir"
		;;

    t) tools_dir="$OPTARG"
    tools_dir=$(readlink -e "$tools_dir")
		echo "Tools directory = $tools_dir"
		;;

    c) threads="$OPTARG"
		echo "Threads = $threads"
		;;

		*) echo "script error: unhandled argument"
		usage
		exit 1
		;;


	esac
done

source $tools_dir/log_eval.sh

#ART="docker run -u 1001:1001 --rm -v $PWD:$PWD -w $PWD vlr37/art_illumina art_illumina"
#ART="art_illumina"
ART="docker run -u 1001:1001 --name art --rm -v $PWD:$PWD -w $PWD vrohnie/art:v2.5.8 /home/art_bin_MountRainier/art_illumina"
SAMBAMBA="docker run -u 1001:1001 --name sambamba --rm -v $PWD:$PWD -w $PWD clinicalgenomics/sambamba:0.8.0"
BWA="docker run -u 1001:1001 --name bwa --rm -v $PWD:$PWD -w $PWD mskcc/bwa_mem:0.7.12 bwa"

#declare -a scripts=("bdmax" "delly" "gridss" "lumpy" "manta" "softsv" "breseq")

declare -a tools=("gridss" "manta" "lumpy" "delly" "bdmax" "softsv" "breseq")
declare -a fractionOfReads=(25 50 75)
seed=87

if [ ! -d $out_dir/bam ]; then
	mkdir $out_dir/bam
fi

if [ ! -d $out_dir/fastq ]; then
	mkdir $out_dir/fastq
fi

if [ ! -d $out_dir/svs ]; then
	mkdir $out_dir/svs
fi


stamp="$(date +'%Y_%d_%m-%H_%M_%S')"
org="CBS7435"

grep -v '^#' $settings | while IFS=$'\t' read -r -a settings_array
do
  settings_string="${settings_array[0]}_f${settings_array[1]}_l${settings_array[2]}_m${settings_array[3]}_s${settings_array[4]}"
  bamdir=$out_dir/bam/$settings_string
  fastqdir=$out_dir/fastq/$settings_string
  svsdir=$out_dir/svs/$settings_string

  echo "THIS IS the DIRECTORY: $bamdir"

  if [ ! -d $bamdir ]; then
    mkdir $bamdir
  fi

  if [ ! -d $fastqdir ]; then
    mkdir $fastqdir
  fi

  if [ ! -d $svsdir ]; then
    mkdir $svsdir
  fi

  for dir in "$base_dir"/*; do

    if [[ -d $dir ]]; then
      base=$(basename "$dir")

      READ1_FILE="${fastqdir}/${base}_1.fq.gz"
      READ2_FILE="${fastqdir}/${base}_2.fq.gz"
      SAM="${bamdir}/${base}.sam"
      BAM_SORTED="${bamdir}/${base}.sorted.bam"

      echo "╔══════════════════════════════════════════════════════════════╗"
      echo "║                   running read generation                    ║"
      echo "╚══════════════════════════════════════════════════════════════╝"

      if [ ! -s "$READ1_FILE" ]; then

        log_eval $PWD "$ART -ss ${settings_array[0]} -i ${dir}/h1.fa -p -na -f ${settings_array[1]} -l ${settings_array[2]} -m ${settings_array[3]} -s ${settings_array[4]} -o ${fastqdir}/${base}_"

        if [ -s "${fastqdir}/${base}_1.fq" -a -s "${fastqdir}/${base}_2.fq" ]; then
          gzip "${fastqdir}/${base}_1.fq" "${fastqdir}/${base}_2.fq"
        fi
      fi

      echo "╔══════════════════════════════════════════════════════════════╗"
      echo "║                      running read mapping                    ║"
      echo "╚══════════════════════════════════════════════════════════════╝"

      if [ ! -s "$BAM_SORTED" ]; then

        log_eval $PWD  "$BWA mem -a -M -R '@RG\tID:${org}\tSM:${base/-/}\tPL:ILLUMINA\tLB:lib' $ref $READ1_FILE $READ2_FILE > $SAM"
        log_eval $PWD  "$SAMBAMBA sambamba view -h -t $threads -f bam -S -o ${SAM/.sam/.bam} $SAM"
        log_eval $PWD  "$SAMBAMBA sambamba sort -t $threads -o $BAM_SORTED ${SAM/.sam/.bam}"
        log_eval $PWD  "$SAMBAMBA sambamba index -t $threads $BAM_SORTED"

        rm $SAM ${SAM/.sam/.bam}
      fi
    fi
  done
done

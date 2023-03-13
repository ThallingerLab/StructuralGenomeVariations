#!/bin/bash

OPTSTRING="hr:s:i:o:t:c:b:"

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

    i) fasta_dir="$OPTARG"
    fasta_dir=$(readlink -e "$fasta_dir")
		echo "Fasta Directory = $fasta_dir"
		;;

#    i) base_dir="$OPTARG"
#    base_dir=$(readlink -e "$base_dir")
#		echo "Basedir = $base_dir"
#		;;

    o) out_dir="$OPTARG"
    out_dir=$(readlink -e "$out_dir")
		echo "Outdir = $out_dir"
		;;

    b) bed_dir="$OPTARG"
    bed_dir=$(readlink -e "bed_dir")
		echo "Bed directory = bed_dir"
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

SAMBAMBA="docker run -u 1001:1001 --name sambamba --rm -v $PWD:$PWD -w $PWD clinicalgenomics/sambamba:0.8.0"
BWA="docker run -u 1001:1001 --name bwa --rm -v $PWD:$PWD -w $PWD mskcc/bwa_mem:0.7.12 bwa"

declare -a tools=("gridss" "manta" "lumpy" "delly" "bdmax" "softsv" "breseq")
declare -a fractionOfReads=(25 50 75)
seed=87

base_in=$(basename "$fasta_dir")

if [ ! -d "$out_dir/${base_in/fasta/bam_plasmids}" ]; then
	mkdir "$out_dir/${base_in/fasta/bam_plasmids}"
fi

if [ ! -d "$out_dir/${base_in/fasta/svs_plasmids}" ]; then
	mkdir "$out_dir/${base_in/fasta/svs_plasmids}"
fi

stamp="$(date +'%Y_%d_%m-%H_%M_%S')"
org="CBS7435"

grep -v '^#' $settings | while IFS=$'\t' read -r -a settings_array
do
#  settings_string="${settings_array[0]}_f${settings_array[1]}_l${settings_array[2]}_m${settings_array[3]}_s${settings_array[4]}"
#  bamdir=$out_dir/bam_plasmids/$settings_string
#  fastqdir=$out_dir/fastq/$settings_string
#  svsdir=$out_dir/svs_plasmids/$settings_string

  settings_string="${settings_array[0]}_f${settings_array[1]}_l${settings_array[2]}_m${settings_array[3]}_s${settings_array[4]}"
  bamdir=$out_dir/${base_in/fasta/bam_plasmids}/$settings_string
  fastqdir=$out_dir/${base_in/fasta/fastq}/$settings_string
  svsdir=$out_dir/${base_in/fasta/svs_plasmids}/$settings_string

  echo "THIS IS the DIRECTORY: $bamdir"

  if [ ! -d $bamdir ]; then
    mkdir $bamdir
  fi

  for dir in "$fasta_dir"/*; do

    if [[ -d $dir ]]; then
      base=$(basename "$dir")

      READ1_FILE="${fastqdir}/${base}_1.fq.gz"
      READ2_FILE="${fastqdir}/${base}_2.fq.gz"

      ref_plasmid=$out_dir/${base_in/fasta/bam_plasmids}/${base}_refPlasmid.fasta

      if [ ! -s "$ref_plasmid" ]; then

        PLASMID=$(awk 'NR==1{print $5}' ${bed_dir}/${base}/${base}_summary.bed)

        printf ">plasmid\n$PLASMID" | fold -80  > "${bed_dir}/${base}/${base}_plasmid.fasta"
        log_eval $PWD "cat $ref ${fastqdir}/${base}/${base}_plasmid.fasta > $ref_plasmid"

        log_eval $PWD  "$BWA index $ref_plasmid"
      fi

      ref="$ref_plasmid"

      SAM="${bamdir}/${base}.sam"
      BAM_SORTED="${bamdir}/${base}.sorted.bam"

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

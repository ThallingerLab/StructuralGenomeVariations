#!/bin/bash

OPTSTRING="hr:s:i:o:t:c:l:b:"

usage()
{
	echo -e "You did it wrong"
}

declare SWITCH
threads=8
bam_base="NONE"

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

    b) bam_base="$OPTARG"
    bam_base=$(readlink -e "$bam_base")
		echo "Bam base Directory = $bam_base"
		;;

    o) out_dir="$OPTARG"
    out_dir=$(readlink -e "$out_dir")
		echo "Outdir = $out_dir"
		;;

    t) tools_dir="$OPTARG"
    tools_dir=$(readlink -e "$tools_dir")
		echo "Tools directory = $tools_dir"
		;;

    l) tools_list="$OPTARG"
    declare -a tools=($tools_list)
		echo "Going to run the following tools = $tools"
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

#declare -a tools=("wham" "svaba" "pindel" "bdmax" "softsv" "manta" "lumpy" "delly" "gridss")
declare -a fractionOfReads=(100 75 50 25)

seed=87

stamp="$(date +'%Y_%d_%m-%H_%M_%S')"
org="CBS7435"

grep -v '^#' $settings | while IFS=$'\t' read -r -a settings_array
do

  base_in=$(basename "$fasta_dir")

  settings_string="${settings_array[0]}_f${settings_array[1]}_l${settings_array[2]}_m${settings_array[3]}_s${settings_array[4]}"

  if [ $bam_base = "NONE" ]; then
    bamdir=$out_dir/${base_in/fasta/bam}/$settings_string
  else
    bamdir=$bam_base/$settings_string
  fi

  fastqdir=$out_dir/${base_in/fasta/fastq}/$settings_string
  svsdir=$out_dir/${base_in/fasta/svs}/$settings_string

  timing=${svsdir}/${stamp}_timing.tsv
  log=${svsdir}/${stamp}_sv_calling.log

  echo "THIS IS the DIRECTORY: $bamdir"

  if [ ! -d $svsdir ]; then
    mkdir $svsdir
  fi
  touch $timing
  touch $log

  for dir in "$fasta_dir"/*; do

    if [[ -d $dir ]]; then
      base=$(basename "$dir")

      if [ ! -d $svsdir/$base ]; then
        mkdir $svsdir/$base
      fi

      READ1_FILE="${fastqdir}/${base}_1.fq.gz"
      READ2_FILE="${fastqdir}/${base}_2.fq.gz"
      SAM="${bamdir}/${base}.sam"
      BAM_SORTED="${bamdir}/${base}.sorted.bam"


      echo "╔══════════════════════════════════════════════════════════════╗"
      echo "║                      starting SV analysis                    ║"
      echo "╚══════════════════════════════════════════════════════════════╝"

      if [ -s "$BAM_SORTED" ]; then

        for fraction in "${fractionOfReads[@]}"; do

          BAM_FRACTION=$BAM_SORTED

          if [ $fraction -ne 100 ]; then
            BAM_FRACTION="${BAM_SORTED/.bam/-${fraction}.bam}"

            log_eval $PWD "$SAMBAMBA sambamba view -h -t $threads -s 0.$fraction -f bam \
              --subsampling-seed=$seed -o $BAM_FRACTION $BAM_SORTED"
          fi

          for tool in "${tools[@]}"; do
            tool_outdir="$svsdir/$base/${tool}_${fraction}"

            if [ ! -d "$tool_outdir" ]; then
              mkdir "$tool_outdir"

              export -f log_eval
              log_eval $PWD "$tools_dir/$tool/${tool}.sh $BAM_FRACTION $ref $READ1_FILE $READ2_FILE $tool_outdir $threads $tools_dir $timing" "$log"
            fi

          done

          if [ $fraction -ne 100 ]; then
            rm ${BAM_FRACTION}*
          fi
        done
      fi
    fi
  done
done

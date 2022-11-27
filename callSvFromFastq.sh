#!/bin/bash

OPTSTRING="hr:s:i:o:t:l:c:"

usage()
{
	echo -e "You did it wrong"
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

    i) fastq_dir="$OPTARG"
    fastq_dir=$(readlink -e "$fastq_dir")
		echo "Fastq Input Directory = $fastq_dir"
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

declare -a tools=("breseq")
declare -a fractionOfReads=(25 50 75 100)
seed=87

stamp="$(date +'%Y_%d_%m-%H_%M_%S')"
org="CBS7435"

grep -v '^#' $settings | while IFS=$'\t' read -r -a settings_array
do
  settings_string="${settings_array[0]}_f${settings_array[1]}_l${settings_array[2]}_m${settings_array[3]}_s${settings_array[4]}"
  svsdir=$out_dir/svs/$settings_string

  echo "THIS IS the DIRECTORY: $fastq_dir"

  if [ ! -d $svsdir ]; then
    mkdir $svsdir
  fi

  for dir in "$fastq_dir"/*; do

    if [[ -d $dir ]]; then
      base=$(basename "$dir")

      if [ ! -d $svsdir/$base ]; then
        mkdir $svsdir/$base
      fi

      READ1_FILE="${fastqdir}/${base}_1.fq.gz"
      READ2_FILE="${fastqdir}/${base}_2.fq.gz"

      echo "╔══════════════════════════════════════════════════════════════╗"
      echo "║                      starting SV analysis                    ║"
      echo "╚══════════════════════════════════════════════════════════════╝"

      if [ -s "$READ1_FILE" -a -s "$READ2_FILE" ]; then

#        for fraction in "${fractionOfReads[@]}"; do

#          READ1_FILE_FRACTION=$READ1_FILE
#          READ2_FILE_FRACTION=$READ2_FILE

#          if [ $fraction -ne 100 ]; then
#            READ1_FILE_FRACTION="${READ1_FILE/.fastq/-${fraction}.fastq}"
#            READ2_FILE_FRACTION="${READ2_FILE/.fastq/-${fraction}.fastq}"

#            log_eval $PWD "$SAMBAMBA sambamba view -h -t $threads -s 0.$fraction -f bam \
#              --subsampling-seed=$seed -o $BAM_FRACTION $BAM_SORTED"
#          fi

          for tool in "${tools[@]}"; do
            tool_outdir="$svsdir/$base/${tool}_${fraction}"

            if [ ! -d "$tool_outdir" ]; then
              mkdir "$tool_outdir"

              export -f log_eval
              log_eval $PWD "$tools_dir/$tool/${tool}.sh $BAM_FRACTION $ref $READ1_FILE $READ2_FILE $tool_outdir $threads $tools_dir $timing" "$log"
            fi

          done

#          if [ $fraction -ne 100 ]; then
#            rm $BAM_FRACTION
#          fi
#        done
      fi
    fi
  done
done

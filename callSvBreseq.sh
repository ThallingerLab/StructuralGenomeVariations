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

    i) fasta_dir="$OPTARG"
    fasta_dir=$(readlink -e "$fasta_dir")
		echo "Fasta Directory = $fasta_dir"
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

tool="breseq"
#declare -a fractionOfReads=(100 75 50 25)
#Calling the other fractions will not produce any junction calls this way
declare -a fractionOfReads=(100)
seed=87

stamp="$(date +'%Y_%d_%m-%H_%M_%S')"
org="CBS7435"

grep -v '^#' $settings | while IFS=$'\t' read -r -a settings_array
do

  base_in=$(basename "$fasta_dir")

  settings_string="${settings_array[0]}_f${settings_array[1]}_l${settings_array[2]}_m${settings_array[3]}_s${settings_array[4]}"
  bamdir=$out_dir/${base_in/fasta/bam}/$settings_string
  fastqdir=$out_dir/${base_in/fasta/fastq}/$settings_string
  svsdir=$out_dir/${base_in/fasta/svs}/$settings_string

  echo "THIS IS the DIRECTORY: $fastqdir"

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

      echo "╔══════════════════════════════════════════════════════════════╗"
      echo "║                      starting SV analysis                    ║"
      echo "╚══════════════════════════════════════════════════════════════╝"

      if [ -s "$READ1_FILE" -a -s "$READ2_FILE" ]; then

        for fraction in "${fractionOfReads[@]}"; do

#          READ1_FILE_FRACTION=$READ1_FILE
#          READ2_FILE_FRACTION=$READ2_FILE

          tool_outdir="$svsdir/$base/${tool}_${fraction}"

          BAM_SORTED="$svsdir/$base/${tool}_100/data/reference.bam"
          BAM_FRACTION="NONE"

          echo "Fraction is: $fraction"

          if [ "$fraction" -ne 100 ]; then
            BAM_FRACTION="${BAM_SORTED/.bam/-${fraction}.bam}"

            echo "Downsampling to: $BAM_FRACTION"

            log_eval $PWD "$SAMBAMBA sambamba view -h -t $threads -s 0.$fraction -f bam \
              --subsampling-seed=$seed -o $BAM_FRACTION $BAM_SORTED"
          fi

          # For the sake of consistency

          if [ ! -d "$tool_outdir" ]; then
            mkdir "$tool_outdir"
            if [ "$BAM_FRACTION" = "NONE" ]; then
              export -f log_eval
              log_eval $PWD "$tools_dir/$tool/${tool}.sh $BAM_FRACTION $ref $READ1_FILE $READ2_FILE $tool_outdir $threads $tools_dir $timing" "$log"
            elif [ -s "$BAM_FRACTION" ]; then
              export -f log_eval
              log_eval $PWD "$tools_dir/$tool/breseq_bam.sh $BAM_FRACTION $ref $READ1_FILE $READ2_FILE $tool_outdir $threads $tools_dir $timing" "$log"
            fi
          fi

#          if [ $fraction -ne 100 ]; then
#            rm $BAM_FRACTION
#          fi
        done
      else
        echo "Could not find files $READ1_FILE, $READ2_FILE"
      fi
    fi
  done
done

#!/bin/bash

OPTSTRING="hr:s:i:o:"

usage()
{
	echo -e  "You did it wrong"
}

declare SWITCH

# Examine individual options
while getopts "$OPTSTRING" SWITCH; do
	case $SWITCH in

		r) ref="$OPTARG"
		echo "Reference = $ref"
		;;

		s) settings="$OPTARG"
		echo "Settings = $settings"
		;;

    i) base_dir="$OPTARG"
		echo "Basedir = $base_dir"
		;;

    o) out_dir="$OPTARG"
		echo "Outdir = out_dir"
		;;

		*) echo "script error: unhandled argument"
		usage
		exit 1
		;;


	esac
done

source log_eval.sh

ART="docker run -u 1001:1001 --rm -v $PWD:$PWD -w $PWD vlr37/art_illumina art_illumina"
SAMBAMBA="docker run -u 1001:1001 --rm -v $PWD:$PWD -w $PWD clinicalgenomics/sambamba:0.8.0"
BWA="docker run -u 1001:1001 --rm -v $PWD:$PWD -w $PWD mskcc/bwa_mem:0.7.12 bwa"

fractionOfReads=(25 50 75)
seed=87

if [ ! -d $out_dir/bam ]; then
	mkdir $out_dir/bam
fi

if [ ! -d $out_dir/fastq ]; then
	mkdir $out_dir/fastq
fi


grep -v '^#' $settings | while IFS=$'\t' read -r -a settings_array

do
  bamdir=bam/${settings_array[0]}_f${settings_array[1]}_l${settings_array[2]}_m${settings_array[3]}_s${settings_array[4]}
  fastqdir=fastq/${settings_array[0]}_f${settings_array[1]}_l${settings_array[2]}_m${settings_array[3]}_s${settings_array[4]}

  echo "THIS IS the DIRECTORY:  cur_outdir"

  if [ ! -d $bamdir ]; then
    mkdir $bamdir
  fi

  if [ ! -d $fastqdir ]; then
    mkdir $fastqdir
  fi

  for dir in "$base_dir"/*; do

    if [[ -d $dir ]]; then
      base=$(basename "$dir")

      READ1_FILE="${fastqdir}/${base}_1.fq.gz"
      READ2_FILE="${fastqdir}/${base}_2.fq.gz"
      SAM="${bamdir}/${base}.sam"
      BAM_SORTED="${bamdir}/${base}.sorted.bam"

      if [ ! -s "$READ1_FILE" ]; then

        log_eval $PWD "$ART -ss ${settings_array[0]} -i ${dir}/h1.fa -p -na -f ${settings_array[1]} -l ${settings_array[2]} -m ${settings_array[3]} -s ${settings_array[4]} -o ${fastqdir}/${base}_"

        if [ -s "${fastqdir}/${base}_1.fq" -a -s "${fastqdir}/${base}_2.fq" ]; then
          gzip "${fastqdir}/${base}_1.fq" "${fastqdir}/${base}_2.fq"

        fi
      fi

      if [ ! -s "$BAM_SORTED" ]; then

        log_eval $PWD  "$BWA mem -a -M -R '@RG\tID:${org}\tSM:${seq}\tPL:ILLUMINA\tLB:lib' $ref $READ1_FILE $READ2_FILE > $SAM"
        log_eval $PWD  "$SAMBAMBA sambamba view -h -t 32 -f bam -S -o ${SAM/.sam/.bam} $SAM"
        log_eval $PWD  "$SAMBAMBA sambamba sort -t 32 -o $BAM_SORTED ${SAM/.sam/.bam}"
        log_eval $PWD  "$SAMBAMBA sambamba index -t 32 $BAM_SORTED"

        rm $SAM ${SAM/.sam/.bam}

        # BAM=${run}.dups_rem.bam
      fi

      if [ ! -s "${BAM_SORTED/.bam/-75.bam}" ]; then

        if [ -s "$BAM_SORTED" ]; then
            echo "yes"
             for fraction in ${fractionOfReads[@]}; do
                    log_eval $PWD "$SAMBAMBA sambamba view -h -t 32 -s 0.$fraction -f bam --subsampling-seed=$seed -o ${BAM_SORTED/.bam/-${fraction}.bam} $BAM_SORTED"
             done
        fi

      fi
    fi
  done
done

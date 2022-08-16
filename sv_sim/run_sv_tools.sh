#!/bin/bash

OPTSTRING="hr:s:b:f:o:"

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

    b) bam_dir="$OPTARG"
		echo "Basedir containing BAM files = $bam_dir"
		;;

    f) fastq_dir="$OPTARG"
		echo "Basedir containing FASTQ files  = $fastq_dir"
		;;

    o) out_dir="$OPTARG"
		echo "Outdir = $out_dir"
		;;

		*) echo "script error: unhandled argument"
		usage
		exit 1
		;;


	esac
done

source log_eval.sh

#!/bin/bash

OPTSTRING="hr:b:c:t:l:p"

usage()
{
	echo -e "You did it wrong"
}

declare SWITCH
plasmid="FALSE"

# Examine individual options
while getopts "$OPTSTRING" SWITCH; do
	case $SWITCH in

		b) basefolder="$OPTARG"
		basefolder=$(readlink -e $basefolder)

		echo "Basefolder = $basefolder"
		;;

		c) callfolder="$OPTARG"
		callfolder=$(readlink -e "$callfolder")

		echo "Callfolder = $callfolder"
		;;

    t) tools_dir="$OPTARG"
    tools_dir=$(readlink -e "$tools_dir")
		echo "Tools directory = $tools_dir"
		;;

    l) tool="$OPTARG"
		echo "Tool = $tool"
		;;

    p) plasmid="TRUE"
		;;

		*) echo "script error: unhandled argument"
		usage
		exit 1
		;;

	esac
done

source "$tools_dir"/log_eval.sh

if [ -s $basefolder -a -s $callfolder ]; then

  for base in  "$basefolder"/*/*summary.bed; do

    type=$(basename $base)
    type=${type/_summary.bed/}

    for call in "${callfolder}/${type}/${tool}"*/${tool}.vcf; do

      call_dir=$(dirname $call)

      if [ ! -d "$call_dir/all" ]; then
        mkdir "$call_dir/all"
      fi

      if [ ! -d "$call_dir/pass" ]; then
        mkdir "$call_dir/pass"
      fi

      if [ -s $call ]; then

        call_pass=${call/.vcf/_pass.vcf}

        if [ ! -s $call_pass ]; then
          log_eval $PWD "egrep '^#|PASS' $call > $call_pass"
        fi

        if [ $plasmid=="TRUE" ]; then
          log_eval $PWD "python3.8 ~/Projects/StructuralGenomeVariations/evaluation/vcfSherlock.py \
          -b $basefolder/${type}/${type}_summary.bed \
          -v $call \
          -f $call_dir/all -d 100 -s $basefolder/summary.tsv -p"

          log_eval $PWD "python3.8 ~/Projects/StructuralGenomeVariations/evaluation/vcfSherlock.py \
          -b $basefolder/${type}/${type}_summary.bed \
          -v $call_pass \
          -f $call_dir/pass -d 100 -s $basefolder/summary.tsv -p"
        else
          log_eval $PWD "python3.8 ~/Projects/StructuralGenomeVariations/evaluation/vcfSherlock.py \
          -b $basefolder/${type}/${type}_summary.bed \
          -v $call \
          -f $call_dir/all -d 100 -s $basefolder/summary.tsv"

          log_eval $PWD "python3.8 ~/Projects/StructuralGenomeVariations/evaluation/vcfSherlock.py \
          -b $basefolder/${type}/${type}_summary.bed \
          -v $call_pass \
          -f $call_dir/pass -d 100 -s $basefolder/summary.tsv"
        fi

      fi
    done
  done
fi

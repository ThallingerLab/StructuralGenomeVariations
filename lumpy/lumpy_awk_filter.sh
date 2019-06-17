#!/bin/bash
#
# Usage:
#  awk_filter.sh <filename> <cutoff value> <output filename>
#
# Examples:
#   lumpy_awk_filter.sh example.vcf 100 cutoff.vcf
#

if [ $# -lt 2 ]; then
  echo -e "Usage:\n  `basename ${0}` <filename> <cutoff value> <output filename>"
  exit 1
fi

fname=${1}
cutoff=${2}
output=${3}
odir="."

awk -v cutoff=${cutoff} -F ";" -e \
'{
   for (i=1;i<=NF;i++) {
       if ($i ~ /SU=/) {
		  tmp=match($i, /[0-9]/);
		  if ((substr(($i),tmp))+0 > cutoff) {
	         print($0);
		  }
	   }
   }
}' ${fname} > ${output}
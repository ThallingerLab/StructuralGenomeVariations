#!/bin/bash
#
# Usage:
#  softsv_awk_filter.sh <directoryname> <cutoff value>
#
# Examples:
#   softsv_awk_filter.sh softsv_directory 100
#

if [ $# -lt 2 ]; then
  echo -e "Usage:\n  `basename ${0}` <directoryname> <cutoff value>"
  exit 1
fi

dname=${1}
cutoff=${2}
odir="."

for file in ${dname}/*.txt ; do 
awk -v cutoff=${cutoff} -F " " -e \
'{
   if ($5+$6 > cutoff) {
      print $0;
   }
}' $file > $file'_cutoff.txt';
done
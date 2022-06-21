#!/bin/bash -x

#Shell script for mapping procedure
#Lea Demelius - 20.2.2019
#Veronika Schusterbauer -08.06.2021

echo "Starting mapping procedure..."

fasta="$1"

fastq1="$2"
fastq2="$3"

org=$4
seq=${fastq1/_*/}

docker pull mskcc/bwa_mem:0.7.12
docker pull staphb/samtools:1.15

#docker pull kfdrc/speedseq
#docker pull rpseq/samblaster
#docker pull erictdawson/lumpy-sv

if [ ! -s "${fasta}.fai" ]; then
  docker run --rm -v $(pwd):/in/ -w /in/ mskcc/bwa_mem:0.7.12 bwa index $fasta
fi

docker run --rm -v $(pwd):/in/ -w /in/ mskcc/bwa_mem:0.7.12 bwa mem -a -M -R '@RG\tID:${org}\tSM:${seq}\tPL:ILLUMINA' $fasta $fastq1 $fastq2 > ${seq}.sam
docker run --rm -v $(pwd):/in/ -w /in/ staphb/samtools:1.15 samtools view -h -b -S ${seq}.sam | samtools sort -o ${seq}_sorted.bam
docker run --rm -v $(pwd):/in/ -w /in/ staphb/samtools:1.15 samtools index ${seq}_sorted.bam"


#docker run --rm -v $(pwd):/in/ -w /in/ rpseq/samblaster -i x.sam --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 > y.sam
#docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools view -S -b y.sam > z.unsorted.bam
#docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools sort z.unsorted.bam -o $name.bam
#docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools view -b -F 1294 $name.bam > z.discordants.unsorted.bam
#docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools view -h $name.bam > z.splitters.before.bam
#docker run --rm -v $(pwd):/in/ -w /in/ erictdawson/lumpy-sv ../app/lumpy-sv/scripts/extractSplitReads_BwaMem -i z.splitters.before.bam > z.splitters.nearlythere.sam
#docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools view -Sb z.splitters.nearlythere.sam > z.splitters.unsorted.bam
#docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools sort z.discordants.unsorted.bam -o $name.discordants.bam
#docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools sort z.splitters.unsorted.bam -o $name.splitters.bam
#docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools index $name.bam
#docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools index $name.discordants.bam
#docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools index $name.splitters.bam
#rm z.splitters.unsorted.bam
#rm z.discordants.unsorted.bam
#rm z.splitters.nearlythere.sam
#rm z.splitters.before.bam
#rm z.unsorted.bam
#rm y.sam
#rm x.sam

echo "End of mapping procedure"

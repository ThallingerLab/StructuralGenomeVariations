#!/bin/bash -x
#Shell script for mapping procedure
#Lea Demelius - 20.2.2019
echo "Starting mapping procedure..."
fasta="$1"
fastq1="$2"
fastq2="$3"
name="$4" 
docker pull lindsayliang/bwa_mem 
docker pull kfdrc/speedseq
docker pull rpseq/samblaster
docker pull jweinstk/samtools
docker pull erictdawson/lumpy-sv
docker run --rm -v $(pwd):/in/ -w /in/ lindsayliang/bwa_mem bwa index $fasta
docker run --rm -v $(pwd):/in/ -w /in/ lindsayliang/bwa_mem bwa mem -R "@RG\tID:id\tSM:100549\tLB:lib" $fasta $fastq1 $fastq2 > x.sam
docker run --rm -v $(pwd):/in/ -w /in/ rpseq/samblaster -i x.sam --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 > y.sam
docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools view -S -b y.sam > z.unsorted.bam
docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools sort z.unsorted.bam -o $name.bam
docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools view -b -F 1294 $name.bam > z.discordants.unsorted.bam
docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools view -h $name.bam > z.splitters.before.bam
docker run --rm -v $(pwd):/in/ -w /in/ erictdawson/lumpy-sv ../app/lumpy-sv/scripts/extractSplitReads_BwaMem -i z.splitters.before.bam > z.splitters.nearlythere.sam
docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools view -Sb z.splitters.nearlythere.sam > z.splitters.unsorted.bam
docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools sort z.discordants.unsorted.bam -o $name.discordants.bam
docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools sort z.splitters.unsorted.bam -o $name.splitters.bam
docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools index $name.bam
docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools index $name.discordants.bam
docker run --rm -v $(pwd):/in/ -w /in/ jweinstk/samtools samtools index $name.splitters.bam
rm z.splitters.unsorted.bam
rm z.discordants.unsorted.bam
rm z.splitters.nearlythere.sam
rm z.splitters.before.bam
rm z.unsorted.bam
rm y.sam
rm x.sam
echo "End of mapping procedure"
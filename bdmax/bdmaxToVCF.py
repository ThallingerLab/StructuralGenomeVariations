import glob
import sys
import argparse
import csv
import os

# adapted from : https://github.com/PapenfussLab/sv_benchmark/blob/master/breakdancer2vcf.py

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='bdmaxToVCF', description='Reformat bdmax output to bedpe')
    parser.add_argument('-i', '--inFile', type=str, nargs='+',
                        help='file containing bdmax tsv file')
    parser.add_argument('-r', '--reference', type=str, nargs='+',
                        help='base reference fasta')
    parser.add_argument('-c', '--contigs', type=str, nargs='+',
                        help='contigs and legths of contigs in form id:length')
    args = parser.parse_args()

in_file = vars(args).get("inFile")[0]

print(args)

if in_file is None:
    parser.print_help()
    sys.exit(1)

vcf_path = in_file.replace(".tsv", ".vcf")

vcf_file = open(vcf_path, "w")

vcf_file.write("##fileformat=VCFv4.2\n")

vcf_file.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
vcf_file.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
vcf_file.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
vcf_file.write("##ALT=<ID=BND,Description=\"Translocation\">\n")
vcf_file.write("##ALT=<ID=INS,Description=\"Insertion\">\n")

vcf_file.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")

vcf_file.write(
    "##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for POS2 coordinate in case of an inter-chromosomal translocation\">\n")
vcf_file.write(
    "##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Genomic position for CHR2 in case of an inter-chromosomal translocation\">\n")
vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n")
vcf_file.write("##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">\n")
vcf_file.write("##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">\n")
vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"sv_type of structural variant\">\n")

vcf_file.write(f"##reference={vars(args).get('reference')[0]}\n")

for cont in vars(args).get("contigs"):
    conts = cont.split(":")
    vcf_file.write(f"##contig=<ID={conts[0]},length={conts[1]}>\n")

file_name = os.path.basename(in_file)

sv_sv_type = "."


def toVcfBreakend(localChr, localPos, localPositive, remoteChr, remotePos, remotePositive):
    if remotePositive:
        remote = "]" + remoteChr + ":" + str(remotePos) + "]"
    else:
        remote = "[" + remoteChr + ":" + str(remotePos) + "["
    if localPositive:
        return "N" + remote
    else:
        return remote + "N"


def isRefToBreakend(orientation):
    str = orientation.replace("-", "+").split("+")
    plusCount = int(str[0])
    minusCount = int(str[1])
    return plusCount >= minusCount


with open(in_file) as file:
    tsv = csv.reader(file, delimiter="\t")
    for idx, line in enumerate(tsv):
        if line[0].startswith("#Chr1"):
            vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCHR1\n")
        elif line[0].startswith("#"):
            # vcf_file.write("#"+'\t'.join(line) + "\n")
            print("Skipping header line:" + '\t'.join(line))
        else:
            # Chr1	Pos1	Orientation1	Chr2	Pos2	Orientation2	sv_type	Size	Score	num_Reads	num_Reads_lib	c1c173e1dd7a0e58ff9bd127c0cdce16.sc.bam
            leftChr = line[0]
            leftPos = int(line[1])
            leftOrientation = line[2]
            rightChr = line[3]
            rightPos = int(line[4])
            rightOrientation = line[5]
            sv_type = line[6]
            id = sv_type + "_" + str(idx)
            size = int(line[7])
            score = int(line[8])
            numreads = int(line[9])

            commonInfo = f";Chr1={leftChr};Pos1={leftPos};Orient1={leftOrientation};Chr2={rightChr};Pos2={rightPos};Orient2={rightOrientation};Type={sv_type};Size={size};Score={score};num_Reads={numreads}"

            if sv_type == "ITX":
                vcf_file.write(
                    f"{leftChr}\t{leftPos}\t{id}a\tN\tN[{rightChr}:{rightPos}[\t{score}\tPASS\tIMPRECISE;UNKNOWN_ORIENTATION;SVLEN={size};SVTYPE=BND;PARID={id}b;EVENT={sv_type}{id}" + commonInfo + "\t.\n" +
                    f"{rightChr}\t{rightPos}\t{id}b\tN\t]{leftChr}:{leftPos + size}]N\t{score}\tPASS\tIMPRECISE;UNKNOWN_ORIENTATION;SVLEN={size};SVTYPE=BND;PARID={id}a;EVENT={sv_type}{id}" + commonInfo + "\t.\n" +
                    f"{rightChr}\t{rightPos}\t{id}c\tN\tN[{leftChr}:{leftPos}[\t{score}\tPASS\tIMPRECISE;UNKNOWN_ORIENTATION;SVLEN={size};SVTYPE=BND;PARID={id}d;EVENT={sv_type}{id}" + commonInfo + "\t.\n" +
                    f"{leftChr}\t{leftPos + size}\t{id}d\tN\t]{rightChr}:{rightPos}]N\t{score}\tPASS\tIMPRECISE;UNKNOWN_ORIENTATION;SVLEN={size};SVTYPE=BND;PARID={id}c;EVENT={sv_type}{id}" + commonInfo + "\t.\n"
                )
            elif sv_type == "CTX":
                vcf_file.write(
                    f"{leftChr}\t{leftPos}\t{id}o\tN\t{toVcfBreakend(leftChr, leftPos, isRefToBreakend(leftOrientation), rightChr, rightPos, isRefToBreakend(rightOrientation))}\t{score}\tPASS\tIMPRECISE;UNKNOWN_ORIENTATION;SVLEN={size};SVTYPE=BND;PARID={id}h;EVENT={sv_type}{id};READS={numreads}" + commonInfo + "\t.\n" +
                    f"{rightChr}\t{rightPos}\t{id}h\tN\t{toVcfBreakend(rightChr, rightPos, isRefToBreakend(rightOrientation), leftChr, leftPos, isRefToBreakend(leftOrientation))}\t{score}\tPASS\tIMPRECISE;UNKNOWN_ORIENTATION;SVLEN={size};SVTYPE=BND;PARID={id}o;EVENT={sv_type}{id};READS={numreads}" + commonInfo + "\t.\n"
                )
            elif sv_type == "INS":
                vcf_file.write(
                    f"{leftChr}\t{leftPos}\t{id}\tN\t<{sv_type}>\t{score}\tPASS\tIMPRECISE;END={leftPos};SVLEN={size};SV_type={sv_type}" + commonInfo + "\t.\n"
                    # ;CIPOS=0,{12}
                )
            elif sv_type == "DEL":
                vcf_file.write(
                    f"{leftChr}\t{leftPos}\t{id}\tN\t<{sv_type}>\t{score}\tPASS\tIMPRECISE;END={leftPos + size};SVLEN={-size};SVTYPE={sv_type}" + commonInfo + "\t.\n"
                    # ;CIPOS={min(0, rightPos - leftPos - size)},{max(0, rightPos - leftPos - size)}
                )
            elif sv_type == "INV":
                # what does rightPos mean for INV? these output lines don't make any sense at all
                vcf_file.write(
                    # chr12	512749	0+2-	chr12	512812	0+2-	INV	74454768	41	2	/home/users/allstaff/cameron.d/i/data.fastcompare/3b57e681753bebdbde7eec2eb9da6510.sc.bam|2	NA
                    f"{leftChr}\t{leftPos}\t{id}\tN\t<{sv_type}>\t{score}\tPASS\tIMPRECISE;END={rightPos};SVLEN={size};SVTYPE={sv_type}" + commonInfo + "\t.\n"
                    # ;CIPOS={12},{13}
                )
            else:
                print("UNHANDLED VARIANT CALL" + line)

vcf_file.close()

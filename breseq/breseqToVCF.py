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

vcf_path = in_file.replace(".gd", ".vcf")

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

with open(in_file) as file:
    tsv = csv.reader(file, delimiter="\t")
    vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCHR1\n")
    for idx, line in enumerate(tsv):
        #if line[0].startswith("#Chr1"):
        #    vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")
        if line[0].startswith("#"):
            # vcf_file.write("#"+'\t'.join(line) + "\n")
            print("Skipping header line:" + '\t'.join(line))
        else:
            sv_type = line[0]
            sv_id = line[1]
            leftChr = line[3]
            leftPos = int(line[4])
            rightChr = leftChr
            rightPos = leftPos
            id = sv_type + "_" + sv_id

            if sv_type == "MC":
                rightPos = int(line[5])
                leftOrientation = True
                rightOrientation = False

                vcf_file.write(f"{leftChr}\t{leftPos}\t{id}a\tN\tN]{leftChr}:{leftPos}]\t.\tPASS\tCHR2={leftChr};END={leftPos};SVLEN=0;SVTYPE=BE;\n")
                vcf_file.write(f"{leftChr}\t{rightPos}\t{id}b\tN\t[{leftChr}:{rightPos}[N\t.\tPASS\tCHR2={leftChr};END={rightPos};SVLEN=0;SVTYPE=BE;\n")

            elif sv_type == "JC":
                leftOrientation = int(line[5]) > 0
                rightChr = line[6]
                rightPos = int(line[7])
                rightOrientation = int(line[9]) > 0

                vcf_file.write(
                    f"{leftChr}\t{leftPos}\t{id}\tN\t{toVcfBreakend(leftChr, leftPos, leftOrientation, rightChr, rightPos, rightOrientation)}\t.\tPASS\tIMPRECISE;END={rightPos};SVLEN=0;SVTYPE=BP"  + "\t.\n"
                    # ;CIPOS={min(0, rightPos - leftPos - size)},{max(0, rightPos - leftPos - size)}
                )
            elif sv_type == "DEL":
                size = int(line[5])
                rightPos = leftPos + size
                leftOrientation = True
                rightOrientation = False

                vcf_file.write(
                    f"{leftChr}\t{leftPos}\t{id}\tN\t<DEL>\t.\tPASS\tIMPRECISE;END={rightPos};SVLEN={size};SVTYPE=DEL" + "\t.\n"
                    # ;CIPOS={min(0, rightPos - leftPos - size)},{max(0, rightPos - leftPos - size)}
                )
            elif sv_type == "INS":
                rightPos = leftPos + 1
                leftOrientation = True
                rightOrientation = False

                insert = line[5]
                size = len(insert)

                vcf_file.write(
                    f"{leftChr}\t{leftPos}\t{id}\tN\t<INS>\t.\tPASS\tIMPRECISE;END={rightPos};SVLEN={size};SVTYPE=INS" + "\t.\n"
                    # ;CIPOS={min(0, rightPos - leftPos - size)},{max(0, rightPos - leftPos - size)}
                )
            elif sv_type == "UN":
                print("UN evidence, not processed")
            else:
                print("UNHANDLED VARIANT CALL")

vcf_file.close()

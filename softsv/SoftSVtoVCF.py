import glob
import sys
import argparse
import csv
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='SoftSVtoBedPE', description='Reformat SoftSV output to bedpe')
    parser.add_argument('-i', '--inFolder', type=str, nargs='+',
                        help='input folder containing SoftSV output')
    parser.add_argument('-r', '--reference', type=str, nargs='+',
                        help='base reference fasta')
    parser.add_argument('-c', '--contigs', type=str, nargs='+',
                        help='contigs and legths of cotigs in form id:length')
    args = parser.parse_args()


in_folder = vars(args).get("inFolder")[0]

print(args)

if in_folder is None:
    parser.print_help()
    sys.exit(1)

path = in_folder + '/*.txt'
files = glob.glob(path, recursive=False)
vcf_path = f"{in_folder}/softsv_parsed.vcf"

vcf_file = open(vcf_path, "w")

vcf_file.write("##fileformat=VCFv4.2\n")

vcf_file.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
vcf_file.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
vcf_file.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
vcf_file.write("##ALT=<ID=BND,Description=\"Translocation\">\n")
vcf_file.write("##ALT=<ID=INS,Description=\"Insertion\">\n")

vcf_file.write("##FILTER=<ID=PASS,Description=\"All filters passed\">\n")

vcf_file.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for POS2 coordinate in case of an inter-chromosomal translocation\">\n")
vcf_file.write("##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Genomic position for CHR2 in case of an inter-chromosomal translocation\">\n")
vcf_file.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the structural variant\">\n")
vcf_file.write("##INFO=<ID=PE,Number=1,Type=Integer,Description=\"Paired-end support of the structural variant\">\n")
vcf_file.write("##INFO=<ID=SR,Number=1,Type=Integer,Description=\"Split-read support\">\n")
vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")

vcf_file.write(f"##reference={vars(args).get('reference')[0]}\n")

for cont in vars(args).get("contigs"):
    conts = cont.split(":")
    vcf_file.write(f"##contig=<ID={conts[0]},length={conts[1]}>\n")

vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\n")

for sv_file in files:

    file_name = os.path.basename(sv_file)

    sv_type = "."

    if file_name.startswith("deletion"):
        sv_type = "DEL"
    elif file_name.startswith("insertion"):
        sv_type = "INS"
    elif file_name.startswith("inversion"):
        sv_type = "INV"
    elif file_name.startswith("tandem"):
        sv_type = "DUP"
    elif file_name.startswith("translocation"):
        sv_type = "BND"

    with open(sv_file) as file:
        tsv = csv.reader(file, delimiter="\t")
        for idx, line in enumerate(tsv):
            if idx != 0:
                chr1 = line[0]
                pos = line[1]
                sv_len = 0
                chr2 = "."
                pe_reads = line[4]
                sr_reads = line[5]
                seq = line[6]

                if sv_type in ("DEL", "INS", "INV", "DUP"):
                    sv_len = line[3]
                    pos2 = -1
                elif sv_type == "BND":
                    pos2 = line[3]
                    chr2 = line[2]

                end = int(pos) + int(sv_len) - 1
                sv_id = sv_type + "_" + str(idx)

                if chr1.startswith("chr"):
                    chr1 = chr1.replace("chr", "")

                if chr2.startswith("chr"):
                    chr2 = chr2.replace("chr", "")

                vcf_file.write(
                    f"{chr1}\t{pos}\t{sv_id}\t.\t<{sv_type}>\t0\tPASS\tSVTYPE={sv_type};END={end};CHR2={chr2};POS2={pos2};PE={pe_reads};SR={sr_reads};CONSENSUS={seq}\t.\n")

vcf_file.close()

from bed_maker_supplementary import *

import sys
import argparse
import csv
from Bio import SeqIO

# adapted from : https://github.com/PapenfussLab/sv_benchmark/blob/master/breakdancer2vcf.py

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='HockBedToVCF', description='Reformat the bedmaker output to VCF')
    parser.add_argument('-i', '--inFile', type=str, nargs='+',
                        help='file containing summary bed file')
    parser.add_argument('-r', '--reference', type=str, nargs='+',
                        help='base reference fasta')
    parser.add_argument('-p', '--hasPlasmid', type=bool, nargs='?',
                        const=True, default=False,
                        help='indicate if ref has plasmids')
    # parser.add_argument('-c', '--contigs', type=str, nargs='+',
    #                     help='contigs and legths of contigs in form id:length')
    args = parser.parse_args()

in_file = vars(args).get("inFile")[0]
reference = vars(args).get("reference")[0]
has_plasmids = vars(args).get("hasPlasmid")
plasmid_length = 0

print(args)

if in_file is None or reference is None:
    parser.print_help()
    sys.exit(1)

import os

vcf_path = in_file.replace(".bed",".vcf")


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
vcf_file.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"BE for breakend and BND for breakpoint \">\n")
vcf_file.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV \">\n")
vcf_file.write("##INFO=<ID=EVENT,Number=1,Type=String,Description=\"EVENT the SV belongs to \">\n")
vcf_file.write("##INFO=<ID=MATEID,Number=1,Type=String,Description=\"ID of mate breakend or breakpoint \">\n")

vcf_file.write(f"##reference={vars(args).get('reference')[0]}\n")

# for cont in vars(args).get("contigs"):
#     conts = cont.split(":")
#     vcf_file.write(f"##contig=<ID={conts[0]},length={conts[1]}>\n")

referenceFile = open(reference, 'r')

for rec in SeqIO.parse(referenceFile, 'fasta'):
    name = rec.id
    seq = rec.seq
    seqLen = len(rec)
    vcf_file.write(f"##contig=<ID={name},length={seqLen}>\n")
    if has_plasmids and name == "plasmid":
        plasmid_length = seqLen
        print(plasmid_length)

vcf_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tCHR1\n")
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
    bed = csv.reader(file, delimiter="\t")
    for idx, line in enumerate(bed):
        chr1 = line[0]
        start = int(line[1])
        end = int(line[2])
        sv_def = line[3]
        insertion = line[4]

        if sv_def == options_dict.get("del"):
            vcf_file.write(f"{chr1}\t{start}\tDEL_{str(idx)}\tN\t<DEL>\t.\tPASS\tCHR2={chr1};END={end};SVLEN={str(end-start)};SVTYPE=DEL;SV_DEF={sv_def};EVENT=DEL_{str(idx)}\n")
        elif sv_def == options_dict.get("ins"):
            if has_plasmids:
                vcf_file.write(f"{chr1}\t{start}\tINS_{str(idx)}a\tN\t{toVcfBreakend(chr1,start,True,'plasmid',0,False)}\t.\tPASS\tCHR2={chr1};END={end};SVLEN={str(end-start)};SVTYPE=BND;SV_DEF={sv_def};EVENT=INS_{str(idx)}\n")
                vcf_file.write(f"{chr1}\t{end}\tINS_{str(idx)}b\tN\t{toVcfBreakend(chr1,end,False,'plasmid',plasmid_length,True)}\t.\tPASS\tCHR2={chr1};END={end};SVLEN={str(end-start)};SVTYPE=BND;SV_DEF={sv_def};EVENT=INS_{str(idx)}\n")
            else:
                vcf_file.write(f"{chr1}\t{start}\tINS_{str(idx)}\tN\t<INS>\t.\tPASS\tCHR2={chr1};END={end};SVLEN={str(end-start)};SVTYPE=INS;SV_DEF={sv_def};EVENT=INS_{str(idx)}\n")
        elif sv_def == options_dict.get("inv"):
            vcf_file.write(f"{chr1}\t{start}\tINV_{str(idx)}\tN\t<INV>\t.\tPASS\tCHR2={chr1};END={end};SVLEN={str(end-start)};SVTYPE=INV;SV_DEF={sv_def};EVENT=INV_{str(idx)}\n")
        elif sv_def == options_dict.get("tdu"):
            vcf_file.write(f"{chr1}\t{start}\tTDU_{str(idx)}\tN\t<DUP>\t.\tPASS\tCHR2={chr1};END={end};SVLEN={str(end-start)};SVTYPE=DUP;SV_DEF={sv_def};EVENT=TDU_{str(idx)}\n")
        elif sv_def == options_dict.get("itd"):
            vcf_file.write(f"{chr1}\t{end}\tITD_{str(idx)}a\tN\t{toVcfBreakend(chr1,end,True,chr1,end,False)}\t.\tPASS\tCHR2={chr1};END={end};SVLEN={str(end-start)};SVTYPE=BND;SV_DEF={sv_def};EVENT=ITD_{str(idx)}\n")
            vcf_file.write(f"{chr1}\t{start}\tITD_{str(idx)}b\tN\t{toVcfBreakend(chr1,start,False,chr1,end+1,True)}\t.\tPASS\tCHR2={chr1};END={end};SVLEN={str(end-start)};SVTYPE=BND;SV_DEF={sv_def};EVENT=ITD_{str(idx)}\n")
        elif sv_def in (options_dict.get("tcp"), options_dict.get("txp"), options_dict.get("tci"), options_dict.get("txi"), options_dict.get("chr")):
            inserts = insertion.split(":")
            chr2 = inserts[1]
            pos2 = int(inserts[2])
            direction = inserts[3]

            if sv_def in {options_dict.get("tcp"), options_dict.get("tci")}:
                if direction == "forward":
                    vcf_file.write(f"{chr1}\t{start}\tTCP_{str(idx)}a\tN\t{toVcfBreakend(chr1,start,False,chr2,pos2,True)}\t.\tPASS\tCHR2={chr2};END={pos2};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TCP_{str(idx)}\n")
                    vcf_file.write(f"{chr1}\t{end}\tTCP_{str(idx)}b\tN\t{toVcfBreakend(chr1,end,True,chr2,pos2+1,False)}\t.\tPASS\tCHR2={chr2};END={pos2+1};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TCP_{str(idx)}\n")
                else:
                    vcf_file.write(f"{chr1}\t{start}\tTCI_{str(idx)}a\tN\t{toVcfBreakend(chr1,start,False,chr2,pos2,False)}\t.\tPASS\tCHR2={chr2};END={pos2};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TCI_{str(idx)}\n")
                    vcf_file.write(f"{chr1}\t{end}\tTCI_{str(idx)}b\tN\t{toVcfBreakend(chr1,end,True,chr2,pos2+1,True)}\t.\tPASS\tCHR2={chr2};END={pos2+1};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TCI_{str(idx)}\n")
            elif sv_def in {options_dict.get("txp"), options_dict.get("txi")}:
                if direction == "forward":
                    vcf_file.write(f"{chr1}\t{start}\tTXP_{str(idx)}a\tN\t{toVcfBreakend(chr1,start,False,chr2,pos2,True)}\t.\tPASS\tCHR2={chr2};END={pos2};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TXP_{str(idx)}\n")
                    vcf_file.write(f"{chr1}\t{end}\tTXP_{str(idx)}b\tN\t{toVcfBreakend(chr1,end,True,chr2,pos2+1,False)}\t.\tPASS\tCHR2={chr2};END={pos2+1};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TXP_{str(idx)}\n")
                    vcf_file.write(f"{chr1}\t{start}\tTXP_{str(idx)}c\tN\t<DEL>\t.\tPASS\tCHR2={chr1};END={end};SVLEN={str(end-start)};SVTYPE=BND;SV_DEF={sv_def};EVENT=TXP_{str(idx)}\n")
                else:
                    vcf_file.write(f"{chr1}\t{start}\tTXI_{str(idx)}a\tN\t{toVcfBreakend(chr1,start,False,chr2,pos2,False)}\t.\tPASS\tCHR2={chr2};END={end};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TXI_{str(idx)}\n")
                    vcf_file.write(f"{chr1}\t{end}\tTXI_{str(idx)}b\tN\t{toVcfBreakend(chr1,end,True,chr2,pos2+1,True)}\t.\tPASS\tCHR2={chr2};END={end};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TXI_{str(idx)}\n")
                    vcf_file.write(f"{chr1}\t{start}\tTXI_{str(idx)}c\tN\t<DEL>\t.\tPASS\tCHR2={chr1};END={end};SVLEN={str(end-start)};SVTYPE=BND;SV_DEF={sv_def};EVENT=TXI_{str(idx)}\n")
            elif sv_def == options_dict.get("chr"):
                if direction == "forward":
                    vcf_file.write(f"{chr1}\t{start}\tCHR_{str(idx)}a\tN\t{toVcfBreakend(chr1,start,True,chr2,pos2+1,False)}\t.\tPASS\tCHR2={chr2};END={pos2+1};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=CHR_{str(idx)}\n")
                    vcf_file.write(f"{chr1}\t{end+1}\tCHR_{str(idx)}b\tN\t{toVcfBreakend(chr1,end+1,False,chr2,pos2,True)}\t.\tPASS\tCHR2={chr2};END={pos2};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=CHR_{str(idx)}\n")
                else:
                    vcf_file.write(f"{chr1}\t{start}\tCHR_{str(idx)}a\tN\t{toVcfBreakend(chr1,start,True,chr2,pos2,True)}\t.\tPASS\tCHR2={chr2};END={pos2};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=CHR_{str(idx)}\n")
                    vcf_file.write(f"{chr1}\t{end+1}\tCHR_{str(idx)}b\tN\t{toVcfBreakend(chr1,end+1,False,chr2,pos2+1,False)}\t.\tPASS\tCHR2={chr2};END={pos2+1};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=CHR_{str(idx)}\n")
        elif sv_def == options_dict.get("sub"):
            if has_plasmids:
                vcf_file.write(f"{chr1}\t{start}\tSUB_{str(idx)}a\tN\t{toVcfBreakend(chr1,start,True,'plasmid',0,False)}\t.\tPASS\tCHR2={chr1};END={start};SVLEN=0;SVTYPE=BE;SV_DEF={sv_def};EVENT=SUB_{str(idx)}\n")
                vcf_file.write(f"{chr1}\t{end}\tSUB_{str(idx)}b\tN\t{toVcfBreakend(chr1,start,False,'plasmid',plasmid_length,True)}\t.\tPASS\tCHR2={chr1};END={end};SVLEN=0;SVTYPE=BE;SV_DEF={sv_def};EVENT=SUB_{str(idx)}\n")
            else:
                vcf_file.write(f"{chr1}\t{start}\tSUB_{str(idx)}a\tN\tN]{chr1}:{start}]\t.\tPASS\tCHR2={chr1};END={start};SVLEN=0;SVTYPE=BE;SV_DEF={sv_def};EVENT=SUB_{str(idx)}\n")
                vcf_file.write(f"{chr1}\t{end}\tSUB_{str(idx)}b\tN\t[{chr1}:{end}[N\t.\tPASS\tCHR2={chr1};END={end};SVLEN=0;SVTYPE=BE;SV_DEF={sv_def};EVENT=SUB_{str(idx)}\n")
        elif sv_def == options_dict.get("tcs"):
            inserts = insertion.split(":")
            chr2 = inserts[1]
            pos1 = int(inserts[0])
            pos2 = int(inserts[2])
            direction = inserts[3]

            if direction == "forward":
                vcf_file.write(f"{chr1}\t{start}\tTCS_{str(idx)}a\tN\t{toVcfBreakend(chr1,start,True,chr2,pos1,False)}\t.\tPASS\tCHR2={chr2};END={pos1};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TCS_{str(idx)}\n")
                vcf_file.write(f"{chr1}\t{end}\tTCS_{str(idx)}b\tN\t{toVcfBreakend(chr1,end,False,chr2,pos2+1,True)}\t.\tPASS\tCHR2={chr2};END={pos2+1};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TCS_{str(idx)}\n")
            else:
                vcf_file.write(f"{chr1}\t{start}\tTCS_{str(idx)}a\tN\t{toVcfBreakend(chr1,start,True,chr2,pos2+1,True)}\t.\tPASS\tCHR2={chr2};END={end};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TSI_{str(idx)}\n")
                vcf_file.write(f"{chr1}\t{end}\tTCS_{str(idx)}b\tN\t{toVcfBreakend(chr1,end,False,chr2,pos1,False)}\t.\tPASS\tCHR2={chr2};END={end};SVLEN=0;SVTYPE=BND;SV_DEF={sv_def};EVENT=TSI_{str(idx)}\n")
        else:
            print("UNHANDLED VARIANT CALL: " + sv_def + " " + options_dict.get("tci"))

vcf_file.close()

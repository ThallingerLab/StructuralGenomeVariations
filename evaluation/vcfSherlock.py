import sys
import io
import getopt
import pandas as pd
import numpy as np

import subprocess
import pysam

# pd.set_option('display.max_rows', 500)
# pd.set_option('display.max_columns', 500)
# pd.set_option('display.width', 1000)

vcf_in = bed_in = out_folder = nuc_in = fasta_in = summary = ""
max_dist = 100
plasmids = False


def input_options(argv):
    global bed_in
    global vcf_in
    global out_folder
    global max_dist
    global summary
    global plasmids

    try:
        opts, args = getopt.getopt(argv, "hb:v:f:d:s:p")
    except getopt.GetoptError:
        print("Usage: -b <bedfile> -v <vcffile> -f <out_folder> -d <max_dist, default:100> -s <summary> ")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("Usage: -b <bedfile> -v <vcffile> -f <out_folder> -d <max_dist, default:100>")
            sys.exit()
        elif opt == "-b":
            bed_in = arg
        elif opt == "-v":
            vcf_in = arg
        elif opt == "-f":
            out_folder = arg
        elif opt == "-s":
            summary = arg
        elif opt == "-d":
            max_dist = int(arg)
        elif opt == "-p":
            plasmids = True

    # print(f"\nComparing:\t\t {bed_in} and {vcf_in}\n")


if __name__ == "__main__":
    input_options(sys.argv[1:])

out_tp = out_folder + "/tp.bed"
out_fn = out_folder + "/fn.bed"
out_fp = out_folder + "/fp.vcf"


def read_vcf(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str, 'FORMAT': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})


def read_bed(path):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    bed = pd.read_csv(
        io.StringIO(''.join(lines)),
        sep='\t', names=["CHROM", "START", "END", "SV_TYPE", "INSERT", "NONE"],
        dtype={'CHROM': str, 'START': int, 'END': int, 'SV_TYPE': str, 'INSERT': str, 'NONE': str})
    try:
        bed[['END2', 'CHR2', 'START2', 'DIR']] = bed.INSERT.str.split(":", expand=True)
    except ValueError:
        bed[['END2', 'CHR2', 'START2', 'DIR']] = [None, None, None, None]
    return bed


def matches(bedNode, vcfNode, max_dist):
    """
    Check if two nodes match

    :param `s1`: range 1 start
    :type `s1`: int
    :param `e1`: range 1 end
    :type `e1`: int
    :param `s2`: range 2 start
    :type `s2`: int
    :param `e2`: range 2 end
    :type `e2`: int

    :return: True if ranges overlap
    :rtype: bool
    """

    bedNode_min = min(int(bedNode[1]), int(bedNode[2]))
    bedNode_max = max(int(bedNode[1]), int(bedNode[2]))
    vcfNode_min = min(int(vcfNode[1]), int(vcfNode[2]))
    vcfNode_max = max(int(vcfNode[1]), int(vcfNode[2]))

    if bedNode[0] == vcfNode[0]:
        if vcfNode_min - max_dist < bedNode_min < vcfNode_min + max_dist and vcfNode_max - max_dist < bedNode_max < vcfNode_max + max_dist:
            return True
        else:
            return False
    else:
        return False


""" Parse the .vcf file and store it in a new list of vcf_records"""

# bed_in = "/Data/Analyses/2023/202308_CloveBiotech/20220819_bed/DEL-1/DEL-1_summary.bed"
# vcf_in = "/Data/Analyses/2023/202308_CloveBiotech/20220819_svs/DEL-1/clovebiotech_100/clovebiotech.vcf"

vcf = read_vcf(vcf_in)
info_strings = '{"' + vcf.INFO.str.replace('CHR2=;', 'CHR2=None;').str.split(';').str.join('","').str.replace('=',
                                                                                                              '":"').str.replace(
    "\"\",", "") + '"}'
info_df = pd.json_normalize(info_strings.apply(eval))

vcfFull = pd.concat([vcf, info_df], axis=1)
bed = read_bed(bed_in)

bed_nodes = pd.DataFrame({'CHR': pd.Series(dtype='str'),
                          'START': pd.Series(dtype='int'),
                          'END': pd.Series(dtype='int'),
                          'INDEX': pd.Series(dtype='int')})

vcf_nodes = pd.DataFrame({'CHR': pd.Series(dtype='str'),
                          'START': pd.Series(dtype='int'),
                          'END': pd.Series(dtype='int'),
                          'INDEX': pd.Series(dtype='int')})

for bedIndex, bedRow in bed.iterrows():

    if pd.isna(bedRow.CHR2):
        if bedRow.SV_TYPE == "substitution" and not plasmids:
            print(bedRow)
            # if its a substitution add two breakends
            bed_nodes.loc[len(bed_nodes)] = {'CHR': bedRow.CHROM, 'START': bedRow.START,
                                             'END': bedRow.START, 'INDEX': bedIndex}
            bed_nodes.loc[len(bed_nodes)] = {'CHR': bedRow.CHROM, 'START': bedRow.END,
                                             'END': bedRow.END, 'INDEX': bedIndex}
        elif bedRow.SV_TYPE in ["substitution", "insertion"] and plasmids:
            bed_nodes.loc[len(bed_nodes)] = {'CHR': bedRow.CHROM, 'START': bedRow.START,
                                             'END': bedRow.END, 'INDEX': bedIndex}
            bed_nodes.loc[len(bed_nodes)] = {'CHR': "plasmid", 'START': 1,
                                             'END': len(bedRow.INSERT), 'INDEX': bedIndex}
        else:
            bed_nodes.loc[len(bed_nodes)] = {'CHR': bedRow.CHROM, 'START': bedRow.START,
                                             'END': bedRow.END, 'INDEX': bedIndex}
    else:
        bed_nodes.loc[len(bed_nodes)] = {'CHR': bedRow.CHROM, 'START': bedRow.START,
                                         'END': bedRow.END, 'INDEX': bedIndex}
        if not bedRow.END2 == "h1":
            bed_nodes.loc[len(bed_nodes)] = {'CHR': bedRow.CHR2, 'START': bedRow.END2,
                                             'END': bedRow.START2, 'INDEX': bedIndex}
        else:
            bed_nodes.loc[len(bed_nodes)] = {'CHR': bedRow.CHR2, 'START': bedRow.START2,
                                             'END': bedRow.START2, 'INDEX': bedIndex}

for vcfIndex, vcfRow in vcfFull.iterrows():

    if vcfRow.CHR2 == vcfRow.CHROM:
        vcf_nodes.loc[len(vcf_nodes)] = {'CHR': vcfRow.CHROM, 'START': vcfRow.POS,
                                         'END': vcfRow.END, 'INDEX': vcfIndex}

        if 'CHR_C' in vcfRow.keys() and not pd.isna(vcfRow.CHR_C):
            if 'END_C' in vcfRow.keys() and not pd.isna(vcfRow.END_C):
                vcf_nodes.loc[len(vcf_nodes)] = {'CHR': vcfRow.CHR_C, 'START': vcfRow.START_C,
                                                 'END': vcfRow.END_C, 'INDEX': vcfIndex}
            else:
                vcf_nodes.loc[len(vcf_nodes)] = {'CHR': vcfRow.CHR_C, 'START': vcfRow.START_C,
                                                 'END': vcfRow.START_C, 'INDEX': vcfIndex}
    elif vcfRow.CHR2 == 'None':
        vcf_nodes.loc[len(vcf_nodes)] = {'CHR': vcfRow.CHROM, 'START': vcfRow.POS,
                                         'END': vcfRow.POS, 'INDEX': vcfIndex}
    else:
        vcf_nodes.loc[len(vcf_nodes)] = {'CHR': vcfRow.CHROM, 'START': vcfRow.POS,
                                         'END': vcfRow.POS, 'INDEX': vcfIndex}
        vcf_nodes.loc[len(vcf_nodes)] = {'CHR': vcfRow.CHR2, 'START': vcfRow.END,
                                         'END': vcfRow.END, 'INDEX': vcfIndex}

hits = pd.DataFrame({'BedEntry': pd.Series(dtype='int'),
                     'BedNode': pd.Series(dtype='int'),
                     'VcfEntry': pd.Series(dtype='int'),
                     'VcfNode': pd.Series(dtype='int')})

for i, bedNode in bed_nodes.iterrows():
    for j, vcfNode in vcf_nodes.iterrows():
        # print(bedNode, vcfNode)
        if matches(bedNode, vcfNode, max_dist):
            hits.loc[len(hits)] = {'BedNode': i, 'BedEntry': int(bedNode[3]),
                                   'VcfNode': j, 'VcfEntry': int(vcfNode[3])}

BedTrueNodes = bed_nodes.INDEX.value_counts()
BedCalledNodes = hits.BedEntry.value_counts()

TP = list()
FN = list()

for key in BedTrueNodes.keys():
    if key in BedCalledNodes.keys():
        if BedTrueNodes[key] <= BedCalledNodes[key]:
            TP.append(key)
        else:
            FN.append(key)
    else:
        FN.append(key)

true_hits = hits[hits['BedEntry'].isin(TP)]

FP = np.unique(vcf_nodes.INDEX)
FP = np.setdiff1d(FP, true_hits.VcfEntry)

tp_bed = bed.iloc[TP]
tp_bed.to_csv(out_tp, index=False, sep="\t")

fn_bed = bed.iloc[FN]
fn_bed.to_csv(out_fn, index=False, sep="\t")

fp_vcf = vcf.iloc[FP]
fp_vcf.to_csv(out_fp, index=False, sep="\t")

if not summary == "":
    outfile = open(summary, "a")
    outfile.write(f"{bed_in}\t{vcf_in}\t{len(TP)}\t{len(FN)}\t{FP.size}\n")

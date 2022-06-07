import random

from Bio import SeqIO
import sys
import getopt
import subprocess
import copy
import numpy as np
import pandas as pd

from bed_maker_supplementary import *

""" STEP 1 - Get input arguments from command line: """

in_genome = excl_bed = outfile = insert = random_dna_length = n_rand_svs_string = rand_sv_type = random_dna_option = svs_string = art_options = ""


def input_options(argv):
    global in_genome
    global excl_bed
    global outfile
    global insert
    global n_rand_svs_string
    global rand_sv_type
    global random_dna_option
    global svs_string
    global art_options
    try:
        opts, args = getopt.getopt(argv, "hg:i:o:n:t:r:s:e:")
    except getopt.GetoptError:
        print(help_dict["input_help"])
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print(help_dict["input_help"])
            sys.exit()
        elif opt == "-g":
            in_genome = arg
        elif opt == "-o":
            outfile = arg
        elif opt == "-i":
            insert = arg
        elif opt == "-n":
            n_rand_svs_string = arg
        elif opt == "-t":
            rand_sv_type = arg
        elif opt == "-r":
            random_dna_option = arg
        elif opt == "-s":
            svs_string = arg
        elif opt == "-e":
            excl_bed = arg
    print(f"\nExecuting on reference genome:\t\t {in_genome}")
    print(f"Using insertion sequences from file:\t {insert}")
    print(f"Writing .bed files to:\t\t\t {outfile}.bed")


if __name__ == "__main__":
    input_options(sys.argv[1:])

""" STEP 2 - Generate Seq list from .fasta inputs: """

input_genome = list(SeqIO.parse(in_genome, "fasta"))
postprocess_chr = False

if insert == "":
    input_insert = [rand_dna(100), rand_dna(50)]
else:
    input_insert = list(SeqIO.parse(insert, "fasta"))

if excl_bed != "":
    excl_df = pd.read_csv(excl_bed, sep='\t', header=None)
    excl_df.columns = ["Chromosome", "Start", "End", "Name"]

""" STEP 3 - Check inputs """

n_rand_svs = 0
try:
    n_rand_svs = int(n_rand_svs_string)
except ValueError:
    n_rand_svs = 0

n_svs = 0
svs_list = [0] * 12
safety_dist = 10000

print(svs_string)

try:
    svs_list = [int(i) for i in svs_string.split(',')]
    if 12 != len(svs_list):
        print("Error: Option -s needs to have 12 integer values seperated by ,")
        sys.exit()
    n_svs = sum(svs_list)
except ValueError:
    n_svs = 0
    svs_list = [0] * 12

if rand_sv_type in abbreviations_list:
    sv_type = options_dict[rand_sv_type]
else:
    sv_type = "rand"

try:
    random_dna_length = int(random_dna_option)
except ValueError:
    # default of random DNA is 5
    random_dna_length = 5

""" STEP 7 - Start generating Files: """

bed_file_visor_path = f"{outfile}/{outfile}_VISOR.bed"
bed_file_visor = open(bed_file_visor_path, "a+")
bed_file_summary_path = f"{outfile}/{outfile}_summary.bed"
bed_file_summary = open(bed_file_summary_path, "a+")

if n_rand_svs > 0:
    if sv_type != "rand":
        print(f"\nAdding {str(n_rand_svs)} SVs of the type: {sv_type} with {random_dna_length} bp of random DNA. \n\n")
        svs_list[options_ind_dict[sv_type]] += n_rand_svs;

    elif sv_type == "rand":
        print(f"\nGenerating {str(n_rand_svs)} random SVs with {random_dna_length} bp of random DNA. \n\n")

        random_keys = ["del", "tdu", "itd", "tcp", "tci", "txp", "txi", "inv"]
        random_ind = [options_ind_dict.get(key) for key in random_keys]

        for i in range(n_rand_svs):
            ind = random.choice(random_ind)
            svs_list[ind] += 1

    n_svs = sum(svs_list)

if n_svs > 0:
    print(f"\nGenerating {str(n_svs)} SVs divided into:\n")

    for key in list(options_ind_dict):
        if svs_list[options_ind_dict[key]] != 0:
            print(f"{options_dict[key]}:\t\t\t{svs_list[options_ind_dict[key]]}")

    print("\n")

    for key in list(options_ind_dict):
        for i in range(svs_list[options_ind_dict[key]]):

            # TO-DO: This should be handled with the excluded regions
            """ Overlap-Checker: """
            all_chroms = list()
            all_starts = list()
            all_ends = list()
            all_breakpoint_chromosomes = list()
            all_breakpoints = list()

            chrom_ranges = list()
            break_point_chrom_ranges = list()

            # Randomly choosing chromosome and start

            random_chromosome_n = random.randrange(len(input_genome))
            random_chromosome = (input_genome[random_chromosome_n])
            random_chromosome_length = len(random_chromosome.seq)

            """ Excluding centromer/rRNA from CBS7435 """
            start = random.randint(safety_dist, random_chromosome_length - safety_dist)

            excl_chrom = excl_df[excl_df['Chromosome'] == random_chromosome.name]

            while any(np.logical_and(excl_chrom['Start'] <= start, excl_chrom['End'] >= start)):
                start = random.randint(safety_dist, random_chromosome_length - safety_dist)

            sv = pd.DataFrame(
                {'Chromosome': random_chromosome, 'Start': start - safety_dist, 'End': start + safety_dist,
                 'Name': 'sv'})
            excl_df = pd.concat([excl_df, sv], ignore_index=True, axis=0)

            # Define insert sequences
            insert_sequence = "None"
            break_point = None
            direction = None

            if key in {"ins", "sub"}:
                end = start + 1
                random_insert_choice = input_insert[random.randrange(len(input_insert))]
                random_insert = (rand_dna(random_dna_length) + random_insert_choice + rand_dna(random_dna_length))
                insert_sequence = str(random_insert.seq)
            elif key in {"del", "inv", "sub"}:
                end = start + random.randint(100, 5000)
            elif key in {"tdu", "itd", "tcp", "tci", "txp", "txi", "tcs","chr"}:
                end = start + random.randint(100, 5000)
                insert_sequence = "2"
                multikeys = {"tcp", "tci", "txp", "txi", "tcs", "chr"}

                if key in multikeys:
                    chrom_range = np.arange(0,len(input_genome)).tolist()

                    if key == "chr":
                        chrom_range.pop(random_chromosome_n)
                        random_chromosome_2_n = random.choice(chrom_range)
                    else:
                        random_chromosome_2_n = random.choice(chrom_range)

                    random_chromosome_2 = (input_genome[random_chromosome_2_n])
                    random_chromosome_2_length = len(random_chromosome_2.seq)

                    # Improve on this
                    break_point = random.randint(safety_dist, random_chromosome_2_length - safety_dist)
                    excl_chrom = excl_df[excl_df['Chromosome'] == random_chromosome_2.name]

                    while any(np.logical_and(excl_chrom['Start'] <= start, excl_chrom['End'] >= start)):
                        break_point = random.randint(safety_dist, random_chromosome_2_length - safety_dist)

                    sv_2 = pd.DataFrame({'Chromosome': random_chromosome_2, 'Start': break_point - safety_dist,
                                         'End': break_point + safety_dist, 'Name': 'sv'})
                    excl_df = pd.concat([excl_df, sv_2], ignore_index=True, axis=0)

                    if key in ("tcp", "txp"):
                        direction = "forward"
                    elif key in ("tci", "txi"):
                        direction = "reverse"
                    elif key == "tcs":
                        tcs_deletion_length = random.randint(500, 5000)
                        tcs_start = break_point - tcs_deletion_length
                        direction = random.choice(["forward", "reverse"])
                    elif key == "chr":
                        direction = random.choice(["forward", "reverse"])
                        postprocess_chr = True

                    insert_sequence = f"h1:{random_chromosome_2.id}:{str(break_point)}:{direction}"

            if key == "sub":
                bed_file_visor.write(
                    f"{random_chromosome.id}\t{str(start)}\t{str(end)}\tdeletion\tNone\t0\n")
                bed_file_visor.write(f"{random_chromosome.id}\t{str(start - 2)}\t{str(start - 1)}\tinsertion\t"
                                     f"{str(insert_sequence)}\t0\n")
                bed_file_summary.write(f"{random_chromosome.id}\t{str(start - 2)}\t{str(end)}\tsubstitution\t"
                                       f"{str(insert_sequence)}\t0\n")
            elif key == "tcs":
                if direction == "forward":
                    base_tcp = "tcp"
                else:
                    base_tcp = "tci"
                bed_file_visor.write(f"{random_chromosome_2.id}\t{str(tcs_start)}\t{str(break_point - 2)}\t"
                                     f"deletion\tNone\t0\n")
                bed_file_visor.write(
                    f"{random_chromosome.id}\t{str(start)}\t{str(end)}\t{options_dict[base_tcp]}\t{insert_sequence}\t0\n")
                bed_file_summary.write(
                    f"{random_chromosome_2.id}\t{str(tcs_start)}\t{str(break_point - 2)}\t{options_dict[key]}\t"
                    f"{start}:{random_chromosome.id}:{end}:{direction}\t0\n")
            elif key == "chr":
                bed_file_summary.write(
                    f"{random_chromosome.id}\t{str(start)}\t{str(start)}\t{options_dict[key]}\t"
                    f"{insert_sequence}\t0\n")
            else:
                bed_file_visor.write(
                    f"{random_chromosome.id}\t{str(start)}\t{str(end)}\t{options_dict[key]}\t{insert_sequence}\t0\n")
                bed_file_summary.write(f"{random_chromosome.id}\t{str(start)}\t{str(end)}\t{options_dict[key]}"
                                       f"\t{insert_sequence}\t0\n")

""" STEP 8 - Bash command to run VISOR HACk: """

bed_file_visor.close()
bed_file_summary.close()

bash_Command_HACK = [f"VISOR HACk -g {in_genome} -b {bed_file_visor_path} -o {outfile}_fa"]
subprocess.run(bash_Command_HACK, shell=True)


if postprocess_chr:
    summary_bed = pd.read_csv(bed_file_summary_path, sep='\t', header=None)
    summary_bed.columns = ["Chromosome", "Start", "End", "Name", "Insert", "Rand"]

    chr_entry = summary_bed[summary_bed['Name'] == options_dict['chr']]
    chr_insert = chr_entry['Insert'].to_string().split(":")

    visor_genome = SeqIO.to_dict(SeqIO.parse(f"{outfile}_fa/h1.fa", "fasta"))

    seq1 = visor_genome[chr_entry['Chromosome'].iloc[0]].seq[0:chr_entry['Start'].iloc[0]]
    seq2 = visor_genome[chr_entry['Chromosome'].iloc[0]].seq[(chr_entry['Start'].iloc[0]):len(visor_genome[chr_entry['Chromosome'].iloc[0]])]

    seq3 = visor_genome[chr_insert[1]].seq[0:int(chr_insert[2])]
    seq4 = visor_genome[chr_insert[1]].seq[int(chr_insert[2]):len(visor_genome[chr_insert[1]])]

    new_genome = copy.deepcopy(visor_genome)

    if chr_insert[3] == "forward":
        new_genome[chr_entry['Chromosome'].iloc[0]].seq = seq1 + seq4
        new_genome[chr_insert[1]].seq = seq3 + seq2
    else:
        new_genome[chr_entry['Chromosome'].iloc[0]].seq = seq1 + seq3.reverse_complement()
        new_genome[chr_insert[1]].seq = seq2.reverse_complement() + seq4

    with open(f"{outfile}_fa/h1_final.fa", 'w') as handle:
        SeqIO.write(new_genome.values(), handle, 'fasta')

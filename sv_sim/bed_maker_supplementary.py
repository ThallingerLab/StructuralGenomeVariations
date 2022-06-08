import random

# Define dictionary for list of options:

#    n_sub = int(string_counter_list[0])
#    n_tcs = int(string_counter_list[1])
#    n_ins = int(string_counter_list[2])
#    n_del = int(string_counter_list[3])
#    n_tdu = int(string_counter_list[4])
#    n_itd = int(string_counter_list[5])
#    n_tcp = int(string_counter_list[6])
#    n_txp = int(string_counter_list[7])
#    n_inv = int(string_counter_list[8])

options_ind_dict = {
    "sub": 0,
    "tcs": 1,
    "ins": 2,
    "del": 3,
    "tdu": 4,
    "itd": 5,
    "tcp": 6,
    "tci": 7,
    "txp": 8,
    "txi": 9,
    "inv": 10,
    "chr": 11
}

options_dict = {
    "sub": "substitution insertion",
    "tcs": "substitution translocation",
    "ins": "insertion",
    "del": "deletion",
    "tdu": "tandem duplication",
    "itd": "inverted tandem duplication",
    "tcp": "translocation copy-paste",
    "tci": "translocation copy-paste",
    "txp": "translocation cut-paste",
    "txi": "translocation cut-paste",
    "inv": "inversion",
    "chr": "chromosomal rearrangement"
}

# List of options for random SVs:
abbreviations_list = options_dict.keys()


options_list = [
    "deletion",
    "insertion",
    "tandem duplication",
    "inverted tandem duplication",
    "translocation copy-paste",
    "translocation cut-paste",
    "substitution",
    "inversion"
]

# Dictionary of prompts:
prompt_dict = {
    "times": "\nPlease define how many SVs should be generated: ",
    "in_option": "\nPlease specify type of desired SV: \n\n"
                 " Deletion [del] \n"
                 " Insertion [ins] \n"
                 " Tandem duplication [tdu] \n"
                 " Inverted tandem duplication [itd] \n"
                 " Translocation copy-paste [tcp] \n"
                 " Inverted Translocation [tci] \n"
                 " Translocation cut-paste [txp] \n"
                 " Inverted Translocation cut-paste [txi] \n"
                 " Substitution [sub] \n"
                 " Inversion [inv] \n\n"
                 "Otherwise, random SVs are generated: ",
    "interchromosomal": "\nAllow for interchromosomal duplications/translocations? [Y/N]: ",
    "random_dna": "\nEnter length (bp) of random DNA to flank insertion points: ",
}

# Dictionary of error notes:
error_dict = {
    "interchromosomal_error": "\n\nError: Please choose valid option! "
                              "Interchromosomal SV default option is: No. ",
    "end_error": "\n\nError: Parameter 'End Coordinate' could not be defined",
    "breakpoint_error": "\n\nError: Parameter 'Breakpoint' could not be defined.",
    "direction_error": "\n\nError: Parameter 'Direction' could not be defined.",
    "random_dna_error": "\n\nError: Parameter 'Random DNA length' has to be an integer."
                        "\nNo random DNA is added to the ends of the insert sequences.",
    "times_error": "\n\nError: Parameter 'Number of SVs' has to be an integer."
                   "\n'Number of SVs' default value is set to 5"
}

# Dictionary of help notes:
help_dict = {
    "input_help": " \n\n Use as follows: \n\n"
                  "\tprogram.py -g <reference_genome.fasta> -o <outfile.bed> -i <insert.fasta> \n "
                  "\t\t   -n <number of SVs> -t <SV type> -r <length of random dna> -s <number of each SV> \n"
                  "\t\t   -a <string with art options> \n\n"
                  "\n Possible SV types: \n\n"
                  "  Deletion [del] \n"
                  "  Insertion [ins] \n"
                  "  Tandem duplication [tdu] \n"
                  "  Inverted tandem duplication [itd] \n"
                  "  Translocation copy-paste [tcp] \n"
                  "  Inverted Translocation [tci] \n"
                  "  Translocation cut-paste [txp] \n"
                  "  Inverted Translocation cut-paste [txi] \n"
                  "  Substitution [sub] \n"
                  "  Inversion [inv] \n"
                  "  Random SVs [rand] \n\n"
                  " The 'number of each SV' is a string of 8 integers divided by a ',' to choose how \n"
                  " many of the following SVs should be randomly generated: \n\n"
                  "   <del,ins,tdu,itd,tcp,txp,sub,inv> \n\n"
                  " (The orientation of Translocations will be chosen at random.) \n\n"
                  " The options for art have to be divided only by a ',' \n\n"
}


# Generate random DNA
def rand_dna(length):
    dna = ""
    for count in range(length):
        dna += random.choice("CGTA")
    return dna

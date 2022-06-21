import random
import sys
import getopt
import subprocess
import os.path

in_genome = plasmids = randDna = exclBed = ""

def input_options(argv):
    global in_genome
    global plasmids
    global randDna
    global exclBed
    try:
        opts, args = getopt.getopt(argv, "hg:p:r:e:")
    except getopt.GetoptError:
        print("You did it wrong")
        sys.exit(2)
    for opt, arg in opts:
        if opt == "-h":
            print("No help here")
            sys.exit()
        elif opt == "-g":
            in_genome = arg
        elif opt == "-p":
            plasmids = arg
        elif opt == "-r":
            randDna = arg
        elif opt == "-e":
            exclBed = arg
    print(f"\nGenerating genomes for reference genome:\t\t {in_genome}")
    print(f"Using insertion sequences from file:\t {plasmids} and  {randDna}")


if __name__ == "__main__":
    input_options(sys.argv[1:])


runBedMakerBase = f"python /scripts/bed_maker.py -g {in_genome}"

# TO-DO: make this more relatable
randSvInd = {
    "del": [3],
    "dup": [4, 5, 6, 7],
    "trans": [8, 9],
    "inv": [10]
}

n = 10
randN = 4


def add_rand_svs(sv_list,ind_dict,n):
    for i in range(n):
        sv = random.choice(list(ind_dict))
        ind = random.choice(ind_dict[sv])
        sv_list[ind] += 1


### DELETIONS

for d in range(n):

    outFolder = "DEL-"+str(d+1)

    if os.path.isdir(outFolder):
        print(f"Directory {outFolder} exist, please remove before rerunning")
    else:
        delArray = [0] * 12
        delArray[3] = 1

        add_rand_svs(delArray,randSvInd,randN)

        delString = ",".join(map(str,delArray))

        if exclBed != "":
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {delString} -e {exclBed}"]
        else:
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {delString}"]

        subprocess.run(f"mkdir {outFolder}", shell=True)
        subprocess.run(runBedMaker, shell=True)

# RANDOM INSERTIONS PLASMIDS

for d in range(n):

    outFolder = "INSP-"+str(d+1)

    if os.path.isdir(outFolder):
        print(f"Directory {outFolder} exist, please remove before rerunning")
    else:
        insArray = [0] * 12
        insArray[2] = 1

        add_rand_svs(insArray,randSvInd,randN)

        insString = ",".join(map(str,insArray))

        if exclBed != "":
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {insString} -i {plasmids} -e {exclBed}"]
        else:
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {insString} -i {plasmids}"]

        subprocess.run(f"mkdir {outFolder}", shell=True)
        subprocess.run(runBedMaker, shell=True)

# RANDOM INSERTIONS E_Coli

for d in range(n):

    outFolder = "INSR-"+str(d+1)

    if os.path.isdir(outFolder):
        print(f"Directory {outFolder} exist, please remove before rerunning")
    else:
        insArray = [0] * 12
        insArray[2] = 1

        add_rand_svs(insArray,randSvInd,randN)

        insString = ",".join(map(str,insArray))

        if exclBed != "":
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {insString} -i {randDna} -e {exclBed}"]
        else:
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {insString} -i {randDna}"]

        subprocess.run(f"mkdir {outFolder}", shell=True)
        subprocess.run(runBedMaker, shell=True)

# SUBSTITUTIONS targeted integration

for d in range(n):

    outFolder = "SUBP-"+str(d+1)

    if os.path.isdir(outFolder):
        print(f"Directory {outFolder} exist, please remove before rerunning")
    else:
        subArray = [0] * 12
        subArray[0] = 1

        add_rand_svs(subArray,randSvInd,randN)

        subString = ",".join(map(str,subArray))

        if exclBed != "":
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {subString} -i {plasmids} -e {exclBed}"]
        else:
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {subString} -i {plasmids}"]

        subprocess.run(f"mkdir {outFolder}", shell=True)
        subprocess.run(runBedMaker, shell=True)

# SUBSTITUTIONS targeted integration E_Coli

for d in range(n):
    outFolder = "SUBR-"+str(d+1)

    if os.path.isdir(outFolder):
        print(f"Directory {outFolder} exist, please remove before rerunning")
    else:
        subArray = [0] * 12
        subArray[0] = 1

        add_rand_svs(subArray,randSvInd,randN)

        subString = ",".join(map(str,subArray))

        if exclBed != "":
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {subString} -i {randDna} -e {exclBed}"]
        else:
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {subString} -i {randDna}"]

        subprocess.run(f"mkdir {outFolder}", shell=True)
        subprocess.run(runBedMaker, shell=True)

# SUBSTITUTIONS KO

for d in range(n):
    outFolder = "SUBKO-"+str(d+1)

    if os.path.isdir(outFolder):
        print(f"Directory {outFolder} exist, please remove before rerunning")
    else:
        koArray = [0] * 12
        koArray[1] = 1

        add_rand_svs(koArray,randSvInd,randN)

        koString = ",".join(map(str,koArray))

        if exclBed != "":
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {koString} -e {exclBed}"]
        else:
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {koString}"]

        subprocess.run(f"mkdir {outFolder}", shell=True)
        subprocess.run(runBedMaker, shell=True)

# Chromosomal Rearrangement

for d in range(n):
    outFolder = "CHR-"+str(d+1)

    if os.path.isdir(outFolder):
        print(f"Directory {outFolder} exist, please remove before rerunning")
    else:
        chrArray = [0] * 12
        chrArray[11] = 1

        add_rand_svs(chrArray, randSvInd, randN)
        chrString = ",".join(map(str, chrArray))

        if exclBed != "":
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {chrString} -e {exclBed}"]
        else:
            runBedMaker = [runBedMakerBase + f" -o {outFolder} -s {chrString}"]

        subprocess.run(f"mkdir {outFolder}", shell=True)
        subprocess.run(runBedMaker, shell=True)



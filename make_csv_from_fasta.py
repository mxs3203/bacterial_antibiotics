import sys
import re
import itertools
from Row import Row
from Bio import SeqIO
import time
# from tqdm import tqdm
bases = ['A','C','T','G']

DEBUG = False

args = sys.argv
fileName = args[1]
outfile = args[2]

def allPossibleBases (size):
    perms = [''.join(p) for p in itertools.product(bases, repeat=size)]
    #print(perms)
    return perms



def compute_border_array(input):
    n = len(input)
    ba = []

    ba.append(0)
    for i in range(1, n):
        b = ba[i-1]
        while b > 0 and input[i] != input[b]:
            b = ba[b-1]

        if input[i] == input[b]:
            ba.append(b+1)
        else:
            ba.append(0)

    return ba


def ba_search(pattern, sequence):
    n = len(sequence)
    m = len(pattern)
    ba = compute_border_array(pattern+"$"+sequence)

    cnt = 0
    b = 0

    for i in range(0, len(ba)):
        if ba[i] == m:
            index = i - m + 1 - (m+1) #or i-2m
            #print("ba_search: Report match on: ", index)
            cnt = cnt + 1

    return (float(cnt)*len(pattern))/len(sequence)

def standardizeHeader(input):
    return input.replace(" ", "|")

def repetitions(s):
    r = re.compile(r"(.+?)\1+")
    for match in r.finditer(str(s)):
        yield (match.group(1), len(match.group(0))/len(match.group(1)))

def calcGC_content(frequencies):
    return (frequencies[3][1] + frequencies[1][1])/(frequencies[0][1] + frequencies[1][1] + frequencies[2][1] + frequencies[3][1])

def GC_AT_ration(frequencies):
    return (frequencies[0][1] + frequencies[2][1])/(frequencies[3][1] + frequencies[1][1])


rows = []
cnt = 0
fastaList = SeqIO.parse(fileName, "fasta")

file = open(outfile, "a+")       # file = open("test.csv", "a+")
for record in fastaList:
    if DEBUG:
        print(record.id)

    parts = record.id.split("|")
    reps = repetitions(record.seq)
    specie = record.description[record.description.find('[')+1:record.description.find(']')]
    totalRep = 0
    for rep in reps:
        if rep[1] >= 3.0 and len(rep[0]) >= 3:
            totalRep = totalRep + rep[1] * len(str(rep[1])) #sumOfAll(number of occurances * the size of repetition)
            if DEBUG:
                print(rep)
    largeRepFreq = totalRep/len(record.seq)
    # print(largeRepFreq)
    all = []
    for i in range(1, 5):
        permutations = allPossibleBases(size=i)
        total = 0
        frequencies = []
        for perm in permutations:
            freq = ba_search(perm, record.seq)
            frequencies.append([perm, freq])
            total = total + freq
            if DEBUG:
                print(perm, " Frequency: ", freq)
                print("This should be 1: ", total)
        all.append(frequencies)

    row = Row(id = parts[1], db_name = parts[0], desc = record.description, one_base_freq = all[0],
              two_base_freq = all[1], three_base_freq = all[2], four_base_freq = all[3],
              specie = specie, strand=parts[2], position_in_gene=parts[3], GC_content=calcGC_content(all[0]),
              AT_GC_ratio = GC_AT_ration(all[0]), record= record, largeRepFreq = largeRepFreq)
    if cnt == 0:
        file.write(row.getCSVHeader() + "\n")
    try:
        file.write(row.toRowCSV() + "\n")
        cnt = cnt + 1
        print("Analyzing gene " + parts[1] + " [ "+str(cnt)+"/4319 ] ") #398

    except IndexError:
        print("Index Error, probably pseudogene", record.id, " :: " , record.description)

file.close()

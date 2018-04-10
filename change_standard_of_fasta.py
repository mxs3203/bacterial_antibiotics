from Bio import SeqIO
import sys
import os

args = sys.argv
fileName = args[1]
outfile = args[2]

# fastaList = SeqIO.parse("Data/NormalGenes/data.fasta.fasta.fasta", "fasta")
# file = open("Data/NormalGenes/standardizedNormalGenes.fasta", "a+")

fastaList = SeqIO.parse(fileName, "fasta")
file = open(outfile, "a+")
for record in fastaList:
    idSplit = record.id.split("|")
    db = idSplit[0].strip()
    ids = idSplit[1].split("cds_")  # Split genome and protein ids
    genome_id = ids[0][0:len(ids[0]) - 1]

    tmp = ids[1].split("_")  # To get rid of final part in protein id
    if len(tmp) > 1:
        id = tmp[0] + "_" + tmp[1]  # Gives protein id
    else:
        id = "Pseudogene"  # When if condition is not met it means it is a Pseudogene
        continue

    descSplit = record.description.split("[")
    locationSplit = next((s for s in descSplit if "location" in s), None)  # Get the split with location information

    locationSplit = locationSplit.split("=")
    if locationSplit[1].find("complement") != -1:
        strand = "-"  # When location=complement(blalblabla) it means is the complementary strand (-)
        location = locationSplit[1].split("(")
        startPosition = str(locationSplit[1])[0:str(locationSplit[1]).find('.')]
        endPosition = str(locationSplit[1])[str(locationSplit[1]).find('.') + 2:str(locationSplit[1]).find(')')]
    else:
        strand = "+"
        startPosition = str(locationSplit[1])[0:str(locationSplit[1]).find('.')]
        endPosition = str(locationSplit[1])[str(locationSplit[1]).find('.') + 2:str(locationSplit[1]).find(']')]


    if startPosition.find("complement") != -1:
        startPosition = str(startPosition)[str(startPosition).find('(') + 1:len(str(startPosition))]

    species = "[Escherichia coli]"

    file.write(">" + db + "|" + id + "|" + strand + "|" + str(startPosition) + "-" + str(
        endPosition) + "|" + genome_id + "|" + "something " + species + "\n")
    file.write(str(record.seq) + "\n")



file.close()

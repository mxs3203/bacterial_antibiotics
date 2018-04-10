from Bio import SeqIO
import sys
from Bio import Entrez
args = sys.argv


def retrieveNCBIrecord(id, db="nucleotide", rettype="gb"):
    """INPUT: id=NCBI Accesion Number. OUTPUT: """
    Entrez.email = "alekss.ro@outlook.es"
    handle = Entrez.efetch(id=id, db=db, rettype=rettype)
    record = SeqIO.read(handle, format=rettype)
    handle.close

    return record


def getFeatures(record):
    try:
        record = retrieveNCBIrecord(getID(record))
        return record.features
    except (Exception, ConnectionError, ConnectionAbortedError, ConnectionRefusedError, RuntimeError):
        print("Something went wrong. Brace yourself")
        pass


def selectFeatures(record, featureNames, connect=True):
    """Give record and string with name of the feature you want to get"""
    if connect:
        out_features = []
        out_features2 = []
        record_features = getFeatures(record)
        for featureName in featureNames:
            if record_features:
                for feature in record_features:
                    tmp_feature = feature.qualifiers.get(featureName, "")
                    if tmp_feature:
                        out_features.append(str(tmp_feature[0]))
                        break
                    else:
                        out_features.append("NA")
                        break
            else:
                return ["Problem conecting with NCBI"]*len(featureNames)

        for out in out_features:
            out_features2.append(out.replace(",", ":"))
        return out_features2

    else:
        return ["NA"]*len(featureNames)

def getID(record):
    parts = record.id.split("|")
    return parts[1]

# For testing
# records = SeqIO.parse("testGenes.fasta", "fasta")
# featuresToGet = ["plasmid", "isolation_source", "host", "country"]
# print(featuresToGet)
# for record in records:
#     print(selectFeatures(record, featuresToGet))

# Some feature names:
# "mol_type", "organism", "strain", "plasmid", "isolation_source", "host", "country"


import FetchFeatures

class Row:
    id = ""
    desc = ""
    one_base_freq = []
    two_base_freq = []
    three_base_freq = []
    four_base_freq = []

    GC_content = -1
    AT_GC_ratio = -1
    large_repetitions = []

    specie = ""
    strand = ""
    db_name = ""
    position_in_gene = ""
    record = None

    largeRepFreq = 0
    def __init__(self, id = "", desc="", db_name= "", one_base_freq = [],
                 two_base_freq = [], three_base_freq = [], four_base_freq = [],
                 specie = "", strand = "",position_in_gene="", GC_content = -1,
                 AT_GC_ratio = -1, record = None, largeRepFreq = 0):
        self.id = id
        self.desc = desc
        self.one_base_freq = one_base_freq
        self.two_base_freq = two_base_freq
        self.three_base_freq = three_base_freq
        self.four_base_freq = four_base_freq
        self.specie = specie
        self.strand = strand
        self.db_name = db_name
        self.position_in_gene = position_in_gene
        self.GC_content = GC_content
        self.AT_GC_ratio = AT_GC_ratio
        self.record = record
        self.largeRepFreq = largeRepFreq


    """
        CSV STUFF
    """
    def listToCSV(self, list, index = 1, lastOne = False):
        output = ""
        for column in list:
            output = output + str(column[index]) + ","
        if lastOne:
            return output[0:len(output)-1]
        else:
            return output

    def getCSVHeader(self):
        return "ID,DESC,DB_name,start_position,end_position,plasmid_notPlasmid, isolation_source, host, country,specie,strand,largeRepFreq,gc_content,AT_GC_ratio," + \
               self.listToCSV(self.one_base_freq, index = 0) + self.listToCSV(self.two_base_freq, index = 0) + \
               self.listToCSV(self.three_base_freq, index = 0) + self.listToCSV(self.four_base_freq, index = 0, lastOne=True)


    def toRowCSV(self):
        positions = self.position_in_gene.split('-')
        return self.id + "," + self.desc + "," + self.db_name + "," + positions[0] + "," + positions[1]+"," +\
                ','.join(FetchFeatures.selectFeatures(self.record, ["plasmid", "isolation_source", "host", "country"], connect=False))+"," + \
                self.specie + "," + self.strand + "," +str(self.largeRepFreq) + ","+\
                str(self.GC_content) + "," + str(self.AT_GC_ratio) + "," + self.listToCSV(self.one_base_freq) + \
                self.listToCSV(self.two_base_freq) + self.listToCSV(self.three_base_freq) + self.listToCSV(self.four_base_freq, lastOne=True)



    """
        Printing stuff
    """
    def toStringDetail(self):
        return "ID: "+self.id + " DESC: " + self.desc + " DB_name: " + self.db_name + " ONEBase: " + \
               str(self.one_base_freq) + " TWOBase: " + str(self.two_base_freq) + " ThreeBase: " + \
               str(self.three_base_freq) + " FourBase: " + str(self.four_base_freq) + " Specie: " + \
               self.specie + " Strand: " + self.strand + " Position in gene: " + self.position_in_gene + " GC_Content: " + \
               str(self.GC_content) + " AT_GC_ratio: " + str(self.AT_GC_ratio)
    def toString(self):
        return "ID: "+self.id + " DESC: " + self.desc + " DB_name: " + self.db_name + " ONEBase: " + \
               str(self.one_base_freq) + " Specie: " + \
               self.specie + " Strand: " + self.strand + " Position in gene: " + self.position_in_gene + " GC_Content: " + \
               str(self.GC_content) + " AT_GC_ratio: " + str(self.AT_GC_ratio)

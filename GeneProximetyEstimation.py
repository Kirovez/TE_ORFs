"""
Script takes two bed files and return table with
closest feature from bed file 2 to the bed file 1 features
"""
import numpy
import sys
from collections import defaultdict

class Helpers():
    @staticmethod
    def gff2bed():
        with open(r"D:\PycharmProjects\sORF\GffParser\Ppatens_318_v3.3.gene_exons.gff3") as gff_file, \
                open("Ppatens_genes.bed", "w") as outfile:
            cnt = 0
            for lines in gff_file:
                if not lines.startswith("##"):
                    sp = lines.rstrip().split("\t")
                    # split gene name (e.g. ID=Pp3c1_10000;Name=Pp3c1_10000;ancestorIdentifier=Pp3c1_10000.v3.1).
                    # !!!Change if it is not necessary
                    feature_name = sp[8].split(";")[0].split("=")[1]
                    chromosome_name = sp[0]
                    strand = sp[6]
                    feature = sp[2]
                    coordinate1 = sp[3]
                    coordinate2 = sp[4]

                    if feature == "gene":
                        cnt += 1
                        outfile.write("\t".join([chromosome_name, coordinate1, coordinate2, strand, feature_name]) + "\n")
        print("number of genes: ", cnt)

    @staticmethod
    def selectTE():
        with open(r"D:\Moss2018\MossRepeatArticle\forIGV\Pp_LTRdigest_modified..bed") as full, \
            open("TE_body.bed", "w") as outFile:
            for lines in full:
                if 'LTR_retrotransposon' in lines:
                    outFile.write(lines)

class FindProximety():
    def __init__(self, bed1, bed2):
        self.bed1, self.bed2 = bed1, bed2
        self.gene_sites = self.generate_coordinates_dic()
        self.run()

    def run(self):
        self.findClosestgene()



    def generate_coordinates_dic(self):
        gene_sites = defaultdict(list)
        with open(self.bed2) as bed2_file:
            for lines in bed2_file:
                sp = lines.rstrip().split("\t")
                chromosome, start, end =  sp[0], sp[1], sp[2]
                gene_sites[chromosome].append(int(start))
                gene_sites[chromosome].append(int(end))
        return gene_sites

    def findClosestgene(self):
        with open(self.bed1) as infile, \
        open("Closest_genes_to_TEs.tab", "w") as outFile:
            for lines in  infile:
                sp = lines.rstrip().split("\t")
                chromosome, start, end =  sp[0], sp[1], sp[2]
                idx_start = numpy.searchsorted(self.gene_sites[chromosome], int(start))
                idx_end = numpy.searchsorted(self.gene_sites[chromosome], int(start))

                if idx_start != 0: # no genes upstream
                    start_distance = abs(self.gene_sites[chromosome][idx_start - 1] - int(start))
                else:
                    start_distance = 1e9 # big number

                if idx_end != len(self.gene_sites[chromosome]):
                    # print(self.gene_sites[chromosome])
                    # print(self.gene_sites[chromosome][idx_end-1])
                    # print(idx_end)

                    # print("*************")
                    end_distance = abs(self.gene_sites[chromosome][idx_end] - int(end))
                else:
                    end_distance = 1e9

                minimum_distance = min(start_distance, end_distance)

                outFile.write(sp[3] + "\t" + str(minimum_distance) + "\n")


#Helpers.gff2bed()
#Helpers.selectTE()
FindProximety(r"D:\PycharmProjects\MossRepeatomArticle\Scripts\TE_body.bed", "Ppatens_genes.bed")

#FindProximety(r"D:\PycharmProjects\MossRepeatomArticle\Scripts\test\test_bed1.txt", r"D:\PycharmProjects\MossRepeatomArticle\Scripts\test\test_bed2.txt")

#print(numpy.searchsorted([0,1,2,3], 4))
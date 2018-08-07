from Bio import SeqIO
from collections import defaultdict
import re
def __parsePEPfile():
    dic_pep_per_TE = defaultdict(list)
    all_pep = defaultdict(list)
    TEs = []
    with open(r"translatable_TE") as infile:
        for lines in infile:
            sp = lines.rstrip().split("\t")

            TE = sp[0]
            pep = sp[5]
            if pep not in all_pep[TE]:
                dic_pep_per_TE[pep].append(TE)
                all_pep[TE].append(pep)
            TEs.append(TE)
    print("Number of unique peptides: ", len(dic_pep_per_TE))
    print("Number of TEs: ", len(TEs))

    return  dic_pep_per_TE


def pepTogenome():
    """
    function takes moss genome and MS peptide table, translates the genome in 6 open reading frames
    :return: table with all possible peptide coordinates in genome and if stop codon is occured before start codon upstream pep then
    in the last colon '+' is written otherwise '-'
    """

    with open("All_alignPEP_TE_backTOgenome.txt", "w") as outfile:
        outfile.write("TE_id" + "\t" +
                      "Peptide" + "\t" +  "Peptide chromosome" + "\t" +
                      "Pep_start" + "\t"+ "Pep_end" + "\t"  + "Strand" + "\t" + "Start.codon" +
                      "\t" + "Stop.codon" +  "\t" + "Comments" + "\t" + "sORF peptide" + "\t" + "sORF length" + "\n")
        dic =  __parsePEPfile()
        for seq in SeqIO.parse("translated_genome.fasta","fasta"):
            strands, frames = seq.description.split(" ")[1:]
            print(seq.id)
            len_chr = len(str(seq.seq))

            tr = str(seq.seq)
            for pep in dic:
                #print("Peptide", pep)
                coords = [[i.start(), i.end()] for i in re.finditer(pep, tr)]
                #print(coords)
                for coordinate in coords:
                    #print(coordinate)
                    pep_start = coordinate[0]
                    pep_end = coordinate[1]



                    # look for start and stop codons upstream
                    upstream_region = tr[0:coordinate[0]]
                    u_stop_codon = upstream_region.rfind("*")

                    if pep.startswith("M"):
                        start_codon = coordinate[0]
                    else:
                        start_codon = upstream_region.rfind("M")

                    stop_codon = tr[coordinate[1]:].find("*") + len(upstream_region) + len(pep)
                    PSC = "-" # comments
                    sORF_pep = ""
                    sORF_size = 0
                    if start_codon == -1:
                        PSC = "SANF" #start not found
                        start_codon = -1
                        stop_codon = -1
                    elif stop_codon == -1:
                        PSC = "SONF" #stop not found
                        start_codon = -1
                        stop_codon = -1
                    elif u_stop_codon > start_codon:
                        PSC = "PSC" # stop codon occuries earlier than start => PSC
                        start_codon = -1
                        stop_codon = -1
                    else:
                        sORF_pep = tr[start_codon:stop_codon + 1]
                        sORF_size = stop_codon - start_codon + 1
                        if strands == -1:
                            start_codon = len_chr - start_codon
                            stop_codon = len_chr - stop_codon


                    if strands == -1:
                        # print(index, stop_codon, start_codon, len_chr)
                        pep_start = len_chr - pep_start
                        pep_end = len_chr - pep_end

                    for TEs in dic[pep]:
                        outfile.write( TEs  + "\t" +
                                       pep + "\t" + seq.id + "\t" +
                                        str(pep_start*3) + "\t" +
                                       str(pep_end*3) + "\t" +
                                        str((int(frames) + 1 )*int(strands)) + "\t" +
                                       str(start_codon*3) + "\t" +
                                       str(stop_codon*3) + "\t" +
                                       str(PSC) + "\t" + sORF_pep + "\t" + str(sORF_size) + "\n")


pepTogenome()

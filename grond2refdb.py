#!/usr/bin/env python

############################################################
# Argument Options
# Header needs to have phylum -> species all separated by '_'

import argparse
parser = argparse.ArgumentParser("Parses UNITE for PIPITS.")
parser.add_argument("--fas",
                    action = "store",
                    dest = "fas",
                    metavar = "fas",
                    help = "[REQUIRED] FASTA input",
                    required = True)
parser.add_argument("--tax",
                    action = "store", 
                    dest = "tax", 
                    metavar = "tax",
                    help = "[REQUIRED] TAXA input", 
                    required = True)
parser.add_argument("-r",
                    action = "store",
                    dest = "out_refseq",
                    metavar = "out_refseq",
                    help = "[REQUIRED] Reference FASTA",
                    required = True)
parser.add_argument("-t",
                    action = "store",
                    dest = "out_reftax",
                    metavar = "out_reftax",
                    help = "[REQUIRED] Reference TAXONOMY",
                    required = True)
options = parser.parse_args()

############################################################


from Bio import SeqIO
import sys

if __name__ == '__main__':

    no_seqs = 0
    with open(options.fas, "r") as infile:
        for line in infile:
            if line[0] == ">":
                no_seqs += 1

    in_fas = open(options.fas, "r")
    in_tax = open(options.tax, "r")
    
    out_refseq = open(options.out_refseq, "w")
    out_reftax = open(options.out_reftax, "w")

    taxDict = {}
    
    counter = 1

    reformatted_lineage = []
    for line in in_tax:

        # C_00000001      d__Bacteria|p__Proteobacteria|c__Gammaproteobacteria|o__Enterobacterales|f__Enterobacteriaceae|g__Enterobacter|s__Enterobacter asburiae_B
        # >Magnaporthaceae_sp|DQ528766|SH174123.07FU|reps|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Magnaporthales;f__Magnaporthaceae;g__unidentified;s__Magnaporthaceae_sp
        SH = line.rstrip().split("\t")[1]
        taxID = line.rstrip().split("\t")[0]
        full_lineage = SH.split("|")

        if len(full_lineage) != 7:
#            print(line)
            continue
        fl = full_lineage
        fl = [level.replace(" ", "_") for level in full_lineage]

        Kingdom =  fl[0].split("__")[1]
        Phylum =   fl[0].split("__")[1] + "|" + fl[1].split("__")[1]
        Class =    fl[0].split("__")[1] + "|" + fl[1].split("__")[1] + "|" + fl[2].split("__")[1]
        Order =    fl[0].split("__")[1] + "|" + fl[1].split("__")[1] + "|" + fl[2].split("__")[1] + "|" + fl[3].split("__")[1]
        Family =   fl[0].split("__")[1] + "|" + fl[1].split("__")[1] + "|" + fl[2].split("__")[1] + "|" + fl[3].split("__")[1] + "|" + fl[4].split("__")[1]
        Genus =    fl[0].split("__")[1] + "|" + fl[1].split("__")[1] + "|" + fl[2].split("__")[1] + "|" + fl[3].split("__")[1] + "|" + fl[4].split("__")[1] + "|" + fl[5].split("__")[1]
        Species =  fl[0].split("__")[1] + "|" + fl[1].split("__")[1] + "|" + fl[2].split("__")[1] + "|" + fl[3].split("__")[1] + "|" + fl[4].split("__")[1] + "|" + fl[5].split("__")[1] + "|" + fl[6].split("__")[1]

        reformatted_lineage = "Root;" + Kingdom + ";" + Phylum + ";" + Class + ";" + Order + ";" + Family + ";" + Genus + ";" + Species

        taxDict[taxID] = reformatted_lineage

#        out_refseq.write(">" + str(counter) + "\t" + reformatted_lineage + "\n")
#        out_refseq.write(str(record.seq) + "\n")
#        out_reftax.write(taxID + "\t" + reformatted_lineage + "\n")


    for record in SeqIO.parse(in_fas, "fasta"):

        seqID = record.description.rstrip()

        if seqID in taxDict.keys():
        
            out_refseq.write(">" + str(counter) + "\t" + taxDict[seqID] + "\n")
            out_refseq.write(str(record.seq) + "\n")

            out_reftax.write(str(counter) + "\t" + taxDict[seqID] + "\n")
        
            counter += 1

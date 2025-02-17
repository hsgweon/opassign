#!/usr/bin/python

import sys

############################################################
# Argument Options
# Header needs to have phylum -> species all separated by '_'

import argparse
parser = argparse.ArgumentParser("Reads and writes each entry as a single file.")
parser.add_argument("-i",
                    action = "store", 
                    dest = "infile", 
                    metavar = "infile",
                    help = "[REQUIRED]", 
                    required = True)
parser.add_argument("-o",
                    action = "store",
                    dest = "outfile",
                    metavar = "outfile",
                    help = "[REQUIRED]",
                    required = True)
parser.add_argument("-l",
                    action = "store",
                    dest = "sampleids",
                    metavar = "sampleids",
                    help = "[REQUIRED]",
                    required = True)
options = parser.parse_args()

############################################################

if __name__ == '__main__':

    infile = open(options.infile, "r")
    sampleids = open(options.sampleids, "r")
    outfile = open(options.outfile, "w")

    line_count = 0

    OTUs = {}

    OTUCount = 0
    
    while True:

        line = infile.readline().rstrip()

        if len(line) == 0:
            break
        
        line_count += 1
        fields = line.split("\t")               

        # Some filters
        if line[0] == '#' or (fields[0] != 'H' and fields[0] != 'S'):
            continue
            
        if len(fields) < 10:
            print >> sys.stderr, "line %d in .uc file has < 10 fields" % line_count
            exit(1)

            
        if fields[0] == "S":
            if fields[8] in OTUs.keys():
                print("something wrong: " + fields[8])
                
            OTUs[fields[8]] = [fields[8].split("_")[0]]
            OTUCount += 1
        elif fields[0] == "H":

      #      print(fields[8].split("_")[0])
            try:
                OTUs[fields[9]].append(fields[8].split("_")[0])
            except:
                OTUs[fields[9]] = [fields[8].split("_")[0]]

    print(len(OTUs))
                
    #for k, v in OTUs.items():
    #    print(k, v)
                
    # Load sample IDs
    sampleids = []
    for line in open(options.sampleids):
        if line[0] != "#":
            sampleids.append(line.split("\t")[0].rstrip())

    
    # header with sampleids
    outfile.write("#OTU_ID\t" + "\t".join(sampleids) + "\n")

    # OTU table
    for OTU, members in OTUs.items():

        # For each OTU
        outfile.write(OTU)

        # For each sample:
        subtotalnumber = 0
        for sampleid in sampleids:
            outfile.write("\t" + str(members.count(sampleid)))
            subtotalnumber += members.count(sampleid)

        if len(members) != subtotalnumber:
            print("Something wrong!")
            exit(1)
        
        outfile.write("\n")

    infile.close()
    outfile.close()

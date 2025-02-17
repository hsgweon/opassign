#!/usr/bin/python

############################################################
# Argument Options

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", 
                    action = "store", 
                    dest = "infile", 
                    metavar = "infile",
                    help = "[REQUIRED] input file", 
                    required = True)
parser.add_argument("-l", 
                    action = "store", 
                    dest = "label", 
                    metavar = "label",
                    help = "[REQUIRED] label", 
                    required = True)
parser.add_argument("-o", 
                    action = "store", 
                    dest = "outfile", 
                    metavar = "outfile",
                    help = "[REQUIRED] output file1", 
                    required = True)
options = parser.parse_args()

############################################################

import sys

in_file = open(options.infile, "r")
out_file = open(options.outfile, "w")

count = 1

for line in in_file:

        if(line[0] == ">"):
                header = ">" + options.label + "_" + line.rstrip().split(">")[1] + "\n"
                out_file.write(header)
                count += 1
        else:
                out_file.write(line)

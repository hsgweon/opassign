#!/usr/bin/env python

import sys

###################################################################
# Processes two input files, a .uc file (from VSEARCH clustering)
# and a taxonomy file to determine and output the 
# Lowest Common Ancestor (LCA) taxonomy for each representative 
# sequence in the .uc file
# 
# Example:
# ./uc2lca.py --uc example.uc --tax example.tsv -o lca.tsv
###################################################################

import argparse
parser = argparse.ArgumentParser("Reads a .uc and a tax file, and performs an LCA.")
parser.add_argument("--uc",
                    action = "store", 
                    dest = "uc", 
                    metavar = "uc",
                    help = "[REQUIRED]", 
                    required = True)
parser.add_argument("--tax",
                    action = "store",
                    dest = "tax",
                    metavar = "tax",
                    help = "[REQUIRED]",
                    required = True)
parser.add_argument("-o",
                    action = "store",
                    dest = "outfile",
                    metavar = "outfile",
                    help = "[REQUIRED]",
                    required = True)
options = parser.parse_args()

############################################################

def read_uc():
    repseq = {}
    with open(options.uc, "r") as f:

        line_count = 0

        for line in f:

            line = line.rstrip()
            line_count += 1

            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")

            if len(fields) < 10:

                print(f"line {line_count} in .uc file has < 10 fields", file = sys.stderr)
                sys.exit(1)

            record_type = fields[0]
            query_id = fields[8]
            target_id = fields[9]

            if record_type == "S":

                if query_id in repseq:
                    print(f"something wrong: {query_id}")

                repseq[query_id] = [tax[query_id]]

            elif record_type == "H":

                # print(query_id)
                repseq.setdefault(target_id, []).append(tax[query_id])

    return repseq

def read_tax():
    tax = {}
    with open(options.tax) as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts:  # Skip empty lines
                key = parts[0]
                taxon = '\t'.join(parts[1:]) if len(parts) > 1 else None
                taxon = taxonomy2list(taxon)
                tax[key] = taxon
    return tax

def taxonomy2list(taxon):
    parts = taxon.split('|')
    result = [part.split('__', 1)[1] if '__' in part else part for part in parts]
    return result

def lca(lists):
    if not lists:
        return []

    min_length = min(len(lst) for lst in lists)
    # print(min_length)
    lca = []

    for i in range(min_length):
        elements = [lst[i] for lst in lists]
        # print(elements)
        if all(element == elements[0] for element in elements):
            lca.append(elements[0])
        else:
            break

    max_length = max(len(lst) for lst in lists)
    lca.extend([''] * (max_length - len(lca)))

    return lca

def add_prefixes(lineage):
    prefixes = ['k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__']
    result = []
    for i, element in enumerate(lineage):
        if element:
            result.append(prefixes[i] + element)
        else:
            result.append(prefixes[i])
    return result

if __name__ == '__main__':

    tax = read_tax()
    repseq = read_uc()

    with open(options.outfile, "w") as outfile:
        for key, value in repseq.items():
            outfile.write(key + "\t" + "|".join(add_prefixes(lca(value))) + "\n")

    exit(0)

#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser("RDP Classifier - multiple CPUs")
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
parser.add_argument("-c",
                    action = "store",
                    dest = "chunksize",
                    metavar = "chunksize",
                    help = "[OPTIONAL] Chunk size",
                    default = "500")
parser.add_argument("-t",
                    action = "store",
                    dest = "threads",
                    metavar = "threads",
                    help = "[REQUIRED]",
                    required = True)
options = parser.parse_args()

import gzip
import multiprocessing
import os, shutil
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
import progressbar
from multiprocessing import Manager

def is_gzipped(filename):
    """Check if a file is gzipped."""
    try:
        with open(filename, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'  # Gzip magic number
    except IOError:
        return False

def open_file(filename, mode='rt'):
    """Open a file in text mode or gzip mode depending on the file type."""
    if is_gzipped(filename):
        return gzip.open(filename, mode)
    else:
        return open(filename, mode)

def count_lines(infile):
    """Counts the total number of sequences in the input file."""
    with open_file(infile, "rt") as f:
        return sum(1 for _ in SeqIO.parse(f, "fasta"))

def divide_file_chunks(infile, chunk_size, temp_dir):
    chunks = []
    with open_file(infile, "rt") as f:
        records = list(SeqIO.parse(f, "fasta"))
        for i in range(0, len(records), chunk_size):
            chunk_filename = os.path.join(temp_dir, f"chunk_{i//chunk_size}.fasta.gz")
            with gzip.open(chunk_filename, "wt") as chunk_file:
                SeqIO.write(records[i:i + chunk_size], chunk_file, "fasta")
            chunks.append(chunk_filename)
    return chunks

def process_chunk(args):
    """Runs a shell command to count lines in a chunk file."""
    chunk_file, progress_counter, lock, temp_file_prefix_cpu, temp_dir, chunk_size = args
    output_file = os.path.join(temp_dir, f"{temp_file_prefix_cpu}_assigned_taxonomy_rdp_raw.txt")
    
    MAX_RAM = "1000g"
    PATH_CLASSIFIER = "/home/xc917132/applications/rdp_classifier_2.14/dist/classifier.jar"
    PATH_RDPCLASSIFIER_DB = "/home/xc917132/applications/rdp_db/GROND_trained_207/GROND_retrained/rRNAClassifier.properties"

    command = f"java -Xms512M -Xmx{MAX_RAM} -jar {PATH_CLASSIFIER} classify -t {PATH_RDPCLASSIFIER_DB} -o {output_file} {chunk_file}"
    subprocess.run(command, shell=True, check=True)

    with lock:
        progress_counter.value += chunk_size  # Increment by the chunk size (number of lines in this chunk)
        # progress_counter.value += 1
        # pbar.update(progress_counter.value)
        pbar.update(min(progress_counter.value, total_lines))
    return output_file

if __name__ == "__main__":
    
    manager = Manager()
    progress_counter = manager.Value('i', 0)
    lock = manager.Lock()

    # Specify the name of the temporary directory
    temp_dir = "temp"
    shutil.rmtree(temp_dir, ignore_errors = True)  # Remove directory if it exists
    os.makedirs(temp_dir)

    # Number of CPUs to use
    num_cpus = multiprocessing.cpu_count()
    num_cpus = int(options.threads)
    chunk_size = int(options.chunksize) # Number of sequences in a chunk
    chunks = divide_file_chunks(options.infile, chunk_size, temp_dir)

    # Count total number of lines (sequences)
    total_lines = count_lines(options.infile)

    print(f"Total number of sequences to process: {total_lines}")

    # Start progressbar with total number of lines
    pbar = progressbar.ProgressBar(max_value = total_lines).start()

    # Pool
    pool_args = [(chunk, progress_counter, lock, f"cpu{i}", temp_dir, chunk_size) for i, chunk in enumerate(chunks)]

    with multiprocessing.Pool(processes = num_cpus) as pool:
        processed_files = pool.map(process_chunk, pool_args)

    # Merge all processed count files into final output
    with open(options.outfile, "w") as outfile:
        subprocess.run(f"cat {' '.join(processed_files)} > {options.outfile}", shell=True, check=True)
    
    shutil.rmtree(temp_dir, ignore_errors=True)
    
    # Finish
    pbar.finish()
    print("All done.")

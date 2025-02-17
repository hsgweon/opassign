#!/usr/bin/env python

import argparse
import gzip
import multiprocessing
import os
import shutil
import subprocess
from Bio import SeqIO
from multiprocessing import Manager, cpu_count
import progressbar

parser = argparse.ArgumentParser("seqdemu: A No-Nonsense Nanopore Demultiplexer.")
parser.add_argument("-i", action="store", dest="infile", metavar="infile", help="[REQUIRED]", required=True)
parser.add_argument("-o", action="store", dest="outfile", metavar="outfile", help="[REQUIRED]", required=True)
parser.add_argument("-c", action = "store", dest = "chunksize", metavar = "chunksize", help = "[OPTIONAL] Chunk size", default = "1000")
parser.add_argument("-t", action="store", dest="threads", metavar="threads", help="[REQUIRED]", required=True)
options = parser.parse_args()

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

def process_chunk(chunk_data, progress_counter, lock, chunk_size, temp_dir):
    """Process a chunk of data."""
    temp_file = os.path.join(temp_dir, f"temp_chunk_{progress_counter.value}.fasta")
    with open_file(temp_file, 'wt') as f:
        SeqIO.write(chunk_data, f, "fasta")

    output_file = os.path.join(temp_dir, f"{progress_counter.value}_assigned_taxonomy_rdp_raw.txt")
    
    MAX_RAM = "1000g"
    PATH_CLASSIFIER = "/home/xc917132/applications/rdp_classifier_2.14/dist/classifier.jar"
    PATH_RDPCLASSIFIER_DB = "/home/xc917132/applications/rdp_db/GROND_trained_207/GROND_retrained/rRNAClassifier.properties"

    command = f"java -XX:ActiveProcessorCount=1 -Xms512M -Xmx{MAX_RAM} -jar {PATH_CLASSIFIER} classify -t {PATH_RDPCLASSIFIER_DB} -o {output_file} {temp_file}"
    subprocess.run(command, shell=True, check=True)

    with lock:
        progress_counter.value += chunk_size  # Increment by the chunk size
        pbar.update(min(progress_counter.value, total_lines))
    return output_file

def chunk_file(infile, chunk_size, temp_dir):
    """Split the file into chunks."""
    chunk_data = []
    with open_file(infile, "rt") as f:
        for record in SeqIO.parse(f, "fasta"):
            chunk_data.append(record)
            if len(chunk_data) >= chunk_size:
                yield chunk_data
                chunk_data = []
        if chunk_data:  # Yield remaining sequences if any
            yield chunk_data

def count_lines(infile):
    """Counts the total number of sequences in the input file."""
    total_lines = 0
    with open_file(infile, "rt") as f:
        for _ in SeqIO.parse(f, "fasta"):
            total_lines += 1
    return total_lines

if __name__ == "__main__":
    manager = Manager()
    progress_counter = manager.Value('i', 0)
    lock = manager.Lock()

    # Specify the name of the temporary directory
    temp_dir = "temp"
    shutil.rmtree(temp_dir, ignore_errors=True)  # Remove directory if it exists
    os.makedirs(temp_dir)

    # Number of CPUs to use
    num_cpus = int(options.threads)
    chunk_size = int(options.chunksize) # Number of sequences in a chunk

    print(f"Counting the total number of sequences to process...")
    total_lines = count_lines(options.infile)
    print(f"Total number of sequences: {total_lines}")

    # Start progressbar with total number of lines
    pbar = progressbar.ProgressBar(max_value=total_lines).start()

    pool_args = [(chunk, progress_counter, lock, chunk_size, temp_dir) for chunk in chunk_file(options.infile, chunk_size, temp_dir)]

    with multiprocessing.Pool(processes=num_cpus) as pool:
        processed_files = pool.starmap(process_chunk, pool_args)

    # Merge all processed count files into final output
    with open(options.outfile, "w") as outfile:
        subprocess.run(f"cat {' '.join(processed_files)} > {options.outfile}", shell=True, check=True)
    
    shutil.rmtree(temp_dir, ignore_errors=True)
    
    # Finish
    pbar.finish()
    print("All done.")

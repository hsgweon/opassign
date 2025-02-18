#!/usr/bin/env python

import argparse
import gzip
import multiprocessing
import os
import shutil
import subprocess
from Bio import SeqIO
from multiprocessing import Manager
import progressbar

def is_gzipped(filename):
    """Check if a file is gzipped."""
    try:
        with open(filename, 'rb') as f:
            return f.read(2) == b'\x1f\x8b'  # Gzip magic number
    except IOError:
        return False

def open_file(filename, mode='rt'):
    """Open a file in text mode or gzip mode depending on the file type."""
    return gzip.open(filename, mode) if is_gzipped(filename) else open(filename, mode)

def count_lines(infile):
    """Counts the total number of sequences in the input file."""
    with open_file(infile, "rt") as f:
        return sum(1 for _ in SeqIO.parse(f, "fasta"))

def process_chunk(chunk_file, temp_file_prefix_cpu, temp_dir, chunk_size, progress_counter, lock):
    """Runs a shell command to classify sequences in a chunk and then deletes the chunk."""
    output_file = os.path.join(temp_dir, f"{temp_file_prefix_cpu}_assigned_taxonomy_rdp_raw.txt")
    
    MAX_RAM = "1000g"
    PATH_CLASSIFIER = "/home/xc917132/applications/rdp_classifier_2.14/dist/classifier.jar"
    #PATH_RDPCLASSIFIER_DB = "/home/xc917132/applications/rdp_db/GROND_trained_207/GROND_retrained/rRNAClassifier.properties"
    PATH_RDPCLASSIFIER_DB = options.rdp_db
    
    command = f"java -Xms512M -Xmx{MAX_RAM} -jar {PATH_CLASSIFIER} classify -t {PATH_RDPCLASSIFIER_DB} -o {output_file} {chunk_file}"
    subprocess.run(command, shell=True, check=True)
    
    with lock:
        progress_counter.value += chunk_size
        pbar.update(min(progress_counter.value, total_lines))
    
    os.remove(chunk_file)  # Delete the chunk after processing
    return output_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser("RDP Classifier - multiple CPUs")
    parser.add_argument("-i", required=True, dest="infile", metavar="infile", help="[REQUIRED]")
    parser.add_argument("-o", required=True, dest="outfile", metavar="outfile", help="[REQUIRED]")
    parser.add_argument("-c", dest="chunksize", metavar="chunksize", help="[OPTIONAL] Chunk size", default="1000")
    parser.add_argument("-t", required=True, dest="threads", metavar="threads", help="[REQUIRED]")
    parser.add_argument("-d", required=True, dest="rdp_db", metavar="rdp_db", help="[REQUIRED]")
    options = parser.parse_args()
    
    manager = Manager()
    progress_counter = manager.Value('i', 0)
    lock = manager.Lock()
    temp_dir = "temp"
    shutil.rmtree(temp_dir, ignore_errors=True)
    os.makedirs(temp_dir)
    
    num_cpus = int(options.threads)
    chunk_size = int(options.chunksize)
    
    print("Counting the total number of sequences to process...")
    total_lines = count_lines(options.infile)
    print(f"Total number of sequences to process: {total_lines}")
    
    pbar = progressbar.ProgressBar(max_value=total_lines).start()
    
    processed_files = []
    with open_file(options.infile, "rt") as f:
        records = list(SeqIO.parse(f, "fasta"))
        
        with multiprocessing.Pool(processes=num_cpus) as pool:
            for i in range(0, len(records), chunk_size * num_cpus):
                chunk_tasks = []
                for j in range(num_cpus):
                    index = i + j * chunk_size
                    if index >= len(records):
                        break
                    chunk_filename = os.path.join(temp_dir, f"chunk_{(index//chunk_size)+1:05d}.fasta.gz")
                    with gzip.open(chunk_filename, "wt") as chunk_file:
                        SeqIO.write(records[index:index + chunk_size], chunk_file, "fasta")
                    
                    chunk_tasks.append(pool.apply_async(
                        process_chunk, (chunk_filename, f"cpu{(index//chunk_size)+1:05d}", temp_dir, chunk_size, progress_counter, lock)
                    ))
                
                processed_files.extend([task.get() for task in chunk_tasks])
    
    with open(options.outfile, "w") as outfile:
        subprocess.run(f"cat {' '.join(processed_files)} > {options.outfile}", shell=True, check=True)
    
    shutil.rmtree(temp_dir, ignore_errors=True)
    pbar.finish()
    print("All done.")

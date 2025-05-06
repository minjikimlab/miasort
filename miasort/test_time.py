import time
from .process_complex import read_complexes
from .process_complex_cython import read_complexes_cython


def time_read_complexes():
    file_path = "../mia-sort_output/GM12878_CTCF-ChIA-PET_LHG0052H.bsorted.ENCLB716IME.hg38.complexes"
    bin_size = 5000000  # Adjust this if needed

    start_time = time.time()
    chroms, ids = read_complexes(file_path, bin_size)
    end_time = time.time()

    print(f"Time taken to read complexes: {end_time - start_time:.2f} seconds")


def time_read_complexes_cython():
    file_path = "../mia-sort_output/GM12878_CTCF-ChIA-PET_LHG0052H.bsorted.ENCLB716IME.hg38.complexes"
    bin_size = 5000000  # Adjust this if needed

    start_time = time.time()
    chroms, ids = read_complexes_cython(file_path, bin_size)
    end_time = time.time()

    print(f"Time taken to read complexes (cython): {end_time - start_time:.2f} seconds")




#!/usr/bin/env python3

import time
from miasort.process_complex import read_complexes
from miasort.process_complex import reduce_search_space
from pybedtools import BedTool


def time_reduce_search_space():
    # File paths
    complexes_file = "../mia-sort_output/GM12878_CTCF-ChIA-PET_LHG0052H.bsorted.ENCLB716IME.hg38.complexes"
    regions_file = "../mia-sort_output/GM12878_cr527_repeated_10.bedte"
    bin_size = 5000000

    # Step 1: Read complexes
    print("Reading complexes...")
    #start_time = time.time()
    chroms, ids = read_complexes_cython(complexes_file, bin_size)
    #end_time = time.time()
    #print(f"Time taken to read complexes: {end_time - start_time:.2f} seconds")

    # Step 2: Read one line from the regions file
    print("Reading one line from regions file...")
    with open(regions_file, 'r') as f:
        for line in f:
            if not line.startswith("#") and line.strip():  # Skip header or empty lines
                region_line = line.strip()
                break  # Use only the first valid line

    # Step 3: Time reduce_search_space
    print("Timing reduce_search_space...")
    start_time = time.time()
    reduced_bedtool = reduce_search_space(region_line, chroms, ids)
    end_time = time.time()
    print(f"Time taken for reduce_search_space: {end_time - start_time:.2f} seconds")

    # Optional: Print the reduced BedTool object (for debugging)
    #print(f"Reduced BedTool: {reduced_bedtool}")

if __name__ == "__main__":
    time_reduce_search_space()


#if __name__ == "__main__":
#    time_read_complexes()
#    time_read_complexes_cython()
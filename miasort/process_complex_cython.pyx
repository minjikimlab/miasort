from libc.stdlib cimport malloc, free
from libc.string cimport memcpy
from cpython cimport array
import csv

def read_complexes_cython(str file_path, int bin_size=5000000):
    """
    Reads the complexes file and builds chroms and ids data structures.
    """
    cdef dict chroms = {}
    cdef dict ids = {}
    cdef dict name_to_id = {}  # Mapping from name (string) to integer ID
    cdef int id_counter = 0    # Counter to assign unique integer IDs
    cdef str line
    cdef list fields
    cdef str chrom, name
    cdef int start, end, num_frags, start_bin, name_id

    with open(file_path, 'r') as f:
        for line in f:
            # Skip header lines or empty lines
            if line.startswith("#") or line.strip() == "":
                continue

            # Parse the line
            fields = line.strip().split()
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            num_frags = int(fields[3])
            name = fields[4]

            # Dynamically assign an integer ID to the name
            if name not in name_to_id:
                name_to_id[name] = id_counter
                id_counter += 1
            name_id = name_to_id[name]  # Get the integer ID for the name

            # Calculate the bin index
            start_bin = start // bin_size

            # Add to chroms
            if chrom not in chroms:
                chroms[chrom] = {}
            if start_bin not in chroms[chrom]:
                chroms[chrom][start_bin] = []
            chroms[chrom][start_bin].append(name_id)

            # Add to ids with (name_id, chrom) as the key
            if (name_id, chrom) not in ids:
                ids[(name_id, chrom)] = []
            ids[(name_id, chrom)].append((start, end, num_frags))

    return chroms, ids
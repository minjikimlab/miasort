# version 0.1: using name as the key
'''def read_complexes(file_path, bin_size=5000000):
    # creates dict for complexes
    chroms = {}
    ids = {}
    with open(file_path, 'r') as f:
            for line in f:
                # Skip header lines if present
                if line.startswith("#") or line.strip() == "":
                    continue

                fields = line.strip().split()
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                num_frags = int(fields[3])
                name = fields[4] 

                # Calculate the bin index
                start_bin = start // bin_size

                if chrom not in chroms:
                    chroms[chrom] = {}
                if start_bin not in chroms[chrom]:
                    chroms[chrom][start_bin] = []
                chroms[chrom][start_bin].append(name)

                if name not in ids:
                    ids[name] = [(chrom, start, end, num_frags)]
                else:
                    ids[name].append((chrom, start, end, num_frags))

    return chroms, ids'''


from collections import defaultdict
# can we change the bin_size for performance improvements?
def read_complexes(file_path, bin_size=5000000):
    # Creates dict for complexes
    chroms = defaultdict(lambda: defaultdict(list))
    ids = defaultdict(list)
    name_to_id = {}
    id_counter = 0

    with open(file_path, 'r') as f:
        for line in f:
            # Skip header lines if present
            if line.startswith("#") or line.strip() == "":
                continue

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
            name_id = name_to_id[name]

            # Calculate the bin index
            start_bin = start // bin_size

            # Populate the chroms dictionary
            chroms[chrom][start_bin].append(name_id)

            # Populate the ids dictionary
            ids[(name_id, chrom)].append((start, end, num_frags))

    return chroms, ids


def get_region_bin(region_line, bin_size=5000000):
    # Extract chromosome and start position from the region line
    start_list = []
    chrom_list = []

    fields = region_line.strip().split()
    num_regs = int((len(fields) - 1)/3)
    for i in range(0, num_regs * 3, 3):  # Process each region in the line
        chrom_list.append(fields[i])
        start_list.append(int(fields[i+1]) // bin_size)

    # get list of (chrom, start_bin) pairs from chrom_list and start_list
    uniq_tuples = list(zip(chrom_list, start_list))

    # add next bins in case of overlap
    # bins are size 5mb so only have to worry about 1 bin overlap
    # Iterate over a copy of the set to avoid modifying it while iterating
    for tup in uniq_tuples[:]:  
        add_tup = (tup[0], tup[1] + 1)
        uniq_tuples.append(add_tup)  # Adding to a set is efficient and avoids duplicates

    # Convert back to a list if needed
    uniq_tuples = list(set(uniq_tuples))
    
    return uniq_tuples



from pybedtools import BedTool
'''
@profile
def reduce_complex_bedtool(chroms, ids, regions):
    """
    For every region in regions, get all complexes with an ID matching the ones in the region.
    Create a BedTool object with chrom, start, and end for the matching complexes.
    """
    # List to store BED lines
    bed_lines = set()

    # Iterate over all regions
    for region in regions:
        chrom = region[0]
        start_bin = region[1]

        # Check if the chromosome and start_bin exist in the chroms dictionary
        if chrom in chroms and start_bin in chroms[chrom]:
            # Collect all IDs (names) for the current region
            matching_ids = set(chroms[chrom][start_bin])

            # Retrieve intervals for matching IDs from the ids dictionary
            for name in matching_ids:
                if name in ids:
                    
                    for chrom, start, end, num_frags in ids[name]:
                        # do we need num frags in here?
                        bed_lines.add(f"{chrom}\t{start}\t{end}\t{num_frags}\t{name}")

    # Create a BedTool object from the BED lines
    bedtool_object = BedTool("\n".join(bed_lines), from_string=True)

    return bedtool_object'''

@profile
def reduce_complex_bedtool(chroms, ids, regions):
    """
    For every region in regions, get all complexes with an ID matching the ones in the region.
    Create a BedTool object with chrom, start, and end for the matching complexes.
    """
    # Set to store BED lines (to avoid duplicates)
    bed_lines = set()

    # Iterate over all regions
    for region in regions:
        chrom = region[0]
        start_bin = region[1]

        # Check if the chromosome and start_bin exist in the chroms dictionary
        if chrom in chroms and start_bin in chroms[chrom]:
            # Collect all IDs (name_ids) for the current region
            matching_ids = set(chroms[chrom][start_bin])

            # Retrieve intervals for matching IDs from the ids dictionary
            for name_id in matching_ids:
                key = (name_id, chrom)  # Use (name_id, chrom) as the key
                if key in ids:
                    for start, end, num_frags in ids[key]:
                        # Add the BED line to the set
                        bed_lines.add(f"{chrom}\t{start}\t{end}\t{num_frags}\t{name_id}")

    # Create a BedTool object from the BED lines
    bedtool_object = BedTool("\n".join(bed_lines), from_string=True)

    return bedtool_object


@profile
def reduce_search_space(filter_regions_line, chroms, ids, bin_size=5000000):

     # read filter file and do process for each line
    fields = filter_regions_line.strip().split("\t")
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    start_bin = start // bin_size
    end_bin = end // bin_size

    tuples = []
    if start_bin == end_bin:
        tuples.append((chrom, start_bin))
        tuples.append((chrom, start_bin+1))
    else:
        for i in range(start_bin, end_bin + 1):
            tuples.append((chrom, i))

    # creates bedtool object for reduced complex space
    return  reduce_complex_bedtool(chroms, ids, tuples)
    

            
def filter_intersections(ChIA_Drop, filter_regions):
    intersected = ChIA_Drop.intersect(filter_regions, wa=True, wb=True)
        # Dictionary to store the intersected regions for each line of b
    filtered_intersections = {}

    for intersection in intersected:
        b_fields = intersection.fields[5:]  # 5 fields in a
        b_fields = ' '.join(b_fields)  # Make the key hashable
        # Check if the key exists, if not, add an empty list
        if b_fields not in filtered_intersections:
            filtered_intersections[b_fields] = []
        # Append the intersection to the list
        filtered_intersections[b_fields].append(intersection)
        # Convert lists to BedTool objects
    for i in filtered_intersections:
        filtered_intersections[i] = BedTool(filtered_intersections[i])

    return filtered_intersections




# test cases for the functions above
'''
import tempfile

def test_read_complexes():
    # Sample input data
    input_data = """\
# Header line
chr1    1000    2000    1    name1
chr1    3000    4000    1    name1
chr1    5000    6000    1    name2
chr2    7000    8000    1    name3
chr2    9000    10000   1    name3
"""

    # Expected output
    expected_output = {
        "chr1": {
            0: {
                "name1": [(1000, 2000), (3000, 4000)],
                "name2": [(5000, 6000)],
            }
        },
        "chr2": {
            0: {
                "name3": [(7000, 8000), (9000, 10000)],
            }
        },
    }

    # Create a temporary file with the input data
    with tempfile.NamedTemporaryFile(mode="w+", delete=False) as temp_file:
        temp_file.write(input_data)
        temp_file.seek(0)  # Go back to the beginning of the file

        # Call the function
        result = read_complexes(temp_file.name, bin_size=5000000)

    # Check if the result matches the expected output
    assert result == expected_output, f"Test failed! Expected {expected_output}, but got {result}"
    print("Test passed!")


def test_get_region_bin():
    # Sample input data
    region_line = "chr1 1000 2000 chr1 3000 4000 chr2 5000 6000 id1"
    bin_size = 5000000

    # Expected output
    expected_output = [
        ("chr1", 0),  # First region
        ("chr1", 1),  # Overlapping bin for first region
        ("chr1", 0),  # Second region
        ("chr1", 1),  # Overlapping bin for second region
        ("chr2", 0),  # Third region
        ("chr2", 1),  # Overlapping bin for third region
    ]

    # Remove duplicates and sort for comparison
    expected_output = sorted(set(expected_output))

    # Call the function
    result = get_region_bin(region_line, bin_size)

    # Remove duplicates and sort the result for comparison
    result = sorted(result)

    # Check if the result matches the expected output
    assert result == expected_output, f"Test failed! Expected {expected_output}, but got {result}"
    print("Test passed!")


def test_reduce_complex_bedtool():
    # Sample input data for chroms and ids
    chroms = {
        "chr1": {
            0: ["name1", "name2"],
        },
        "chr2": {
            1: ["name3", "name1"],
        },
    }
    ids = {
        "name1": [("chr1", 1000, 2000), ("chr1", 3000, 4000), ("chr2", 11000, 12000)],
        "name2": [("chr1", 5000, 6000)],
        "name3": [("chr2", 7000, 8000), ("chr2", 9000, 10000)],
    }

    # Sample regions
    regions = [("chr1", 0), ("chr2", 1)]

    # Expected output
    expected_output = """\
chr1    1000    2000    name1
chr1    3000    4000    name1
chr2    11000   12000   name1
chr1    5000    6000    name2
chr2    7000    8000    name3
chr2    9000    10000   name3
"""

    # Call the function
    bedtool_object = reduce_complex_bedtool(chroms, ids, regions)

    # Check if the result matches the expected output
    assert str(bedtool_object) == expected_output.strip(), f"Test failed! Expected:\n{expected_output}\nBut got:\n{bedtool_object}"
    print("Test passed!")

# Run the test
test_reduce_complex_bedtool()


# Run the test
test_get_region_bin()


# Run the test
test_read_complexes()
'''
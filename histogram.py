def generate_file(ranked_gems, output_file):
    histogram = {}
    num_gems = len(ranked_gems)
    for gem_tuple in ranked_gems:
        num_fragments = len(gem_tuple[1])
        if num_fragments not in histogram.keys():
            histogram[num_fragments] = 1
        else:
            histogram[num_fragments] += 1

    # write to a file
    with open(f"{output_file}.txt", "w") as file:
        # write the header
        file.write("num_frag_per_GEM\tnum_comp\n")
        # write the histogram data
        for key, value in histogram.items():
            file.write(f"{key}\t{value}\n")
        file.write(f"total\t{num_gems}")
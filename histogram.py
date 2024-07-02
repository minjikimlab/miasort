import os

def generate_file(ranked_gems, output_file, out_dir):
    histogram = {}
    num_gems = len(ranked_gems)
    for gem_tuple in ranked_gems:
        num_fragments = len(gem_tuple[1])
        if num_fragments not in histogram.keys():
            histogram[num_fragments] = 1
        else:
            histogram[num_fragments] += 1

    # determine the output path
    if out_dir != "/":
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        output_path = os.path.join(out_dir, f"{output_file}.txt")
    else:
        output_path = f"{output_file}.txt"

    # write to a file
    with open(output_path, "w") as file:
        # write the header
        file.write("num_frag_per_GEM\tnum_comp\n")
        # write the histogram data
        for key, value in histogram.items():
            file.write(f"{key}\t{value}\n")
        file.write(f"total\t{num_gems}")

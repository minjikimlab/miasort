import pybedtools
from pybedtools import BedTool
import argparse
import time
import plot
import sort


def main(path1, path2, type, anchor_line, output_file):
    pybedtools.helpers.cleanup()

    # Specify for the path of ___ALL.region.PEanno and import it (GEMs)
    ChIA_Drop = BedTool(path1)

    # Specify for the path of ____PETcnt_G9.motifannot and import it (anchors, regions)
    Region = BedTool(path2)

    # Use only the first line of the region file for the left and right references
    first_region = list(Region)[anchor_line]
    left_anchor = f"{first_region.chrom}\t{first_region.start}\t{first_region.end}"
    right_anchor = f"{first_region.fields[3]}\t{first_region.fields[4]}\t{first_region.fields[5]}"
    region = f"{first_region.chrom}\t{first_region.start}\t{first_region.fields[5]}"

    if type == "left":
        sort_start_time = time.time()
        ranked_gems = sort.process_left(ChIA_Drop, left_anchor, right_anchor, region)
        print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")

        plot_start_time = time.time()
        plot.plot_ranked_gems(ranked_gems, output_file, left_anchor, right_anchor)
        print(f"It took {time.time() - plot_start_time} secs in total to plot the GEMs")

    elif type == "right":
        sort_start_time = time.time()
        ranked_gems = sort.process_right(ChIA_Drop, left_anchor, right_anchor, region)
        print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")

        plot_start_time = time.time()
        plot.plot_ranked_gems(ranked_gems, output_file, left_anchor, right_anchor)
        print(f"It took {time.time() - plot_start_time} secs in total to plot the GEMs")

    elif type == "both":
        sort_start_time = time.time()
        ranked_gems = sort.process_both(ChIA_Drop, left_anchor, right_anchor, region)
        print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")

        plot_start_time = time.time()
        plot.plot_ranked_gems(ranked_gems, output_file, left_anchor, right_anchor)
        print(f"It took {time.time() - plot_start_time} secs in total to plot the GEMs")

    else:
        print(f"Type '{type}' is not supported.")


if __name__ == '__main__':
    start_time = time.time()  # record execution time

    parser = argparse.ArgumentParser()
    parser.add_argument('--path1', type=str, required=True, help='Path to the GEMs file (.PEanno)')
    parser.add_argument('--path2', type=str, required=True, help='Path to the regions file (.domains)')
    parser.add_argument('--type', type=str, required=True, help='Type of processing: left, right, or both')
    parser.add_argument('--anchor', type=str, required=True, help='The index of anchor line (0-based)')
    parser.add_argument('--output_file', type=str, required=True, help='Output file for the plot')
    args = parser.parse_args()
    path1 = args.path1
    path2 = args.path2
    type = args.type
    anchor_line = int(args.anchor)
    output_file = args.output_file

    main(path1, path2, type, anchor_line, output_file)

    print(f"It took {time.time() - start_time} secs in total to finish this program")

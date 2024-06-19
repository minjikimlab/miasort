import pybedtools
from pybedtools import BedTool
import argparse
import time
import plot
import sort


class ConditionalArgument(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values == 'middle' or values == "only-middle" or values == "only-middle-1frag":
            setattr(namespace, 'region_required', True)
        else:
            setattr(namespace, 'region_required', False)
        setattr(namespace, self.dest, values)


def main(path1, path2, type, anchor_line, output_file, middle_region):
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

    elif type == "right":
        sort_start_time = time.time()
        ranked_gems = sort.process_right(ChIA_Drop, left_anchor, right_anchor, region)
        print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")

    elif type == "both":
        sort_start_time = time.time()
        ranked_gems = sort.process_both(ChIA_Drop, left_anchor, right_anchor, region)
        print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")

    elif type == "middle":
        sort_start_time = time.time()
        ranked_gems = sort.process_middle(ChIA_Drop, left_anchor, right_anchor, region, middle_region)
        print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")

    elif type == "only-middle":
        sort_start_time = time.time()
        ranked_gems = sort.process_only_middle(ChIA_Drop, middle_region)
        print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")

    elif type == "only-middle-1frag":
        sort_start_time = time.time()
        ranked_gems = sort.process_only_middle_1frag(ChIA_Drop, middle_region)
        print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")

    else:
        print(f"Type '{type}' is not supported.")
        return

    plot_start_time = time.time()
    plot.plot_ranked_gems(ranked_gems, output_file, left_anchor, right_anchor, middle_region, type)
    print(f"It took {time.time() - plot_start_time} secs in total to plot the GEMs")


if __name__ == '__main__':
    start_time = time.time()  # record execution time

    parser = argparse.ArgumentParser()
    parser.add_argument('--path1', type=str, required=True, help='Path to the GEMs file (.PEanno)')
    parser.add_argument('--path2', type=str, required=True, help='Path to the regions file (.domains)')
    parser.add_argument('--type', type=str, required=True, help='Type of processing: left, right, middle, only-middle, only-middle-1frag or both', action=ConditionalArgument)
    parser.add_argument('--anchor', type=str, required=True, help='The index of anchor line (0-based)')
    parser.add_argument('--output_file', type=str, required=True, help='Output file for the plot')
    parser.add_argument('--region', type=str, help='Specific region to process when type is middle, only-middle')

    args = parser.parse_args()

    if args.region_required and not args.region:
        parser.error("--region is required when --type is 'middle'")

    path1 = args.path1
    path2 = args.path2
    processing_type = args.type
    anchor_line = int(args.anchor)
    output_file = args.output_file
    region = args.region if processing_type == 'middle' or processing_type == 'only-middle' or processing_type == 'only-middle-1frag' else None

    main(path1, path2, processing_type, anchor_line, output_file, region)

    print(f"It took {time.time() - start_time} secs in total to finish this program")

import pybedtools
from pybedtools import BedTool
import argparse
import time
import plot
import histogram
import sort


def main(start_time, path1, path2, type, num_fragments, anchor_line, output_file, middle_region):
    pybedtools.helpers.cleanup()

    ChIA_Drop = BedTool(path1)
    Region = BedTool(path2)

    # Use only the first line of the region file for the left and right references
    first_region = list(Region)[anchor_line]
    left_anchor = f"{first_region.chrom}\t{first_region.start}\t{first_region.end}"
    right_anchor = f"{first_region.fields[3]}\t{first_region.fields[4]}\t{first_region.fields[5]}"
    region = f"{first_region.chrom}\t{first_region.start}\t{first_region.fields[5]}"

    sort_start_time = time.time()

    if type == "left":
        ranked_gems = sort.process_left(ChIA_Drop, num_fragments, left_anchor, right_anchor, region)

    elif type == "right":
        ranked_gems = sort.process_right(ChIA_Drop, num_fragments, left_anchor, right_anchor, region)

    elif type == "both":
        ranked_gems = sort.process_both(ChIA_Drop, left_anchor, right_anchor, region)

    elif type == "middle":
        ranked_gems = sort.process_middle(ChIA_Drop, num_fragments, left_anchor, right_anchor, region, middle_region)

    elif type == "only-middle":
        ranked_gems = sort.process_only_middle(ChIA_Drop, middle_region)

    elif type == "only-middle-1frag":
        ranked_gems = sort.process_only_middle_1frag(ChIA_Drop, middle_region)

    else:
        print(f"Type '{type}' is not supported.")
        return

    print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")

    plot.plot_ranked_gems(start_time, ranked_gems, output_file, left_anchor, right_anchor, middle_region, type)

    histogram.generate_file(ranked_gems, output_file)


class ConditionalArgument(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values == 'middle' or values == "only-middle" or values == "only-middle-1frag":
            setattr(namespace, 'region_required', True)
        else:
            setattr(namespace, 'region_required', False)
        setattr(namespace, self.dest, values)


if __name__ == '__main__':
    start_time = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--path1', type=str, required=True,
                        help='Path to the GEMs file (.PEanno)')
    parser.add_argument('--path2', type=str, required=True,
                        help='Path to the regions file (.domains)')
    parser.add_argument('--type', type=str, required=True,
                        help='Type of processing: left, right, middle, only-middle, only-middle-1frag or both',
                        action=ConditionalArgument)
    parser.add_argument('--numfrag', type=str, required=True,
                        help='Number of fragments allowed in the region')
    parser.add_argument('--anchor', type=str, required=True,
                        help='The index of anchor line (0-based)')
    parser.add_argument('--output_file', type=str, required=True,
                        help='Output file for the plot')
    parser.add_argument('--region', type=str,
                        help='Specific region to process when type is middle, only-middle, or only-middle-1frag')

    args = parser.parse_args()

    if args.region_required and not args.region:
        parser.error("--region is required when --type is 'middle, only-middle or only-middle-1frag'")

    path1 = args.path1
    path2 = args.path2
    processing_type = args.type
    num_fragments = int(args.numfrag)
    anchor_line = int(args.anchor)
    output_file = args.output_file
    region = args.region \
        if processing_type == 'middle' \
        or processing_type == 'only-middle' \
        or processing_type == 'only-middle-1frag' \
        else None

    main(start_time, path1, path2, processing_type, num_fragments, anchor_line, output_file, region)

import pybedtools
from pybedtools import BedTool
import argparse
import time
import plot
import histogram
import sort
from helper import process_multiple_regions


def main(start_time, path1, path2, processing_type, num_fragments, anchor_line, output_file, region, operation):
    pybedtools.helpers.cleanup()

    ChIA_Drop = BedTool(path1)

    if processing_type != "multiple":
        Region = BedTool(path2)
        first_region = list(Region)[anchor_line]
        left_anchor = f"{first_region.chrom}\t{first_region.start}\t{first_region.end}"
        right_anchor = f"{first_region.fields[3]}\t{first_region.fields[4]}\t{first_region.fields[5]}"
        filter_region = f"{first_region.chrom}\t{first_region.start}\t{first_region.fields[5]}"

        sort_start_time = time.time()

        if processing_type == "left":
            ranked_gems = sort.process_left(ChIA_Drop, num_fragments, left_anchor, right_anchor, filter_region)
        elif processing_type == "right":
            ranked_gems = sort.process_right(ChIA_Drop, num_fragments, left_anchor, right_anchor, filter_region)
        elif processing_type == "both":
            ranked_gems = sort.process_both(ChIA_Drop, left_anchor, right_anchor, filter_region)
        elif processing_type == "middle":
            ranked_gems = sort.process_middle(ChIA_Drop, num_fragments,
                                              left_anchor, right_anchor, filter_region, region)
        elif processing_type == "only-middle":
            ranked_gems = sort.process_only_middle(ChIA_Drop, region)
        elif processing_type == "only-middle-1frag":
            ranked_gems = sort.process_only_middle_1frag(ChIA_Drop, region)
        else:
            print(f"Type '{processing_type}' is not supported.")
            return

        print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")
        plot.plot_ranked_gems_scaled(start_time, ranked_gems, output_file,
                              left_anchor, right_anchor, region, processing_type)
        histogram.generate_file(ranked_gems, output_file)

    else:
        if operation:
            sort_start_time = time.time()
            yes_chroms, no_chroms = process_multiple_regions(region, operation)
            # Process with the specified operation for multiple
            ranked_gems = sort.process_multiple(ChIA_Drop, num_fragments, yes_chroms, no_chroms)
            print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")
            plot.plot_ranked_gems_multiple_regions(start_time, ranked_gems, output_file, yes_chroms+no_chroms)
            histogram.generate_file(ranked_gems, output_file)

        else:
            print("Operation is required for type 'multiple'.")
            return


class ConditionalArgument(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        if values in ['middle', 'only-middle', 'only-middle-1frag', 'multiple']:
            setattr(namespace, 'region_required', True)
        else:
            setattr(namespace, 'region_required', False)
        if values == 'multiple':
            setattr(namespace, 'multiple_type', True)
        else:
            setattr(namespace, 'multiple_type', False)
        setattr(namespace, self.dest, values)


if __name__ == '__main__':
    start_time = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--path1', type=str, required=True,
                        help='Path to the GEMs file (.PEanno)')
    parser.add_argument('--path2', type=str,
                        help='Path to the regions file (.domains)')
    parser.add_argument('--type', type=str, required=True,
                        help='Type of processing: left, right, middle, only-middle, only-middle-1frag, both, or multiple',
                        action=ConditionalArgument)
    parser.add_argument('--numfrag', type=str,
                        help='Number of fragments allowed in the region')
    parser.add_argument('--anchor', type=str,
                        help='The index of anchor line (0-based)')
    parser.add_argument('--output_file', type=str,
                        help='Output file for the plot')
    parser.add_argument('--region', type=str,
                        help='Specific region to process when type is middle, only-middle, only-middle-1frag, or multiple')
    parser.add_argument('--operation', type=str,
                        help='Operation to perform when type is multiple')

    args = parser.parse_args()

    if args.multiple_type:
        if not args.operation or not args.region:
            parser.error("--operation and --region are required when --type is 'multiple'")
        if any([args.path2, args.anchor]):
            parser.error("--path2 and --anchor are invalid options when --type is 'multiple'")
    else:
        if args.region_required and not args.region:
            parser.error("--region is required when --type is 'middle, only-middle, only-middle-1frag, or multiple'")
        if args.type != "multiple" and args.operation:
            parser.error("--operation is not valid unless --type is 'multiple'")

    path1 = args.path1
    path2 = args.path2
    processing_type = args.type
    num_fragments = int(args.numfrag) if args.numfrag else None
    anchor_line = int(args.anchor) if args.anchor else None
    output_file = args.output_file
    region = args.region
    operation = args.operation

    main(start_time, path1, path2, processing_type, num_fragments, anchor_line, output_file, region, operation)

    """Redesign the command line."""



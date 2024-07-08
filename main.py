import pybedtools
from pybedtools import BedTool
import argparse
import time
import plot
import histogram
import sort
from helper import process_multiple_regions, process_graphs_arg, \
    create_filename, process_color_arg


def main(start_time, path1, path2, processing_type, graphs,
         num_fragments_min, num_fragments_max, region, operation,
         dataset, out_dir, colors, anchor_options):
    pybedtools.helpers.cleanup()

    ChIA_Drop = BedTool(path1)

    colors_flags = process_color_arg(colors)

    if processing_type != "multiple":  # abc processing
        Region = BedTool(path2)
        graphs_flags = process_graphs_arg(graphs)

        for anchors in Region:
            anchors = anchors.fields
            id = anchors[9]
            A = f"{anchors[0]}\t{anchors[1]}\t{anchors[2]}"
            B = f"{anchors[3]}\t{anchors[4]}\t{anchors[5]}"
            C = f"{anchors[6]}\t{anchors[7]}\t{anchors[8]}"

            # error check
            if anchors[1] >= anchors[2] or anchors[4] >= anchors[5] or anchors[7] >= anchors[8] \
            or anchors[2] >= anchors[4] or anchors[5] >= anchors[7]:
                print(f"Error for {id}: left is larger than right, please check the input file")
                continue

            filter_region = f"{anchors[0]}\t{anchors[1]}\t{anchors[8]}"
            # only do this once and then save it
            region_bed = BedTool(filter_region, from_string=True)
            ChIA_Drop_anchor = ChIA_Drop.intersect(region_bed, wa=True, wb=True)

            filter = f"{anchors[0]}\t{anchors[1]}\t{anchors[5]}"
            region_bed = BedTool(filter, from_string=True)
            ChIA_Drop_ab = ChIA_Drop_anchor.intersect(region_bed, wa=True, wb=True)

            filter = f"{anchors[0]}\t{anchors[4]}\t{anchors[8]}"
            region_bed = BedTool(filter, from_string=True)
            ChIA_Drop_bc = ChIA_Drop_anchor.intersect(region_bed, wa=True, wb=True)

            if graphs_flags["AtoB"]:
                ranked_gems = sort.process_left(ChIA_Drop_ab, num_fragments_min, num_fragments_max, A, B, filter_region)
                output_file = create_filename(dataset, id, "AtoB", num_fragments_min, num_fragments_max, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, B, C, out_dir,
                                             colors_flags, anchor_options, id, path1, "AtoB")
                histogram.generate_file(ranked_gems, output_file, out_dir)

            if graphs_flags["AtoC"]:
                ranked_gems = sort.process_left(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, A, C, filter_region)
                output_file = create_filename(dataset, id, "AtoC", num_fragments_min, num_fragments_max, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, C, B, out_dir,
                                             colors_flags, anchor_options, id, path1, "AtoC")
                histogram.generate_file(ranked_gems, output_file, out_dir)

            if graphs_flags["BtoA"]:
                ranked_gems = sort.process_right(ChIA_Drop_ab, num_fragments_min, num_fragments_max, A, B, filter_region)
                output_file = create_filename(dataset, id, "BtoA", num_fragments_min, num_fragments_max, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, B, C, out_dir,
                                             colors_flags, anchor_options, id, path1, "BtoA")
                histogram.generate_file(ranked_gems, output_file, out_dir)

            if graphs_flags["BtoC"]:
                ranked_gems = sort.process_left(ChIA_Drop_bc, num_fragments_min, num_fragments_max, B, C, filter_region)
                output_file = create_filename(dataset, id, "BtoC", num_fragments_min, num_fragments_max, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, B, C, A, out_dir,
                                             colors_flags, anchor_options, id, path1, "BtoC")
                histogram.generate_file(ranked_gems, output_file, out_dir)

            if graphs_flags["CtoA"]:
                ranked_gems = sort.process_right(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, A, C, filter_region)
                output_file = create_filename(dataset, id, "CtoA", num_fragments_min, num_fragments_max, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, C, B, out_dir,
                                             colors_flags, anchor_options, id, path1, "CtoA")
                histogram.generate_file(ranked_gems, output_file, out_dir)

            if graphs_flags["CtoB"]:
                ranked_gems = sort.process_right(ChIA_Drop_bc, num_fragments_min, num_fragments_max, B, C, filter_region)
                output_file = create_filename(dataset, id, "CtoB", num_fragments_min, num_fragments_max, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, B, C, A, out_dir,
                                             colors_flags, anchor_options, id, path1, "CtoB")
                histogram.generate_file(ranked_gems, output_file, out_dir)

            if graphs_flags["AandC"]:
                region = f"{anchors[0]}:{anchors[1]}-{anchors[2]};{anchors[6]}:{anchors[7]}-{anchors[8]}"
                yes_chroms, no_chroms = process_multiple_regions(region, "yes;yes")
                ranked_gems = sort.process_multiple(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, yes_chroms, no_chroms)
                output_file = create_filename(dataset, id, "AandC", num_fragments_min, num_fragments_max, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, C, B, out_dir,
                                             colors_flags, anchor_options, id, path1, "AandC")
                histogram.generate_file(ranked_gems, output_file, out_dir)

            if graphs_flags["Bcentered"]:
                ranked_gems = sort.process_middle(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, A, C, filter_region, B)
                output_file = create_filename(dataset, id, "Bcentered", num_fragments_min, num_fragments_max, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, C, B, out_dir,
                                             colors_flags, anchor_options, id, path1, "Bcentered")
                histogram.generate_file(ranked_gems, output_file, out_dir)

    else:
        if operation:
            sort_start_time = time.time()
            yes_chroms, no_chroms = process_multiple_regions(region, operation)
            # Process with the specified operation for multiple
            ranked_gems = sort.process_multiple(ChIA_Drop, num_fragments_min, num_fragments_max, yes_chroms, no_chroms)
            print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")
            plot.plot_ranked_gems_multiple_regions(start_time, ranked_gems, "output_file", yes_chroms+no_chroms)  # TODO: revise file name
            histogram.generate_file(ranked_gems, "output_file")  # TODO: revise file name

        else:
            print("Operation is required for type 'multiple'.")
            return


if __name__ == '__main__':
    start_time = time.time()

    parser = argparse.ArgumentParser()
    parser.add_argument('--path1', type=str, required=True,
                        help='Path to the GEMs file (.PEanno)')
    parser.add_argument('--path2', type=str,
                        help='Path to the regions file')
    parser.add_argument('--type', type=str, required=True,
                        help='Type of processing: abc or multiple')
    parser.add_argument('--graphs', type=str, required=True,
                        help='Graphs to genereate when type is abc')
    parser.add_argument('--numfrag_min', type=str, default="2",
                        help='Minimum number of fragments allowed in the region')
    parser.add_argument('--numfrag_max', type=str, default="1000",
                        help='Maximum number of fragments allowed in the region')
    parser.add_argument('--region', type=str,
                        help='Specific region to process when type is middle, only-middle, only-middle-1frag, or multiple')
    parser.add_argument('--operation', type=str,
                        help='Operation to perform when type is multiple')
    parser.add_argument('--out_dir', type=str, default="/",
                        help='The output directory name that the output files will be put in')
    parser.add_argument('--colors', type=str, default="red;green;#525252",
                        help='The color of the anchors, fragments and lines (in order, seperated by semicolons)')
    parser.add_argument('--anchor_options', type=str, default="no",
                        help='Three options: yes_complete (draw anchors complete on the graph), yes_top (draw anchors on the top of the graph) and no (do not draw anchors)')

    args = parser.parse_args()

    path1 = args.path1
    dataset = path1.split("-")[0]
    path2 = args.path2
    processing_type = args.type
    graphs = args.graphs
    num_fragments_min = int(args.numfrag_min)
    num_fragments_max = int(args.numfrag_max)
    region = args.region
    operation = args.operation
    out_dir = args.out_dir
    colors = args.colors
    anchor_options = args.anchor_options

    main(start_time, path1, path2, processing_type, graphs, num_fragments_min, num_fragments_max,
         region, operation, dataset, out_dir, colors, anchor_options)

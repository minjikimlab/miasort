import pybedtools
from pybedtools import BedTool
import argparse
import time
import plot
import histogram
import sort
from helper import process_multiple_regions, process_graphs_arg, \
    create_filename


def main(start_time, path1, path2, processing_type, graphs,
         num_fragments, region, operation, dataset):
    pybedtools.helpers.cleanup()

    ChIA_Drop = BedTool(path1)

    if processing_type != "multiple":  # abc processing
        Region = BedTool(path2)
        graphs_flags = process_graphs_arg(graphs)

        for anchors in Region:
            anchors = anchors.fields
            id = anchors[9]
            A = f"{anchors[0]}\t{anchors[1]}\t{anchors[2]}"
            B = f"{anchors[3]}\t{anchors[4]}\t{anchors[5]}"
            C = f"{anchors[6]}\t{anchors[7]}\t{anchors[8]}"
            filter_region = f"{anchors[0]}\t{anchors[1]}\t{anchors[8]}"

            # only do this once and then save it
            region_bed = BedTool(filter_region, from_string=True)
            ChIA_Drop = ChIA_Drop.intersect(region_bed, wa=True, wb=True)

            if graphs_flags["AtoB"]:
                ranked_gems = sort.process_left(ChIA_Drop, num_fragments, A, B, filter_region)
                output_file = create_filename(dataset, id, "AtoB", num_fragments, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, B, C)
                histogram.generate_file(ranked_gems, output_file)

            if graphs_flags["AtoC"]:
                ranked_gems = sort.process_left(ChIA_Drop, num_fragments, A, C, filter_region)
                output_file = create_filename(dataset, id, "AtoC", num_fragments, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, C, B)
                histogram.generate_file(ranked_gems, output_file)

            if graphs_flags["BtoA"]:
                ranked_gems = sort.process_right(ChIA_Drop, num_fragments, A, B, filter_region)
                output_file = create_filename(dataset, id, "BtoA", num_fragments, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, B, C)
                histogram.generate_file(ranked_gems, output_file)

            if graphs_flags["BtoC"]:
                ranked_gems = sort.process_left(ChIA_Drop, num_fragments, B, C, filter_region)
                output_file = create_filename(dataset, id, "BtoC", num_fragments, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, B, C, A)
                histogram.generate_file(ranked_gems, output_file)

            if graphs_flags["CtoA"]:
                ranked_gems = sort.process_right(ChIA_Drop, num_fragments, A, C, filter_region)
                output_file = create_filename(dataset, id, "CtoA", num_fragments, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, C, B)
                histogram.generate_file(ranked_gems, output_file)

            if graphs_flags["CtoB"]:
                ranked_gems = sort.process_right(ChIA_Drop, num_fragments, B, C, filter_region)
                output_file = create_filename(dataset, id, "CtoB", num_fragments, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, B, C, A)
                histogram.generate_file(ranked_gems, output_file)

            if graphs_flags["AandC"]:
                region = f"{anchors[0]}:{anchors[1]}-{anchors[2]};{anchors[6]}:{anchors[7]}-{anchors[8]}"
                yes_chroms, no_chroms = process_multiple_regions(region, "yes;yes")
                ranked_gems = sort.process_multiple(ChIA_Drop, num_fragments, yes_chroms, no_chroms)
                output_file = create_filename(dataset, id, "AandC", num_fragments, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, C, B)
                histogram.generate_file(ranked_gems, output_file)

            if graphs_flags["Bcentered"]:
                ranked_gems = sort.process_middle(ChIA_Drop, num_fragments,
                                                    A, C, filter_region, B)
                output_file = create_filename(dataset, id, "Bcentered", num_fragments, len(ranked_gems))
                plot.plot_ranked_gems_scaled(ranked_gems, output_file, A, C, B)
                histogram.generate_file(ranked_gems, output_file)


    else:
        if operation:
            sort_start_time = time.time()
            yes_chroms, no_chroms = process_multiple_regions(region, operation)
            # Process with the specified operation for multiple
            ranked_gems = sort.process_multiple(ChIA_Drop, num_fragments, yes_chroms, no_chroms)
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
    parser.add_argument('--numfrag', type=str,
                        help='Number of fragments allowed in the region')
    parser.add_argument('--region', type=str,
                        help='Specific region to process when type is middle, only-middle, only-middle-1frag, or multiple')
    parser.add_argument('--operation', type=str,
                        help='Operation to perform when type is multiple')

    args = parser.parse_args()

    path1 = args.path1
    dataset = path1.split("-")[0]
    path2 = args.path2
    processing_type = args.type
    graphs = args.graphs
    num_fragments = int(args.numfrag) if args.numfrag else None
    region = args.region
    operation = args.operation

    main(start_time, path1, path2, processing_type, graphs, num_fragments, region, operation, dataset)


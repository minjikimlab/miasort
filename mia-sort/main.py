import pybedtools
from pybedtools import BedTool
import pyranges as pr
import argparse
import time
import plot
import histogram
import sort
import records
import os
import shutil
import csv
from helper import process_multiple_regions, process_graphs_arg, \
    create_plot_filename, process_color_arg, create_histogram_filename, \
    create_csv_filename, generate_filter_regions


def main(start_time, path1, path2, processing_type, graphs,
         num_fragments_min, num_fragments_max, region, operation,
         dataset, out_dir, colors, anchor_options, intersection_options,
         graph_flag, extension, histogram_options):
    pybedtools.helpers.cleanup()

    ChIA_Drop = BedTool(path1)

    colors_flags = process_color_arg(colors)

    # delete the out_dir folder if it exists
    if out_dir != "/" and os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    if processing_type != "multiple":  # abc processing
        filter_regions_filename = "filter_regions.bed"
        generate_filter_regions(path2, filter_regions_filename)
        filter_regions = BedTool(filter_regions_filename)

        start = time.time()
        intersected = ChIA_Drop.intersect(filter_regions, wa=True, wb=True)
        print(f"Intersect a and b: {time.time() - start} secs")

        # Dictionary to store the intersected regions for each line of b
        filtered_intersections = {}

        start = time.time()
        # Iterate over each intersection
        for intersection in intersected:
            # Extract the fields of the intersected line from b
            if path1[:3] == 'LHG' or "SPRITE" in path1 or "4DNFIACOTIGL" in path1 or "ChIA-Drop" in path1:
                # print("Assume 5 fields in region file")
                b_fields = intersection.fields[5:]  # 5 fields in a
            elif "PoreC" in path1:
                # print("Assume 4 fields in region file")
                b_fields = intersection.fields[4:]  # 4 fields in a
            else:
                # print("Assume 6 fields in region file")
                b_fields = intersection.fields[6:]  # 6 fields in a
            b_fields = ' '.join(b_fields)  # Make the key hashable
            # Check if the key exists, if not, add an empty list
            if b_fields not in filtered_intersections:
                filtered_intersections[b_fields] = []
            # Append the intersection to the list
            filtered_intersections[b_fields].append(intersection)
        # Convert lists to BedTool objects
        for i in filtered_intersections:
            filtered_intersections[i] = BedTool(filtered_intersections[i])
        print(f"Process intersected regions: {time.time() - start} secs\n-------------------------------------\n")

        graphs_flags = process_graphs_arg(graphs)

        csv_file = create_csv_filename(dataset, path2)
        if out_dir != "/":
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            path = os.path.join(out_dir, csv_file)
        else:
            path = csv_file
        # Write the header of the comp records file
        with open(path, 'a', newline='') as file:
            writer = csv.writer(file)
            field = ["Region ID", "A", "B", "C", "Region", "Sort Scheme",
                     "num_complexes", "num_1frag", "num_2frag", "num_3frag", "num_4frag", "num>=5frag"]
            writer.writerow(field)

        for key, ChIA_Drop_anchor in filtered_intersections.items():
            anchors = key.split(" ")[3:]
            id = anchors[9]

            # Error check
            if anchors[1] >= anchors[2] or anchors[4] >= anchors[5] or anchors[7] >= anchors[8] \
            or anchors[2] >= anchors[4] or anchors[5] >= anchors[7]:
                print(f"Error for {id}: left is larger than right, please check the input file")
                continue

            A = f"{anchors[0]}\t{anchors[1]}\t{anchors[2]}"
            B = f"{anchors[3]}\t{anchors[4]}\t{anchors[5]}"
            C = f"{anchors[6]}\t{anchors[7]}\t{anchors[8]}"
            filter_region = f"{anchors[0]}\t{anchors[1]}\t{anchors[8]}"

            filter = f"{anchors[0]}\t{anchors[1]}\t{anchors[5]}"
            region_bed = BedTool(filter, from_string=True)
            ChIA_Drop_ab = ChIA_Drop_anchor.intersect(region_bed, wa=True, wb=True)

            filter = f"{anchors[0]}\t{anchors[4]}\t{anchors[8]}"
            region_bed = BedTool(filter, from_string=True)
            ChIA_Drop_bc = ChIA_Drop_anchor.intersect(region_bed, wa=True, wb=True)

            # ---------------- First plot ----------------

            ranked_gems_list = []
            left_anchor_list = []
            right_anchor_list = []
            middle_anchor_list = []
            commands_list = ["AtoC", "CtoA", "AandC"]

            if graphs_flags["AtoC"]:
                start = time.time()
                ranked_gems = sort.process_left(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, A, C, filter_region)
                output_file = create_plot_filename(dataset, id, "AtoC", num_fragments_min, num_fragments_max, len(ranked_gems))
                ranked_gems_list.append(ranked_gems)
                left_anchor_list.append(A)
                right_anchor_list.append(C)
                middle_anchor_list.append(B)
                if histogram_options == "yes":
                    histogram.generate_file(ranked_gems, output_file, out_dir)
                records.generate_file(id, A, B, C, "AtoC", len(ranked_gems), csv_file, out_dir, ranked_gems)
                print(f"AtoC: {time.time() - start}")

            if graphs_flags["CtoA"]:
                start = time.time()
                ranked_gems = sort.process_right(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, A, C, filter_region)
                output_file = create_plot_filename(dataset, id, "CtoA", num_fragments_min, num_fragments_max, len(ranked_gems))
                ranked_gems_list.append(ranked_gems)
                left_anchor_list.append(A)
                right_anchor_list.append(C)
                middle_anchor_list.append(B)
                if histogram_options == "yes":
                    histogram.generate_file(ranked_gems, output_file, out_dir)
                records.generate_file(id, A, B, C, "CtoA", len(ranked_gems), csv_file, out_dir, ranked_gems)
                print(f"CtoA: {time.time() - start}")

            if graphs_flags["AandC"]:
                start = time.time()
                region = f"{anchors[0]}:{anchors[1]}-{anchors[2]};{anchors[6]}:{anchors[7]}-{anchors[8]}"
                yes_chroms, no_chroms = process_multiple_regions(region, "yes;yes")
                ranked_gems = sort.process_multiple(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, yes_chroms, no_chroms)
                output_file = create_plot_filename(dataset, id, "AandC", num_fragments_min, num_fragments_max, len(ranked_gems))
                ranked_gems_list.append(ranked_gems)
                left_anchor_list.append(A)
                right_anchor_list.append(C)
                middle_anchor_list.append(B)
                if histogram_options == "yes":
                    histogram.generate_file(ranked_gems, output_file, out_dir)
                records.generate_file(id, A, B, C, "AandC", len(ranked_gems), csv_file, out_dir, ranked_gems)
                print(f"AandC: {time.time() - start}")

            if graph_flag == "yes":
                output_file = create_plot_filename(dataset, id, "stripes", num_fragments_min, num_fragments_max, len(ranked_gems))
                plot.plot_ranked_gems(ranked_gems_list, output_file, left_anchor_list,
                                            right_anchor_list, middle_anchor_list, out_dir,
                                            colors_flags, anchor_options, id, dataset, commands_list, extension)

            # ---------------- Second plot ----------------

            ranked_gems_list = []
            left_anchor_list = []
            right_anchor_list = []
            middle_anchor_list = []
            commands_list = ["BtoAC", "BtoA", "BtoC"]

            if graphs_flags["Bcentered"]:  # B centered to A & C
                start = time.time()
                ranked_gems = sort.process_middle(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, A, C, filter_region, B)
                output_file = create_plot_filename(dataset, id, "Bcentered", num_fragments_min, num_fragments_max, len(ranked_gems))
                ranked_gems_list.append(ranked_gems)
                left_anchor_list.append(A)
                right_anchor_list.append(C)
                middle_anchor_list.append(B)
                if histogram_options == "yes":
                    histogram.generate_file(ranked_gems, output_file, out_dir)
                records.generate_file(id, A, B, C, "BtoAC", len(ranked_gems), csv_file, out_dir, ranked_gems)
                print(f"Bcentered: {time.time() - start}")

            if graphs_flags["BtoA"]:
                start = time.time()
                ranked_gems = sort.process_right(ChIA_Drop_ab, num_fragments_min, num_fragments_max, A, B, filter_region)
                output_file = create_plot_filename(dataset, id, "BtoA", num_fragments_min, num_fragments_max, len(ranked_gems))
                ranked_gems_list.append(ranked_gems)
                left_anchor_list.append(A)
                right_anchor_list.append(B)
                middle_anchor_list.append(C)
                if histogram_options == "yes":
                    histogram.generate_file(ranked_gems, output_file, out_dir)
                records.generate_file(id, A, B, C, "BtoA", len(ranked_gems), csv_file, out_dir, ranked_gems)
                print(f"BtoA: {time.time() - start}")

            if graphs_flags["BtoC"]:
                start = time.time()
                ranked_gems = sort.process_left(ChIA_Drop_bc, num_fragments_min, num_fragments_max, B, C, filter_region)
                output_file = create_plot_filename(dataset, id, "BtoC", num_fragments_min, num_fragments_max, len(ranked_gems))
                ranked_gems_list.append(ranked_gems)
                left_anchor_list.append(B)
                right_anchor_list.append(C)
                middle_anchor_list.append(A)
                if histogram_options == "yes":
                    histogram.generate_file(ranked_gems, output_file, out_dir)
                records.generate_file(id, A, B, C, "BtoC", len(ranked_gems), csv_file, out_dir, ranked_gems)
                print(f"BtoC: {time.time() - start}")

            if graph_flag == "yes":
                output_file = create_plot_filename(dataset, id, "jets", num_fragments_min, num_fragments_max, len(ranked_gems))
                plot.plot_ranked_gems(ranked_gems_list, output_file, left_anchor_list,
                                            right_anchor_list, middle_anchor_list, out_dir,
                                            colors_flags, anchor_options, id, dataset, commands_list, extension)

            print("-------------------------------------")

    else:
        if out_dir != "/" and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        yes_chroms, no_chroms = process_multiple_regions(region, operation)
        # Process with the specified operation for multiple
        ranked_gems = sort.process_multiple(ChIA_Drop, num_fragments_min, num_fragments_max, yes_chroms, no_chroms)
        output_file = create_plot_filename(dataset, None, "multiple", num_fragments_min,
                                           num_fragments_max, len(ranked_gems))
        plot.plot_ranked_gems([ranked_gems], output_file, [""], [""], [""], out_dir,
                                    colors_flags, anchor_options, 0, dataset, ["multiple"],
                                    extension, flag="multiple", regions=yes_chroms+no_chroms)
        if histogram_options == "yes":
            histogram.generate_file(ranked_gems, "output_file", out_dir)  # TODO: revise file name


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
    parser.add_argument('--intersection_options', type=str, default="bedtools",
                        help='Two options: bedtools or rust')
    parser.add_argument('--graph_flag', type=str, default="yes",
                        help='Special graph flag: yes or no - whether or not the toolkit should generate plots')
    parser.add_argument('--extension', type=str, default="6000",
                        help='Three options: default width=6000b, user-specified width, or natural width')
    parser.add_argument('--historgram_options', type=str, default="no",
                        help='Whether or not the program should draw histograms seperately: yes or no')

    args = parser.parse_args()

    path1 = args.path1
    dataset = path1.split("/")[-1].split(".region")[0]
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
    intersection_options = args.intersection_options
    graph_flag = args.graph_flag
    extension = args.extension
    histogram_options = args.historgram_options

    print("Finish command line parsing.")

    main(start_time, path1, path2, processing_type, graphs, num_fragments_min, num_fragments_max,
         region, operation, dataset, out_dir, colors, anchor_options, intersection_options, graph_flag,
         extension, histogram_options)

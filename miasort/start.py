import pybedtools
from pybedtools import BedTool
import os
import shutil
import csv

from .plot import plot_ranked_gems
from .histogram import generate_file
from .sort import *
from .records import write_to_csv_file, write_to_csv_file_multiple
from .helper import process_multiple_regions, process_graphs_arg, \
    create_plot_filename, process_color_arg, \
    create_csv_filename, generate_filter_regions

def start(path1, path2, processing_type, graphs,
         num_fragments_min, num_fragments_max, region, operation,
         dataset, out_dir, colors, anchor_options,
         graph_flag, extension, histogram_options, frag_height,
         line_width, plot_width, subplots_margins):
    pybedtools.helpers.cleanup()

    ChIA_Drop = BedTool(path1)

    colors_flags = process_color_arg(colors)

    frag_description = ""
    if extension == "natural":
        frag_description = "fragexact"
    else:
        frag_description = f"frag{extension}bp"

    # delete the out_dir folder if it exists
    if out_dir != "/" and os.path.exists(out_dir):
        shutil.rmtree(out_dir)

    if processing_type == "abc":
        filter_regions_filename = "filter_regions.bed"
        generate_filter_regions(path2, filter_regions_filename)
        filter_regions = BedTool(filter_regions_filename)

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

        graphs_flags = process_graphs_arg(graphs)

        csv_file = create_csv_filename(dataset, path2, frag_description)
        if out_dir != "/":
            if not os.path.exists(out_dir):
                os.makedirs(out_dir)
            path = os.path.join(out_dir, csv_file)
        else:
            path = csv_file
        # Write the header of the comp records file
        with open(path, 'a', newline='') as file:
            writer = csv.writer(file)
            field = ["Region ID", "Site A", "Site B", "Site C", "Region", "Sort Scheme",
                     "num_complexes", "num_1frag", "num_2frag", "num_3frag",
                     "num_4frag", "num>=5frag"]
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

            if graphs_flags["AtoC"] or graphs_flags["CtoA"] or graphs_flags["AandC"]:
                if graphs_flags["AtoC"]:
                    ranked_gems = process_left(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, A, C, filter_region)
                    output_file = create_plot_filename(dataset, id, "AtoC", num_fragments_min, num_fragments_max, len(ranked_gems), frag_description)
                    ranked_gems_list.append(ranked_gems)
                    left_anchor_list.append(A)
                    right_anchor_list.append(C)
                    middle_anchor_list.append(B)
                    if histogram_options == "yes":
                        generate_file(ranked_gems, output_file, out_dir)
                    write_to_csv_file(id, A, B, C, "AtoC", len(ranked_gems), csv_file, out_dir, ranked_gems)

                if graphs_flags["CtoA"]:
                    ranked_gems = process_right(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, A, C, filter_region)
                    output_file = create_plot_filename(dataset, id, "CtoA", num_fragments_min, num_fragments_max, len(ranked_gems), frag_description)
                    ranked_gems_list.append(ranked_gems)
                    left_anchor_list.append(A)
                    right_anchor_list.append(C)
                    middle_anchor_list.append(B)
                    if histogram_options == "yes":
                        generate_file(ranked_gems, output_file, out_dir)
                    write_to_csv_file(id, A, B, C, "CtoA", len(ranked_gems), csv_file, out_dir, ranked_gems)

                if graphs_flags["AandC"]:
                    region = f"{anchors[0]}:{anchors[1]}-{anchors[2]};{anchors[6]}:{anchors[7]}-{anchors[8]}"
                    yes_chroms, no_chroms = process_multiple_regions(region, "yes;yes")
                    ranked_gems = process_multiple(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, yes_chroms, no_chroms)
                    output_file = create_plot_filename(dataset, id, "AandC", num_fragments_min, num_fragments_max, len(ranked_gems), frag_description)
                    ranked_gems_list.append(ranked_gems)
                    left_anchor_list.append(A)
                    right_anchor_list.append(C)
                    middle_anchor_list.append(B)
                    if histogram_options == "yes":
                        generate_file(ranked_gems, output_file, out_dir)
                    write_to_csv_file(id, A, B, C, "AandC", len(ranked_gems), csv_file, out_dir, ranked_gems)

                if graph_flag == "yes":
                    output_file = create_plot_filename(dataset, id, "stripes", num_fragments_min, num_fragments_max,
                                                    len(ranked_gems), frag_description)
                    plot_ranked_gems(ranked_gems_list, output_file, left_anchor_list,
                                                right_anchor_list, middle_anchor_list, out_dir,
                                                colors_flags, anchor_options, id, dataset, commands_list, extension,
                                                frag_height, line_width, plot_width, subplots_margins, frag_description)

            # ---------------- Second plot ----------------

            ranked_gems_list = []
            left_anchor_list = []
            right_anchor_list = []
            middle_anchor_list = []
            commands_list = ["BtoAC", "BtoA", "BtoC"]

            if graphs_flags["Bcentered"] or graphs_flags["BtoA"] or graphs_flags["BtoC"]:
                if graphs_flags["Bcentered"]:  # B centered to A & C
                    ranked_gems = process_middle(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, A, C, filter_region, B)
                    output_file = create_plot_filename(dataset, id, "Bcentered", num_fragments_min, num_fragments_max, len(ranked_gems), frag_description)
                    ranked_gems_list.append(ranked_gems)
                    left_anchor_list.append(A)
                    right_anchor_list.append(C)
                    middle_anchor_list.append(B)
                    if histogram_options == "yes":
                        generate_file(ranked_gems, output_file, out_dir)
                    write_to_csv_file(id, A, B, C, "BtoAC", len(ranked_gems), csv_file, out_dir, ranked_gems)

                if graphs_flags["BtoA"]:
                    ranked_gems = process_right(ChIA_Drop_ab, num_fragments_min, num_fragments_max, A, B, filter_region)
                    output_file = create_plot_filename(dataset, id, "BtoA", num_fragments_min, num_fragments_max, len(ranked_gems), frag_description)
                    ranked_gems_list.append(ranked_gems)
                    left_anchor_list.append(A)
                    right_anchor_list.append(B)
                    middle_anchor_list.append(C)
                    if histogram_options == "yes":
                        generate_file(ranked_gems, output_file, out_dir)
                    write_to_csv_file(id, A, B, C, "BtoA", len(ranked_gems), csv_file, out_dir, ranked_gems)

                if graphs_flags["BtoC"]:
                    ranked_gems = process_left(ChIA_Drop_bc, num_fragments_min, num_fragments_max, B, C, filter_region)
                    output_file = create_plot_filename(dataset, id, "BtoC", num_fragments_min, num_fragments_max, len(ranked_gems), frag_description)
                    ranked_gems_list.append(ranked_gems)
                    left_anchor_list.append(B)
                    right_anchor_list.append(C)
                    middle_anchor_list.append(A)
                    if histogram_options == "yes":
                        generate_file(ranked_gems, output_file, out_dir)
                    write_to_csv_file(id, A, B, C, "BtoC", len(ranked_gems), csv_file, out_dir, ranked_gems)

                if graph_flag == "yes":
                    output_file = create_plot_filename(dataset, id, "jets", num_fragments_min, num_fragments_max,
                                                    len(ranked_gems), frag_description)
                    plot_ranked_gems(ranked_gems_list, output_file, left_anchor_list,
                                                right_anchor_list, middle_anchor_list, out_dir,
                                                colors_flags, anchor_options, id, dataset, commands_list, extension,
                                                frag_height, line_width, plot_width, subplots_margins, frag_description)

    elif processing_type == "AandBandC":
        filter_regions_filename = "filter_regions.bed"
        generate_filter_regions(path2, filter_regions_filename)
        filter_regions = BedTool(filter_regions_filename)

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

        if out_dir != "/" and not os.path.exists(out_dir):
            os.makedirs(out_dir)

        csv_file = create_csv_filename(dataset, path2, frag_description)
        path = os.path.join(out_dir, csv_file)
        # Write the header of the comp records file
        with open(path, 'a', newline='') as file:
            writer = csv.writer(file)
            field = ["Region ID", "Site A", "Site B", "Site C", "Region", "Sort Scheme",
                     "num_complexes", "num_1frag", "num_2frag", "num_3frag",
                     "num_4frag", "num>=5frag"]
            writer.writerow(field)

        for key, ChIA_Drop_anchor in filtered_intersections.items():
            region = key.split(" ")[3:]
            id = region[9]

            # Error check
            if region[1] >= region[2] or region[4] >= region[5] or region[7] >= region[8] \
            or region[2] >= region[4] or region[5] >= region[7]:
                print(f"Error for {id}: left is larger than right, please check the input file")
                continue

            A = f"{region[0]}:{region[1]}-{region[2]}"
            B = f"{region[3]}:{region[4]}-{region[5]}"
            C = f"{region[6]}:{region[7]}-{region[8]}"
            r = f"{region[0]}:{region[1]}-{region[8]}"
            region = f"{A};{B};{C}"

            if operation != "yes;yes;yes":
                print("The operation is automatically set as `yes;yes;yes`")
            operation = "yes;yes;yes"

            yes_chroms, no_chroms = process_multiple_regions(region, operation)
            ranked_gems = process_multiple(ChIA_Drop_anchor, num_fragments_min, num_fragments_max, yes_chroms, no_chroms)
            output_file = create_plot_filename(dataset, id, "AandBandC", num_fragments_min,
                                            num_fragments_max, len(ranked_gems), frag_description)
            plot_ranked_gems([ranked_gems], output_file, [""], [""], [""], out_dir,
                                        colors_flags, anchor_options, id, dataset, ["multiple"],
                                        extension, frag_height, line_width, plot_width, subplots_margins, frag_description,
                                        flag="multiple_abc", regions=yes_chroms+no_chroms)
            if histogram_options == "yes":
                generate_file(ranked_gems, "output_file", out_dir)  # TODO: revise file name
            write_to_csv_file_multiple(id, A, B, C, r, "AandBandC", len(ranked_gems), csv_file, out_dir, ranked_gems)

    else:
        if out_dir != "/" and not os.path.exists(out_dir):
            os.makedirs(out_dir)
        yes_chroms, no_chroms = process_multiple_regions(region, operation)
        ranked_gems = process_multiple(ChIA_Drop, num_fragments_min, num_fragments_max, yes_chroms, no_chroms)
        output_file = create_plot_filename(dataset, None, "multiple", num_fragments_min,
                                           num_fragments_max, len(ranked_gems), frag_description)
        plot_ranked_gems([ranked_gems], output_file, [""], [""], [""], out_dir,
                                    colors_flags, anchor_options, 0, dataset, ["multiple"],
                                    extension, frag_height, line_width, plot_width, subplots_margins, frag_description,
                                    flag="multiple", regions=yes_chroms+no_chroms)
        if histogram_options == "yes":
            generate_file(ranked_gems, "output_file", out_dir)  # TODO: revise file name

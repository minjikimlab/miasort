from .start import start


def abc_sort(path1, path2, graphs, out_dir='/', plot=True, histogram=False, anchor_option='no',
             colors='red;green;#525252', num_frag_min=2, num_frag_max=1000, extension='6000'):
    """Sort Three Regions."""
    dataset = path1.split("/")[-1].split(".region")[0]

    if plot:
        graph_flag = "yes"
    else:
        graph_flag = "no"

    if histogram:
        histogram_options = "yes"
    else:
        histogram_options = "no"

    start(path1, path2, "abc", graphs, num_frag_min, num_frag_max,
         "", "", dataset, out_dir, colors, anchor_option, graph_flag,
         extension, histogram_options)


def multiple_sort(path1, path2, out_dir='/', plot=True, histogram=False, anchor_option='no',
                    colors='red;green;#525252', num_frag_min=2, num_frag_max=1000, extension='6000'):
    """Sort with A and B and C."""
    dataset = path1.split("/")[-1].split(".region")[0]

    if plot:
        graph_flag = "yes"
    else:
        graph_flag = "no"

    if histogram:
        histogram_options = "yes"
    else:
        histogram_options = "no"

    start(path1, path2, "AandBandC", "", num_frag_min, num_frag_max,
         "", "", dataset, out_dir, colors, anchor_option, graph_flag,
         extension, histogram_options)


def unlimited_multiple_sort(path1, regions, operations, out_dir='/', plot=True, histogram=False, anchor_option='no',
                            colors='red;green;#525252', num_frag_min=2, num_frag_max=1000, extension='6000'):
    """Sort with A and B and C."""
    dataset = path1.split("/")[-1].split(".region")[0]

    if plot:
        graph_flag = "yes"
    else:
        graph_flag = "no"

    if histogram:
        histogram_options = "yes"
    else:
        histogram_options = "no"

    start(path1, "", "multiple", "", num_frag_min, num_frag_max,
         regions, operations, dataset, out_dir, colors, anchor_option, graph_flag,
         extension, histogram_options)

# if __name__ == '__main__':
#     start_time = time.time()

#     parser = argparse.ArgumentParser()
#     parser.add_argument('--path1', type=str, required=True,
#                         help='Path to the GEMs file (.PEanno)')
#     parser.add_argument('--path2', type=str,
#                         help='Path to the regions file')
#     parser.add_argument('--type', type=str, required=True,
#                         help='Type of processing: abc, AandBandC, or random_multiple')
#     parser.add_argument('--graphs', type=str, required=True,
#                         help='Graphs to genereate when type is abc')
#     parser.add_argument('--numfrag_min', type=str, default="2",
#                         help='Minimum number of fragments allowed in the region')
#     parser.add_argument('--numfrag_max', type=str, default="1000",
#                         help='Maximum number of fragments allowed in the region')
#     parser.add_argument('--region', type=str,
#                         help='Specific region to process when type is middle, only-middle, only-middle-1frag, or multiple')
#     parser.add_argument('--operation', type=str,
#                         help='Operation to perform when type is multiple')
#     parser.add_argument('--out_dir', type=str, default="/",
#                         help='The output directory name that the output files will be put in')
#     parser.add_argument('--colors', type=str, default="red;green;#525252",
#                         help='The color of the anchors, fragments and lines (in order, seperated by semicolons)')
#     parser.add_argument('--anchor_options', type=str, default="no",
#                         help='Three options: yes_complete (draw anchors complete on the graph), yes_top (draw anchors on the top of the graph) and no (do not draw anchors)')
#     parser.add_argument('--intersection_options', type=str, default="bedtools",
#                         help='Two options: bedtools or rust')
#     parser.add_argument('--graph_flag', type=str, default="yes",
#                         help='Special graph flag: yes or no - whether or not the toolkit should generate plots')
#     parser.add_argument('--extension', type=str, default="6000",
#                         help='Three options: default width=6000b, user-specified width, or natural width')
#     parser.add_argument('--historgram_options', type=str, default="no",
#                         help='Whether or not the program should draw histograms seperately: yes or no')

#     args = parser.parse_args()

#     path1 = args.path1
#     dataset = path1.split("/")[-1].split(".region")[0]
#     path2 = args.path2
#     processing_type = args.type
#     graphs = args.graphs
#     num_fragments_min = int(args.numfrag_min)
#     num_fragments_max = int(args.numfrag_max)
#     region = args.region
#     operation = args.operation
#     out_dir = args.out_dir
#     colors = args.colors
#     anchor_options = args.anchor_options
#     graph_flag = args.graph_flag
#     extension = args.extension
#     histogram_options = args.historgram_options

#     print("Finish command line parsing.")

#     main(path1, path2, processing_type, graphs, num_fragments_min, num_fragments_max,
#          region, operation, dataset, out_dir, colors, anchor_options, graph_flag,
#          extension, histogram_options)

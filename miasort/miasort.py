from .start import start

def abc_sort(path1, path2, graphs, out_dir='/', plot=True, histogram=False, anchor_option='no',
             colors='red;green;#525252', num_frag_min=2, num_frag_max=1000, extension='6000',
             frag_height=0.6, line_width=2.0, plot_width=50, subplots_margins=(0.9, 0.05, 0.9)):
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
         extension, histogram_options, frag_height,
         line_width, plot_width, subplots_margins)


def multiple_sort(path1, path2, out_dir='/', plot=True, histogram=False, anchor_option='no',
                    colors='red;green;#525252', num_frag_min=2, num_frag_max=1000, extension='6000',
                    frag_height=0.6, line_width=1.5, plot_width=50, subplots_margins=(0.9, 0.05, 0.9)):
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
         extension, histogram_options, frag_height,
         line_width, plot_width, subplots_margins)


def unlimited_multiple_sort(path1, regions, operations, out_dir='/', plot=True, histogram=False, anchor_option='no',
                            colors='red;green;#525252', num_frag_min=2, num_frag_max=1000, extension='6000',
                            frag_height=0.6, line_width=1.5, plot_width=50, subplots_margins=(0.9, 0.05, 0.9)):
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
         extension, histogram_options, frag_height,
         line_width, plot_width, subplots_margins)

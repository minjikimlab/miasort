#!/usr/bin/env python3

import os
import miasort
from miasort import start  # Explicitly import the start function
#from line_profiler import LineProfiler
from memory_profiler import profile
import matplotlib.pyplot as plt
#from memory_profiler import memory_usage


curr_out_dir = "./cr527_profiling_tests/current/"
test_out_dir = "./cr527_profiling_tests/test10_cython/"

# Start iostat to monitor disk I/O
#iostat_process = subprocess.Popen(["iostat", "-dx", "1"], stdout=open("iostat_output.txt", "w"))
#@profile
def main():
    # running from /mia-sort
    # cd ../../mia-sort_output/
    # can change to any repeat
    miasort.abc_sort("../mia-sort_output/GM12878_CTCF-ChIA-PET_LHG0052H.bsorted.ENCLB716IME.hg38.complexes",
                 "../mia-sort_output/GM12878_cr527_repeated_10.bedte",
                 "AtoC;CtoA;AandC",
                 #out_dir=curr_out_dir,
                 out_dir=test_out_dir,
                 plot=False,
                 colors="red;#FF0000;#525252",
                 anchor_option="yes_top")
    """miasort.multiple_sort("./data/100lines_GM12878_SPRITE_4DNFIBEVVTN5.hg38.complexes",
                 "./data/100lines_GM12878-conv-loops-loading-regions_uniqanchors.bedte",
                 out_dir="./GM12878_SPRITE_no_plot_profiling_conv-loops_stripes",
                 plot=False,
                 colors="red;#FF0000;#525252",
                 anchor_option="yes_top")"""

if __name__ == "__main__":
    main()

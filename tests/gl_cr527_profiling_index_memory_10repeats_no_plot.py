#!/usr/bin/env python3

import miasort

#from line_profiler import LineProfiler
from memory_profiler import profile


test_out_dir = "./cr527_profiling_tests/test10_memory/"

@profile
def main():
    # running from /mia-sort
    # cd ../../mia-sort_output/
    # can change to any repeat
    miasort.abc_sort("../mia-sort_output/GM12878_CTCF-ChIA-PET_LHG0052H.bsorted.ENCLB716IME.hg38.complexes",
                 "../mia-sort_output/GM12878_cr527_repeated_10.bedte",
                 "AtoC;CtoA;AandC",
                 out_dir=test_out_dir,
                 plot=False,
                 colors="red;#FF0000;#525252",
                 anchor_option="yes_top")

if __name__ == "__main__":
    main()

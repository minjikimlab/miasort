#!/usr/bin/env python3

#!/usr/bin/env python3

import os
import miasort
from miasort import start  # Explicitly import the start function
from line_profiler import LineProfiler
#from memory_profiler import profile
import matplotlib.pyplot as plt
from memory_profiler import memory_usage

#@profile
def main():
    print("test abc_sort")
    miasort.abc_sort("./data/test_input.region",
                     "./data/test_input.domains",
                     "AtoC;CtoA;AandC;Bcentered;BtoA;BtoC",
                     out_dir="./test_folder_syn_6000",
                     anchor_option="yes_complete")

if __name__ == "__main__":
    # memory profiler
    main()
    """mem_usage = memory_usage((main,))
    plt.plot(mem_usage)
    plt.xlabel('Time (s)')
    plt.ylabel('Memory Usage (MiB)')
    plt.title('Memory Usage Over Time')
    plt.savefig('memory_usage_plot.png')  # Save the plot to a file
    plt.show()  # Display the plot"""
    
    # line profiler stuff
    """profiler = LineProfiler()
    profiler.add_function(miasort.abc_sort)
    profiler.add_function(start)  # Add the start function to the profiler
    profiler_wrapper = profiler(main)
    profiler_wrapper()
    profiler.print_stats()"""
    
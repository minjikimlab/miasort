import os
import miasort
from miasort import start  # Explicitly import the start function
from line_profiler import LineProfiler
#from memory_profiler import profile
import matplotlib.pyplot as plt
#from memory_profiler import memory_usage

import subprocess


"""
input_file = "/nfs/turbo/umms-minjilab/mia-sort_output/GM12878_SPRITE_4DNFIBEVVTN5.hg38.complexes"
output_file = "/nfs/turbo/umms-minjilab/mia-sort_output/100lines_GM12878_SPRITE_4DNFIBEVVTN5.hg38.complexes"
# Use the head command to get the first 10 rows
subprocess.run(["head", "-n", "100", input_file], stdout=open(output_file, "w"))


input_file2 = "/nfs/turbo/umms-minjilab/mia-sort_output//GM12878-conv-loops-loading-regions_uniqanchors.bedte"
output_file2 = "/nfs/turbo/umms-minjilab/mia-sort_output/100lines_GM12878-conv-loops-loading-regions_uniqanchors.bedte"
# Use the head command to get the first 10 rows
subprocess.run(["head", "-n", "100", input_file2], stdout=open(output_file2, "w"))
"""
# Define the output directory as a subfolder in the current directory
#current_directory = os.getcwd()
#out_dir = os.path.join(current_directory, "GM12878_SPRITE_profiling_conv-loops_stripes")

# Ensure the output directory exists
#os.makedirs(out_dir, exist_ok=True)

# Start iostat to monitor disk I/O
#iostat_process = subprocess.Popen(["iostat", "-dx", "1"], stdout=open("iostat_output.txt", "w"))
@profile
def main():
    miasort.abc_sort("./data/100lines_GM12878_SPRITE_4DNFIBEVVTN5.hg38.complexes",
                 "./data/100lines_GM12878-conv-loops-loading-regions_uniqanchors.bedte",
                 "AtoC;CtoA;AandC",
                 out_dir="./GM12878_SPRITE_profiling_conv-loops_stripes",
                 colors="red;#FF0000;#525252",
                 anchor_option="yes_top")

if __name__ == "__main__":
    # memory profiler
    #main()
    """mem_usage = memory_usage((main,))
    plt.plot(mem_usage)
    plt.xlabel('Time (s)')
    plt.ylabel('Memory Usage (MiB)')
    plt.title('Memory Usage Over Time')
    plt.savefig('memory_usage_plot.png')  # Save the plot to a file
    plt.show()  # Display the plot"""

    #disk I/O profiler
    main()
    # Stop iostat
    #iostat_process.terminate()

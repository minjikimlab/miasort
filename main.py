import pybedtools
from pybedtools import BedTool
import argparse
import matplotlib.pyplot as plt

def process_left(ChIA_Drop, Region):
    results = []
    for region in Region:
        left_anchor = f"{region.chrom}\t{region.start}\t{region.end}"
        right_anchor = f"{region.fields[3]}\t{region.fields[4]}\t{region.fields[5]}"

        left_anchor_bed = BedTool(left_anchor, from_string=True)
        right_anchor_bed = BedTool(right_anchor, from_string=True)

        # Filter GEMs that meet the left anchor condition
        filtered_gems = ChIA_Drop.intersect(left_anchor_bed, wa=True)

        valid_gems = []
        for gem in filtered_gems:
            gem_start = int(gem[1])
            gem_end = int(gem[2])
            right_anchor_end = int(right_anchor.split('\t')[2])

            if gem_end <= right_anchor_end:
                valid_gems.append(gem)

        # Convert valid GEMs to BedTool object
        if valid_gems:
            valid_gems_bed = BedTool(valid_gems)

            # Sort the valid GEMs by length using pybedtools sort
            sorted_gems = valid_gems_bed.sort(sizeA=True)

            results.extend(sorted_gems)

    return results

def plot_ranked_gems(ranked_gems, output_file):
    fig, ax = plt.subplots(figsize=(10, 5))

    # Assuming each GEM is a bed-like object with chrom, start, end fields
    for gem in ranked_gems:
        chrom, start, end = gem.chrom, int(gem.start), int(gem.end)
        ax.plot([start, end], [1, 1], marker='|', color='blue')

    # Customize the plot with labels, titles, etc.
    ax.set_title("Ranked GEMs Plot")
    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("GEMs")

    # Save plot to file
    plt.savefig(output_file)
    plt.show()

def main(path1, path2, type, output_file):
    pybedtools.helpers.cleanup()

    # Specify for the path of ___ALL.region.PEanno and import it (GEMs)
    ChIA_Drop = BedTool(path1)

    # Specify for the path of ____PETcnt_G9.motifannot and import it (anchors, regions)
    Region = BedTool(path2)

    if type == "left":
        ranked_gems = process_left(ChIA_Drop, Region)
    else:
        print(f"Type '{type}' is not supported yet.")
        return

    plot_ranked_gems(ranked_gems, output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--path1', type=str, required=True, help='Path to the GEMs file (.PEanno)')
    parser.add_argument('--path2', type=str, required=True, help='Path to the regions file (.domains)')
    parser.add_argument('--type', type=str, required=True, help='Type of processing: left, right, or middle')
    parser.add_argument('--output_file', type=str, required=True, help='Output file for the plot')
    args = parser.parse_args()

    path1 = args.path1
    path2 = args.path2
    type = args.type
    output_file = args.output_file

    main(path1, path2, type, output_file)

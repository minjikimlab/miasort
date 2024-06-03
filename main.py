import pybedtools
from pybedtools import BedTool
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import sys

def process_left(ChIA_Drop, left_anchor, right_anchor):
    left_anchor_bed = BedTool(left_anchor, from_string=True)
    right_anchor_end = int(right_anchor.split('\t')[2])

    # Group fragments by GEM ID
    gem_fragments = {}
    for fragment in ChIA_Drop:
        gem_id = fragment[4]  # Correct column for the unique ID (5th column, 0-indexed 4)
        if gem_id not in gem_fragments:
            gem_fragments[gem_id] = []
        gem_fragments[gem_id].append(fragment)

    valid_gems = []
    for gem_id, fragments in gem_fragments.items():
        # Check if any fragment meets the left anchor condition
        meets_left_anchor = False
        rightmost_end = 0

        fragment_start = sys.maxsize
        for fragment in fragments:
            if int(fragment[1]) < fragment_start:
                fragment_start = int(fragment[1])
                fragment_end = int(fragment[2])
            rightmost_end = max(rightmost_end, int(fragment[2]))

        if fragment_start >= int(left_anchor.split('\t')[1]) and fragment_end <= int(left_anchor.split('\t')[2]):
            meets_left_anchor = True

        if meets_left_anchor and rightmost_end <= right_anchor_end:
            valid_gems.append((gem_id, fragments))

    # Sort GEMs by their length
    valid_gems.sort(key=lambda x: max(int(frag[2]) for frag in x[1]) - min(int(frag[1]) for frag in x[1]))

    return valid_gems

def plot_ranked_gems(ranked_gems, output_file, left_anchor, right_anchor):
    fig, ax = plt.subplots(figsize=(15, 20))

    # Plotting the GEMs
    current_y = 0
    gem_positions = {}
    gem_fragments = {}

    for gem_id, fragments in ranked_gems:
        current_y += 1
        gem_positions[gem_id] = current_y
        gem_fragments[gem_id] = fragments

        for fragment in fragments:
            chrom, start, end = fragment.chrom, int(fragment.start), int(fragment.end)
            y_pos = current_y
            ax.plot([start, end], [y_pos, y_pos], marker='|', color='blue')

    # Connect fragments of the same GEM with lines
    for gem_id, fragments in gem_fragments.items():
        fragments.sort(key=lambda x: int(x.start))  # Sort fragments by their start position
        for i in range(len(fragments) - 1):
            start1, end1 = int(fragments[i].start), int(fragments[i].end)
            start2, end2 = int(fragments[i + 1].start), int(fragments[i + 1].end)
            y_pos = gem_positions[gem_id]
            line = Line2D([end1, start2], [y_pos, y_pos], color='gray', linestyle='--')
            ax.add_line(line)

    # Adding the left and right anchor regions
    left_start, left_end = int(left_anchor.split('\t')[1]), int(left_anchor.split('\t')[2])
    right_start, right_end = int(right_anchor.split('\t')[1]), int(right_anchor.split('\t')[2])

    rect_left = patches.Rectangle((left_start, 0), left_end - left_start, current_y + 1, linewidth=1, edgecolor='r', facecolor='r', alpha=0.2)
    rect_right = patches.Rectangle((right_start, 0), right_end - right_start, current_y + 1, linewidth=1, edgecolor='r', facecolor='r', alpha=0.2)

    ax.add_patch(rect_left)
    ax.add_patch(rect_right)

    ax.set_title(f"Ranked GEMs Plot - {left_anchor.split('\t')[0]}")
    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("GEMs")
    ax.set_yticks(range(1, current_y + 1))
    ax.set_yticklabels(range(1, current_y + 1))
    ax.set_xlim(left_start - 1000, right_end + 1000)

    plt.savefig(output_file)
    plt.show()

def main(path1, path2, type, output_file):
    pybedtools.helpers.cleanup()

    # Specify for the path of ___ALL.region.PEanno and import it (GEMs)
    ChIA_Drop = BedTool(path1)

    # Specify for the path of ____PETcnt_G9.motifannot and import it (anchors, regions)
    Region = BedTool(path2)

    # Use only the first line of the region file for the left and right references
    first_region = list(Region)[0]
    left_anchor = f"{first_region.chrom}\t{first_region.start}\t{first_region.end}"
    right_anchor = f"{first_region.fields[3]}\t{first_region.fields[4]}\t{first_region.fields[5]}"

    if type == "left":
        ranked_gems = process_left(ChIA_Drop, left_anchor, right_anchor)
        plot_ranked_gems(ranked_gems, output_file, left_anchor, right_anchor)
    else:
        print(f"Type '{type}' is not supported yet.")

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

import pybedtools
from pybedtools import BedTool
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import pandas as pd
import time

def process_left(ChIA_Drop_old, left_anchor, right_anchor, region):
    left_anchor_chrom = left_anchor.split('\t')[0]
    left_anchor_start = int(left_anchor.split('\t')[1])
    left_anchor_end = int(left_anchor.split('\t')[2])
    right_anchor_end = int(right_anchor.split('\t')[2])

    print("Filter GEMs to only include those in the region")
    ChIA_Drop = ChIA_Drop_old.intersect(BedTool(region, from_string=True), wa=True, wb=True)

    print("Create BedTool object for the left anchor")
    left_anchor_bed = BedTool(left_anchor, from_string=True)

    print("Get the GEM IDs that intersect with the left anchor")
    intersecting_fragments = ChIA_Drop.intersect(left_anchor_bed, wa=True, wb=True)
    intersecting_gem_ids = set(fragment.fields[4] for fragment in intersecting_fragments)

    print("Filter ChIA_Drop to only include GEMs with intersecting IDs")
    ChIA_Drop = ChIA_Drop.filter(lambda x: x.fields[4] in intersecting_gem_ids)

    ChIA_Drop = ChIA_Drop.filter(lambda x: x.chrom == left_anchor_chrom)

    print("Group GEM fragments by their GEM ID and get the min start and max end positions")
    grouped_gems = {}

    for fragment_interval in ChIA_Drop:
        fragment = fragment_interval.fields
        gem_id = fragment[4]
        start = int(fragment[1])
        end = int(fragment[2])

        if gem_id not in grouped_gems:
            grouped_gems[gem_id] = {
                'min_start': start,
                'max_end': end,
                'fragments': [fragment]
            }
        else:
            grouped_gems[gem_id]['min_start'] = min(grouped_gems[gem_id]['min_start'], start)
            grouped_gems[gem_id]['max_end'] = max(grouped_gems[gem_id]['max_end'], end)
            grouped_gems[gem_id]['fragments'].append(fragment)

    print("Add valid gems")
    valid_gems = []
    for gem_id, gem_info in grouped_gems.items():
        leftmost_fragment_start = gem_info['min_start']
        rightmost_fragment_end = gem_info['max_end']
        gem_length = rightmost_fragment_end - leftmost_fragment_start

        if (
            leftmost_fragment_start >= left_anchor_start
            and leftmost_fragment_start <= left_anchor_end
            and rightmost_fragment_end <= right_anchor_end
        ):
            # Get all fragments of the valid GEM
            fragments = [
                pybedtools.create_interval_from_list(fragment)
                for fragment in gem_info['fragments']
            ]
            valid_gems.append((gem_id, fragments, gem_length))

    print("Sort the valid GEMs by their length")
    valid_gems.sort(key=lambda x: x[2])

    return valid_gems

def plot_ranked_gems(ranked_gems, output_file, left_anchor, right_anchor):
    fig, ax = plt.subplots(figsize=(15, 14))

    # Plotting the GEMs
    gem_positions = {}
    gem_fragments = {}

    for i, (gem_id, fragments, length) in enumerate(ranked_gems):
        y_pos = i + 1
        gem_positions[gem_id] = y_pos
        gem_fragments[gem_id] = fragments

        for fragment in fragments:
            chrom, start, end = fragment.chrom, fragment.start, fragment.end
            ax.plot([start, end], [y_pos, y_pos], marker='|', color='blue')

    # Connect fragments of the same GEM with lines
    for gem_id, fragments in gem_fragments.items():
        for i in range(len(fragments) - 1):
            start1, end1 = fragments[i].start, fragments[i].end
            start2, end2 = fragments[i + 1].start, fragments[i + 1].end
            y_pos = gem_positions[gem_id]
            line = Line2D([end1, start2], [y_pos, y_pos], color='gray', linestyle='--')
            ax.add_line(line)

    # Adding the left and right anchor regions
    left_start, left_end = int(left_anchor.split('\t')[1]), int(left_anchor.split('\t')[2])
    right_start, right_end = int(right_anchor.split('\t')[1]), int(right_anchor.split('\t')[2])

    rect_left = patches.Rectangle((left_start, 0), left_end - left_start, len(ranked_gems) + 1, linewidth=1, edgecolor='r', facecolor='r', alpha=0.2)
    rect_right = patches.Rectangle((right_start, 0), right_end - right_start, len(ranked_gems) + 1, linewidth=1, edgecolor='r', facecolor='r', alpha=0.2)

    ax.add_patch(rect_left)
    ax.add_patch(rect_right)

    ax.set_title(f"Ranked GEMs Plot - {left_anchor.split('\t')[0]}")
    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("GEMs")
    ax.set_yticks(range(1, len(ranked_gems) + 1))
    ax.set_yticklabels(range(1, len(ranked_gems) + 1))
    ax.set_xlim(left_start - 1000, right_end + 1000)
    ax.invert_yaxis()

    plt.savefig(output_file)

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
    region = f"{first_region.chrom}\t{first_region.start}\t{first_region.fields[5]}"

    if type == "left":
        sort_start_time = time.time()
        ranked_gems = process_left(ChIA_Drop, left_anchor, right_anchor, region)
        print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")

        plot_start_time = time.time()
        plot_ranked_gems(ranked_gems, output_file, left_anchor, right_anchor)
        print(f"It took {time.time() - plot_start_time} secs in total to plot the GEMs")

    elif type == "right":
        sort_start_time = time.time()
        ranked_gems = process_right(ChIA_Drop, left_anchor, right_anchor)
        print(f"It took {time.time() - sort_start_time} secs in total to sort the GEMs")

        plot_start_time = time.time()
        plot_ranked_gems(ranked_gems, output_file, left_anchor, right_anchor)
        print(f"It took {time.time() - plot_start_time} secs in total to plot the GEMs")

    else:
        print(f"Type '{type}' is not supported.")

if __name__ == '__main__':
    start_time = time.time()  # record execution time

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

    print(f"It took {time.time() - start_time} secs in total to finish this program")
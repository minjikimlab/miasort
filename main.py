import pybedtools
from pybedtools import BedTool
import argparse
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import time

def process_left(ChIA_Drop, left_anchor, right_anchor):
    left_anchor_chrom = left_anchor.split('\t')[0]
    left_anchor_start = int(left_anchor.split('\t')[1])
    left_anchor_end = int(left_anchor.split('\t')[2])
    right_anchor_end = int(right_anchor.split('\t')[2])

    # Create BedTool object for the left anchor
    left_anchor_bed = BedTool(left_anchor, from_string=True)

    # Intersect ChIA_Drop with left_anchor to get candidate GEMs
    candidate_gems = ChIA_Drop.intersect(left_anchor_bed, wa=True)

    # Filter candidate GEMs to only include those on the same chromosome as the left anchor
    candidate_gems = candidate_gems.filter(lambda x: x.chrom == left_anchor_chrom)

    # Convert candidate_gems to a list of lists
    candidate_gems_list = [gem.fields for gem in candidate_gems]

    # Group GEM fragments by their GEM ID and get the minimum start and maximum end positions
    grouped_gems = {}
    for gem in candidate_gems_list:
        gem_id = gem[4]
        start = int(gem[1])
        end = int(gem[2])

        if gem_id not in grouped_gems:
            grouped_gems[gem_id] = [start, end]
        else:
            grouped_gems[gem_id][0] = min(grouped_gems[gem_id][0], start)
            grouped_gems[gem_id][1] = max(grouped_gems[gem_id][1], end)

    valid_gems = []
    for gem_id, (leftmost_fragment_start, rightmost_fragment_end) in grouped_gems.items():
        gem_length = rightmost_fragment_end - leftmost_fragment_start

        if (
            leftmost_fragment_start >= left_anchor_start
            and leftmost_fragment_start <= left_anchor_end
            and rightmost_fragment_end <= right_anchor_end
        ):
            # Get the fragments of the valid GEM
            fragments = [
                pybedtools.create_interval_from_list(gem)
                for gem in candidate_gems_list
                if gem[4] == gem_id
            ]
            valid_gems.append((gem_id, fragments, gem_length))

    # Sort the valid GEMs by their length
    valid_gems.sort(key=lambda x: x[2])
    print(valid_gems)

    return valid_gems

def plot_ranked_gems(ranked_gems, output_file, left_anchor, right_anchor):
    fig, ax = plt.subplots(figsize=(15, 14))

    # Plotting the GEMs
    current_y = 0
    gem_positions = {}
    gem_fragments = {}

    for gem_id, fragments, length in ranked_gems:
        current_y += 1
        gem_positions[gem_id] = current_y
        gem_fragments[gem_id] = fragments

        for fragment in fragments:
            chrom, start, end = fragment.chrom, fragment.start, fragment.end
            y_pos = current_y
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
        sort_start_time = time.time()
        ranked_gems = process_left(ChIA_Drop, left_anchor, right_anchor)
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

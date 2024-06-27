import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import time

def plot_ranked_gems(start_time, ranked_gems, output_file, left_anchor,
                     right_anchor, middle_anchor=None, flag="normal"):
    fig, ax = plt.subplots(figsize=(15, 18))

    # Define the vertical spacing
    vertical_spacing = 0.5  # Smaller value for less spacing

    # Plotting the GEMs
    gem_positions = {}
    gem_fragments = {}

    for i, (gem_id, fragments, length) in enumerate(ranked_gems):
        y_pos = i * vertical_spacing + 1
        gem_positions[gem_id] = y_pos
        gem_fragments[gem_id] = fragments

        for fragment in fragments:
            chrom, start, end = fragment.chrom, fragment.start, fragment.end
            rect = patches.Rectangle((start, i - 0.05),
                                     end - start, 0.1, linewidth=1, edgecolor='g', facecolor='g')
            ax.add_patch(rect)

    count = 0
    # Connect fragments of the same GEM with solid lines
    for gem_id, fragments in gem_fragments.items():
        for i in range(len(fragments) - 1):
            end1 = fragments[i].end
            start2 = fragments[i + 1].start
            y_pos = gem_positions[gem_id]
            line = Line2D([end1, start2], [count, count], color='black', linestyle='-')
            ax.add_line(line)
        count += 1

    left_start, left_end = int(left_anchor.split('\t')[1]), int(left_anchor.split('\t')[2])
    right_start, right_end = int(right_anchor.split('\t')[1]), int(right_anchor.split('\t')[2])
    rect_left = patches.Rectangle((left_start, -1), left_end - left_start,
                                  (len(ranked_gems) + 1), linewidth=1,
                                  edgecolor='r', facecolor='r', alpha=0.2)
    rect_right = patches.Rectangle((right_start, -1), right_end - right_start,
                                   (len(ranked_gems) + 1), linewidth=1,
                                   edgecolor='r', facecolor='r', alpha=0.2)

    if flag != "only-middle" and flag != "only-middle-1frag":
        # Adding the left and right anchor regions
        ax.add_patch(rect_left)
        ax.add_patch(rect_right)

    if (flag == "middle" or flag == "only-middle" \
    or flag == "only-middle-1frag") and middle_anchor:
        _, positions = middle_anchor.split(':')
        middle_start, middle_end = positions.split('-')
        rect_middle = patches.Rectangle(
            (int(middle_start), 0),
            int(middle_end) - int(middle_start),
            (len(ranked_gems) + 1) * vertical_spacing,
            linewidth=1,
            edgecolor='red',
            facecolor='red',
            alpha=0.2
        )
        ax.add_patch(rect_middle)

    ax.set_title(f"Ranked GEMs Plot - {left_anchor.split('\t')[0]}")
    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("GEMs")
    ax.set_yticks([i for i in range(len(ranked_gems))], labels=range(1, len(ranked_gems) + 1))
    ax.set_xlim(left_start - 1000, right_end + 1000)
    ax.set_ylim(-1, len(ranked_gems) + 2)
    ax.invert_yaxis()  # labels read top-to-bottom

    print(f"It took {time.time() - start_time} secs in total to finish this program")

    plt.savefig(output_file)

    # for displaying the plot in a complete way,
    # delete this in the case of running it on GreatLakes servers
    # plt.show()


def plot_ranked_gems_multiple_regions(start_time, ranked_gems, output_file, regions):
    fig, ax = plt.subplots(figsize=(15, 18))

    # Define the vertical spacing
    vertical_spacing = 0.5  # Smaller value for less spacing

    # Plotting the GEMs
    gem_positions = {}
    gem_fragments = {}

    for i, (gem_id, fragments, length) in enumerate(ranked_gems):
        y_pos = i * vertical_spacing + 1
        gem_positions[gem_id] = y_pos
        gem_fragments[gem_id] = fragments

        for fragment in fragments:
            chrom, start, end = fragment.chrom, fragment.start, fragment.end
            rect = patches.Rectangle((start, i - 0.05),
                                     end - start, 0.1, linewidth=1, edgecolor='g', facecolor='g')
            ax.add_patch(rect)

    count = 0
    # Connect fragments of the same GEM with solid lines
    for gem_id, fragments in gem_fragments.items():
        for i in range(len(fragments) - 1):
            end1 = fragments[i].end
            start2 = fragments[i + 1].start
            y_pos = gem_positions[gem_id]
            line = Line2D([end1, start2], [count, count], color='black', linestyle='-')
            ax.add_line(line)
        count += 1

    left_start = regions[0][1]
    right_end = regions[0][2]
    for region in regions:
        _, start, end = region
        left_start = min(left_start, start)
        right_end = max(right_end, end)
        rect = patches.Rectangle((start, -1), end - start,
                                    (len(ranked_gems) + 1), linewidth=1,
                                    edgecolor='r', facecolor='r', alpha=0.2)
        ax.add_patch(rect)

    ax.set_title(f"Ranked GEMs Plot")
    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("GEMs")
    ax.set_yticks([i for i in range(len(ranked_gems))], labels=range(1, len(ranked_gems) + 1))
    ax.set_xlim(left_start - 1000, right_end + 1000)
    ax.set_ylim(-1, len(ranked_gems) + 2)
    ax.invert_yaxis()  # labels read top-to-bottom

    print(f"It took {time.time() - start_time} secs in total to finish this program")

    plt.savefig(output_file)

    # for displaying the plot in a complete way,
    # delete this in the case of running it on GreatLakes servers
    # plt.show()
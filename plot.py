import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import time
import os

# def plot_ranked_gems(ranked_gems, output_file, left_anchor,
#                      right_anchor, middle_anchor):
#     """Plot accurate plot."""
#     fig, ax = plt.subplots(figsize=(30, 30))

#     # Plotting the GEMs
#     gem_positions = {}
#     gem_fragments = {}

#     for i, (gem_id, fragments, length) in enumerate(ranked_gems):
#         gem_fragments[gem_id] = fragments

#         for fragment in fragments:
#             chrom, start, end = fragment.chrom, fragment.start, fragment.end
#             rect = patches.Rectangle((start, i - 0.25),
#                                      end - start, 0.5, linewidth=1, edgecolor='g', facecolor='g')
#             ax.add_patch(rect)

#     count = 0
#     # Connect fragments of the same GEM with solid lines
#     for gem_id, fragments in gem_fragments.items():
#         for i in range(len(fragments) - 1):
#             end1 = fragments[i].end
#             start2 = fragments[i + 1].start
#             line = Line2D([end1, start2], [count, count], color='grey', linestyle='-', linewidth=0.1)
#             ax.add_line(line)
#         count += 1

#     left_start, left_end = int(left_anchor.split('\t')[1]), int(left_anchor.split('\t')[2])
#     right_start, right_end = int(right_anchor.split('\t')[1]), int(right_anchor.split('\t')[2])
#     middle_start, middle_end = int(middle_anchor.split('\t')[1]), int(middle_anchor.split('\t')[2])

#     rect_left = patches.Rectangle((left_start, -1), left_end - left_start,
#                                   (len(ranked_gems) + 1), linewidth=1,
#                                   edgecolor='r', facecolor='r', alpha=0.2)
#     rect_right = patches.Rectangle((right_start, -1), right_end - right_start,
#                                    (len(ranked_gems) + 1), linewidth=1,
#                                    edgecolor='r', facecolor='r', alpha=0.2)
#     rect_middle = patches.Rectangle((middle_start, -1), middle_end - middle_start,
#                                 (len(ranked_gems) + 1), linewidth=1,
#                                 edgecolor='r', facecolor='r', alpha=0.2)

#     # Adding the left and right anchor regions
#     ax.add_patch(rect_left)
#     ax.add_patch(rect_right)
#     ax.add_patch(rect_middle)

#     ax.set_title(f"Ranked GEMs Plot - {left_anchor.split('\t')[0]}")
#     ax.set_xlabel("Genomic Position")
#     ax.set_ylabel("GEMs")
#     ax.set_yticks([i for i in range(len(ranked_gems))], labels=range(1, len(ranked_gems) + 1))
#     ax.set_xlim(left_start - 1000, right_end + 1000)
#     ax.set_ylim(-1, len(ranked_gems) + 2)
#     ax.invert_yaxis()  # labels read top-to-bottom

#     plt.savefig(output_file)

#     # for displaying the plot in a complete way,
#     # delete this in the case of running it on GreatLakes servers
#     # plt.show()


def plot_ranked_gems_scaled(ranked_gems, output_file, left_anchor,
                        right_anchor, middle_anchor, out_dir, colors_flags, anchor_options):
    directory_str = output_file
    if out_dir != "/":
        # create directory if it doesn't exist
        directory_str = f"./{out_dir}/{output_file}"
        directory = os.path.dirname(directory_str)
        if not os.path.exists(directory):
            os.makedirs(directory)

    fig, ax = plt.subplots(figsize=(50, 50))

    gem_fragments = {}

    for i, (gem_id, fragments, _) in enumerate(ranked_gems):
        gem_fragments[gem_id] = fragments
        for fragment in fragments:
            _, start, end = fragment.chrom, fragment.start, fragment.end
            rect = patches.Rectangle((start, i - 0.2),
                                     (end - start) * 3, 0.4, linewidth=1,
                                     edgecolor=colors_flags["fragments"],
                                     facecolor=colors_flags["fragments"])
            ax.add_patch(rect)

    count = 0
    # connect fragments of the same GEM with solid lines
    for gem_id, fragments in gem_fragments.items():
        start = fragments[0].start
        end = fragments[-1].end
        line = Line2D([start, end], [count, count],
                      color=colors_flags["lines"],
                      linestyle='-', linewidth=0.4)
        ax.add_line(line)
        count += 1

    left_start, left_end = int(left_anchor.split('\t')[1]), int(left_anchor.split('\t')[2])
    right_start, right_end = int(right_anchor.split('\t')[1]), int(right_anchor.split('\t')[2])
    middle_start, middle_end = int(middle_anchor.split('\t')[1]), int(middle_anchor.split('\t')[2])

    if anchor_options == "yes_complete":
        rect_left = patches.Rectangle((left_start, -1), left_end - left_start,
                                    (len(ranked_gems) + 1), linewidth=1,
                                    edgecolor=colors_flags["anchors"],
                                    facecolor=colors_flags["anchors"],
                                    alpha=0.2)
        rect_right = patches.Rectangle((right_start, -1), right_end - right_start,
                                    (len(ranked_gems) + 1), linewidth=1,
                                    edgecolor=colors_flags["anchors"],
                                    facecolor=colors_flags["anchors"],
                                    alpha=0.2)
        rect_middle = patches.Rectangle((middle_start, -1), middle_end - middle_start,
                                        (len(ranked_gems) + 1), linewidth=1,
                                        edgecolor=colors_flags["anchors"],
                                        facecolor=colors_flags["anchors"],
                                        alpha=0.2)
        # adding the left and right anchor regions
        ax.add_patch(rect_left)
        ax.add_patch(rect_right)
        ax.add_patch(rect_middle)

    elif anchor_options == "yes_top":
        print("Waiting to be implemented")

    ax.set_title(f"Ranked GEMs Plot - {left_anchor.split('\t')[0]}")
    ax.set_xlabel("Genomic Position")
    ax.set_ylabel("GEMs")
    ax.set_yticks([i for i in range(len(ranked_gems))], labels=range(1, len(ranked_gems) + 1))
    ax.set_xlim(min(left_start, right_start, middle_start) - 1000,
                max(left_end, right_end, middle_end) + 1000)
    ax.set_ylim(-1, len(ranked_gems) + 2)
    ax.invert_yaxis()  # labels read top-to-bottom

    plt.savefig(directory_str)

    plt.close(fig)  # close the figure to free up memory

    # for displaying the plot in a complete way,
    # delete this in the case of running it on GreatLakes servers
    # plt.show()


def plot_ranked_gems_multiple_regions(start_time, ranked_gems, output_file, regions):
    fig, ax = plt.subplots(figsize=(30, 30))

    # Plotting the GEMs
    gem_positions = {}
    gem_fragments = {}

    for i, (gem_id, fragments, length) in enumerate(ranked_gems):
        gem_fragments[gem_id] = fragments

        for fragment in fragments:
            chrom, start, end = fragment.chrom, fragment.start, fragment.end
            rect = patches.Rectangle((start, i - 0.25),
                                     end - start, 0.5, linewidth=1, edgecolor='g', facecolor='g')
            ax.add_patch(rect)

    count = 0
    # Connect fragments of the same GEM with solid lines
    for gem_id, fragments in gem_fragments.items():
        for i in range(len(fragments) - 1):
            end1 = fragments[i].end
            start2 = fragments[i + 1].start
            line = Line2D([end1, start2], [count, count], color='black', linestyle='-', linewidth=0.5)
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

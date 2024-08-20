import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec
from helper import figsize_height_scaler, kb_format, create_plot_title


def plot_three_ranked_gems(ranked_gems_list, output_file, left_anchor_list,
                           right_anchor_list, middle_anchor_list, out_dir, colors_flags,
                           anchor_options, id, path1, commands_list, extension, flag="abc", regions=[]):
    directory_str = output_file
    if out_dir != "/":
        directory_str = f"./{out_dir}/{output_file}"

    heights = [round(figsize_height_scaler(len(ranked_gems_list[0]))),
                round(figsize_height_scaler(len(ranked_gems_list[1]))),
                round(figsize_height_scaler(len(ranked_gems_list[2])))]

    fig = plt.figure(figsize=(50, sum(heights)*1.5))

    # Create GridSpec with custom heights
    gs = GridSpec(3, 1, height_ratios=heights, figure=fig)

    for idx in range(3):
        ax = fig.add_subplot(gs[idx])

        ranked_gems = ranked_gems_list[idx]
        left_anchor = left_anchor_list[idx]
        right_anchor = right_anchor_list[idx]
        middle_anchor = middle_anchor_list[idx]

        gem_fragments = {}

        for i, (gem_id, fragments, _) in enumerate(ranked_gems):
            gem_fragments[gem_id] = fragments
            for fragment in fragments:
                _, start, end = fragment.chrom, fragment.start, fragment.end
                if extension == "natural":
                    width = end - start
                else:
                    width = int(extension)
                rect = patches.Rectangle((start - width / 2, i - 0.225),
                                         width, 0.6, linewidth=2,
                                         edgecolor="black",
                                         facecolor=colors_flags["fragments"],
                                         zorder=2)
                ax.add_patch(rect)

        count = 0
        for gem_id, fragments in gem_fragments.items():
            start = fragments[0].start
            end = fragments[-1].start
            line = Line2D([start, end], [count, count],
                          color=colors_flags["lines"],
                          linestyle='-', linewidth=0.4,
                          zorder=1)
            ax.add_line(line)
            count += 1

        if flag == "abc":
            left_start, left_end = int(left_anchor.split('\t')[1]), int(left_anchor.split('\t')[2])
            right_start, right_end = int(right_anchor.split('\t')[1]), int(right_anchor.split('\t')[2])
            middle_start, middle_end = int(middle_anchor.split('\t')[1]), int(middle_anchor.split('\t')[2])

        if anchor_options == "yes_complete":
            if flag == "abc":
                rect_left = patches.Rectangle((left_start, -1), left_end - left_start,
                                              (len(ranked_gems) + 2), linewidth=1,
                                              edgecolor=colors_flags["anchors"],
                                              facecolor=colors_flags["anchors"],
                                              alpha=0.2)
                rect_right = patches.Rectangle((right_start, -1), right_end - right_start,
                                               (len(ranked_gems) + 2), linewidth=1,
                                               edgecolor=colors_flags["anchors"],
                                               facecolor=colors_flags["anchors"],
                                               alpha=0.2)
                rect_middle = patches.Rectangle((middle_start, -1), middle_end - middle_start,
                                                (len(ranked_gems) + 2), linewidth=1,
                                                edgecolor=colors_flags["anchors"],
                                                facecolor=colors_flags["anchors"],
                                                alpha=0.2)
                ax.add_patch(rect_left)
                ax.add_patch(rect_right)
                ax.add_patch(rect_middle)

            else:
                left_start = regions[0][1]
                right_end = regions[0][2]
                for region in regions:
                    _, start, end = region
                    left_end = min(left_start, start)
                    right_end = max(right_end, end)
                    rect = patches.Rectangle((start, -1), end - start,
                                             (len(ranked_gems) + 2), linewidth=1,
                                             edgecolor=colors_flags["anchors"],
                                             facecolor=colors_flags["anchors"],
                                             alpha=0.2)
                    ax.add_patch(rect)

            ax.set_ylim(-1, len(ranked_gems) + 1)

        elif anchor_options == "yes_top":
            if flag == "abc":
                rect_left = patches.Rectangle((left_start, -3), left_end - left_start,
                                              1, linewidth=1, edgecolor="black",
                                              facecolor="black")
                rect_right = patches.Rectangle((right_start, -3), right_end - right_start,
                                               1, linewidth=1, edgecolor="black",
                                               facecolor="black")
                rect_middle = patches.Rectangle((middle_start, -3), middle_end - middle_start,
                                                1, linewidth=1, edgecolor="black",
                                                facecolor="black")
                ax.add_patch(rect_left)
                ax.add_patch(rect_right)
                ax.add_patch(rect_middle)
            else:
                left_start = regions[0][1]
                right_end = regions[0][2]
                for region in regions:
                    _, start, end = region
                    left_end = min(left_start, start)
                    right_end = max(right_end, end)
                    rect = patches.Rectangle((start, -3), end - start,
                                             1, linewidth=1,
                                             edgecolor="black",
                                             facecolor="black")
                    ax.add_patch(rect)

            ax.set_ylim(-3, len(ranked_gems) + 1)

        else:  # anchor_options == "no"
            ax.set_ylim(-1, len(ranked_gems) + 1)

        title_font = {'fontsize': 30, 'fontweight': 'bold'}
        label_font = {'fontsize': 30}
        tick_font_size = 23

        if flag == "abc":
            anchors = [left_anchor, middle_anchor, right_anchor]
            anchors.sort(key=lambda x: int(x.split('\t')[1]))
            ax.set_title(create_plot_title(id, path1, commands_list[idx], anchors),
                         fontdict=title_font)

            left_end = min(left_start, right_start, middle_start)
            right_end = max(left_end, right_end, middle_end)
        else:
            ax.set_title("Multiple")

        distance = right_end - left_end
        margin = round(0.05 * distance)
        ax.set_xlim(left_end - margin, right_end + margin)

        ax.set_xlabel("Genomic Position", fontdict=label_font)
        ax.set_ylabel("Chromatin Complexes", fontdict=label_font)
        ax.set_yticks([i for i in range(len(ranked_gems))], labels=range(1, len(ranked_gems) + 1))

        ax.invert_yaxis()  # labels read top-to-bottom
        ax.tick_params(axis='x', labelsize=tick_font_size)
        ax.tick_params(axis='y', labelsize=tick_font_size)

        ax.xaxis.set_major_formatter(plt.FuncFormatter(kb_format))

    plt.subplots_adjust(top=0.9, bottom=0.05, hspace=0.9)
    plt.savefig(directory_str)
    plt.close(fig)

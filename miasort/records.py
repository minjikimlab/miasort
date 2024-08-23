import os
import csv

def write_to_csv_file(id, A, B, C, command, num_complexes, csv_file, out_dir, ranked_gems):
    # Determine the output path
    if out_dir != "/":
        output_path = os.path.join(out_dir, csv_file)
    else:
        output_path = csv_file

    histogram = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
    for gem_tuple in ranked_gems:
        num_fragments = int(len(gem_tuple[1]))
        if num_fragments >= 5:
            histogram[5] += 1
        else:
            histogram[num_fragments] += 1

    l = A.split('\t')
    anchor_a = f"{l[0]}:{l[1]}-{l[2]}"

    m = B.split('\t')
    anchor_b = f"{m[0]}:{m[1]}-{m[2]}"

    r = C.split('\t')
    anchor_c = f"{r[0]}:{r[1]}-{r[2]}"

    region = f"{l[0]}:{l[1]}-{r[2]}"

    with open(output_path, "a") as file:
        writer = csv.writer(file)
        field = [id, anchor_a, anchor_b, anchor_c, region, command, num_complexes,
                 histogram[1], histogram[2], histogram[3], histogram[4], histogram[5]]
        writer.writerow(field)


def write_to_csv_file_multiple(id, A, B, C, region, command, num_complexes, csv_file, out_dir, ranked_gems):
    # Determine the output path
    if out_dir != "/":
        output_path = os.path.join(out_dir, csv_file)
    else:
        output_path = csv_file

    histogram = {1: 0, 2: 0, 3: 0, 4: 0, 5: 0}
    for gem_tuple in ranked_gems:
        num_fragments = int(len(gem_tuple[1]))
        if num_fragments >= 5:
            histogram[5] += 1
        else:
            histogram[num_fragments] += 1

    with open(output_path, "a") as file:
        writer = csv.writer(file)
        field = [id, A, B, C, region, command, num_complexes,
                 histogram[1], histogram[2], histogram[3], histogram[4], histogram[5]]
        writer.writerow(field)

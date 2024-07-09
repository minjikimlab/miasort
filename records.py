import os
import csv

def generate_file(id, A, B, C, command, num_complexes, csv_file, out_dir):
    # determine the output path
    if out_dir != "/":
        output_path = os.path.join(out_dir, csv_file)
    else:
        output_path = csv_file

    l = A.split('\t')
    anchor_a = f"{l[0]}:{l[1]}-{l[2]}"

    m = B.split('\t')
    anchor_b = f"{m[0]}:{m[1]}-{m[2]}"

    r = C.split('\t')
    anchor_c = f"{r[0]}:{r[1]}-{r[2]}"

    region = f"{l[0]}:{l[1]}-{r[2]}"

    with open(output_path, "a") as file:
        writer = csv.writer(file)
        field = [id, anchor_a, anchor_b, anchor_c, region, command, num_complexes]
        writer.writerow(field)

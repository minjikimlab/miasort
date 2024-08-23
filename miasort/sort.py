import pybedtools
from pybedtools import BedTool

def process_left(ChIA_Drop, num_fragments_min, num_fragments_max, left_anchor, right_anchor, region):
    right_anchor_start = int(right_anchor.split('\t')[1])
    right_anchor_end = int(right_anchor.split('\t')[2])

    left_anchor_bed = BedTool(left_anchor, from_string=True)

    intersecting_fragments = ChIA_Drop.intersect(left_anchor_bed, wa=True, wb=True)
    intersecting_gem_ids = set(fragment.fields[4] for fragment in intersecting_fragments)

    ChIA_Drop = ChIA_Drop.filter(lambda x: x.fields[4] in intersecting_gem_ids)

    grouped_gems = {}
    bad_gem_ids = set()
    for fragment_interval in ChIA_Drop:
        fragment = fragment_interval.fields
        gem_id = fragment[4]
        gem_size = int(fragment[3])
        start = int(fragment[1])
        end = int(fragment[2])

        if (start >= right_anchor_start and end <= right_anchor_end) or end >= right_anchor_start:
            bad_gem_ids.add(gem_id)

        if gem_id not in grouped_gems.keys():
            grouped_gems[gem_id] = {
                'min_start': start,
                'max_end': end,
                'fragments': [fragment],
                'gem_size': gem_size
            }
        else:
            grouped_gems[gem_id]['min_start'] = min(grouped_gems[gem_id]['min_start'], start)
            grouped_gems[gem_id]['max_end'] = max(grouped_gems[gem_id]['max_end'], end)
            grouped_gems[gem_id]['fragments'].append(fragment)

    for bad_gem_id in bad_gem_ids:
        del grouped_gems[bad_gem_id]

    valid_gems = []
    for gem_id, gem_info in grouped_gems.items():
        leftmost_fragment_start = gem_info['min_start']
        rightmost_fragment_end = gem_info['max_end']
        gem_size = gem_info['gem_size']
        gem_length = rightmost_fragment_end - leftmost_fragment_start

        num_frgaments = len(gem_info['fragments'])
        if num_frgaments >= num_fragments_min and num_frgaments <= num_fragments_max:
            fragments = [
                pybedtools.create_interval_from_list(fragment)
                for fragment in gem_info['fragments']
            ]
            valid_gems.append((gem_id, fragments, gem_length))

    valid_gems.sort(key=lambda x: x[2])

    # test_gem_id_difference(valid_gems, "cr1491_SE_Left.bed")

    return valid_gems


def process_right(ChIA_Drop, num_fragments_min, num_fragments_max, left_anchor, right_anchor, region):
    left_anchor_start = int(left_anchor.split('\t')[1])
    left_anchor_end = int(left_anchor.split('\t')[2])

    right_anchor_bed = BedTool(right_anchor, from_string=True)
    intersecting_fragments = ChIA_Drop.intersect(right_anchor_bed, wa=True, wb=True)
    intersecting_gem_ids = set(fragment.fields[4] for fragment in intersecting_fragments)

    ChIA_Drop = ChIA_Drop.filter(lambda x: x.fields[4] in intersecting_gem_ids)

    grouped_gems = {}
    bad_gem_ids = set()
    for fragment_interval in ChIA_Drop:
        fragment = fragment_interval.fields
        gem_id = fragment[4]
        gem_size = int(fragment[3])
        start = int(fragment[1])
        end = int(fragment[2])

        if (start >= left_anchor_start and end <= left_anchor_end) or start <= left_anchor_end:
            bad_gem_ids.add(gem_id)

        if gem_id not in grouped_gems.keys():
            grouped_gems[gem_id] = {
                'min_start': start,
                'max_end': end,
                'fragments': [fragment],
                'gem_size': gem_size
            }
        else:
            grouped_gems[gem_id]['min_start'] = min(grouped_gems[gem_id]['min_start'], start)
            grouped_gems[gem_id]['max_end'] = max(grouped_gems[gem_id]['max_end'], end)
            grouped_gems[gem_id]['fragments'].append(fragment)

    for bad_gem_id in bad_gem_ids:
        del grouped_gems[bad_gem_id]

    valid_gems = []
    for gem_id, gem_info in grouped_gems.items():
        leftmost_fragment_start = gem_info['min_start']
        rightmost_fragment_end = gem_info['max_end']
        gem_size = gem_info['gem_size']
        gem_length = rightmost_fragment_end - leftmost_fragment_start

        num_frgaments = len(gem_info['fragments'])
        if num_frgaments >= num_fragments_min and num_frgaments <= num_fragments_max:
            fragments = [
                pybedtools.create_interval_from_list(fragment)
                for fragment in gem_info['fragments']
            ]
            valid_gems.append((gem_id, fragments, gem_length))

    valid_gems.sort(key=lambda x: x[2])

    # test_gem_id_difference(valid_gems, "cr1491_SE_Right.bed")

    return valid_gems


def process_middle(ChIA_Drop, num_fragments_min, num_fragments_max, left_anchor, right_anchor, region, middle_anchor):
    middle_anchor_chrom, middle_anchor_start, middle_anchor_end = middle_anchor.split('\t')

    left_anchor_start = int(left_anchor.split('\t')[1])
    left_anchor_end = int(left_anchor.split('\t')[2])
    right_anchor_start = int(right_anchor.split('\t')[1])
    right_anchor_end = int(right_anchor.split('\t')[2])

    # Get unique GEM IDs from candidate GEMs
    candidate_gem_ids = set(fragment[4] for fragment in ChIA_Drop)

    # Make sure the GEM has at least 1 fragment in left area (area_1)
    # and at least 1 fragment in right area (area_2)
    area_1 = BedTool(f"{middle_anchor_chrom}\t{left_anchor_end}\t{middle_anchor_start}",
                     from_string=True)
    area_2 = BedTool(f"{middle_anchor_chrom}\t{middle_anchor_end}\t{right_anchor_start}",
                     from_string=True)

    area_1_frags = ChIA_Drop.intersect(area_1, wa=True, wb=True)
    area_1_gem_ids = set(area_1_frag[4] for area_1_frag in area_1_frags)

    area_2_frags = ChIA_Drop.intersect(area_2, wa=True, wb=True)
    area_2_gem_ids = set(area_2_frag[4] for area_2_frag in area_2_frags)

    in_area_gem_ids = area_1_gem_ids.intersection(area_2_gem_ids)

    # Filter ChIA_Drop to retain only GEMs with IDs in candidate_gem_ids
    valid_gems = []
    gem_fragments = {}
    gem_lengths = {}
    bad_gem_ids = set()

    for fragment in ChIA_Drop:
        gem_id = fragment[4]
        if gem_id in candidate_gem_ids and gem_id in in_area_gem_ids:
            if gem_id not in gem_fragments:
                gem_fragments[gem_id] = []
            gem_fragments[gem_id].append(fragment)

            start = int(fragment[1])
            end = int(fragment[2])

            if (start >= left_anchor_start and end <= left_anchor_end) \
            or (start >= right_anchor_start and end <= right_anchor_end):
                bad_gem_ids.add(gem_id)

            if gem_id in gem_lengths:
                gem_lengths[gem_id] = (min(gem_lengths[gem_id][0], start),
                                       max(gem_lengths[gem_id][1], end))
            else:
                gem_lengths[gem_id] = (start, end)

    for bad_gem_id in bad_gem_ids:
        del gem_fragments[bad_gem_id]

    # Further filter valid GEMs based on the leftmost fragment and right anchor
    for gem_id, fragments in gem_fragments.items():
        # leftmost_fragment_start = int(fragments[0][1])
        # leftmost_fragment_end = int(fragments[0][2])
        start, end = gem_lengths[gem_id]

        if start > left_anchor_end and end < right_anchor_start \
        and len(fragments) >= num_fragments_min and len(fragments) <= num_fragments_max:
            valid_gems.append((gem_id, fragments, end - start, start, end))

    # Sort GEMs by their length
    if len(valid_gems) > 0:
        valid_gems.sort(key=lambda x: x[2])

        prev_start = valid_gems[0][3]
        prev_end = valid_gems[0][4]
        new_valid_gems = [valid_gems[0][:3]]

        for valid_gem in valid_gems[1:]:
            curr_start = valid_gem[3]
            curr_end = valid_gem[4]

            if not (prev_start <= curr_start or prev_end >= curr_end):
                new_valid_gems.append(valid_gem[:3])
                prev_start = curr_start
                prev_end = curr_end
    else:
        new_valid_gems = valid_gems

    return new_valid_gems


def process_multiple(ChIA_Drop, num_fragments_min, num_fragments_max, yes_chroms, no_chroms):
    # reduce search space
    chr_id = yes_chroms[0][0]
    if not len(yes_chroms):
        left_most_end = no_chroms[0][1]
        right_most_end = no_chroms[-1][2]
    elif not len(no_chroms):
        left_most_end = yes_chroms[0][1]
        right_most_end = yes_chroms[-1][2]
    else:
        left_most_end = min(yes_chroms[0][1], no_chroms[0][1])
        right_most_end = max(yes_chroms[-1][2], no_chroms[-1][2])
    region = BedTool(f"{chr_id}\t{left_most_end}\t{right_most_end}", from_string=True)

    # Process the first chromosome to initialize the valid_gem_ids
    chr_id, left, right = yes_chroms[0]
    frags = ChIA_Drop.intersect(
        BedTool(f"{chr_id}\t{left}\t{right}", from_string=True),
        wa=True,
        wb=True,
    )
    valid_gem_ids = set(frag[4] for frag in frags)

    # Intersect with GEM IDs from the remaining chromosomes
    for yes_chrom in yes_chroms[1:]:
        chr_id, left, right = yes_chrom[:3]
        frags = ChIA_Drop.intersect(
            BedTool(f"{chr_id}\t{left}\t{right}", from_string=True),
            wa=True,
            wb=True,
        )
        candidate_gem_ids = set(frag[4] for frag in frags)
        valid_gem_ids.intersection_update(candidate_gem_ids)

    # Exclude GEM IDs from the no_chroms regions
    for no_chrom in no_chroms:
        chr_id, left, right = no_chrom[:3]
        frags = ChIA_Drop.intersect(
            BedTool(f"{chr_id}\t{left}\t{right}", from_string=True),
            wa=True,
            wb=True,
        )
        candidate_gem_ids = set(frag[4] for frag in frags)
        valid_gem_ids.difference_update(candidate_gem_ids)

    valid_gems = []
    gem_fragments = {}
    gem_lengths = {}

    for fragment in ChIA_Drop:
        gem_id = fragment[4]
        start = int(fragment[1])
        end = int(fragment[2])
        if gem_id in valid_gem_ids and start >= left_most_end and end <= right_most_end:
            if gem_id not in gem_fragments:
                gem_fragments[gem_id] = []
            gem_fragments[gem_id].append(fragment)

            if gem_id in gem_lengths:
                gem_lengths[gem_id] = (min(gem_lengths[gem_id][0], start),
                                       max(gem_lengths[gem_id][1], end))
            else:
                gem_lengths[gem_id] = (start, end)

    # Further filter valid GEMs based on the leftmost fragment and right anchor
    for gem_id, fragments in gem_fragments.items():
        # leftmost_fragment_start = int(fragments[0][1])
        # leftmost_fragment_end = int(fragments[0][2])
        start, end = gem_lengths[gem_id]
        fragments.sort(key=lambda x: x[1])
        if len(fragments) >= num_fragments_min and len(fragments) <= num_fragments_max:
            valid_gems.append((gem_id, fragments, end - start))

    return valid_gems

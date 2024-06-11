import pybedtools
from pybedtools import BedTool


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

    print("Filter ChIA_Drop to only include GEMs with the correct chrom id")
    ChIA_Drop = ChIA_Drop.filter(lambda x: x.chrom == left_anchor_chrom)

    print("Group GEM fragments by their GEM ID and get the min start and max end positions")
    grouped_gems = {}
    for fragment_interval in ChIA_Drop:
        fragment = fragment_interval.fields
        gem_id = fragment[4]
        gem_size = int(fragment[3])
        start = int(fragment[1])
        end = int(fragment[2])

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

    print("Add valid gems")
    valid_gems = []
    for gem_id, gem_info in grouped_gems.items():
        leftmost_fragment_start = gem_info['min_start']
        rightmost_fragment_end = gem_info['max_end']
        gem_size = gem_info['gem_size']
        gem_length = rightmost_fragment_end - leftmost_fragment_start

        if (
            leftmost_fragment_start >= left_anchor_start
            and leftmost_fragment_start <= left_anchor_end
            and rightmost_fragment_end <= right_anchor_end
            and len(gem_info['fragments']) == gem_size
        ):
            # Get all fragments of the valid GEM
            fragments = [
                pybedtools.create_interval_from_list(fragment)
                for fragment in gem_info['fragments']
            ]
            valid_gems.append((gem_id, fragments, gem_length))

    print("Sort the valid GEMs by their length")
    valid_gems.sort(key=lambda x: x[2])
    print(valid_gems)
    return valid_gems


def process_both(ChIA_Drop_old, left_anchor, right_anchor, region):
    left_anchor_chrom = left_anchor.split('\t')[0]
    left_anchor_start = int(left_anchor.split('\t')[1])
    left_anchor_end = int(left_anchor.split('\t')[2])
    right_anchor_start = int(right_anchor.split('\t')[1])
    right_anchor_end = int(right_anchor.split('\t')[2])

    print("Filter GEMs to only include those in the region")
    ChIA_Drop = ChIA_Drop_old.intersect(BedTool(region, from_string=True), wa=True, wb=True)

    print("Filter ChIA_Drop to only include GEMs with the correct chrom id")
    ChIA_Drop = ChIA_Drop.filter(lambda x: x.chrom == left_anchor_chrom)

    print("Group GEM fragments by their GEM ID and get the min start and max end positions")
    grouped_gems = {}
    for fragment_interval in ChIA_Drop:
        fragment = fragment_interval.fields
        gem_id = fragment[4]
        gem_size = int(fragment[3])
        start = int(fragment[1])
        end = int(fragment[2])

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

    print("Add valid gems")
    valid_gems = []
    for gem_id, gem_info in grouped_gems.items():
        leftmost_fragment_start = gem_info['min_start']
        rightmost_fragment_end = gem_info['max_end']
        gem_size = gem_info['gem_size']
        gem_length = rightmost_fragment_end - leftmost_fragment_start

        if (
            leftmost_fragment_start <= left_anchor_end and
            leftmost_fragment_start >= left_anchor_start and
            rightmost_fragment_end <= right_anchor_end and
            rightmost_fragment_end >= right_anchor_start and
            len(gem_info['fragments']) == gem_size
        ):
            # Get all fragments of the valid GEM
            fragments = [
                pybedtools.create_interval_from_list(fragment)
                for fragment in gem_info['fragments']
            ]
            valid_gems.append((gem_id, fragments, gem_length))

    print("Sort the valid GEMs by their length")
    valid_gems.sort(key=lambda x: x[2])
    print(valid_gems)
    return valid_gems


def process_right(ChIA_Drop_old, left_anchor, right_anchor, region):
    left_anchor_chrom = left_anchor.split('\t')[0]
    left_anchor_start = int(left_anchor.split('\t')[1])
    left_anchor_end = int(left_anchor.split('\t')[2])
    right_anchor_end = int(right_anchor.split('\t')[2])

    print("Filter GEMs to only include those in the region")
    ChIA_Drop = ChIA_Drop_old.intersect(BedTool(region, from_string=True), wa=True, wb=True)

    print("Create BedTool object for the right anchor")
    right_anchor_bed = BedTool(right_anchor, from_string=True)

    print("Get the GEM IDs that intersect with the right anchor")
    intersecting_fragments = ChIA_Drop.intersect(right_anchor_bed, wa=True, wb=True)
    intersecting_gem_ids = set(fragment.fields[4] for fragment in intersecting_fragments)
    print("Filter ChIA_Drop to only include GEMs with intersecting IDs")
    ChIA_Drop = ChIA_Drop.filter(lambda x: x.fields[4] in intersecting_gem_ids)

    print("Filter ChIA_Drop to only include GEMs with the correct chrom id")
    ChIA_Drop = ChIA_Drop.filter(lambda x: x.chrom == left_anchor_chrom)

    print("Group GEM fragments by their GEM ID and get the min start and max end positions")
    grouped_gems = {}
    for fragment_interval in ChIA_Drop:
        fragment = fragment_interval.fields
        gem_id = fragment[4]
        gem_size = int(fragment[3])
        start = int(fragment[1])
        end = int(fragment[2])

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

    print("Add valid gems")
    valid_gems = []
    for gem_id, gem_info in grouped_gems.items():
        leftmost_fragment_start = gem_info['min_start']
        rightmost_fragment_end = gem_info['max_end']
        gem_size = gem_info['gem_size']
        gem_length = rightmost_fragment_end - leftmost_fragment_start

        if (
            rightmost_fragment_end <= right_anchor_end
            and leftmost_fragment_start >= left_anchor_start
            and len(gem_info['fragments']) == gem_size
        ):
            # Get all fragments of the valid GEM
            fragments = [
                pybedtools.create_interval_from_list(fragment)
                for fragment in gem_info['fragments']
            ]
            valid_gems.append((gem_id, fragments, gem_length))

    print("Sort the valid GEMs by their length")
    valid_gems.sort(key=lambda x: x[2])
    print(valid_gems)
    return valid_gems
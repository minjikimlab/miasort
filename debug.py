import pybedtools
from pybedtools import BedTool


def test_gem_id_difference(ranked_gems, correct_input_name):
    answer = BedTool(correct_input_name)
    correct_gem_ids = [gem.fields[3] for gem in answer]
    curr_gem_ids = [t[0] for t in ranked_gems]
    show_disparity(correct_gem_ids, curr_gem_ids)



def show_disparity(a, b):
    # Convert lists to sets
    set_a = set(a)
    set_b = set(b)

    # Find elements in a but not in b
    diff_a_b = set_a - set_b

    # Find elements in b but not in a
    diff_b_a = set_b - set_a

    # Combine both differences
    difference = list(diff_a_b | diff_b_a)

    print("Elements in a but not in b:", diff_a_b)
    print("Elements in b but not in a:", diff_b_a)
    print("All differing elements:", difference)

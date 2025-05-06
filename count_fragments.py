import os
import mmap
from collections import Counter
from concurrent.futures import ProcessPoolExecutor

# same problem as count_complexes
# trying to get histogram of number of fragments per complex

def generate_histogram(file_path):
    """
    Generate a histogram of the number of fragments per complex (column 4),
    using mmap for efficient large file reading.
    """
    fragment_counts = Counter()

    with open(file_path, 'r') as f:
        with mmap.mmap(f.fileno(), length=0, access=mmap.ACCESS_READ) as mm:
            for line in iter(mm.readline, b''):
                line = line.decode()
                if line.startswith("#") or not line.strip():
                    continue
                fields = line.strip().split()
                try:
                    num_fragments = int(fields[3])  # Column 4 (index 3)
                    fragment_counts[num_fragments] += 1
                except (IndexError, ValueError):
                    continue  # Handle malformed lines gracefully

    normalized_counts = {k: v / k for k, v in fragment_counts.items()}
    return os.path.basename(file_path), fragment_counts, normalized_counts


def process_directory_parallel(directory):
    """
    Process all .complexes files in the directory in parallel.
    """
    files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".complexes")]
    results = {}

    with ProcessPoolExecutor(max_workers=4) as executor:
        for file_name, raw_counts, norm_counts in executor.map(generate_histogram, files):
            results[file_name] = {
                "fragment_counts": raw_counts,
                "normalized_counts": norm_counts
            }

    return results


# Example usage
if __name__ == "__main__":
    directory = "../mia-sort_output/"
    results = process_directory_parallel(directory)
    for file_name, counts in results.items():
        print(f"{file_name}: {counts['normalized_counts']}")

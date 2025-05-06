from concurrent.futures import ProcessPoolExecutor
import os

# trying to count unique complexes in multiple files in parallel
# problem is that in situ hi-c is 271g and run out of memory
# maybe dont need to do it in parallel
# keep upping the memory limit until it works?

def count_unique_complexes(file_path):
    unique_complexes = set()
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split()
            if len(fields) > 4:
                unique_complexes.add(fields[4])
    return file_path, len(unique_complexes)

def count_in_directory_parallel(directory):
    results = {}
    files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".complexes")]
    with ProcessPoolExecutor(max_workers=4) as executor:
        for file_path, count in executor.map(count_unique_complexes, files):
            results[os.path.basename(file_path)] = count
    return results

# Example usage
directory = "../mia-sort_output/"
results = count_in_directory_parallel(directory)
for file_name, count in results.items():
    print(f"{file_name}: {count} unique complexes")

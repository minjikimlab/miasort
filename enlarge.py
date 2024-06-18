import sys

def copy_text(source_file, target_file, num_times):
    try:
        # Open the source file in read mode
        with open(source_file, 'r') as src:
            # Read the content of the source file
            content = src.read()

        # Open the target file in write mode
        with open(target_file, 'w') as tgt:
            for _ in range(num_times):
                # Write the content to the target file
                tgt.write(content)

        print(f"Content from {source_file} has been copied to {target_file}.")

    except FileNotFoundError:
        print(f"Error: The file {source_file} was not found.")
    except IOError as e:
        print(f"Error: {e}")

if __name__ == "__main__":
    num_times = 2  # user defined
    if len(sys.argv) != 3:
        print("Usage: python copy_text_twice.py <source_file> <target_file>")
    else:
        source_file = sys.argv[1]
        target_file = sys.argv[2]
        copy_text(source_file, target_file, num_times)
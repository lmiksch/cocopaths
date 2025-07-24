import ast

def load_paths_from_file(filename):
    """
    Reads paths from a text file.
    
    Each line in the file should be a valid Python literal (e.g. a list of strings)
    representing one path.
    
    Returns:
        A list of paths (each path is a list).
    """
    paths = []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line:
                try:
                    # Safely evaluate the string to a Python object.
                    path = ast.literal_eval(line)
                    if isinstance(path, list):
                        paths.append(path)
                    else:
                        print(f"Warning: Skipping line (not a list): {line}")
                except Exception as e:
                    print(f"Error parsing line: {line}\n{e}")
    return paths

def compare_path_sets_from_files(file1, file2, output_shared_file="shared_paths.txt"):
    """
    Loads paths from two files, compares them, prints the comparison,
    and writes the shared paths to a file.
    
    Args:
        file1 (str): Path to the first text file.
        file2 (str): Path to the second text file.
        output_shared_file (str): Output file to save shared paths.
    """
    paths1 = load_paths_from_file(file1)
    paths2 = load_paths_from_file(file2)
    
    set1 = {tuple(path) for path in paths1}
    set2 = {tuple(path) for path in paths2}
    
    common = set1.intersection(set2)
    only_in_file1 = set1 - set2
    only_in_file2 = set2 - set1
    
    print("Common paths (appear in both files):")
    if common:
        for path in sorted(common):
            print(list(path))
    else:
        print("None")
    
    print("\nPaths only in", file1, ":")
    if only_in_file1:
        for path in sorted(only_in_file1):
            print(list(path))
    else:
        print("None")
        
    print("\nPaths only in", file2, ":")
    if only_in_file2:
        for path in sorted(only_in_file2):
            print(list(path))
    else:
        print("None")




# Example usage:
if __name__ == "__main__":
    file1 = "6_steps_fixed.txt"          # Replace with your actual file name/path.
    file2 = "6_FPs_new_consistency.txt"       # Replace with your actual file name/path.
    compare_path_sets_from_files(file1, file2)

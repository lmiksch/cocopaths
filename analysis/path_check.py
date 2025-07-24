#!/usr/bin/env python3
import multiprocessing, copy
from cocopaths.utils import path_to_pairtablepath
from cocopaths.cocopath import *  # This brings in build_graph, graph, etc.
import os
import io
import contextlib,math
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

# Try to import upsetplot and pandas for an additional UpSet plot.
try:
    import pandas as pd
    from upsetplot import UpSet, from_indicators
    HAVE_UPSETPLOT = True
except ImportError:
    HAVE_UPSETPLOT = False

# -------------------------------------------------------------------
# Helper functions (unchanged or adapted)
# -------------------------------------------------------------------

def generate_structures_with_unpaired(n):
    """
    Generate all possible dot‐bracket structures for a given length n,
    allowing unpaired positions.
    """
    structures = []
    def helper(open_count, close_count, dot_count, current_structure):
        if open_count + close_count + dot_count == n and open_count == close_count:
            structures.append(current_structure)
            return

        if open_count < n // 2:
            helper(open_count + 1, close_count, dot_count, current_structure + '(')

        if close_count < open_count:
            helper(open_count, close_count + 1, dot_count, current_structure + ')')

        if open_count + close_count + dot_count < n:
            helper(open_count, close_count, dot_count + 1, current_structure + '.')
    helper(0, 0, 0, '')
    return structures

def coco_legal(path):
    """
    Checks if a given abstract cotranscriptional folding path (aCFP) is legal.
    Returns a dictionary with keys:
      - 'consistent'
      - 'assignable'
      - 'favorable'
    """
    checks = dict()
    acfp = path

    # Build the pairtable if necessary.
    if acfp[0][0] == ".":
        pairtable_acfp = path_to_pairtablepath(acfp)
    else:
        pairtable_acfp = acfp

    acfp_graph = graph()
    acfp_graph.create_nodes_from_pairtable(pairtable_acfp)
    # Set the global node_list so that static methods (e.g. check_complementarity) work.
    graph.node_list = sorted(list(acfp_graph.graph.keys()), key=lambda node: node.name)

    nodes = list(acfp_graph.graph.keys())
    nodes.insert(0, 0)

    acfp_graph.create_edges(pairtable_acfp, nodes)
    acfp_graph.get_edges()
    acfp_graph.set_edge_neighbors()

    connected_components = find_connected_modules(pairtable_acfp)
    for x, component in enumerate(connected_components):
        for node_name in component:
            acfp_graph.get_node_by_name(node_name).connected = x
    
    # Assignability check.
    if acfp_graph.bipartite_check(connected_components):
        checks['assignable'] = True
    else:
        checks['assignable'] = False

    # Consistency check.
    if acfp_graph.consistency_check(acfp):
        checks['consistent'] = True
    else:
        checks['consistent'] = False


    # Favorability check.
    if acfp_graph.is_favorable(acfp):
        checks['favorable'] = True
    else:
        checks['favorable'] = False

    return checks



def count_paths(n):
    """
    Return the number of distinct paths of length n
    you'd get from your generate_paths routine,
    without ever enumerating them.
    """

    possible_structs = [generate_structures_with_unpaired(x + 2) for x in range(n - 1)]

    # if n == 1, your only path is ["."], so count = 1
    if n <= 1:
        return 1

    # multiply together the sizes of each step's choice-list
    # possible_structs[0] is used to go from length 1 to 2, etc.
    return math.prod(len(possible_structs[i]) for i in range(n-1))

def generate_paths(n, possible_structs, current_path=None):
    """
    Recursively generate all paths of length n using the precomputed
    possible structures for each step and yield one complete path at a time.
    """
    if current_path is None:
        current_path = ["."]
    if len(current_path) == n:
        yield current_path
    else:
        step = len(current_path)
        for pos_struc in possible_structs[step - 1]:
            new_path = current_path + [pos_struc]
            yield from generate_paths(n, possible_structs, new_path)


def process_path(path):
    with io.StringIO() as buf, contextlib.redirect_stdout(buf):
        # Reinitialize the global node_list explicitly:
        if path[0][0] == ".":
            pairtable_acfp = path_to_pairtablepath(path)
        else:
            pairtable_acfp = path
        acfp_graph = graph()
        acfp_graph.create_nodes_from_pairtable(pairtable_acfp)
        # Explicitly set the global node_list for use in static methods:
        graph.node_list = sorted(list(acfp_graph.graph.keys()), key=lambda node: node.name)
        
        # Now run the legal checks
        try:
            check_results = coco_legal(path)
        except SystemExit:
            check_results = {'consistent': False, 'assignable': False, 'favorable': False}
        try:
            translate_acfp(path)
            sls_translatable = True
        except (SystemExit, Exception) as e:
            print(f"translate_acfp raised {type(e).__name__}: {e}", file=sys.stderr)
            sls_translatable = False
        check_results['sls_translatable'] = sls_translatable
    return (path, check_results)


def generate_coco_legal_paths(n):
    """
    Generates all possible paths of length n, processes them in parallel,
    writes translatable paths to a file, and yields each processed result.
    """
    possible_structs = [generate_structures_with_unpaired(x + 2) for x in range(n - 1)]
    path_generator = generate_paths(n, possible_structs)
    output_file = f"{n}_steps_all.txt"
    with open(output_file, "a") as f_all:
        with multiprocessing.Pool() as pool:
            for path, check_results in pool.imap_unordered(process_path, path_generator, chunksize=100):
                if check_results.get('sls_translatable'):
                    f_all.write(f"{str(path)}\n")
                yield (path, check_results)

# -------------------------------------------------------------------
# Main analysis: Process lengths 3-7, create figures, and compare counts
# -------------------------------------------------------------------
def main():
    results_folder = "analysis_results"
    os.makedirs(results_folder, exist_ok=True)
    
    overall_summary_lines = []
    # Dictionary to store paths that are NOT favorable or NOT SLS-translatable, keyed by length.
    not_ok_paths_all = {}
    
    # Dictionary to store all paths by their 4-tuple combination key (C, A, F, S)
    overall_combo_paths = {}
    
    # Dictionary to store counts for each length for later comparison.
    counts_by_length = {}
    
    # Loop over lengths 3 to 7.
    for n in range(4, 7):
        print(f"Analyzing paths for length {n} ...")
        total_paths_processed = 0
        
        # Individual counts.
        consistent_count = 0
        assignable_count = 0
        favorable_count = 0
        sls_count = 0  # SLS-translatable count
        
        # Combination counts: key is a 4-tuple: (consistent, assignable, favorable, sls_translatable)
        combo_counts = {}
        # Dictionary mapping combination key to list of paths (for this length).
        combo_paths = {}
        
        # List to store paths that are NOT favorable or NOT SLS-translatable.
        # (Here we record paths that satisfy Consistent and Assignable but fail either favorable or SLS-translatable.)
        not_ok_paths = []
        
        output_file = os.path.join(results_folder, f"{n}_steps_all.txt")
        
        for path, results in generate_coco_legal_paths(n):
            total_paths_processed += 1
            print(path)
            # Get booleans from results; default to False if missing.
            c = bool(results.get('consistent', False))
            a = bool(results.get('assignable', False))
            f = bool(results.get('favorable', False))
            s = bool(results.get('sls_translatable', False))
            
            # Update individual counts.
            if c: 
                consistent_count += 1
            if a:
                assignable_count += 1
            if f:
                favorable_count += 1
            if s:
                sls_count += 1

            # Update combination counts and store corresponding paths.
            combo = (c, a, f, s)
            combo_counts[combo] = combo_counts.get(combo, 0) + 1
            if combo not in combo_paths:
                combo_paths[combo] = []
            combo_paths[combo].append(path)
            
            # If the path is both consistent and assignable but is not favorable or not SLS-translatable, record it.
            if c and a and not (f and s):
                not_ok_paths.append(path)

            if total_paths_processed % 100000 == 0:
                print(f"Processed {total_paths_processed} paths for length {n}.")

        # Save the not_ok_paths for this length.
        not_ok_paths_all[n] = not_ok_paths
        # Merge combo_paths into overall_combo_paths.
        for combo, paths in combo_paths.items():
            overall_combo_paths[combo] = overall_combo_paths.get(combo, []) + paths

        # Save counts for this length for later comparison.
        counts_by_length[n] = {
            'total': total_paths_processed,
            'consistent': consistent_count,
            'assignable': assignable_count,
            'favorable': favorable_count,
            'sls': sls_count
        }
        
        # Build summary for this length.
        summary_lines = [
            f"--- Summary for length {n} ---",
            f"Total paths processed: {total_paths_processed}",
            f"Consistent paths: {consistent_count}",
            f"Assignable paths: {assignable_count}",
            f"Favorable paths: {favorable_count}",
            f"SLS-translatable paths: {sls_count}",
            "Combination counts (C, A, F, S):"
        ]
        # Print all 16 possible combinations.
        all_combos = [(c, a, f, s) for c in (False, True)
                               for a in (False, True)
                               for f in (False, True)
                               for s in (False, True)]
        for combo in sorted(all_combos):
            count = combo_counts.get(combo, 0)
            summary_lines.append(f"  {combo}: {count}")
        summary_text = "\n".join(summary_lines)
        print(summary_text)
        overall_summary_lines.append(summary_text)
        
        # --- Venn diagram using venn3 ---
        # We create a Venn diagram for three sets: Assignable, Consistent, and Favorable.
        # Since we want to ignore SLS in the diagram, we sum over both SLS values.
        def sum_combo(combo_dict, c, a, f):
            total = 0
            for s in (False, True):
                total += combo_dict.get((c, a, f, s), 0)
            return total

        only_assignable = sum_combo(combo_counts, False, True, False)
        only_consistent = sum_combo(combo_counts, True, False, False)
        only_favorable  = sum_combo(combo_counts, False, False, True)
        assignable_consistent_only = sum_combo(combo_counts, True, True, False)
        assignable_favorable_only = sum_combo(combo_counts, False, True, True)
        consistent_favorable_only = sum_combo(combo_counts, True, False, True)
        triple = sum_combo(combo_counts, True, True, True)
        

        print(f"{only_assignable = }")
        print(f"{only_consistent = }")
        print(f"{only_favorable = }")
        print(f"{assignable_consistent_only = }")
        print(f"{assignable_favorable_only = }")
        print(f"{consistent_favorable_only = }")
        print(f"{triple = }")

        # venn3 expects the subsets in the order:
        # (only A, only B, only C, A∩B only, A∩C only, B∩C only, A∩B∩C)
        venn_data = (
            only_assignable,            # Only in Assignable
            only_consistent,            # Only in Consistent
            only_favorable,             # Only in Favorable
            assignable_consistent_only, # In A∩B only (Assignable and Consistent, but not Favorable)
            assignable_favorable_only,  # In A∩F only (Assignable and Favorable, but not Consistent)
            consistent_favorable_only,  # In C∩F only (Consistent and Favorable, but not Assignable)
            triple                      # In all three sets
        )
        plt.figure(figsize=(6,6))
        v = venn3(subsets=venn_data, set_labels=('Assignable', 'Consistent', 'Favorable'))
        for text in v.set_labels:
            if text:
                text.set_fontsize(12)
        for text in v.subset_labels:
            if text:
                text.set_fontsize(10)
        plt.tight_layout()
        plt.title(f"Venn Diagram for Length {n}\nSLS-translatable count: {sls_count}")
        venn_filename = os.path.join(results_folder, f"length_{n}_venn.pdf")
        plt.savefig(venn_filename)
        plt.clf()

        # --- Bar Chart ---
        categories = ['Total', 'Consistent', 'Assignable', 'Favorable', 'SLS-translatable', 'None']
        none_count_val = combo_counts.get((False, False, False, False), 0)
        counts = [
            total_paths_processed,
            consistent_count,
            assignable_count,
            favorable_count,
            sls_count,
            none_count_val
        ]
        plt.figure(figsize=(8,5))
        bars = plt.bar(categories, counts, color='skyblue')
        plt.title(f"Path Count Summary for Length {n}")
        plt.ylabel("Count")
        for bar in bars:
            height = bar.get_height()
            plt.annotate(f'{height}',
                         xy=(bar.get_x() + bar.get_width() / 2, height),
                         xytext=(0, 3),
                         textcoords="offset points",
                         ha='center', va='bottom')
        bar_filename = os.path.join(results_folder, f"length_{n}_bar.pdf")
        plt.savefig(bar_filename)
        plt.clf()

    # Write overall summary to a text file.
    summary_file = os.path.join(results_folder, "summary.txt")
    with open(summary_file, "w") as f:
        f.write("\n\n".join(overall_summary_lines))
    print(f"Analysis complete. Results and figures are saved in the folder '{results_folder}'.")
        
    # --- Compare Counts Across Lengths ---
    comparison_lines = []
    header = "Length\tTotal\tConsistent\tAssignable\tFavorable\tSLS-translatable\tFavorable==SLS?"
    comparison_lines.append(header)
    for n in sorted(counts_by_length.keys()):
        data = counts_by_length[n]
        equal_flag = "Yes" if data['favorable'] == data['sls'] else "No"
        line = f"{n}\t{data['total']}\t{data['consistent']}\t{data['assignable']}\t{data['favorable']}\t{data['sls']}\t{equal_flag}"
        comparison_lines.append(line)
        # Append the combination counts for this length:
        comparison_lines.append("Combination counts (C, A, F, S):")
        for combo in sorted(counts_by_length[n]):
            comparison_lines.append(f"  {combo}: {counts_by_length[n][combo]}")
        comparison_lines.append("")  # Blank line for separation
    comparison_text = "\n".join(comparison_lines)
    comparison_file = os.path.join(results_folder, "counts_comparison.txt")
    with open(comparison_file, "w") as cf:
        cf.write(comparison_text)
    print("\nCounts comparison:\n", comparison_text)
    
    #Print out the paths for the combination that are consistent, assignable, SLS-translatable but NOT favorable.
    target_combo = (True, True, False, False)
    print(f"\nPaths for combination {target_combo} (count: {len(overall_combo_paths.get(target_combo, []))}):")
    for path in overall_combo_paths.get(target_combo, []):
        print(path)
    
    ## --- EXTRA ---
    #Print out the paths that are consistent, assignable, and favorable but NOT SLS-translatable.
    target_combo_extra = (False, True, True, False)
    print(f"\nPaths for combination {target_combo_extra} (count: {len(overall_combo_paths.get(target_combo_extra, []))}):")
    for path in overall_combo_paths.get(target_combo_extra, []):
            process_path(path)
            print(path)
##
#
    target_combo_extra = (True, True, True, True)
    print(f"\nPaths for combination {target_combo_extra} (count: {len(overall_combo_paths.get(target_combo_extra, []))}):")
    for path in overall_combo_paths.get(target_combo_extra, []):
        if len(path) == 6:
            process_path(path)
            print(path)
    #target_combo_extra = (True, True, False, False)
    #print(f"\nPaths for combination {target_combo_extra} (count: {len(overall_combo_paths.get(target_combo_extra, []))}):")
    #for path in overall_combo_paths.get(target_combo_extra, []):
    #    if len(path) == 6:
    #        process_path(path)
    #        print(path)
            
if __name__ == "__main__":
    # Uncomment to run the full analysis:
    main()

    #for n in range(1,9):
    #    print(f'For length {n}')
    #    print(count_paths(n))
        
        
"""
    test_acfps = [
        ['.', '()', '().', '....', '(()).', '(()())'],
        ['.', '..', '.()'],
        ['.', '()', '.()', '()()', '(()).', '(())()'],
        ['.', '()', '().', '.().', '(()).', '(())()'],
        ['.', '()', '().', '(())', '(()).', '(())..'],
        ['.', '()', '().', '(())', '(()).', '(()..)'],
    ]

    print("Running favorability check on the test set...\n")
    for i, acfp in enumerate(test_acfps, start=1):
        print(f"\n=== Test path {i} ===")
        try:
            result = build_graph(acfp)
            print(f"\nTest path {i}: {acfp} --> {'FAVORABLE' if result else 'NOT favorable'}")
        except SystemExit as e:
            print(f"\nTest path {i}: {acfp} --> SystemExit caught: {e}")
        except Exception as err:
            print(f"\nTest path {i}: {acfp} --> ERROR: {err}")
"""
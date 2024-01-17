"""
This module creates all possible folding paths given a domain sequence. This is mainly used to test our path_to_seq algortihm.
"""
import random

import math

from cocopaths.utils import path_to_pairtablepath

def generate_structures_with_unpaired(n):
    structures = [ ]
    def helper(open_count, close_count, dot_count, current_structure):
        if open_count + close_count + dot_count == n and open_count == close_count:
            structures.append(current_structure)
            return current_structure

        if open_count < n // 2:
            helper(open_count + 1, close_count, dot_count, current_structure + '(')

        if close_count < open_count:
            helper(open_count, close_count + 1, dot_count, current_structure + ')')

        if open_count + close_count + dot_count < n:
            helper(open_count, close_count, dot_count + 1, current_structure + '.')

    helper(0, 0, 0, '')

    return structures


def is_legal_fp(struct,fp):

    """
    Checks if the structure is legal for addition.

    Args:
        struct(str): Dot-bracket annotation of the to be added structure
        fp(list): List where each entry represents a folding step 
    """

    """" .() .(). 
     
         ()..
       
         """
    #print("\n\nBegin is legal FP ",struct,fp)
    defined_pairs = []
    
    for i,step in enumerate(fp[1:],start=2): 
        if hamming_distance(step,struct[:i]) < 2 or step == struct or step[-2:] == "()" and struct[-3:-1] == ".(":
            continue
        else:
            return False
        
    return True


def hamming_distance(str1, str2):
    #print("STr1 and 2",str1,str2)
    if len(str1) != len(str2):
        raise ValueError("Input strings must have the same length")

    return sum(c1 != c2 for c1, c2 in zip(str1, str2))






def generate_path(n):

    final_paths = [['.']]
    for _ in range(n-1):
        final_paths.append([])



    possible_structs = [generate_structures_with_unpaired(x + 2) for x in range(n-1)]
    #print("Possible Sturcts", possible_structs)
    
    for x in range(1, n):
        #print("P_structs", possible_structs[x - 1])
        new_paths = []
        for curr_path in final_paths[x - 1]:
            #print("\n\nCurrent Path: ", curr_path, "X = ",x)
            for pos_struc in possible_structs[x - 1]:
                #print('\npossible struct', pos_struc)
                #if is_legal_fp(pos_struc, curr_path):
                new_path = list(curr_path).copy()
                new_path.append(pos_struc)
                final_paths[x].append(new_path)
            #print(f"New paths{new_paths}")
        final_paths.extend(new_paths)
                        
    print("\n\nDone\n\n")
    for step in final_paths:
        for path in step:
            print(path)     

    return final_paths        
    
    




if __name__ == "__main__":
    # Example: generate structures for length 4
    generate_path(4)


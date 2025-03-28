import argparse
from cocopaths.cocopath import translate_acfp
from cocopaths.cocodesign import rna_design,domain_path_to_nt_path,acfp_to_domainfp,constrained_efe,objective_function,call_findpath,extend_domain_seq,score_sequence,toehold_structures
from cocopaths.utils import only_logic_domain_struct
import RNA
import numpy as np 
import subprocess,os,re
import pandas as pd
import ast
from filelock import FileLock
import string


def get_default_parameters():

    # Combine lowercase and uppercase letters
    letters = string.ascii_lowercase 

    # Extend each letter with a '*' without space between the letter and '*'
    extended_letters = ' '.join(f"{letter}*" for letter in letters)

    # Generate L0 to L100 and S0 to S100
    l_values = ' '.join(f"L{num}" for num in range(101))
    l_star_values = ' '.join(f"L{num}*" for num in range(101))

    s_values = ' '.join(f"Z{num}" for num in range(101))

    # Combine the original letters, extended letters, and L/S sequences
    result = ' '.join(letters) + ' ' + extended_letters + ' ' + l_values + ' ' + s_values + ' ' + l_star_values

    domain_seq_fp2 = 'a b c d e f g h i j k l m n o p q r s t u v w x y z a* b* c* d* e* f* g* h* i* j* k* l* m* n* o* p* q* r* s* t* u* v* w* x* y* z* L0 L1 L2 L3 L4 L5 L6 L7 L8 L9 L10 L11 L12 L13 L14 L15 L16 L17 L18 L19 L20 L21 L22 L23 L24 Z0 Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20 S21 S22 S23 S24 S25'

    domain_seq = result
    d_length = {}
    
    for domain in domain_seq.split():
        if domain[0] == "L":
            d_length[domain] = 8
        elif domain[0] == 'Z':
            d_length[domain] = round(int(domain[1]) * 4)  
        else: 
            d_length[domain] = 3
    parameters = {"k_slow": 0.00001,'k_fast': 20, "cutoff": 0.05,"d_length":d_length,"d_seq":domain_seq,"logic":True}

    return parameters

def evaluate_drt_output(dr_out_file, extended_fp):
    """
    Evaluate the DrTransformer output file and calculate the populations for each step
    in the nucleotide folding path. For eachs step corresponding to the step in the domain 
    level sequence the max occupancy of all structures satisfying the structural constraints is taken. 
    Since multiple time points for the specific lengths are available. 

    Args:
        dr_out_file (str): Path to the DrTransformer output file (e.g., *.drf).
        extended_fp (list): The nucleotide folding path, where each step represents
                            the target structure for comparison.

    Returns:
        populations (list): Populations for each folding step.
    """
    # Read the DrTransformer output
    with open(dr_out_file, 'r') as file:
        dr_out = [line.rstrip().split() for line in file]

    # Initialize data structures for tracking populations and seen structures
    populations = [0 for _ in range(len(extended_fp))]
    seen_drt_out = [[] for _ in range(len(extended_fp))]
    new_time = True
    # Process DrTransformer output
    i,j = 0,0  # Index for the current folding step
    cur_time = 0
    curr_pop = [0]
    for line in dr_out[1:]:  # Skip the header if present


        structure = line[3]  # Extract structure
        probability = float(line[2])  # Extract probability
        time = line[1]

        
        # Check if the structure belongs to the current folding step
        if len(structure) > len(extended_fp[i]):
            populations[i] += max(curr_pop)
            curr_pop = [0]
            i += 1  # Move to the next folding step

        if len(structure) == len(extended_fp[i]): 
            if time != cur_time:
    
                cur_time = time
                curr_pop.append(0)
        
        if drt_match(structure, extended_fp[i]):# and structure not in seen_drt_out[i]:
            #populations[i] += probability
            curr_pop[-1] += probability
            seen_drt_out[i].append(structure)
            if populations[i] > 1:
                raise SystemExit(f'{structure = }  {extended_fp[i] = }')

    #Add population of last step 
    populations[i] += max(curr_pop)

    # Print results for debugging
    print("Populations", "Structure")
    for step, population in zip(extended_fp, populations):
        print(population, step)

    return populations



def drt_match(drt_out_struct,fp_struct):
    if len(drt_out_struct) != len(fp_struct)    :
        return False
    for drt_char,fp_char in zip(drt_out_struct,fp_struct):
        if fp_char == ".":
            continue
        elif drt_char != fp_char:
            return False 
    return True

'''
def drt_match(drt_out_struct,fp_struct,parameters):

    #
    if 'x' in fp_struct: 
        if len(drt_out_struct) != len(fp_struct)    :
            return False
        for drt_char,fp_char in zip(drt_out_struct,fp_struct):
            if fp_char == ".":
                continue
            elif drt_char != fp_char:
                return False 
        return True
    else:    
        split_fp_struct = split_extend_structure(domain_seq,parameters['domain_lengths'],fp_struct)
        split_drt_struct = split_extend_structure(domain_seq,parameters['domain_lengths'],drt_out_struct)

        for fp_sub_struct, drt_sub_struct in zip(split_fp_struct,split_drt_struct):
            if '.' in fp_sub_struct and not is_legal_dotbracket(drt_sub_struct):
                return False
            else:
                if fp_sub_struct != drt_sub_struct: 
                    return False
        return True

'''
 
if __name__ == '__main__': 
    '''
    acfp = ['.', '..', '...', '....', '..(.)', '....()']
    d_seq_old = 'L0* Z0 L1* Z1 a* L2* b* Z2 L3* Z3 c* b L2 a d* Z4 d a* L2* b* c Z5'
    d_seq_new = '''

    acfp = ['.', '..', '...', '..()', '..(.)', '..(..)']
    d_seq = 'L0*  Z0  L1*  Z1 a* b* L2* c* d* Z2  L2  Z3 c L2 b Z4 d c L2 b a Z5'

    parameters = get_default_parameters()

    domain_fp = acfp_to_domainfp(acfp,d_seq)
    print(domain_fp)
    extended_fp = domain_path_to_nt_path(domain_fp,d_seq,parameters)
    extended_fp = [x.replace('x','.') for x in extended_fp]
    print(extended_fp)
    evaluate_drt_output('/home/mescalin/miksch/Documents/test/drtrafo/NoName.drf', extended_fp)
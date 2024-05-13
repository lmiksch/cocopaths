# -*- coding: utf-8 -*-
"""
This module contains all the different util functions. 
"""

def is_balanced_structure(s):
    stack = []
    for char in s:
        if char == '(':
            stack.append('(')
        elif char == ')':
            if not stack:
                return False  # Unmatched closing parenthesis
            stack.pop()
    
    return len(stack) == 0

def path_to_pairtablepath(path):
    """ Function which turns a structural path in dot bracket annotation into a pairtable path

    args: 
        path(list): List of the folding path in dot bracket
    output: 
        pairtable_path(List): List of the folding path using the pairtable annotation

    """
    pairtable_path = [[] for x in path]
    pairtable_path[0] = [1,0]
    for x in range(1,len(path)):
        pairtable_path[x] = [0 for i in range(len(path[x])+1)]
        pairtable_path[x][0] = len(path[x]) 
        for z in range(len(path[x])):
            if path[x] == ".":
                pairtable_path[x][z] = 0
            o = 1
            i = 0
            match = False
            if path[x][z] == "(":
                while match == False:
                    i += 1
                    if path[x][z+i] == "(":
                        o += 1
                    elif path[x][z+i] == ")" and o != 0:
                        o -= 1

                    if path[x][z+i] == ")" and o == 0:
                        
                        pairtable_path[x][z+1] = z + i + 1
                        pairtable_path[x][z+i+1] = z + 1
                        match = True
                    
    return pairtable_path


def find_connected_modules(afp):
    """ Identifies connected Modules and returns the indices of the connected modules in a list. Each sublist corresponds to connected modules. 


    Args: 
        afp(list): abstract folding path in the form of a pairtable
    
    Returns: 
        connected_modules(list): sublist corresponds to connected modules with indices 
    """

    all_module = list(range(1,len(afp)+1))
    connected_modules = []
    curr_module = 1
    i,j = 0,0 
    while all_module:
        all_list = [all_module[0]]
        connected_modules.append(all_list)
        all_module.remove(all_module[0])
        
        while curr_module <= len(afp):
            for step in range(1,len(afp)): 
                
                try:
                    if afp[step][curr_module] != 0 and afp[step][curr_module] in all_module:
                        connected_modules[i].append(afp[step][curr_module])
                        all_module.remove(afp[step][curr_module])
                except:
                    continue
            
            j += 1 

            try:
                curr_module = connected_modules[i][j]
                

            except:
                if len(all_module) == 0:    
                    return connected_modules
                i += 1 
                j = 0 
                all_list = [all_module[0]]
                connected_modules.append(all_list)
                all_module.remove(all_module[0])
                curr_module = connected_modules[i][j]
                
                
    print("Error Something went wrong please check input")


def cv_db2kernel(domain_seq,dot_bracket_seq):
    """
    Takes a dot-bracket sequence and a domain level sequence and converts it into a kernel annotated sequence.

    Args:
        domain_seq(str): domain level seq with space annotation
        dot_bracket_seq(str): dot bracket seq 

    Returns: 
        kernel_seq(str): kernel annotated string with space annotation
    """
    

    
    split_dot_bracket = list(dot_bracket_seq)

    bps = dot_bracket_seq.count('(')
    kernel_annotation = [[] for x in split_dot_bracket]

    open_domains = [] 

    for i, i_domain in enumerate(domain_seq):

        if split_dot_bracket[i] == "(":
            kernel_annotation[i].append(str(i_domain + "("))

        if split_dot_bracket[i] == ")":
            kernel_annotation[i].append(")")
    
        if split_dot_bracket[i] == "." or split_dot_bracket[i] == "x":
            kernel_annotation[i].append(domain_seq[i])
        
    joined_kernel_annotation = []

    


    for x in kernel_annotation:
        joined_kernel_annotation.append(x[0])

    

    return " ".join(joined_kernel_annotation)


def kernel_to_dot_bracket(kernel_seq):
    """
    Takes a kernel annotated sequence and converts it back to a dot-bracket sequence.

    Args:
        kernel_seq(str): kernel annotated string with space annotation

    Returns:
        dot_bracket_seq(str): dot bracket sequence
    """
    kernel_annotations = kernel_seq.split()
    dot_bracket_seq = ""

    for annotation in kernel_annotations:
        if "(" in annotation:
            dot_bracket_seq += "("
        elif ")" in annotation:
            dot_bracket_seq += ")"
        else:
            dot_bracket_seq += "."


    
    return dot_bracket_seq



def is_star_pair(a,b):
    """Checks if two domain pair in Star Annotation"""
    if a[0] == b[0]:
        if len(a) < len(b) or len(a) > len(b):
            return True
        
def only_logic_domain_struct(seq,path):
    """Takes a seq and the extended folding path only returns the folding path of only the b domains

    Args: 
        seq(list): sequence in form of a list where each entry in the list corresponds to one module
        path(str): module folding path in form of a list where each entry corresponds to one module being transcribed

    """
    logic_domain_struct = []
    for x in range(len(path)):
        
        if seq[x].startswith("L"):
            
            logic_domain_struct.append(path[x])

    
    logic_domain_struct = "".join(logic_domain_struct)
    return logic_domain_struct

def afp_terminal_input():

    '''Function to let user put in an abstract folding path in the terminal while checking if it's viable
    '''
    afp = []
    while True:
        print("\n")
        print(f"Current Input: {afp}")
        print("Please input a folding path in dot-bracket annotation or use '$' to exit input and continue use 'r' to reset input:")
        user_input = input()
        # Check for exit conditions
        if user_input == "$":
            print(f"\n\nFinal Input:\n{afp}\n\n")
            break
        elif user_input == "r" or user_input == "R":
            afp = []
            print("Input cleared")
            continue
        
        if is_balanced_structure(user_input):

            # Check if the user input contains only ".", "(", and ")"
            
            if all(char == "." or char in ("(", ")") for char in user_input):
                if len(user_input) == len(afp) + 1:
                    afp.append(user_input)
                else:
                    print("Please add 1 character per step")
            else:
                print("Error: Invalid character in the folding path. Only '.', '(', and ')' are allowed.")
        else:
            print("Structure is not balanced -> closing/opening brackets don't match")

    return afp
if __name__=="__main__":
    print(' ')
# -*- coding: utf-8 -*-
"""
This module contains all the different convert functions which can be used to convert outputs into different annotations. 
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
        
        while curr_module <= len(afp) + 1:
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
            
            



if __name__=="__main__":
    print(' ')
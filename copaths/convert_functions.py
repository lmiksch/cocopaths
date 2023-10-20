# -*- coding: utf-8 -*-
"""
This module contains all the different convert functions which can be used to convert outputs into different annotations. 
"""


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

def pairtable_to_path(pairtable):
    """Gives out a dot-bracket path 

        Args:
            pairtable (list): list of lists of a pairtable resembling a folding path
        
        Returns:
            path (list): sublits correspond to the given path but in pairtable format
    """
    path = []
    for x in range(len(pairtable)):
        struct = []
        
        for i in range(1,pairtable[x][0]+1):
            if pairtable[x][i] > i:
                struct.append("(")
            elif pairtable[x][i] < i and pairtable[x][i] != 0:
                struct.append(")")
            elif pairtable[x][i] == 0:
                struct.append(".")
        
        path.append("".join(struct))

    return path

def pairtable_to_struct(pairtable):
    """ Given a pairtable as input returns the structure in dot-bracket annotation

        Args:
            pairtable (list): pairtable
        Returns:
            struct (list): dot-bracket annotated structure in a list

    """

    struct = []
    for i in range(1,pairtable[0]+1):
        if pairtable[i] > i:
            struct.append("(")
        elif pairtable[i] < i and pairtable[i] != 0:
            struct.append(")")
        elif pairtable[i] == 0:
            struct.append(".")
    for x in range(1,len(struct)-1):
        if struct[x] == ".":
            return    
    struct = ("".join(struct))
    
    return struct

def convert_pts_out_to_nussi_in(string):
    """Removes spaces from string
    
    """


    new_string = string.replace(" ","")
    return new_string

def convert_to_UL(string):
    """Converts sting annotated domain sequence into a upper lower case annotated string

    Args: 
        string (string): space annotated domain level sequence

    Returns: 
        UL_seq (str): Upper lower case annotated domain level sequence
    
    """


    counts = string.count("*")
    c_string = list(string)
    
    for x in range(1,len(c_string)):
        if string[x] == "*":
            c_string[x-1] = c_string[x-1].upper()
    for x in range(counts):		
        c_string.remove("*")

    c_string = "".join(c_string)
    return c_string.replace(" ","")


def module_folding_path(structures,seq):
    """Takes output of nussinov algorithm and converts it just to the folding path of the b domain for comparison against the input path. 
        It will only consider adding the strucure after a module has fully transcribed. 

        Args: 
            Structure(list): Output of nussinov algorithm 
            seq(str): sequences which was used in the nussinov algorithm with space and star annotation
        
        Returns: 
            module_structure(list): 
    """
    
    modules = seq.split("l")
    for x in range(len(modules)):
        modules[x] = modules[x].split() 
        
    
    seq = convert(seq)
    module_path = []
    

    lengths = [0 for x in modules]

    for x in range(len(modules)):
        for z in range(x+1):
            lengths[x] = len(modules[z]) + lengths[x] 
    t = 1 
    for x in range(len(seq)):
        if seq[x] == "l":
            lengths[t] += t 
            t += 1
    


    for x in range(len(lengths)):
        module_path.append(structures[lengths[x]])
    liste = [[] for x in module_path]

    print(module_path)
    for x in range(len(module_path)):
        for z in range(len(module_path[x])):
            if seq[z] == "b" or seq[z] == "B":
                liste[x].append(module_path[x][z])
        
    for x in range(len(liste)):
        liste[x] = "".join(liste[x])
    return(liste)	
  
def only_b_domainfp(seq,path):
    """Takes a seq and the extended folding path only returns the folding path of only the b domains

    Args: 
        seq(list): sequence in form of a list where each entry in the list corresponds to one module
        path(list): module folding path in form of a list where each entry corresponds to one module being transcribed

    """
    b_domainfp = [[] for x in path]

    
    for x in range(len(path)):

        for y in range(len(path[x][0])):
           
            if seq[y] == "b" or seq[y] == "B":
                
                b_domainfp[x].append(path[x][0][y])

    for x in range(len(b_domainfp)):
        b_domainfp[x] = "".join(b_domainfp[x])
    return b_domainfp

def d_length(domain):
    if domain[0] == "b" or domain[0] == "B":
        return 5
    elif domain[0] == "l":
        return 5 # + round(int(domain[1])) 
    return 4 


def convert_UL_list(seq):
    """takes a domain level sequence in form of a list and converts it into UL annotation
    """    
    for x in range(len(seq)):
        if seq[x][-1] == "*":
            seq[x] = seq[x][:-1].upper()
    return seq
  
def extended_domain_path(domain_path):
    """ Extends domain level path to the corresponding nt path but with domains: abc -> aaa bbbbb ccc

    Args:
        domain_path (list): sublist correspond to  path sequence
    
    Returns: 
        full_path (list): extended path 

    """

    
    UL_domain_path = convert_UL_list(domain_path)
    full_path = []
        
    for z in range(len(UL_domain_path)):
        full_path.append(UL_domain_path[z] * d_length(UL_domain_path[z]))
            


    return "".join(full_path)

def split_ntseq_to_domainfp(nt_seq,domain_seq):
    """Takes a nucleotide sequence and splits it up in subsequences where each sequence i corresponds to the sequence at transcription step i 

        Args: 
            nt_seq (string): nucleotide sequence
            domain_seq (string: domain level sequence 

        
        Returns:   
            nt_path (list): list where each sublist corresponds to sequence at transcription step
    
    """

    split_seq = domain_seq.split()
        
    split_nt_sequence = []
    UL_seq = UL_list(split_seq)
        
    l_pointer = 0
    for z in split_seq:
        r_pointer = l_pointer + d_length(z)
            
        split_nt_sequence.append(nt_seq[l_pointer:r_pointer])
        l_pointer = r_pointer
        
    
        
    nt_path = []

    for  x in range(len(UL_seq)):
        if UL_seq[x][0] == "l":
                
            nt_path.append("".join(split_nt_sequence[:x+1]))
    nt_path.append(nt_seq)

    return nt_path


def is_star_pair(a,b):
    """Checks if two domain pair in Star Annotation"""
    if a[0] == b[0]:
        if len(a) < len(b) or len(a) > len(b):
            return True
        


def UL_list(list):
    """ Converts a list of domains into UL list
    """
    UL_liste = []
    for x in range(len(list)):
        if list[x][-1] == "*":
            
            UL_liste.append(list[x][:-1].upper())
        else:
            UL_liste.append(list[x])

    return(UL_liste)

def convert(string):
    """converts string in star annotation to UL_seq and keeps the spaces


    """
    counts = string.count("*")
    c_string = list(string)
    
    for x in range(1,len(c_string)):
        if string[x] == "*":
            c_string[x-1] = c_string[x-1].upper()
    for x in range(counts):		
        c_string.remove("*")

    c_string = "".join(c_string)
    return c_string


def differnt_l_domains(domain_seq):
    
    split_seq = domain_seq.split()

    count = 0
    
    for x,domain in enumerate(split_seq):
        
        if domain == "l":
            split_seq[x] = "l" + str(count)
            count += 1 
        
    return " ".join(split_seq)


"""
def afp_to_domainfp(afp,domain_seq):
    Takes an abstract folding path and a domain level sequence and converts it into a domain level path 
        Example: 
            afp: [["."],["()"],[".()"],["()()"]],domain_seq="b l b* a* l a b c d l d* c* b*"

            returns [[".."],["(.).."],["..((.))..."],["(.)...(((.)))"]]
    
    domain_fp = [[] for _ in afp]
    domain_seq = domain_seq.split()

    for x,cur_path in enumerate(afp):

        print(x)
        print(cur_path[0])

        module_index = -1
        path = ""
        
        for i,domain in enumerate(domain_seq):
            print("\n")
            #print("module index", module_index)
            print("domain",domain)
            if domain == "l":
                path += "."
                module_index += 1

                if module_index == x:
                    print("done with folding step")
                    domain_fp[x].append(path)
                    break 
                else:
                    print("test1:",module_index,x)
                
                    

            print(cur_path[0][module_index+1])
            if cur_path[0][module_index+1] == ".":
                path += "."


            elif cur_path[0][module_index + 1] == "(" and domain != "l": 
                print("(--------------------")
                
                j = i
                k = x -1
                while k <= len(afp[x]):
                    print(k)
                    print("her:",cur_path[0][k])
                    if cur_path[0][k] == ")":
                        print(domain_seq[i],domain_seq[j])
                        if is_star_pair(domain_seq[i],domain_seq[j]):
                            print("success")
                            path += "("
                            print(path)
                            break
                    j += 1
                    if domain_seq[j] == "l":

                        k += 1 
                


            elif cur_path[0][module_index + 1] == ")" and domain != "l":
                print("xxxxxxxxxx)")
                
                j = i - 1
                k = x 
                while  0 <= k:
                    print("k = ",k)
                    print("her:",cur_path[0][k])

                    if cur_path[0][k] == ")":
                        print(domain_seq[i],domain_seq[j])
                        if is_star_pair(domain_seq[i],domain_seq[j]):
                            print("success",path,"end path")
                            path += ")"
                            print("success",path,"end path")
                            break
                    
                    j -= 1
                    if domain_seq[j] == "l":

                        k -= 1 
                path += "."
        print("---------------------------------------------")
        print("Path:",path)
        domain_fp[x].append(path)

    return domain_fp
"""


def afp_to_domainfp(afp,domain_seq):
    """Takes an abstract folding path and a domain level sequence and converts it into a domain level path 
        Example: 
            afp: [["."],["()"],[".()"],["()()"]],domain_seq="b l b* a* l a b c d l d* c* b*"

            returns [[".."],["(.).."],["..((.))..."],["(.)...(((.)))"]]
    """
    domain_fp = [[] for _ in afp]
    domain_seq = domain_seq.split()

    for x,cur_path in enumerate(afp):
        module_index = 0
        path = ""

        for i,domain in enumerate(domain_seq):
            #print("module index", module_index)

            if domain == "l":    
                module_index += 1

                if module_index == x + 1:
                    path += '.'
                    break 

                else:
                    path += "."
                    

            if cur_path[0][module_index] == ".":
                path += "."
                
            elif cur_path[0][module_index] == "(" and domain != "l": 
                j = i
                k = module_index
                added_notation = False
                
                while k <= len(afp[x][0]) and j <= len(domain_seq) -1 :
                    
                    if cur_path[0][k] == ")":
                        
                        if is_star_pair(domain_seq[i],domain_seq[j]):
                            path += "("
                            added_notation = True
                            break
                    
                    if domain_seq[j] == "l":
                        k += 1 
                    j += 1    
                if added_notation != True:
                    path += "."
                    

            elif cur_path[0][module_index] == ")" and domain != "l":
                j = i 
                k = module_index
                added_notation = False 
                while  0 <= k and 0 <= j:
                    if cur_path[0][k] == "(":
                        
                        if is_star_pair(domain_seq[i],domain_seq[j]):
                            path += ")"     
                            added_notation = True
                            break

                    j -= 1

                    if domain_seq[j] == "l":
                        k -= 1 

                if added_notation != True:
                    path += "."


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

    #print("convert_functions")
    #print(path_to_pairtablepath(['.', '()', '.()', '(())']))


    #get_module_fp_sequences("AAABBBBBCCCLLLCCCBBBBBAAALLLBBBBBBBBBB")
    #print(extended_domain_path("vbulj*d*a*b*c*e*k*ltmifcbaghnslr*o*h*g*a*b*c*f*i*p*q*lzwqpifcbaghorxyly*x*r*o*h*g*a*b*c*f*i*p*q*w*z*lkecbadjlu*b*v*lblb*"))

    #print(UL_list("a aa* c av* a d e b b* bbb*".split()))
    #only_b_domainfp("b   l  A B C  l  c b a".split(),[['..'], ['(..)..'], ['..(((.)))']])
    #print(differnt_l_domains("l l l l l l l l l l")),, "
    #print(afp_to_domainfp([["."],["()"]],domain_seq="b l b* a* l"))
    #print(afp_to_domainfp([["."],["()"],[".()"],["()()"],[".()()"],["()()()"],[".()()()"],["()()()()"],[".()()()()"],["()()()()()"],[".()()()()()"],["()()()()()()"]],domain_seq="b   l   b* a*  l  a b c  l  c* b* a* d*  l  d a b c e  l  e* c* b* a* d* f*  l  f d a b c   e g  l  g* e* c* b* a* d* f* h*  l  h f d a b c   e   g i  l  i* g* e* c* b* a* d* f* h* j*  l  j h f d a b c   e   g   i  l   b* ".split())
    #print(extended_domain_path(("b   l   b* a*  l  a b c  l  c* b* a* d*  l  d a b c e  l  e* c* b* a* d* f*  l  f d a b c   e g  l  g* e* c* b* a* d* f* h*  l  h f d a b c   e   g i  l  i* g* e* c* b* a* d* f* h* j*  l  j h f d a b c   e   g   i  l   b* ").split()))
    #print(afp_to_domainfp([["."],["()"],["()."],["(())"]],"a b c l b* a* l b l c* b* a* "))


    print(find_connected_modules([[0, 0], [2, 2, 1], [3, 2, 1, 0], [4, 0, 0, 4, 3],[5,0,0,0,0,0],[6,0,0,0,0,6,5]]))
    print(extended_domain_path(convert_UL_list("b l0 b* a* l1 a b c l2 c* b* a* d* l3 d a b c e l4 e* c* b* a* d* f* l5 f d a b c e g l6 g* e* c* b* a* d* f* h* l7 h f d a b c e g i l8".split())))
                
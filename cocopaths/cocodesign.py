import infrared as ir 
from infrared import rna
import RNA 


def domain_path_to_nt_path(path,UL_seq,parameters):
    """ Domain_path_to_nt_path

    Takes the domain level path and extends the number of the base pairings corresponding to their length: ( with length 5 --> (((((

    args:
        path (list): domain level path 
        UL_seq (list): upper lower case sequence 
    
    Returns:
        ext_path (list): where each sublist corresponds to the extended domain path
    """
    UL_seq = UL_seq.split()
    ext_path = [[] for x in path]
    
    for x in range(len(path)):
        
        z_found = False 
        for z in range(len(path[x][0])):
            
            if UL_seq[z][0] == "z" and path[x][0][z] == "(":
                ext_path[x].append("(((")
                z_found = True
                
            elif UL_seq[z][0] == "b" and z_found and path[x][0][z] == ")":
                ext_path[x].append(")))...")
                z_found = False

            else:    
                ext_path[x].append(path[x][0][z] * parameters["d_length"][UL_seq[z]])

        ext_path[x] = "".join(ext_path[x])

    return ext_path

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

def identical_domains_constraint(domain,UL_seq,model):
    """Identical_domains_constraint
    Applies the constraint to our model, that the value of indices of the same domain must be equal.

    Args:
        domain: domain name
        UL_seq: upper lower case domain sequence
        model: model created in the infrared sampler
    """ 
    UL_seq = UL_seq.split()   
    for x in range(len(UL_seq)):
        if UL_seq[x] == domain:
            i_pointer = 0
            

            for u in UL_seq[:x]:
                
                i_pointer += cv.d_length(u)

            j_pointer = i_pointer    

            for q in range(len(UL_seq[:x]),len(UL_seq)):
           
                if domain not in UL_seq[x+1:]:
                    return    
                if UL_seq[q] == domain and j_pointer != i_pointer:
                   
                    break
                j_pointer += cv.d_length(UL_seq[q])

        
            if i_pointer > j_pointer:
                return
            for z in range(cv.d_length(domain)):
                if i_pointer != j_pointer:
                    
                    model.add_constraints(IdenticalDomains(i_pointer,j_pointer))
                    
                i_pointer += 1

def rna_design(seq,path,parameters):
    print("Given sequence: ",seq)
    print("Given path: ",path)

    split_seq = seq.split()


    d_seq_len = sum([parameters["d_seq"][x] for x in split_seq])
    
    
    # define constraints
    #Identical Domains
    ir.def_constraint_class(
            "IdenticalDomains",
            lambda i,j: [i,j],
            lambda x,y: x == y,
            module  = __name__ 
        )

    model = ir.Model(d_seq_len,4)

    UL_seq = convert_to_UL(seq)

    ext_path = domain_path_to_nt_path(path,UL_seq,parameters)

    # Folding path constraint
    for x in ext_path: 
        print(x)
        cons = []
        bps = rna.parse(x)
        cons = [rna.BPComp(i,j) for (i,j) in bps]
        

        model.add_constraints(cons)


    #Identical domain constrain
    unique_domains = "".join(set(UL_seq))
    for domain in  unique_domains:

        if domain != "l":
            identical_domains_constraint(domain,UL_seq,model)
    


    extended_fp = domain_path_to_nt_path(path,UL_seq)



def main():
    pass 


if __name__ == "__main__":
    main()
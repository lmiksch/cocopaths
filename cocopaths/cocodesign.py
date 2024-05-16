import infrared as ir 
from infrared import rna
import RNA, random, math, argparse, logging, os,sys
from .utils import(is_balanced_structure,afp_terminal_input,afp_to_domainfp)
from cocopaths import __version__
from cocopaths.cocosim import verify_domain_foldingpath
from peppercornenumerator.input import read_pil


#______define_logger_____#
cocodesign_logger = logging.getLogger('cocodesign')
console_handler = logging.StreamHandler()
formatter = logging.Formatter('# %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
cocodesign_logger.addHandler(console_handler)

def valid_cutoff(value):
    value = float(value)
    if 0 <= value < 1:
        return value
    else:
        raise argparse.ArgumentTypeError(f"Cutoff must be between 0 and 1, but got {value}")

def objective_function(fe,efe,mse,barrier):
    """
    Function exists to only change one thing when changing the obj function
    """
    
    obj_fun = "abs(fe - efe) + barrier*0.1"
    score = eval(obj_fun.format(fe=fe,efe=efe,barrier=barrier))
    return score,obj_fun

def couple(pair):
	if pair[0][0].upper() == pair[1][0].upper() and pair[0] != pair[1]:
		return True
     


def domain_path_to_nt_path(path,d_seq,parameters):
    """ Domain_path_to_nt_path

    Takes the domain level path and extends the number of the base pairings corresponding to their length: ( with length 5 --> (((((

    args:
        path (list): domain level path 
        UL_seq (string): upper lower case sequence 
    
    Returns:
        ext_path (list): where each sublist corresponds to the extended domain path
    """

    split_seq = d_seq.split()
    ext_path = []
    

    for step_index,step in enumerate(path): 
        ex_step_path = ""
        for nt,domain in zip(list(step),split_seq):
            extendo = nt * parameters["d_length"][domain]
            ex_step_path += extendo
        ext_path.append(ex_step_path)


    return ext_path



def is_star_pair(a,b):
    """Checks if two domain pair in Star Annotation"""
    if a[0] == b[0]:
        if len(a) < len(b) or len(a) > len(b):
            return True
        





def convert_to_UL(string):
    """Converts sting annotated domain sequence into a upper lower case annotated string

    Args: 
        string (string): space annotated domain level sequence

    Returns: 
        UL_seq (str): Upper lower case annotated domain level sequence
    
    """

    split_seq = string.split()
    counts = string.count("*")
    c_string = split_seq
    
    for x in range(len(split_seq)):
        
        if split_seq[x][0] == "L" and split_seq[x][-1] != "*":
            c_string[x] = split_seq[x].lower()

        elif split_seq[x][0] == "L" and split_seq[x][-1] == "*":
            c_string[x] = c_string[x][:-1]

        elif split_seq[x][-1] == "*":
            c_string[x] = c_string[x].upper()
            c_string[x] = c_string[x][:-1]
    
    c_string = " ".join(c_string)
    
    return c_string

def identical_domains_constraint(domain,split_seq,model,parameters):
    """Identical_domains_constraint
    Applies the constraint to our model, that the value of indices of the same domain must be equal.

    Args:
        domain: domain name
        UL_seq: upper lower case domain sequence
        model: model created in the infrared sampler
    """ 
    for x in range(len(split_seq)):
        if split_seq[x] == domain:
            i_pointer = 0
            

            for u in split_seq[:x]:
                
                i_pointer += parameters["d_length"][u]

            j_pointer = i_pointer    

            for q in range(len(split_seq[:x]),len(split_seq)):
           
                if domain not in split_seq[x+1:]:
                    return    
                if split_seq[q] == domain and j_pointer != i_pointer:
                   
                    break
                j_pointer += parameters["d_length"][split_seq[q]]

        
            if i_pointer > j_pointer:
                return
            for z in range(parameters["d_length"][domain]):
                if i_pointer != j_pointer:
                    
                    model.add_constraints(IdenticalDomains(i_pointer,j_pointer))
                    
                i_pointer += 1
                j_pointer += 1


def extend_domain_seq(d_seq,parameters):

    UL_domain_seq = convert_to_UL(d_seq)


    split_domain_seq = d_seq.split()
    UL_split_seq = UL_domain_seq.split()
    extended_domain_seq = ""

    for UL_domain,domain in zip(UL_split_seq,split_domain_seq):
        if UL_domain[0] == "L":
            extended_domain_seq += "L" * int(parameters["d_length"][domain])
        elif UL_domain[0] == "l":
            extended_domain_seq += "l" * int(parameters["d_length"][domain])
        elif UL_domain[0] == "S":
            extended_domain_seq += "S" * int(parameters["d_length"][domain])
        else:
            extended_domain_seq += UL_domain * int(parameters["d_length"][domain])

    return extended_domain_seq


def call_findpath(seq, ss1, ss2, md, fpw, mxb = float('inf')):
    """ Call ViennaRNA findpath. Modified from DrTransformerd_length
    """
    fc = RNA.fold_compound(seq)
    
    if mxb == float('inf'):
        path = fc.path_findpath(ss1, ss2, width = fpw)
        
    else:
        e1 = round(fc.eval_structure(ss1), 4)
        dcal_bound = int(round((mxb + e1)))
        path = fc.path_findpath(ss1, ss2, maxE = dcal_bound, width = fpw)
        
    del fc

    if len(path):
        mypath = []
        barrier = None
        for step in path:
            struct = step.s
            energy = float(round(step.en,4))
            mypath.append((struct, energy))
            if barrier is None or barrier < energy:
                barrier = energy
        barrier -= float(round(path[0].en,4))
        del step, path # potentially a good idea.
        return mypath, barrier
    return None, None

def constrained_efe(sequence,c):
        """Calculates the ensemble free energy of sequence
        
        """
        fc = RNA.fold_compound(sequence)
        # fc.hc_add_from_db(c) # not sure if it's necessary to add constraint
        return fc.pf()[1]

def mc_optimize(model, objective, steps, temp, start=None):
        sampler = ir.Sampler(model)
        cur = sampler.sample() if start is None else start
        curval = objective(cur)
        best, bestval = cur, curval
        
        ccs = model.connected_components()
        weights = [1/len(cc) for cc in ccs]
        
        for i in range(steps):
            cc = random.choices(ccs,weights)[0]
            new = sampler.resample(cc, cur)
            newval = objective(new)
            print("\rCurrent Score: ",newval, end="")
            if (newval >= curval
                or random.random() <= math.exp((newval-curval)/temp)):
                cur, curval = new, newval
                if curval > bestval:
                    best, bestval = cur, curval
        
        return (best, bestval)

def rna_design(seq,path,parameters):


    split_seq = seq.split()

    d_seq_len = sum([int(parameters["d_length"][x]) for x in split_seq])
    
    
    # define constraints
    #Identical Domains
    ir.def_constraint_class(
            "IdenticalDomains",
            lambda i,j: [i,j],
            lambda x,y: x == y,
            module  = __name__ 
        )

    model = ir.Model(d_seq_len,4)

    

    ext_path = path
    
    
    # Folding path constraint
    for x in ext_path: 
        cons = []
        bps = rna.parse(x)
        cons = [rna.BPComp(i,j) for (i,j) in bps]
        

        model.add_constraints(cons)


    

    def rstd_objective(sequence,score_list = False):
        
        split_nt_sequence = []
        #creates a list of list where each sublist i corresponds to the sequence at transcription step i 
        l_pointer = 0
        for z in split_seq:
            
            r_pointer = l_pointer + parameters["d_length"][z]
            
            split_nt_sequence.append(sequence[l_pointer:r_pointer])
            l_pointer = r_pointer
        
        
        nt_path = []
        for x in range(len(split_seq)):
            if split_seq[x][0][0] == "S":
            
                
                nt_path.append("".join(split_nt_sequence[:x+1]))
        nt_path.append(sequence)
        
       
        total = []
        # add score for first folding step
        fc = RNA.fold_compound(nt_path[0])
        fe = fc.eval_structure(ext_path[0].replace(".","x"))
        efe = constrained_efe(nt_path[0],ext_path[0])
        mse = 0
        barrier = 0
        obj_score = objective_function(fe,efe,mse,barrier)[0]
        total.append(obj_score)

        for x in range(1,len(ext_path)):
            #prepare input for finpath 
    

            ss1 = ext_path[x-1] + ("." * (len(nt_path[x])-len(nt_path[x-1])))
         
            efe = constrained_efe(nt_path[x],ext_path[x])
            fc = RNA.fold_compound(nt_path[x])
            
            fe = fc.eval_structure(ext_path[x].replace(".","x"))
            mypath, barrier = call_findpath(nt_path[x],ss1,ext_path[x],0,30)

           
            if mypath != None:
                deltaE = abs(mypath[-1][1]) - abs(mypath[0][1])
                
            else: 
                deltaE = 99
                barrier = 99


            mse = 0
            global factor


            #print("fe", fe)


            
            #print("efe",efe)
            #print("barrierIdenticalDomains(i_pointer,j_pointer)",barrier)
            obj_score = objective_function(fe,efe,mse,barrier)[0]
            
            total.append(obj_score) 

        #calculate MSE of scores to keep all scores equal 


        total_mean = sum(total)/ len(total)

        squared_error = [(x - total_mean) ** 2 for x in total]
        #print("total after", total)   
        for i,score in enumerate(total):
            total[i] = score + squared_error[i]

        if score_list:
            return total
        #print("total after", total)    
        return sum(total)

    objective = lambda x: -rstd_objective(rna.ass_to_seq(x))

    best, best_val = mc_optimize(model, objective,steps = parameters['steps'], temp = 0.04)


    print("Calculated NT sequence:")
    print(rna.ass_to_seq(best), -best_val)

    return rna.ass_to_seq(best), -best_val


def set_verbosity(cocodesign_logger, console_handler, verbosity):
    if verbosity == 0:
        console_handler.setLevel(logging.CRITICAL)
        cocodesign_logger.setLevel(logging.CRITICAL)
    elif verbosity == 1:
        console_handler.setLevel(logging.WARNING)
        cocodesign_logger.setLevel(logging.WARNING)
    elif verbosity == 2:
        console_handler.setLevel(logging.INFO)
        cocodesign_logger.setLevel(logging.INFO)
    elif verbosity >= 3:
        console_handler.setLevel(logging.DEBUG)
        cocodesign_logger.setLevel(logging.DEBUG)
    else:
        console_handler.setLevel(logging.CRITICAL)
        cocodesign_logger.setLevel(logging.CRITICAL)

def main():
    

    #_________________Argparse_________________#
    parser = argparse.ArgumentParser(
        formatter_class = argparse.RawDescriptionHelpFormatter,
        description=f"""Cocodesign version {__version__}: generates a nucleotide level sequence based on a domain level sequence 
)""")

    parser.add_argument("-i","--input" , nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="Reads file in the form of a pil  format as input if not specified user can input via console.")
    parser.add_argument("-v", "--verbose",action="count", default = 0,
        help = "Track process by writing verbose output to STDOUT during calculations.")
    parser.add_argument("-s","--steps",help="Number of steps in the optimization",default=2000,type=int)
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.00001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-l", "--logic", action="store_true", default=False,help="Visualizes Logic domain pairings. Best used when analyzing cocopaths generated sequences. (default = False)")
    parser.add_argument("-f", "--force", action="store_true", default=False,help="Forces the design disregarding the outcome of the simulation (default = False)")
    parser.add_argument("-cutoff", "--cutoff", action="store", type=valid_cutoff, default=float('-inf'),help="Cutoff value at which structures won't get accepted (default: -inf, valid range: 0 to 1)")
    




    args = parser.parse_args()

    #_________________Set_logger_verbosity________________#

    set_verbosity(cocodesign_logger, console_handler, args.verbose)
    d_length = {}

    if args.input.isatty():
            
            print("\n")
            print("Please input a domain level sequence:")
            d_seq = input()
    else:

        file_extension = os.path.splitext(args.input.name)[1]  # Get the file extension

        if file_extension == ".pil":
            afp = None
            pil_input = read_pil(args.input,True)
            print(pil_input)
            if len(pil_input[0]) == 1:
                input_complex = next(iter(pil_input[0].values()))
                print("Value of the entry:", input_complex)
                print(f"{input_complex._sequence = }  ")
                print(f"{input_complex.kernel_string = }  ")
                

                d_seq = input_complex.kernel_string

                if all(char.isalpha() or char == '*' for char in d_seq):
                    raise SystemExit("SystemExit: Only a domain level sequence is accepted. (No structural information)")

                for domain in input_complex._sequence:
                    length = domain._length
                    name = domain._name
                    if name not in d_length:
                        d_length[name] = length

                print(f"{d_length = }")

                print(f'\nInput domain level sequence: {d_seq}\n')
                        
            else:
                raise SystemExit("SystemExit:More than one kernel sequence in input. We can only simulate one kernel string at a time. ")              
    
        else: 
            #Handle input from cocosim here
            raise SystemExit("Only data in the .pil format is currently accepted.")

    print("Please input the afp by which the domain level sequence was designed:")
    afp = afp_terminal_input()  


    print(f"{d_seq = }")

    if len(d_length) == 0:    
        for domain in d_seq.split():
            if domain[0] == "L":
                d_length[domain] = 8
            elif domain[0] == 'S':
                d_length[domain] =  round(int(domain[1]) * 4)  
            else: 
                d_length[domain] = 3 

   
    parameters = {"k_slow": args.k_slow,'k_fast': args.k_fast, "cutoff": args.cutoff,"d_length":d_length,"d_seq":d_seq,"logic":args.logic,'steps':args.steps}


    # add check to see if number of folding steps = number of spacer domains

    #___simulate_domain-level-foldingpath______# 

    if not args.force:
        dominant_fp = verify_domain_foldingpath(afp,d_seq,parameters)

        if not dominant_fp:
            print(f"Simulated folding path does not match afp")
            exit()


    #How to get from AFP to domain fp 

    domain_fp = afp_to_domainfp(afp,d_seq)

    ext_folding_path = domain_path_to_nt_path(domain_fp,d_seq,parameters)


    print(f"Input AFP: {afp}\n\n")
    print(f"Input Domain level sequence: {d_seq}\n\n")    
    nt_seq, score = rna_design(d_seq,ext_folding_path,parameters)



    print(f"\n\nRNA design done\nSequence = {nt_seq} \n {score = }\n{extend_domain_seq(d_seq,parameters)}")


    #print output and extended nucleotide sequence for easier analysis 



if __name__ == "__main__":
    
    main()
import infrared as ir 
from infrared import rna
import RNA, random, math, argparse, logging, os,sys
from .utils import(is_balanced_structure,afp_terminal_input,afp_to_domainfp,path_to_pairtablepath,split_into_modules)
from cocopaths import __version__
from cocopaths.cocosim import verify_domain_foldingpath
from peppercornenumerator.input import read_pil
import numpy as np

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

def objective_function(cefe,fe,efe,c_barrier,barrier,ensemble_defect):
    """
    Function exists to only change one thing when changing the obj function
    """
    
    #obj_fun =  "abs(cefe - fe) + 0.1*barrier + ensemble_defect" 
    obj_fun = "abs(efe - cefe) + 0.1*c_barrier"
    score = eval(obj_fun.format(cefe=cefe,fe=fe,barrier=barrier,ensemble_defect=ensemble_defect))
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
        
def find_toeholdstruct(domain_seq,domain_fp,step_index,pair_index):
    '''Finds the indices of the toehold domains for the step. 
    
        Args:
            domain_seq(str): 
            domain_fp(list):
        Returns:
            toehold_struct(tuple): structure where toeholds are formed
    '''

    module_list = split_into_modules(domain_seq)
    
    #print(f'\n Finding toehold structures')
    #print(f'{domain_seq = }\n{domain_fp = }\n{step_index = }\n{pair_index}')
    #print(f'to go  step = {domain_fp[step_index]}\npairing step {domain_fp[pair_index]}\n{module_list = }')

    #attacker is the first if no bp between them 
    #attacker second unit if bp between them 



    attack_module = module_list[pair_index]
    pair_module = module_list[step_index]

    #print(f'\n{attack_module = }')
    #print(f'\n{pair_module = }')

    outside_toehold_indices = []
    inside_toehold_indices = []

    inside_domain = str(-1)
    outside_domain = str(-1)


    after_Logic = False
    for d_index_1,domain_1 in enumerate(attack_module.split()):
        for d_index_2,domain_2 in enumerate(pair_module.split()): 
            if is_star_pair(domain_1,domain_2) and domain_1[0] != 'L': 
                #print(f'star pair {domain_1}, {domain_2}')


                if after_Logic:
                    outside_domain = domain_1
                        
                else: 
                    inside_domain = domain_1

        if domain_1[0] == 'L':
            after_Logic = True

    #print(f'\nOutside {outside_domain}')

     
    #print(f'Inside {inside_domain}')

    s_index = 0
    between_struct = False
    between = False
    #print(domain_fp[step_index])
    #print(domain_seq.split())

    for char,domain in zip(domain_fp[step_index-1],domain_seq.split()):
        #print(char,between)
        if domain[0] == 'S':
            s_index = domain[1]
        
        if domain[0] == inside_domain[0]:
            if between:
                between = False  
            else:
                between = True


        if between and char != '.':
            between_struct = True


    if between_struct:
        toehold = outside_domain
        for i,ind in zip(attack_module.split(),range(len(attack_module))):
            #print(i,'index = ',ind)

            if i[0] == toehold[0]:
                toehold_struct = domain_fp[step_index-1] + '.' * ind

    else:
        toehold = inside_domain
        for i,ind in zip(attack_module.split(),range(len(pair_module))):
            #print(i,'i = ')
            if i[0] == toehold[0]:
                toehold_struct = domain_fp[step_index-1] + '.' * ind


    


    s_index = 0

    for char,domain,index in zip(domain_fp[step_index],domain_seq.split(),range(len(domain_fp[step_index]))):
        if domain[0] == 'S':
            s_index += 1
        
        if domain[0] == toehold[0] and (s_index == step_index or s_index == pair_index):
            toehold_struct = toehold_struct[:index] + domain_fp[step_index][index] + toehold_struct[index + 1:]

    return toehold_struct
        

    exit()

def toehold_structures(domain_fp,domain_seq,afp):
    """Calculates the structures after the toeholds have formed in the folding path.

    Args:
        domain_fp(list): each entry corresponds to a step in folding path
        domain_seq(str): domain_level sequence 
    Returns:
    """

    split_d_seq = domain_seq.split()

    pt_path = path_to_pairtablepath(domain_fp)
    afp_pt = path_to_pairtablepath(afp)
   

    module_seq = split_into_modules(domain_seq)


    toehold_path = []
    for i,step in enumerate(domain_fp[:-1]):
        toehold_path.append(step)
        if step != domain_fp[i+1][:len(step)]:#checks for refolding event if structure has changed
            #print(f'{step = }, {domain_fp[i+1][:len(step)] = }')

            p_module_index = afp_pt[i+1][-1] -1

            #print(f'{p_module_index = }')

            #print(f'{module_seq = }')

            if len(module_seq[i+1].split()) > 1 and len(module_seq[p_module_index].split()) > 1: #skips steps where no domains are employed
                toehold_struct = find_toeholdstruct(domain_seq,domain_fp,i+1,p_module_index)

                toehold_path.append(toehold_struct)

    toehold_path.append(domain_fp[-1])


    return toehold_path

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
        split_seq: upper lower case domain sequence
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
            #print("\rCurrent Score: ",newval, end="")
            if (newval >= curval
                or random.random() <= math.exp((newval-curval)/temp)):
                cur, curval = new, newval
                if curval > bestval:
                    best, bestval = cur, curval
        
        return (best, bestval)


def score_sequence(seq,d_seq,parameters,afp,domain_fp):
    '''Scores the sequence with the current objective function from cocodesign
    
    Args: 
        seq(str): nucleotide sequence 
        parameters(dict): dict of parameters important is that d_length dict exists
        afp(list): each entry corresponds to a step in the abstract folding path

    Returns: 
        output(str): Formatted output for writing in file or displaying in terminal
        score(float): max score of the sequence

    '''
    split_seq = d_seq.split()

    split_nt_sequence = []
        #creates a list of list where each sublist i corresponds to the sequence at transcription step i 
    l_pointer = 0
    for z in split_seq:
        
        r_pointer = l_pointer + parameters["d_length"][z]
        
        split_nt_sequence.append(seq[l_pointer:r_pointer])
        l_pointer = r_pointer
    
    ext_path = domain_path_to_nt_path(domain_fp,d_seq,parameters)
    nt_path = []
    
    for x in ext_path:    
        nt_path.append("".join(seq[:len(x)]))
    nt_path.append(seq)




    total = []
    output_matrix = []

    for x in range(0,len(ext_path)):
        #prepare input for finpath 
        mse = 0
        if x == 0: 
            ss1 = ext_path[x] 
        else:
            ss1 = ext_path[x-1] + ("." * (len(nt_path[x])-len(nt_path[x-1])))
         
        fc = RNA.fold_compound(nt_path[x])
        
        efe = fc.pf()[1]

        fe = fc.eval_structure(ext_path[x])

        fc.hc_add_from_db(ext_path[x])
        
        tfe = fc.eval_structure(ext_path[x])
        
        cefe = fc.pf()[1]
        

        

        
        mypath, barrier = call_findpath(nt_path[x],ss1,ext_path[x],0,30)
        if mypath != None:
            deltaE = abs(mypath[-1][1]) - abs(mypath[0][1])
        else: 
            deltaE = 99
            barrier = 99
        global factor
        


        #constraint barrier
        
        #print(f"{ss1 = }")
        #print(f"{ext_path[x] = }")


        fc_1 = RNA.fold_compound(nt_path[x])
        fc_1.hc_add_from_db(ss1)
        c_mfe_1 = fc_1.mfe()[0]
        #print(f"{ext_path[x] = }")

        fc_2 = RNA.fold_compound(nt_path[x])
        fc_2.hc_add_from_db(ext_path[x].replace('.','x'))
        c_mfe_2 = fc_2.mfe()[0]

        #print(f"{c_mfe_1 = }")
        #print(f"{c_mfe_2 = }")


        mypath, c_barrier = call_findpath(nt_path[x],c_mfe_1,c_mfe_2,0,30)
        if mypath != None:
            deltaE = abs(mypath[-1][1]) - abs(mypath[0][1])
        else: 
            deltaE = 99
            c_barrier = 99
        global factor

        
        #print(f'{efe = }, {cefe = }, {efe = }, {c_barrier = }, {barrier = }')

        ensemble_defect = fc.ensemble_defect(ext_path[x])
        obj_score,obj_function = objective_function(cefe,fe,efe,c_barrier,barrier,ensemble_defect)
        total.append(obj_score) 
        prob = fc.pr_structure(ext_path[x].replace(".","x"))
        output_matrix.append([obj_score,prob,abs(efe - cefe),c_barrier,len(mypath),ensemble_defect])

        #output += f"{obj_score:<12.6g}\t{prob:<8.6g}\t{(efe - fe) ** 2:<11.6g}\t{barrier:<8.6g}\t{len(mypath):<9}\t{ensemble_defect:<9}\n"

    mean_efe_fe_squared = sum(row[2] for row in output_matrix) / len(output_matrix)
    mean_ensemble_defect = sum(row[5] for row in output_matrix) / len(output_matrix)

    #Option to include MSE 

    for i, row in enumerate(output_matrix):
        efe_fe_squared_error = abs(row[2] - mean_efe_fe_squared) 
        ensemble_defect_error = abs(row[5] - mean_ensemble_defect) 
        
        squared_error = efe_fe_squared_error + ensemble_defect_error
        #total[i] += squared_error
        
        output_matrix[i].append(efe_fe_squared_error)
        output_matrix[i].append(ensemble_defect_error)
        output_matrix[i].append(squared_error)
        output_matrix[i].append(total[i])
        mean = sum(total)/len(total)

     
    #format output 
    output = f"{obj_function = }\nObj Score\tProb    \tEfe-Fe^2   \tBarrier \tPath Length\tEnsemble Defect\tEfe-Fe^2 Error\tEnsemble Defect Error\tSquared Error\tTotal Score\n"
    output += "\n".join(
        f"{row[0]:<12.6g}\t{row[1]:<8.6g}\t{row[2]:<11.6g}\t{row[3]:<8.6g}\t{row[4]:<9}\t{row[5]:<9.5}\t{row[6]:<15.6g}\t{row[7]:<21.6g}\t{row[8]:<13.6g}\t{row[9]:<11.6g}" 
        for row in output_matrix)

    #Return maximum score   

    return max(total),output

def nt_path_to_afp(nt_path,domain_seq,parameters):
    afp = []
    for i,step in enumerate(nt_path):
        cur_struct = []
        seen_l = 0            
        index = 0
        for domain in domain_seq.split():
            if domain.endswith(str(i)):
                break
            elif not domain.startswith('L') and not domain.startswith('S'):
                index += parameters['d_length'][domain[0]]
            elif domain.startswith('S'):
                index += parameters['d_length'][domain]
            elif domain.startswith('L'):
                cur_struct.append(step[index + 1])
                index += parameters['d_length'][domain]

        afp.append(''.join(cur_struct))

    return afp 


def rna_design(seq,path,parameters,afp,domain_fp):

    print(f"{seq =}")
    print("hello there")
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


    for domain in set(parameters["d_length"]):
        if domain != "l":
            identical_domains_constraint(domain,split_seq,model,parameters)
    

    ext_path = path
    
    # Folding path constraint
    for x in ext_path: 
        cons = []
        bps = rna.parse(x)
        cons = [rna.BPComp(i,j) for (i,j) in bps]
        model.add_constraints(cons)

    objective = lambda x: -score_sequence(rna.ass_to_seq(x),seq,parameters,afp,domain_fp)[0]


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
    parser.add_argument("-aCFP", "--aCFP", action="store", type=str, default=None,help="aCFP where each step is seperated by a comma. If not specified the user needs to use the Terminal as input.")





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
    if not args.aCFP:
        afp = afp_terminal_input()
    else:
        afp = args.aCFP.split(',')


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

    #domain_fp = toehold_structures(domain_fp,d_seq,afp)

    for step in domain_fp:
        print(step)

    ext_folding_path = domain_path_to_nt_path(domain_fp,d_seq,parameters)

    for path in ext_folding_path:
        print(path)

    print(f"Input AFP: {afp}\n")
    print(f"Input Domain level sequence: {d_seq}\n")    
    nt_seq, score = rna_design(d_seq,ext_folding_path,parameters,afp,domain_fp)

    print(f"\n\nRNA design done\nSequence = {nt_seq} \n {score = }\n{extend_domain_seq(d_seq,parameters)}")


    #print output and extended nucleotide sequence for easier analysis 
    return nt_seq


if __name__ == "__main__":
    
    main()

    
    afp = ['.','()','.()','()()']
    domain_seq = ' L0*  S0 a L0 b S1 c* b* L0* a* d* S2 d a L0 b c S3'

    domain_fp =  afp_to_domainfp(afp,domain_seq)

    #print((domain_fp,domain_seq,afp))

    '''

    afp = ['.','()','().','(())']
    domain_seq = ' a* b* L0* c* d* S0 c L0 b S1  L0*  S2 d c L0 b a S3'

    domain_fp =  afp_to_domainfp(afp,domain_seq)

    toehold_structures(domain_fp,domain_seq,afp)
    '''
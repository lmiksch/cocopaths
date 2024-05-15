"""Compiler for cotranscriptional folding simulator

Simulates cotranscriptional folding of a domain level sequence.

uses Function from Peppercornenumerator to find reaction types/ possible products of reactions
"""
import sys,argparse,logging,copy,os,random,inspect
from peppercornenumerator import peppercorn 
from dsdobjects.objectio import read_pil as dsd_read_pil
from peppercornenumerator.input import read_pil
from peppercornenumerator import Enumerator 
from peppercornenumerator.enumerator import BI_REACTIONS
from peppercornenumerator.reactions import bind21
from peppercornenumerator.objects import PepperComplex
from peppercornenumerator.condense import is_outgoing , SetOfFates

from crnsimulator import ReactionGraph, get_integrator
from crnsimulator.odelib_template import add_integrator_args

import numpy as np
from io import StringIO
from natsort import natsorted
from .utils import cv_db2kernel, kernel_to_dot_bracket, only_logic_domain_struct, afp_terminal_input, afp_to_domainfp, domainfp_to_afp
import numpy as np 


logger = logging.getLogger('cocosim')
console_handler = logging.StreamHandler()
formatter = logging.Formatter('# %(levelname)s \n - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)

file_handler = logging.FileHandler('cocosim.log')
file_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
file_handler.setFormatter(file_formatter)
logger.addHandler(file_handler)

def run_sim(d_seq, parameters):
    """
    Runs a cotranscriptional simulation of the given d_seq using peppercorn from the peppercornenumerator

    Args:
        d_seq(str): domain level sequence 
        parameters(dict)
    """    

    # Initiate datastructures
    used_structure_names = 1
    all_complexes = {}
    folding_step_complexes = []
    d_seq_split = d_seq.split()
    
    print(f'\n\n\nCocosim\nDomain Level seq: {d_seq}\n\n')
    


    # before input parsing define name of complexes 
    
    complexes, reactions = input_parsing(d_seq_split[0:2], {"E0" : [cv_db2kernel(d_seq_split[0:2], ".."),1]},parameters)
    

    resulting_complexes,transient_complexes,enum = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters, all_complexes=all_complexes)
    
    
    
    map_transient_states(resulting_complexes,transient_complexes,all_complexes,enum,parameters)
    
    
    
    # all complex keeps track of complexes throughtout the sim
    

    
    #only in first step this should only result in 1 structure
    for complex in resulting_complexes:
        all_complexes["Id_" + str(str(len(all_complexes) + 1))] = [complex,complex.occupancy]
        complex.occupancy = 1
        
    

    folding_step_complexes.append(resulting_complexes)	
   



    #Print continous output
    if parameters["logic"]:
        for x, complex in complexes.items():
            print(f'Transcription Step | Occupancy  |  Logic domain pairing | Structure	 \n')
            kernel_string = kernel_to_dot_bracket(complex.kernel_string)
            db_struct = (only_logic_domain_struct(d_seq.split(),kernel_string))
            #occupancy = np.float128(complex.occupancy)
            if round(complex.occupancy,4) > 1e-4:
                print(f"{2:3}\t|\t{complex.occupancy:^8.4f}\t|\t{db_struct:16}|\t{complex.kernel_string}")
        print()    

    else:
        for x, complex in complexes.items():
            print(f'Transcription Step | Occupancy  | Structure	 \n')
            kernel_string = kernel_to_dot_bracket(complex.kernel_string)
            db_struct = (only_logic_domain_struct(d_seq.split(),kernel_string))
            #occupancy = np.float128(complex.occupancy)
            if round(complex.occupancy,4) > 1e-4:
                print(f"{2:3}\t|\t{complex.occupancy:^8.4f}\t|\t{complex.kernel_string}")
        print()  
        

    for step in range(2, len(d_seq_split)):
        logger.info("\n\n\n______________________________________________")
        logger.info(f"Step: {step} Current Complexes: {parameters}")

        next_complexes = []
        
        
        #put here add one domain to struct and give new name
    
        new_complexes = {}
        
        occ_sum = 0 
        for complex in folding_step_complexes[-1]:
            occ_sum += complex.occupancy
        assert abs(occ_sum - 1 ) < 1e-10,f"Difference at the beginning = {abs(occ_sum - 1 )} -> previous step was fucked"

        #______Function_for_extending_and_updating_names__ 

        old_names = []
        for complex in folding_step_complexes[-1]:
            
            c_name = complex._name
            next_struct = complex.kernel_string + " " + d_seq_split[step]

            new_complexes["E" + str(used_structure_names)] = [next_struct,complex.occupancy]
            
            used_structure_names += 1
            if complex.occupancy > 0: 
                old_names.append(c_name)

        
        
        complexes, reactions = input_parsing(d_seq_split[0:step + 1], new_complexes,parameters)

        #Update all_complexes with new complexes 

        
        new_complexes = complexes.copy() 
        new_complex_items = []
      
        for x, c_name in enumerate(old_names):
            for cid, complex_obj in all_complexes.items():
                
                if c_name == complex_obj[0]._name:

                    first_key = next(iter(complexes))
                    next_complex = complexes.pop(first_key)
                    next_complex.occupancy = np.float128(complex_obj[0].occupancy)
                    new_complex_items.append((cid, next_complex))
                
        # Add the new items to the dictionary
        for cid, complex_item in new_complex_items:
            all_complexes[cid] = [complex_item,complex_item.occupancy]	


        #Add remaining ones to all_complexes should incluce concentration aswell
        if len(complexes) > 0:
            for key,complex in complexes.items():
                all_complexes["Id_" + str(str(len(all_complexes) + 1))] = [complex,complex.occupancy] 
                

        

        for key,complex in all_complexes.items():
            if complex[0].occupancy == None or complex[1] == None:
                raise SystemExit("Unoccupied Complex in all complex after updating")



        complexes = new_complexes

        
        
        total_occupancy = 0

        #Print continous output
        if parameters["logic"]:
            for x, complex in complexes.items():
                try:
                    occupancy = np.float128(complex.occupancy)
                    total_occupancy += occupancy
                    kernel_string = kernel_to_dot_bracket(complex.kernel_string)
                    db_struct = (only_logic_domain_struct(d_seq.split(),kernel_string))
                    if occupancy != 0:
                        print(f"{step + 1:3}\t|\t{occupancy:^8.4f}\t|\t{db_struct:16}|\t{complex.kernel_string}\t")
                except: 
                    pass
                    #print(f"{step + 1:3}   |	{0:8.4f}      | {complex.kernel_string} \n")
            print()
        else:
            for x, complex in complexes.items():
                try:
                    occupancy = np.float128(complex.occupancy)
                    total_occupancy += occupancy
                    kernel_string = kernel_to_dot_bracket(complex.kernel_string)
                    if occupancy != 0:
                        print(f"{step + 1:3}\t|\t{occupancy:^8.4f}\t|\t{complex.kernel_string}\t")
                except: 
                    pass
        


        logger.debug(f"Step:{step} Total Occupancy:{total_occupancy:20.20f}\n")
        assert abs(total_occupancy - 1) < 1e-8,f"Occupancy is not equal 1 difference is {abs(total_occupancy - 1)}"         
        resting_complexes,transient_complexes,enum = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters,all_complexes=all_complexes)
        
        
        
        #checks if there is a transient state which needs to be mapped
        
        map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)

        resting_complexes = list(set(resting_complexes))
        oc_sum = 0
        for complex in resting_complexes:
            try: 
                oc_sum += complex.occupancy 
            except: 
                pass
        oc_sum = round(oc_sum,8) #can be later adjusted       
        if abs(oc_sum - 1) >= 1e-6:	

            raise SystemExit(f"SystemExit: Occupancies summ up to {oc_sum}, Line 205")
        
        #Additionally maps resting states 
        calc_macro_pop(enum,all_complexes,resting_complexes,parameters)




        #simulate the condensed reactions
        resting_complexes = simulate_system(enum,parameters,resting_complexes,all_complexes,parameters["d_length"][d_seq_split[step]])



        
        oc_sum = 0
        macro_sum = 0

        for complex in resting_complexes:
            oc_sum += complex.occupancy
        for macro in enum._resting_macrostates:
            macro_sum += macro.occupancy

        logger.info(f"Sum of macrostate occupancies before cutoff: {macro_sum:.20f}")
        logger.info(f"Sum of occupancy before cutoff: {oc_sum:.20f}")
        oc_sum = round(oc_sum,10)
        if abs(oc_sum - 1) <= 1e-10:	
            for complex in resting_complexes:
                if not is_complex_in_all_complexes(complex,all_complexes):
                    raise SystemExit("resulting resting complex not in all complexes")
        else:
            raise SystemExit(f"SystemExit: Before cutoff Occupancies summ up to {oc_sum}")
        
        # apply cutoffs 
        
        enum = apply_cutoff(enum,parameters,all_complexes)

        
        resting_complexes = list(set(enum._resting_complexes))
        oc_sum = 0
        
        for complex in enum._resting_complexes:
            oc_sum += complex.occupancy
        oc_sum = round(oc_sum,10) #can be later adjusted

        # checks to see of occupancies sum up to 1
        if abs(oc_sum - 1) <= 1e-10:	
            folding_step_complexes.append(resting_complexes)
            for complex in resting_complexes:
                if not is_complex_in_all_complexes(complex,all_complexes):
                    raise SystemExit("resulting resting complex not in all complexes")
        else:
            raise SystemExit(f"SystemExit: Occupancies summ up to {oc_sum}")
        

    return folding_step_complexes


def simulate_system(enum,parameters,resting_complexes,all_complexes,new_domain_length):
    
    
    condensed_reactions = enum.condensation._condensed_reactions
    logger.info("\n\n\nBegin with simulating system")
    reactions = []

    oc_vector = []
    seen = set()
    if condensed_reactions and new_domain_length > 0:
        #create reactions list and occupancy_vector needed for crn simulations
        for reaction in condensed_reactions:
            #if reaction._reactants[0].name not in seen: 
                seen.add(reaction._reactants[0].name)
                
                reactants = [reaction._reactants[0].name]
                products = [reaction._products[0].name]

                occupancy1,occupancy2 = np.float128(0),np.float128(0) 

                logger.info("\n\n calc macro occupancies\n")
                for r_complex in reaction._reactants[0]._complexes: 
                        
                    #set conc to 0 for undefined complexes 
                    if r_complex.occupancy == None:
                        r_complex.occupancy = np.float128(0)
                    
                    
                    
                    
                    occupancy1 += r_complex.occupancy
                    
                    
                for p_complex in reaction._products[0]._complexes: 
                    occupancy2 += np.float128(p_complex.occupancy)


                rate = reaction._const
                reactions.append([reactants,products,[rate]])
                oc_vector.append(float(occupancy1))
                if occupancy2:
                    oc_vector.append(float(occupancy2))
                else:	
                    oc_vector.append(0)




        resulting_occupancies =	sim_condensed_rates(reactions,oc_vector,parameters,new_domain_length)
        logger.debug(f"resulting occupancies: {resulting_occupancies}")
        #Update occupancies



        resting_complexes = update_macrostates(resulting_occupancies,all_complexes = all_complexes,enum= enum,resting_complexes=resting_complexes,parameters=parameters)
        logger.debug(f"Numger of Resulting resting complexes {len(resting_complexes)}")
        
    return resting_complexes

def update_macrostates(result_occ,all_complexes,enum,resting_complexes,parameters):
    
    """
    parameters: 
        result_occ(dict): Dictionary, key = complex._name, value float
        all_complexes(dict): Id: [Complex,conc]
        enum(object): enumerated object from peppercorn
        resting_complexes(list)
    """

    logger.info("\n\nUpdating Macrostates\n\n")
    macro_seen = set()
    macro_sum = 0
    for reaction in result_occ.values():
        for key,value in reaction.items():
            for macrostate in enum._resting_macrostates:
                if any(complex_obj._name == key for complex_obj in macrostate.complexes):
                    macrostate.occupancy = np.float128(value) 
                    for complex in macrostate._complexes:
                        complex.occupancy = value/len(macrostate._complexes)
                    if macrostate not in macro_seen:
                        macro_sum += value
                        macro_seen.add(macrostate)

    occ_sum = round(sum([complex.occupancy for complex in resting_complexes]),4) 

    assert abs(occ_sum -1) < 1e-10 , f"Sum of occupancies is {occ_sum} not 1"
    calc_macro_pop(enum,all_complexes,resting_complexes,parameters)

    
    logger.info(f"Summe over resting complexes:{sum([complex.occupancy for complex in resting_complexes])}")
    return resting_complexes

def is_complex_in_all_complexes(complex,all_complexes):

    for id,value in all_complexes.items():
        if complex == value[0]:	
            return True
    
    return False


def sim_condensed_rates(reactants,concvect,parameters,d_length):
    """Script is based on Pillsimulator script from Peppercorn and uses CRNsimulator to simulate a chemical reaction network. 
    In this case it calculates the occupancies after accounting for the condensed rates.
    """
    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    
    parser.add_argument("-i", "--input_file", nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="Input file. If not provided, reads from stdin.")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity level. -v only shows peppercorn output")
    parser.add_argument("-cutoff", "--cutoff", action="store", type=valid_cutoff, default=float('-inf'),help="Cutoff value at which structures won't get accepted (default: -inf, valid range: 0 to 1)")
    parser.add_argument("-l", "--logic", action="store_true", default=False,help="Visualizes Logic domain pairings. Best used when analyzing cocopaths generated sequences. (default = False)")


    s_args = parser.parse_args()
    s_args.cutoff = float(parameters["cutoff"])

    #need to hardcode all args:
    logger.info("\n\n\nBegin simulation setup\n\n")

    

    s_args.pyplot_labels = None
    
    s_args.list_labels = None

    s_args.nxy = False

    s_args.t8 = 10 * d_length




    s_args.t0 = 0

    s_args.t_log = False

    s_args.t_lin = 2

    s_args.atol = None
    s_args.rtol = None
    s_args.mxstep = 0 
    s_args.labels_strict = False
    s_args.header = False

    s_args.pyplot = False
    s_args.labels = []
    s_args.pyplot_xlim = None 
    s_args.pyplot_ylim = None

    

    
    
    
    filename = "filename" + str(reactants[0])

    odename = "odename"


    V = []
    for reaction in reactants:
        V.append(reaction[0][0])
        V.append(reaction[1][0])	

    

    # Split CRN into irreversible reactions
    crn = reactants
    new = []
    for [r, p, k] in crn:
        assert not (None in k)
        if len(k) == 2:
            new.append([r, p, k[0]])
            new.append([p, r, k[1]])
        else:
            new.append([r, p, k[0]])
    crn = new

    RG = ReactionGraph(crn)
    

    C = []
    Vars = []
    seen = set()
    for s,c in zip(V,concvect):
        if s not in seen:
            seen.add(s)
            Vars.append(s)
            C.append(c)


    
    # set p0 
    p0 = []
    for i,s in enumerate(C,start=1):
        p0.append(str(i) + "=" + str(s))


    s_args.p0 = p0

    #if os.path.exists(filename):
    #    os.remove(filename)
    
    
    logger.debug(f"{s_args.p0 = }")
    logger.debug(f"{Vars = },{C = }")

    filename, odename = RG.write_ODE_lib(sorted_vars = Vars, concvect = C,
                                              jacobian= True,
                                              const = [False for x in Vars],
                                             filename = filename,
                                             odename = "odename")
        
        



    integrate = get_integrator(filename)




    #sys redirection necessary to redirect output of integrate
    sys.stderr = open(os.devnull, 'w')
    logger.info("\n\n\nBegin simulation\n\n")
    try:
        time,occupancies = integrate(s_args) 
    finally:
        sys.stderr = sys.__stderr__
    

    os.remove(filename)

    end_conc = [oc for oc in occupancies[1:]]
    resulting_concentrations = {}
    for i in range(len(reactants)):
        resulting_concentrations[i] = {}
        for complex, conc in zip(Vars,end_conc):
            resulting_concentrations[i][complex] = conc
    return resulting_concentrations



def calc_macro_pop(enum,all_complexes,resulting_complexes,parameters):
    """Function to calculate the distribution of occupancies in each Macrostate. 
    Uses the stationary distribution to distribute the occpancies. 

    Args: 
        enum(object): Enumerator object from peppercorn
        all_complexes(dict): Dictionary of all seen complexes in the simulation
        resulting_complexes(dict): Dictionary of the resulting complexes in this step. 

    
    """


    logger.info("\n______________\nCalc Macropop\n\n")
    resting_macrostates = enum._resting_macrostates

    stat_dist = enum.condensation.stationary_dist
    summe = 0 
    for macrostate in resting_macrostates:
        stat_dist_copy = dict(stat_dist[macrostate])
        
        macro_pop = np.float128(0)
        for complex in macrostate._complexes:
                
                # add to macrostate population if no occupancy is defined -> complex new 
                try:
                    macro_pop += complex.occupancy
                except:
                    macro_pop += 0

        macrostate.occupancy = macro_pop            
        
        for stat_complex,pop in stat_dist_copy.items():
            new_dist = pop * macro_pop
        
            if is_complex_in_all_complexes(stat_complex,all_complexes):
                
                for c_id, all_complex in all_complexes.items():
                    if stat_complex == all_complex[0]:
                        all_complexes[c_id] = [stat_complex,new_dist]
                
            else:
                
                all_complexes["Id_" + str(str(len(all_complexes) + 1))] = [stat_complex,new_dist]
                
        
            for r_complex in resulting_complexes:
                if stat_complex == r_complex: 
                    r_complex.occupancy = np.float128(new_dist)

    if len(resulting_complexes) == 1:
        resulting_complexes[0].occupancy = 1
  
    for r_complex in resulting_complexes:

            summe += np.float128(r_complex.occupancy)
    

    logger.debug(f"Sum over all Macrostates after cal macropop {summe:.20f} {type(summe)}")
    for complex in resulting_complexes:
        logger.info(f"{complex},{complex.occupancy}")
    
    assert round(summe,10) == 1, f'Occupancies sum up to {summe} and not 1' 

    return resulting_complexes
            
def apply_cutoff(enum,parameters,all_complexes):
    """Applies the user defined cutoff on the system. 
    First looks at the each macrostate if the whole macrostate is under the treshhold. 

    If macrostate < threshold: (Solution not fixed below is a possible solution)
        - look if there are outgoing condensed rates if yes readjust macrostate occupancies based on the outgoing rates of the removed macrostate
        - if no outgoing condensed rate, readjust whole system (don't know if this is best solution)

    If only a complex in a macrostate is under the threshold, the complex gets removed and the whole macrostate gets readjusted. 
    
    """
    logger.info("\n\nEnforcing Cutoffs\n\n")
    
    stat_dist = enum.condensation.stationary_dist

    #necessary to not reassign occupancies to cut complexes later on 
    cut_macrostates = set()
    macrostates = enum._resting_macrostates 
    


    rm_occ = 0 # keeps track of the allready removed occupancy of the system

    occupancy_groups = {}
    for macrostate in enum._resting_macrostates:
        occupancy = macrostate.occupancy
        if occupancy not in occupancy_groups:
            occupancy_groups[occupancy] = []
        occupancy_groups[occupancy].append(macrostate)

    # Iterate through sorted groups and shuffle within each group
    for occupancy, group in sorted(occupancy_groups.items()):
        random.shuffle(group)

        # Iterate through shuffled group
        for macro in group:

            if rm_occ + macro.occupancy > parameters["cutoff"]:
                break

            #cut_macrostates.add(macro)

            #enforce_cutoff_macrostate(macro, enum, all_complexes, cut_macrostates)


            enum = enforce_via_transient(macro,enum,all_complexes,cut_macrostates,parameters)



    oc_sum = 0
    """
    for macro in macrostates:
        for complex in macro._complexes:
            oc_sum += complex.occupancy
    if abs(oc_sum - 1) > 1e-10:
        raise SystemExit(f"SystemExit: After macro cutoff Occupancies summ up to {oc_sum}")
        """

    return enum
        

def enforce_via_transient(cut_macrostate,enum,all_complexes,cut_macrostates,parameter):

    """Idea: By setting all outgoing reactions of the macrostate to > k_fast resting_complex -> transient complex. 
    """
    
    #exiting the function if non cutable macrostates
    #no outgoing reaction
    out_rxns = []
    
    for reaction in enum.condensation._condensed_reactions:
        if reaction._reactants[0] == cut_macrostate:
            out_rxns.append(reaction)


    if len(out_rxns) == 0:
        
        
        logger.info(f"\n\n\n______No outgoing reaction of cut macrostate --> can't enforce cutoff for {cut_macrostate} _________")
        
        return enum

    cut_macrostates.add(cut_macrostate)

    output = enum.to_pil(condensed=True,detailed = False) 

            
    logger.debug(f"\n\nCut Macrostate: {cut_macrostate} {cut_macrostate.occupancy}  \n{cut_macrostates}\nBefore pruning: \n {output} \n\n\n")

    
    mult_factor = 10e6

    #safe for old_condensed for redistribution later
    old_condensed_rxn = [reaction for reaction in enum.condensation._condensed_reactions]

    #Condense network again
    cut_complexes = cut_macrostate._complexes


    #now problem with double condensation 
    for reaction in enum._reactions:
        
        # Maybe change to just outgoing reaction
        if reaction._reactants[0] in cut_complexes: 
            reaction._const = reaction._const * mult_factor
            #modified_reactions.append(reaction)
            #fin_complex = reaction._products[0]
    #enum.condensation._condensed_reactions = None #solution? 
            
    for complex in cut_complexes:
        if complex in enum._resting_complexes:
            enum._resting_complexes.remove(complex)
            enum._transient_complexes.append(complex)

    enum._resting_macrostates = None


    enum.condense()
    
    
    output = enum.to_pil(condensed=True,detailed = False) 

    logger.info(f"\n\n\nAfter Pruning: \n {output} \n\n\n")

    #readjust occupancies 
    cut_occ = cut_macrostate.occupancy
    
    logger.debug(f"Cutocc:{cut_occ}")
    tot_out = 0
    for reaction in old_condensed_rxn:
        if reaction._reactants[0].name == cut_macrostate.name:
            tot_out += reaction._const


    for reaction in old_condensed_rxn:
            if reaction._reactants[0].name == cut_macrostate.name:
                reaction._products[0].occupancy += cut_occ * (reaction._const/tot_out)

    
        
    
    sum = 0
    for complex in enum._resting_complexes:
        sum += complex.occupancy


    logger.debug(f"Sum: {sum}")

    stat_dist = enum.condensation.stationary_dist



    #redistribute occupancies in macrostate 
    
    for macro in enum.condensation.stationary_dist.keys(): 
        stat_dist = enum.condensation.stationary_dist
        if macro not in cut_macrostates:
            logger.debug(f"Macro:{macro}")
            macro_dist = dict(stat_dist[macro])

            for complex,value in macro_dist.items(): 
                new_occ = value * macro.occupancy
                
                complex.occupancy = new_occ


    enum._resting_macrostates = enum.condensation.resting_macrostates
   

    summe = 0
    for complex in enum._resting_complexes: 
        summe += complex.occupancy


    #update all_complexes
    for complex in enum._resting_complexes: 
        logger.debug(f"{complex},{complex.occupancy}")
        update_complex_in_all_complexes(complex,all_complexes)

    

    assert abs(summe - 1) <= 1e-10,f"Redistribution was not succesfull summe = {summe}"

    assert not any([cut_complex in enum._resting_complexes for cut_complex in cut_complexes ])

    return enum

def flatten(lst):
    result = []
    for item in lst:
        if isinstance(item, list):
            result.extend(flatten(item))
        else:
            result.append(item)
    return result

def enforce_cutoff_macrostate(macrostate,enum,all_complexes,cut_macrostates):
        
    logger.info(f'\n\n\nEnforcing Macrostate Cutoff for {macrostate} with occupancy {macrostate.occupancy}\n\n')

    cut_complexes = macrostate._complexes
    cut_occ = sum([complex.occupancy for complex in cut_complexes])

    #set occupancies of macrostate and complexes to 0 

    for complex in cut_complexes:
        complex.occupancy = 0 
    macrostate.occupancy = 0

    



    if len(enum.condensation._condensed_reactions) == 0:
        #redistribute removed occupancy to remaining macrostates 
        for macro in enum._resting_macrostates:
            if macro != macrostate:
                macro.occupancy += cut_occ*(1/(len(enum._resting_macrostates) - len(cut_macrostates)))

    else: 
        
        
        sum_rates = sum([reaction._const for reaction in enum.condensation.condensed_reactions if all([macrostate == reaction._reactants[0],  macrostate in cut_macrostates])])
        
        if sum_rates == 0: #no outgoing reactions same fate as above

            for macro in enum._resting_macrostates: 
                if macro != macrostate and macro not in cut_macrostates:
                    macro.occupancy += cut_occ*(1/(len(enum._resting_macrostates) - len(cut_macrostates)))
                  

            assert abs(1 - sum([macro.occupancy for macro in enum._resting_macrostates])) < 1e-10,f"Macro occupancies sum up to {sum([macro.occupancy for macro in enum._resting_macrostates])}"
            
            #raise SystemExit("Cut complex has no outgoing reaction") 
        else:

            for reaction in enum.condensation.condensed_reactions:
                if macrostate == reaction._reactants[0] and macrostate in cut_macrostates: #check if the cut macrostate has an outgoing reaction
                    reaction._products[0].occupancy += cut_occ * (reaction._const/sum_rates)
                    

    
    stat_dist = enum.condensation.stationary_dist
    
    #redistribute occupancies in macrostate 

    macro_dist = dict(stat_dist[macrostate])
    
    for macro in enum._resting_macrostates: 
        
        
        if macro not in cut_macrostates:
            macro_dist = dict(stat_dist[macro])

            for complex,value in macro_dist.items(): 
                new_occ = value * macro.occupancy
                
                complex.occupancy = new_occ
            
    summe = 0
    for macro in enum._resting_macrostates: 
        for complex in macro._complexes: 
            summe += complex.occupancy

    assert abs(summe - 1) <= 1e-10,f"Redistribution was not succesfull summe = {summe}"


    #raise SystemExit('Whole macrostate is under threshold')


def enforce_cutoff_complex(enum,macrostate,cut_complex,parameters,all_complexes,cut_complexes):

    logger.info(f'\n\n\nEnforcing Complex Cutoff for {cut_complex} with occupancy {cut_complex.occupancy}\n\n')
    
    #check if total macrostate gets cut
    macro_occ = np.float128(0)

    for complex in macrostate._complexes:
        macro_occ += complex.occupancy
   

    free_occ = cut_complex.occupancy
    stat_dist = enum.condensation.stationary_dist


    stat_dist_copy = dict(stat_dist[macrostate])
    #macro_pop = np.float128(macrostate.occupancy)

    sum_remaining_occs = macro_occ - free_occ

    logger.debug(f"Sum remaining occs: {sum_remaining_occs}") 
    check_num = 0 
    for complex in macrostate._complexes:
            if complex != cut_complex and complex not in cut_complexes:
                new_dist = np.float128(complex.occupancy) / sum_remaining_occs  
                new_dist = new_dist * macro_occ
                
                if is_complex_in_all_complexes(complex,all_complexes):
                    
                    for c_id, all_complex in all_complexes.items():
                        if complex == all_complex[0]:
                            all_complexes[c_id][1] = new_dist
                            complex.occupancy = new_dist
                            check_num = new_dist


    cut_complex.occupancy = 0 

    update_complex_in_all_complexes(cut_complex,all_complexes)

    updated_macro_occupancy = 0 

    for complex in macrostate._complexes:
        updated_macro_occupancy += complex.occupancy
    
    

    assert abs(updated_macro_occupancy - macro_occ) <= 1e-10,f'Updated Macro Occupancy {updated_macro_occupancy:.10f} initial macro_occ {macro_occ:.10f}' 

def enumerate_step(complexes, reactions, parameter, all_complexes):
    """Takes complexes in form of PepperComplexes and uses Peppercorn Enumerator to find possible structures and their respective occupancies. 
    
    Args:
        complexes(dict): {Name:PepperComplex}
        reactions(class): PepperReaction not really used in our program
    Returns:
        resulting_complexes(list): Resulting PepperComplexes
    """
    k_slow = parameter["k_slow"]
    k_fast = parameter['k_fast']
    
    if bind21 in BI_REACTIONS:
        BI_REACTIONS.remove(bind21)

    logger.info(f"\n\n\nBeginning of Enumerate step:\nNumber of Complexes to be enumerated:{len(complexes)}\nReactions:{reactions}\nParameters:{parameter}\n\n")
    init_cplxs = [x for x in complexes.values() if x.concentration is None or x.concentration[1] != 0 and x.occupancy != 0]
    name_cplxs = list(complexes.values())
    enum = Enumerator(init_cplxs, reactions, named_complexes=name_cplxs)

    if k_slow: 
        enum.k_slow = k_slow



    #need to find suitable k_fast
    if k_fast: 
        enum.k_fast= k_fast
    
    
    # Start to enumerate
    enum.enumerate()
    
    
 
    enum.condense()

    logger.info("\n\nDone Enumerating\n\n")
    resulting_complexes = [cplx for cplx in natsorted(enum.resting_complexes)] 
    
    transient_complexes = [cplx for cplx in natsorted(enum.transient_complexes)]
    
    
    output = enum.to_pil(condensed=True,detailed = True) 

        
    logger.warning(f"\n\n\n\nOutput: \n {output} \n\n\n")
    return resulting_complexes, transient_complexes, enum

    



def map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters):

    
    logger.info(f"\n\nMap Transient states\n\n")
    
    
    new_complexes = {}

    for key, complex in copy.deepcopy(all_complexes).items():
        for t_complex in transient_complexes:
            # Check if the current complex matches a transient complex
            if complex[0] == t_complex:
                for tup, va in enum.condensation.cplx_decay_prob.items():
                    
                    if tup[0] == t_complex:
                        macrostate = tup[1][0]
                        stat_dist = enum.condensation.stationary_dist[macrostate]

                        
                        for stat_complex, value in stat_dist.items():
                            ##new_conc = curr_conc + (va(exit prob of transient) * value(stationary of resting state in macrostate) * t_c(concentration of transient state)) 
                            try:
                                if stat_complex.occupancy:
                                    new_conc = np.float128(stat_complex.occupancy) + (va * value * np.float128(t_complex.occupancy)) 
                                else:
                                    new_conc = va * value * np.float128(t_complex.occupancy)
                            except: 
                                    new_conc = va * value * np.float128(t_complex.occupancy)
                            stat_complex.occupancy = np.float128(new_conc)


                            # Check if the new complex is not already in all_complexes
                            if not is_complex_in_all_complexes(stat_complex,all_complexes) and not is_complex_in_all_complexes(stat_complex,new_complexes):
                                all_complexes["Id_" + str(len(all_complexes) + 1)] = [stat_complex, new_conc]
                                
                            elif is_complex_in_all_complexes(stat_complex,all_complexes):
                                update_complex_in_all_complexes(stat_complex,new_complexes)	

    logger.info("\n\nEnd Mapping transient states")
    
    

def update_complex_in_all_complexes(complex,all_complexes):
    
    for key,a_complex in all_complexes.items():
        if a_complex[0] == complex:
            all_complexes[key] = [complex,complex.occupancy]


    

def write_output(final_structures,d_seq,parameters = None):

    data_output = ""
    

    spacer_indices = [index for index, entry in enumerate(d_seq.split()) if entry.startswith('S')]
    ts = 0
    data_output += ("\nResting Complexes after each Spacer:\n\n")
    data_output += "Transcription Step |  Occupancy   |    Structure  \n"
    for x in final_structures:
        if x and x[0].kernel_string[-2] == "S": 
            for complex in x:
                if complex.occupancy >= 0.001:
                    data_output += f"{spacer_indices[ts]+ 1:3}  |	{complex.occupancy:^7.5f}	|	{complex.kernel_string} \n"
            ts += 1 

            data_output += "\n"
    
    
    ts = 0
    data_output += ("\n\nOnly Logic Domain pairings:\n\n")
    data_output += "Transcription Step | Occupancy  |  Structure	 \n"



    struct_list = []
    for x in final_structures:
        struct_dict = {}
        if x and x[0].kernel_string[-2] == "S": 
            data_output += "\n"
            for complex in x:
                if complex.occupancy >= 0.0001:
                    kernel_string = kernel_to_dot_bracket(complex.kernel_string)
                    db_struct = (only_logic_domain_struct(d_seq.split(),kernel_string))
                    data_output += f"{spacer_indices[ts]+ 1:3}  |	{np.float128(complex.occupancy):^8.5f}   |   {db_struct} 	\n"
                    if db_struct in struct_dict: 
                        struct_dict[db_struct] += complex.occupancy
                    else: 	
                        struct_dict[db_struct] = complex.occupancy	
            ts += 1
            struct_list.append(struct_dict)

    data_output += ("\n\nOnly Logic Domain pairings condensed:\n\n")
    data_output += "Transcription Step | Occupancy  |  Structure	 \n"
    t_s = 0
    for s_dict in sorted(struct_list, key=lambda d: len(list(d.keys())[0]), reverse=False):
        data_output += "\n"
        for struct, occ in sorted(s_dict.items(), key=lambda item: item[1], reverse=True):
            data_output += f"{spacer_indices[t_s]+ 1:3}  |  {float(occ):8.5f}   |   {struct}  \n"
        t_s += 1


    return data_output

    

def input_parsing(d_seq, complexes,parameters):
    """Takes domain level sequence and structures and formats those into an input format for peppercorn

    Args:
        d_seq(str) : domain level sequence
        complexes(dict)	: dict of the complexes in the form of {Name:[kernel structure,population]}
    """
    
    logger.info(f"Input Parsing: Structure input \n {complexes}")
    unique_domains = set([domain for domain in d_seq])


    system_input = f""
    for unique_domain in unique_domains:
        system_input += f"length {unique_domain} = {parameters['d_length'][unique_domain]} \n"
      
    
    system_input += "\n"
    for name, lst in complexes.items():
            if lst[1] > 0:
                system_input += f"{name} = {lst[0]}\n"
        

    complexes, reactions = read_pil(system_input)
    
    for key,complex in complexes.items():
        complex.occupancy = np.float128(None)
        complex.id = None

    return complexes, reactions


def extract_domain_sequence(input_lines):
    try:
    
        for line in input_lines:
            if not line.startswith("#") and not line.startswith("[") and not line.startswith("I") and len(line) > 1 and not line.startswith("Resulting Domain Level sequence:"):
                result_sequence = line
                return result_sequence
            elif line.startswith("Resulting Domain Level sequence:"):
                result_sequence = line.split(":",1)[1].strip()
                return result_sequence
    except:
        raise ImportError("Check Input: can't find anything")
        


def extract_afp(input_lines):
    try: 
        if input_lines[0] == "\n":
            return eval(next((line for line in input_lines if line.startswith("[")),None))
            
    except:
        return None 


def set_verbosity(console_handler, verbosity):
    if verbosity == 0:
        console_handler.setLevel(logging.CRITICAL)
        file_handler.setLevel(logging.CRITICAL)

    elif verbosity == 1:
        console_handler.setLevel(logging.WARNING)
        file_handler.setLevel(logging.WARNING)

    elif verbosity == 2:
        console_handler.setLevel(logging.INFO)
        file_handler.setLevel(logging.INFO)

    elif verbosity >= 3:
        console_handler.setLevel(logging.DEBUG)
        file_handler.setLevel(logging.DEBUG)


def valid_cutoff(value):
    value = float(value)
    if 0 <= value < 1:
        return value
    else:
        raise argparse.ArgumentTypeError(f"Cutoff must be between 0 and 1, but got {value}")

def verify_domain_foldingpath(afp,domain_seq,parameters,simulated_structures = None):

    """Function to verfiy that a domain level sequence folds according to the afp, by which it was generated. 

    Args: 
        afp(list): Each entry corresponds to a step in the folding path. 
        domain_seq(string): domain level sequence

    Returns: 
        if correct folding path:
            domain_path(list): domain level folding path 

        False folding path: 
            False 
    """

    #Minimum occupancy of desired structure for succesfull sim
    threshold = 0.5

    #converting afp to domain level afp for easier comparison

    domain_afp = afp_to_domainfp(afp,domain_seq)

    if simulated_structures == None:
        simulated_structures = run_sim(domain_seq,parameters)


    sim_domain_fp = [[] for _ in simulated_structures]

    for i,step in enumerate(simulated_structures):
        for complex in step: 
            sim_domain_fp[i].append(''.join(complex._structure))
    
    spacer_indices = [index - 1 for index, entry in enumerate(domain_seq.split()) if entry.startswith('S')]
    sim_unit_path = [simulated_structures[i] for i in spacer_indices]



    dominant_structures_fp = []

    logic_indices = [index for index,entry in enumerate(domain_seq.split()) if entry.startswith("L") ]
    for x in range(len(sim_unit_path)):
        target_occ = 0
        max_occ = 0

        for struct in sim_unit_path[x]: 
            
            #convert struct here to only logic domains and compare to afp

            logic_struct = ''.join([struct._structure[i] for i in logic_indices[0:x]])


            
            if "".join(struct._structure) == domain_afp[x]: 
                target_occ  += struct.occupancy
                if struct.occupancy >= max_occ:
                    max_occ = struct.occupancy
                    max_occ_struct = "".join(struct._structure)
        
        if target_occ >= threshold:
            dominant_structures_fp.append(max_occ_struct)
            
            
        else:
            return False


    assert len(dominant_structures_fp) == len(afp),f"Some structures are missing {dominant_structures_fp = }"


    return dominant_structures_fp
                





def main():
    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    
    parser.add_argument("-i", "--input_file", nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="Input file. If not provided, reads from stdin.")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.00001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity level. -v only shows peppercorn output")
    parser.add_argument("-l", "--logic", action="store_true", default=False,help="Visualizes Logic domain pairings. Best used when analyzing cocopaths generated sequences. (default = False)")
    parser.add_argument("-cutoff", "--cutoff", action="store", type=valid_cutoff, default=float('-inf'),help="Cutoff value at which structures won't get accepted (default: -inf, valid range: 0 to 1)")
    
    

    args = parser.parse_args()
    args.cutoff = float(args.cutoff)
    set_verbosity(logger,args.verbose)
    console_handler.setLevel(logging.DEBUG)
    
    

    # Read input     
    d_length = {} #dict where domain lengths get saved

    if args.input_file.isatty():
        print("Please enter a domain level sequence:")
        d_seq = input()
        if len(d_seq) == 0:
            raise SystemExit("No Input given")
        afp = None
    else:
        file_extension = os.path.splitext(args.input_file.name)[1]  # Get the file extension

        if file_extension == ".pil":
            afp = None

            pil_input = read_pil(args.input_file,True)
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

                



            else:
                raise SystemExit("SystemExit:More than one kernel sequence in input. We can only simulate one kernel string at a time. ")

            
            
        else:
            try:
                input_lines = args.input_file.readlines()
                d_seq = extract_domain_sequence(input_lines)
                afp = extract_afp(input_lines)
                args.input_file.close()
            except: 
                raise SystemExit('No valid input was given. Please make sure it is in a valid pil format.')

    # create dictionary for domain lengths 
    if not d_length:

        for domain in d_seq.split():
            if domain[0] == "L":
                d_length[domain] = 8
            elif domain[0] == 'S':
                d_length[domain] = round(int(domain[1]) * 4)  
            else: 
                d_length[domain] = 3 
        

    parameters = {"k_slow": args.k_slow,'k_fast': args.k_fast, "cutoff": args.cutoff,"d_length":d_length,"d_seq":d_seq,"logic":args.logic}


    print(parameters)
    
    logger.warning(parameters)


    print("Given Domain Sequence:", d_seq)

   
    if afp == None:
        afp = afp_terminal_input()  


    #______Running_Simulation___________# 
    simulated_structures = run_sim(d_seq, parameters)

    
    last_step = simulated_structures[-1]
    if args.logic: 
        #Print continous output
        for complex in last_step:
            kernel_string = kernel_to_dot_bracket(complex.kernel_string)
            db_struct = (only_logic_domain_struct(d_seq.split(),kernel_string))
            #occupancy = np.float128(complex.occupancy)
            if round(complex.occupancy,4) > 1e-4:
                print(f"END\t|\t{complex.occupancy:^8.4f}\t|\t{db_struct:16}|\t{complex.kernel_string}")
        print()    
    else: 
        #Print continous output
        for complex in last_step:
            kernel_string = kernel_to_dot_bracket(complex.kernel_string)
            #occupancy = np.float128(complex.occupancy)
            if round(complex.occupancy,4) > 1e-4:
                print(f"END\t|\t{complex.occupancy:^8.4f}\t|\t{complex.kernel_string}")
        print()  

    
    #_____Writing_and_printing_output___# 
    if args.logic:
        output = ""

        d_seq
        ts = 2

        
        if 'S' in d_seq:
            output += write_output(simulated_structures,d_seq)
        
        output += f"\n\nFollowing sequence was simulated:\n{d_seq}"

        if afp: 
            output += f"\n\nFollowing AFP was given for Sequence design"
            for step in afp:
                
                output += f"{step}\n"

        print(output)

    #______Verify_Simulation
    print("---------------")
    dominant_path = verify_domain_foldingpath(afp,d_seq,parameters,simulated_structures)

    if dominant_path:
        print("\n\nSimulation verified that the domain sequences folds according to the input abstract folding path.")

    else: 
        print("\n\nSimulated path differs from the abstract folding path. Please adjust domain lengths.\n\n")    
    logger.removeHandler(file_handler)  # Remove the file handler from the logger
    file_handler.close()  # Close the file handler to release the file
    os.remove("cocosim.log")

if __name__ == "__main__":
    main()
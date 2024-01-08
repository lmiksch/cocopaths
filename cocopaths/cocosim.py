"""Compiler for cotranscriptional folding simulator

Simulates cotranscriptional folding of a domain level sequence.

uses Function from Peppercornenumerator to find reaction types/ possible products of reactions
"""
import sys
from peppercornenumerator import peppercorn 
from dsdobjects.objectio import read_pil as dsd_read_pil
from peppercornenumerator.input import read_pil
from peppercornenumerator import Enumerator 
from peppercornenumerator.enumerator import BI_REACTIONS
from peppercornenumerator.reactions import bind21
from peppercornenumerator.objects import PepperComplex
from crnsimulator import ReactionGraph, get_integrator
from crnsimulator.odelib_template import add_integrator_args
import numpy as np
from io import StringIO
from natsort import natsorted
from .utils import cv_db2kernel, kernel_to_dot_bracket, only_logic_domain_struct
import argparse
import logging 
import numpy as np 
import copy
import os


logger = logging.getLogger('cocosim')
console_handler = logging.StreamHandler()
formatter = logging.Formatter('# %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)




def run_sim(d_seq, parameters,args):
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
    
    complexes, reactions = input_parsing(d_seq_split[0:2], {"E0" : [cv_db2kernel(d_seq_split[0:2], "..")]},parameters)
    

    resulting_complexes,transient_complexes,enum = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters, all_complexes=all_complexes)
    
    
    
    map_transient_states(resulting_complexes,transient_complexes,all_complexes,enum,args)
    
    
    
    # all complex keeps track of complexes throughtout the sim
    

    
    #only in first step this should only result in 1 structure
    for complex in resulting_complexes:
        all_complexes["Id_" + str(str(len(all_complexes) + 1))] = [complex,complex.occupancy]
        complex.occupancy = 1
        
    

    folding_step_complexes.append(resulting_complexes)	
   
    
    #Print continous output
    for x, complex in complexes.items():
        
        occupancy = np.float128(complex.occupancy)
        
        print(f"{2:3}   |	      {occupancy:8.4f}   | {complex.kernel_string}")
    print()    
    

    for step in range(2, len(d_seq_split)):
        logger.warning("\n\n\n______________________________________________")
        logger.warning(f"Step: {step} Current Complexes: {parameters}")

        next_complexes = []
        
        
        #put here add one domain to struct and give new name
    
        new_complexes = {}
        
        occ_sum = 0 
        for complex in folding_step_complexes[-1]:
            occ_sum += complex.occupancy

        assert abs(occ_sum - 1 ) < 1e-10,"Difference at the beginning -> previous step was fucked"

        #______Function_for_extending_and_updating_names__ 

        old_names = []
        for complex in folding_step_complexes[-1]:
            
            c_name = complex._name
            next_struct = complex.kernel_string + " " + d_seq_split[step]

            new_complexes["E" + str(used_structure_names)] = [next_struct,complex.occupancy]
            
            used_structure_names += 1
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
        for x, complex in complexes.items():
            try:
                occupancy = np.float128(complex.occupancy)
                total_occupancy += occupancy
                print(f"{step + 1:3}   |	      {occupancy:8.4f}   | {complex.kernel_string}")
            except: 
                pass
                print(f"{step + 1:3}   |	{0:8.4f}      | {complex.kernel_string} \n")
        print()


        logger.debug(f"Step:{step} Total Occupancy:{total_occupancy:20.20f}\n")
        assert abs(total_occupancy - 1) < 1e-10,f"Occupancy is not equal 1 difference is {abs(total_occupancy - 1)}"         
        resting_complexes,transient_complexes,enum = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters,all_complexes=all_complexes)
        
        
        
        #checks if there is a transient state which needs to be mapped
        
        map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,args)
        resting_complexes = list(set(resting_complexes))
        oc_sum = 0
        for complex in resting_complexes:
            try: 
                oc_sum += complex.occupancy 
            except: 
                pass
        oc_sum = round(oc_sum,8) #can be later adjusted       
        if abs(oc_sum - 1) >= 1e-6:	

            raise SystemExit(f"SystemExit: Occupancies summ up to {oc_sum}")
        
        #Additionally maps resting states 
        calc_macro_pop(enum,all_complexes,resting_complexes,args)


        #simulate the condensed reactions
        resting_complexes = simulate_system(enum,args,resting_complexes,all_complexes,parameters["d_length"][d_seq_split[step]])


        resting_complexes = list(set(resting_complexes))
        oc_sum = 0
        for complex in resting_complexes:
            oc_sum += complex.occupancy 
        oc_sum = round(oc_sum,6) #can be later adjusted


        # checks to see of occupancies sum up to 1
        if abs(oc_sum - 1) <= 1e-6:	
            folding_step_complexes.append(resting_complexes)
            for complex in resting_complexes:
                if not is_complex_in_all_complexes(complex,all_complexes):
                    raise SystemExit("resulting resting complex not in all complexes")
        else:
            raise SystemExit(f"SystemExit: Occupancies summ up to {oc_sum}")
        
        # apply cutoffs 

        apply_cutoff(enum,parameters,all_complexes)



    

    return folding_step_complexes


def simulate_system(enum,args,resting_complexes,all_complexes,new_domain_length):
    
    
    condensed_reactions = enum.condensation._condensed_reactions
    logger.info("Begin with simulating system")
    reactions = []

    oc_vector = []
    seen = set()
    if condensed_reactions:
        #create reactions list and occupancy_vector needed for crn simulations
        for reaction in condensed_reactions:
            if reaction._reactants[0].name not in seen: 
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




        resulting_occupancies =	sim_condensed_rates(reactions,oc_vector,args,new_domain_length)
        #Update occupancies

        resting_complexes = update_macrostates(resulting_occupancies,all_complexes = all_complexes,enum= enum,resting_complexes=resting_complexes,args=args)

    return resting_complexes

def update_macrostates(result_occ,all_complexes,enum,resting_complexes,args):
    
    """
    Args: 
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
                if any(complex_obj._name == key for complex_obj in macrostate.complexes) and value > args.cutoff:
                    macrostate.occupancy = np.float128(value) 
                    if macrostate not in macro_seen:
                        macro_sum += value
                        macro_seen.add(macrostate)

    assert round(sum([complex.occupancy for complex in resting_complexes]),4) == 1
    calc_macro_pop(enum,all_complexes,resting_complexes,args)


    
    logger.info(f"Summe over resting complexes:{sum([complex.occupancy for complex in resting_complexes])}")
    return resting_complexes

def is_complex_in_all_complexes(complex,all_complexes):

    for id,value in all_complexes.items():
        if complex == value[0]:	
            return True
    
    return False


def sim_condensed_rates(reactants,concvect,args,d_length):
    """Script is based on Pillsimulator script from Peppercorn and uses CRNsimulator to simulate a chemical reaction network. 
    In this case it calculates the occupancies after accounting for the condensed rates.
    """
    

    #need to hardcode all args:
    logger.info("\n\n\nBegin simulation setup\n\n")

    args.pyplot_labels = None
    
    args.list_labels = None

    args.nxy = False

    args.t8 = 0.02 * d_length

    

    args.t0 = 0

    args.t_log = False

    args.t_lin = 2

    args.atol = None
    args.rtol = None
    args.mxstep = 0 
    args.labels_strict = False
    args.header = True

    args.pyplot = False
    args.labels = []
    args.pyplot_xlim = None 
    args.pyplot_ylim = None

    

    p0 = []
    for i,s in enumerate(concvect,start=1):
        p0.append(str(i) + "=" + str(s))


    
    args.p0 = p0
    
    filename = "filename"

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


    args.p0 = p0
    
    filename, odename = RG.write_ODE_lib(sorted_vars = Vars, concvect = C,
                                              jacobian= True,
                                              const = [False for x in Vars],
                                             filename = filename,
                                             odename = odename)
        

    integrate = get_integrator(filename)


    #sys redirection necessary to redirect output of integrate
    sys.stderr = open(os.devnull, 'w')
    logger.info("\n\n\nBegin simulation\n\n")

    try:
        time,occupancies = integrate(args) 
    finally:
        sys.stderr = sys.__stderr__
    


    

    end_conc = [oc for oc in occupancies[1:]]
    resulting_concentrations = {}
    for i in range(len(reactants)):
        resulting_concentrations[i] = {}
        for complex, conc in zip(Vars,end_conc):
            if float(conc) >= args.cutoff:
                resulting_concentrations[i][complex] = conc
            else: 
                resulting_concentrations[i][complex] = 0


    return resulting_concentrations



def calc_macro_pop(enum,all_complexes,resulting_complexes,args):
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
        #print('macrostate',macrostate)
        stat_dist_copy = dict(stat_dist[macrostate])
        
        macro_pop = np.float128(0)
        for complex in macrostate._complexes:
                #print("complex",complex)
                # add to macrostate population if no occupancy is defined -> complex new 
                try:
                    macro_pop += complex.occupancy
                    #print(complex.occupancy)
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
                #logger.debug(f"New_dist: {new_dist}")
                if stat_complex == r_complex: 
                    r_complex.occupancy = np.float128(new_dist)

  
    logger.debug("endo of update macrostates")
    for r_complex in resulting_complexes:
            #logger.debug(f"New_dist: {new_dist}")
            logger.debug(r_complex)
            summe += np.float128(r_complex.occupancy)

    logger.debug(f"Sum over all Macrostates {summe:.20f} {type(summe)}")


    
    assert round(summe,8) == 1, f'Occupancies sum up to {summe} and not 1' 
            
def apply_cutoff(enum,parameters,all_complexes):
    """Applies the user defined cutoff on the system. 
    First looks at the each macrostate if the whole macrostate is under the treshhold. 

    If macrostate < threshold: (Solution not fixed below is a possible solution)
        - look if there are outgoing condensed rates if yes radjust macrostate occupancies based on the outgoing rates of the removed macrostate
        - if no outgoing condensed rate, readjust whole system (don't know if this is best solution)

    If only a complex in a macrostate is under the threshold, the complex gets removed and the whole macrostate gets readjusted. 
    
    
    """

    macrostates = enum._resting_macrostates 

    for macrostate in macrostates:
        macro_occ = sum([occ.occupancy for occ in macrostate._complexes]) 
        if macro_occ <= parameters['cutoff']:
            print("Macrostate to be cut and occupancy",macrostate,macro_occ)
            enforce_cutoff_macrostate(macrostate)
        
        for complex in macrostate._complexes: 
            if complex.occupancy <= parameters['cutoff']:
                enforce_cutoff_complex(enum,macrostate,complex,parameters,all_complexes)

        
def enforce_cutoff_macrostate(macrostate):


    print("Macrostate to be cut ",macrostate,macrostate.occupancy)
    raise SystemExit('Whole macrostate is under threshold')


def enforce_cutoff_complex(enum,macrostate,cut_complex,args,all_complexes):

    logger.debug(f'\n\n\nEnforcing Cutoff for {cut_complex} with occupancy {cut_complex.occupancy}\n\n')
    
    #check if total macrostate gets cut
    macro_occ = np.float128(0)

    for complex in macrostate._complexes:
        macro_occ += complex.occupancy
   

    free_occ = cut_complex.occupancy
    stat_dist = enum.condensation.stationary_dist


    stat_dist_copy = dict(stat_dist[macrostate])
    #macro_pop = np.float128(macrostate.occupancy)

    sum_remaining_occs = macro_occ - free_occ

    logger.debug('Sum remaining occs: ',sum_remaining_occs) 
    check_num = 0 
    for complex in macrostate._complexes:
            if complex != cut_complex:
                logger.debug("1: ",complex,complex.occupancy,type(complex.occupancy),type(free_occ))
                new_dist = np.float128(complex.occupancy) / sum_remaining_occs  
                logger.debug("New dist",new_dist)
                new_dist = new_dist * macro_occ
                logger.debug("after * macro_occ New dist",new_dist)
                
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
        logger.debug(complex,complex.occupancy)
        updated_macro_occupancy += complex.occupancy
        
    logger.debug(free_occ, type(free_occ),type(check_num),check_num)

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
    init_cplxs = [x for x in complexes.values() if x.concentration is None or x.concentration[1] != 0]
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
    
    
    output = enum.to_pil(condensed=True,detailed = True) # should be not condensed for same output format like peppercorn
    

    

    macro_states = enum._resting_macrostates
    
    logger.warning(f"\n\n\n\nOutput: \n {output} \n\n\n")
    return resulting_complexes, transient_complexes, enum

    



def map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,args):

    
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

                            if new_conc > args.cutoff:#exclude structures below the cutoff 
            
                                # Check if the new complex is not already in all_complexes
                                if not is_complex_in_all_complexes(stat_complex,all_complexes) and not is_complex_in_all_complexes(stat_complex,new_complexes):
                                    all_complexes["Id_" + str(len(all_complexes) + 1)] = [stat_complex, new_conc]
                                    
                                elif is_complex_in_all_complexes(stat_complex,all_complexes):
                                    update_complex_in_all_complexes(stat_complex,new_complexes)	

                                
    #all_complexes.update(new_complexes)
    logger.info("\n\nEnd Mapping transient states")

    

def update_complex_in_all_complexes(complex,all_complexes):
    
    for key,a_complex in all_complexes.items():
        if a_complex[0] == complex:
            all_complexes[key] = [complex,complex.occupancy]


    

def write_output(final_structures,d_seq,parameters = None):

    data_output = ""
    

    ts = 1
    data_output += ("\nResting Complexes after each Spacer:\n\n")
    data_output += "Transcription Step |  Occupancy   |    Structure  \n"
    for x in final_structures:
        if x and x[0].kernel_string[-2] == "S": 
            for complex in x:
                if complex.occupancy >= 0.001:
                    data_output += f"{ts:3}  |	{complex.occupancy:7.5f}	|	{complex.kernel_string} \n"
            ts += 1 

            data_output += "\n"
    
    
    ts = 1
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
                    
                    data_output += f"{ts:3}  |	{np.float128(complex.occupancy):8.5f}   |   {db_struct} 	\n"
                    if db_struct in struct_dict: 
                        struct_dict[db_struct] += complex.occupancy
                    else: 	
                        struct_dict[db_struct] = complex.occupancy	
            ts += 1
            struct_list.append(struct_dict)

    data_output += ("\n\nOnly Logic Domain pairings condensed:\n\n")
    data_output += "Transcription Step | Occupancy  |  Structure	 \n"
    t_s = 1
    for s_dict in sorted(struct_list, key=lambda d: len(list(d.keys())[0]), reverse=False):
        data_output += "\n"
        for struct, occ in sorted(s_dict.items(), key=lambda item: item[1], reverse=True):
            data_output += f"{t_s:3}  |  {float(occ):8.5f}   |   {struct}  \n"
        t_s += 1


    return data_output

    

def input_parsing(d_seq, complexes,parameters):
    """Takes domain level sequence and structures and formats those into an input format for peppercorn

    Args:
        d_seq(str) : domain level sequence
        complexes(dict)	: dict of the complexes in the form of {Name:[kernel structure,population]}
    """
    toehold_length = 5
    logic_length = 8
    space_length = 8
    logger.info(f"Input Parsing: Structure input \n {complexes}")
    unique_domains = set([domain.replace('*', '') for domain in d_seq])

    system_input = f""
    for unique_domain in unique_domains:
        if unique_domain.islower():
            system_input += f"length {unique_domain} = {toehold_length} \n"
        elif unique_domain.startswith("L"):
            system_input += f"length {unique_domain} = {logic_length} \n"
        elif unique_domain.startswith("S"):
            system_input += f"length {unique_domain} = {space_length} \n"

    system_input += "\n"

    for name, lst in complexes.items():
        system_input += f"{name} = {lst[0]}\n"
        


    logger.warning(f"\n\n\nSystem Input \n\n{system_input}\n\n")
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
    elif verbosity == 1:
        console_handler.setLevel(logging.WARNING)
    elif verbosity == 2:
        console_handler.setLevel(logging.INFO)
    elif verbosity == 3:
        console_handler.setLevel(logging.DEBUG)


def main():
    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    
    parser.add_argument("-i", "--input_file", nargs='?', type=argparse.FileType('r'), default=sys.stdin,
                        help="Input file. If not provided, reads from stdin.")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity level. -v only shows peppercorn output")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    


    args = parser.parse_args()
    args.cutoff = float(args.cutoff)
    set_verbosity(logger,args.verbose)

    if args.input_file.isatty():
        print("Please enter a domain level sequence:")
        d_seq = input()
        if len(d_seq) == 0:
            raise SystemExit("No Input given")
        afp = None
    else:
        input_lines = args.input_file.readlines()
        d_seq = extract_domain_sequence(input_lines)
        afp = extract_afp(input_lines)
        args.input_file.close()
        
        print("Following AFP was given for Sequence design:")
        for step in afp:
            print(step)
        print("\n")


    # create dictionary for domain lengths 
    d_length = {}

    for domain in d_seq.split():
        if domain[0] == "L" or domain[0] == "S":
            d_length[domain] = 8
        else: 
            d_length[domain] = 3 


    parameters = {"k_slow": args.k_slow,'k_fast': args.k_fast, "cutoff": args.cutoff,"d_length":d_length}



    
    logger.warning(parameters)


    print("Given Domain Sequence:", d_seq)

    




    #______Running_Simulation___________# 
    simulated_structures = run_sim(d_seq, parameters,args)
    

    #_____Writing_and_printing_output___# 
    output = ""

    
    ts = 1
    output += ("\nResting Complexes after each transcription step:\n\n")
    output += "Transcription Step |    Occupancy | Structure \n"
    for x in simulated_structures: 
        for complex in x:
            output += f"{ts:3}   |	{complex.occupancy:7.5f}      |   {complex.kernel_string}\n"
        ts += 1 
        output += "\n"
    
    if 'S' in d_seq:
        output += write_output(simulated_structures,d_seq)
    
    output += f"\n\nFollowing sequence was simulated:\n{d_seq}"

    if afp: 
        output += f"\n\nFollowing AFP was given for Sequence design"
        for step in afp:
            
            output += f"{step}\n"

    print(output)




if __name__ == "__main__":
    main()
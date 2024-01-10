"""pytest moduel to test cocosim
"""
import pytest
import logging, os, re, subprocess
from unittest.mock import patch
from cocopaths.cocosim import map_transient_states,calc_macro_pop,enumerate_step,run_sim, enforce_cutoff_macrostate,apply_cutoff
from peppercornenumerator.input import read_pil 
from peppercornenumerator.enumerator import Enumerator
from peppercornenumerator.reactions import bind21
from peppercornenumerator.enumerator import BI_REACTIONS
import argparse
from natsort import natsorted
from decimal import Decimal,getcontext
import numpy as np 
#setup logger                                      
@pytest.fixture(scope="function")
def configure_logger():
    logger = logging.getLogger("copaths")  # Use the same logger name as in your main code
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter("# %(levelname)s - %(message)s")
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)
    logger.setLevel(logging.DEBUG)  # Set the desired logging level (e.g., DEBUG, INFO, etc.)
    yield
    logger.handlers = []

@pytest.fixture
def input_path():
    return os.path.join(os.path.dirname(__file__), 'test_cocosim_inputs')

def test_enumerate_step():

    pass 


def create_input(filename,occupancies,input_path):
    
    """Creates complexes and structures for further testing
    

        Args():
            filename(str): filename of input file (.pil format)
            occupancies(list): list of occupancies of the corresponding complexes as in the input file 
    """
    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    parser.add_argument("-v","--verbosity",action="count")
    args = args = parser.parse_args()
    args.cutoff = float(args.cutoff)
    file_path = os.path.join(input_path,filename)
    with open(file_path, 'r') as file:
        file_content = file.read()
    complexes, reactions = read_pil(file_content)

    assert len(complexes) == len(occupancies) 
    complex_count = 1 

    all_complexes = {}

    

    for complex,occ in zip(complexes.values(),occupancies):
        
        complex.occupancy = np.float128(occ)
        all_complexes["Id_" + str(complex_count)] = [complex,occ]
        complex_count += 1

    init_cplxs = [x for x in complexes.values() if x.concentration is None or x.concentration[1] != 0]
    name_cplxs = list(complexes.values())
    enum = Enumerator(init_cplxs, reactions, named_complexes=name_cplxs)

    enum.k_slow = 1
    if bind21 in BI_REACTIONS:
        BI_REACTIONS.remove(bind21)
    enum.k_fast = 20 
    
    enum.enumerate()
    enum.condense()

    print(enum.to_pil(condensed=True,detailed = True))

    resting_complexes = [cplx for cplx in natsorted(enum.resting_complexes)] 
    
    transient_complexes = [cplx for cplx in natsorted(enum.transient_complexes)]


    return enum,resting_complexes,transient_complexes,all_complexes,complex_count,args

def test_map_transient_states_1(configure_logger,input_path):
    #2 transient to 1 resting state

    
    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args = create_input("map_transient_1.txt",[0.5,0.5],input_path)

    print(transient_complexes)
    print(resting_complexes)
    old_all_complexes = all_complexes.copy()
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,args)
    
    for id,complex in all_complexes.items():
        print(id,complex,complex[0].occupancy)

    
    new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
    
    assert len(new_entries) == 1 

    for complex in new_entries.values(): 
        assert round(complex[0].occupancy,6) == 1 

     

def test_map_transient_states_2(configure_logger,input_path):
    #2 transient to 1 resting state
    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args = create_input("map_transient_2.txt",[1],input_path)
    print(transient_complexes)
    print(resting_complexes)
    old_all_complexes = all_complexes.copy()
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,args)
    
    for id,complex in all_complexes.items():
        print(id,complex,complex[0].occupancy)

    
    new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
    
    assert len(new_entries) == 1 

    for complex in new_entries.values(): 
        assert complex[0].occupancy == 1 


def test_map_transient_states_3(configure_logger,input_path):
    #2 transient to 1 resting state
    enum,resting_complexes,transient_complexes,all_complexes, complex_count ,args = create_input("map_transient_3.txt",[0.3,0.3,0.3],input_path)

    print(transient_complexes)
    print(resting_complexes)
    old_all_complexes = all_complexes.copy()
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,args)
    
    for id,complex in all_complexes.items():
        print(id,complex,complex[0].occupancy)

    
    new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
    
    assert len(new_entries) == 5 
    
    solution_occs = [0.1280320558,0.0006771507,0.1280320558,0.6331525165,0.0101062211]
    
    print(sum(solution_occs))
    for complex,occ in zip(new_entries.values(),solution_occs): 
        assert float(round(complex[0].occupancy,10)) in solution_occs
        solution_occs.remove(float(round(complex[0].occupancy,10)))

def test_map_transient_states_4(configure_logger,input_path):
    #2 transient to 1 resting state
    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args = create_input("map_transient_4.txt",[0.2,0.2,0.2,0.2,0.2],input_path)

    print(transient_complexes)
    print(resting_complexes)
    old_all_complexes = all_complexes.copy()
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,args)
    
    for id,complex in all_complexes.items():
        print(id,complex,complex[0].occupancy)

    
    new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
    
    assert len(new_entries) == 1 
    solution_occs = [0.2,0.2,0.2,0.2,0.4]
    for complex,occ in zip(new_entries.values(),solution_occs): 
        assert float(complex[0].occupancy) in solution_occs
        solution_occs.remove(float(complex[0].occupancy))

def test_calc_macrostate_oc_1(configure_logger,input_path):


    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args = create_input("calc_macrostate1.txt",[0.5,0.5],input_path)
    old_all_complexes = all_complexes.copy()

    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,args)


    calc_macro_pop(enum,all_complexes,resting_complexes,args)

    print(all_complexes)
    for key,complex in all_complexes.items():
        print(complex)
    
    new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
    
    assert len(new_entries) == 1 

    for complex in new_entries.values(): 
        assert round(complex[0].occupancy,6) == 1

    pass 



def test_calc_macrostate_oc_2(configure_logger,input_path):

    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args = create_input("calc_macrostate2.txt",[1],input_path)
    old_all_complexes = all_complexes.copy()

    
    print('\n\n\n\n complex count ',complex_count)
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,args)

    print("\n\nResting Complexes before\n")
    for complex in resting_complexes:
        print(complex)
    calc_macro_pop(enum,all_complexes,resting_complexes,args)
    print("\n\nResting Complexes after\n")
    for complex in resting_complexes:
        print(complex)


    print(all_complexes)
    for key,complex in all_complexes.items():
        print(complex)
    
    
    assert len(resting_complexes) == 1 

    for complex in resting_complexes: 
        assert complex.occupancy == 1 

    

def test_calc_macrostate_oc_3(configure_logger,input_path):

    
    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args = create_input("calc_macrostate3.txt",[0.25,0.25,0.25,0.25],input_path)
    old_all_complexes = all_complexes.copy()
    print('\n\n\n\n complex count ',complex_count)
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,args)


    print("\n\nResting Complexes before\n")
    for complex in resting_complexes:
        print(complex,complex.occupancy)
    print("\n\nMAcrocalc")
    calc_macro_pop(enum,all_complexes,resting_complexes,args)
    print("\n\nResting Complexes after\n")
    for complex in resting_complexes:
        print(complex,complex.occupancy)

    
    for key,complex in all_complexes.items():
        print(complex,complex)
    
    
    assert len(resting_complexes) == 4 

    solution_occs = [0.3740109438,0.0019781123,0.3740109438,0.25]
    for complex,occ in zip(resting_complexes,solution_occs): 
        assert float(round(complex.occupancy,10)) in solution_occs
        solution_occs.remove(float(round(complex.occupancy,10)))

def test_enforce_cutoff_complex(input_path):

    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v","--verbosity",action="count")

    args  = parser.parse_args()
    args.cutoff = float(args.cutoff)
    args.condensed = True


    
    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args = create_input("calc_macrostate3.txt",[np.float128(0.2),np.float128(0.3),np.float128(0.25),np.float128(0.25)],input_path)
    
    
    for complex1,complex2 in zip(resting_complexes,transient_complexes):
        
        try:
            complex1.occupancy = round(complex1.occupancy,15)
            complex2.occupancy = round(complex2.occupancy,15)
        except:
            pass

    parameters = {"k_slow": None,'k_fast': None, "cutoff": 0.22,"d_length":None}

    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,args)
    macrostates = enum._resting_macrostates
    calc_macro_pop(enum,all_complexes,resting_complexes,args)



    below_t_count_all = 0
    below_t_count_resting = 0 
   
        
    apply_cutoff(enum,parameters,all_complexes)

    for complex in all_complexes.values():
        print(complex)
        if complex[1] == 0: 
            below_t_count_all += 1
    for complex in resting_complexes:
        if complex.occupancy == 0: 
                below_t_count_resting += 1

    
    assert below_t_count_all == below_t_count_resting
    
    assert below_t_count_all == 1

   

def test_enforce_cutoff_macrostate(input_path):
 
    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v","--verbosity",action="count")

    args  = parser.parse_args()
    args.cutoff = float(args.cutoff)
    args.condensed = True


    
    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args = create_input("calc_macrostate3.txt",[np.float128(0.25),np.float128(0.3),np.float128(0.2),np.float128(0.25)],input_path)
    
    
    for complex1,complex2 in zip(resting_complexes,transient_complexes):
        
        try:
            complex1.occupancy = round(complex1.occupancy,15)
            complex2.occupancy = round(complex2.occupancy,15)
        except:
            pass

    parameters = {"k_slow": None,'k_fast': None, "cutoff": 0.22,"d_length":None}

    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,args)
    macrostates = enum._resting_macrostates
    calc_macro_pop(enum,all_complexes,resting_complexes,args)



    below_t_count_all = 0
    below_t_count_resting = 0 
   
        
    apply_cutoff(enum,parameters,all_complexes)

    calc_macro_pop(enum,all_complexes,resting_complexes,args)
    occ_sum_all = 0
    for complex in all_complexes.values():
        if complex[0] in resting_complexes:
            print(complex)
            occ_sum_all += complex[1]
            if complex[1] == 0: 
                below_t_count_all += 1
            


    occ_sum_resting = 0
    for complex in resting_complexes:
        print(complex,complex.occupancy)
        occ_sum_resting += complex.occupancy
        if complex.occupancy == 0: 
                below_t_count_resting += 1
                

    
    assert below_t_count_all == below_t_count_resting, f"Below t count all {below_t_count_all} == {below_t_count_resting} below t count resting "
    
    assert occ_sum_all == occ_sum_resting 

    assert abs(1 - occ_sum_all) <= 1e-10



   





def test_run_sim():
    


    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v","--verbosity",action="count")

    args = args = parser.parse_args()
    args.cutoff = float(args.cutoff)
    args.condensed = True

    # create dictionary for domain lengths 
    d_length = {}

    for domain in 'a b a*'.split():
        if domain[0] == "L" or domain[0] == "S":
            d_length[domain] = 8
        else: 
            d_length[domain] = 3 

    parameters = {"k_slow": args.k_slow, 'k_fast': args.k_fast,"condensed":args.condensed, "cutoff": args.cutoff,'complexes': {},"d_length":d_length}

    simulated_structures = run_sim('a b a*',parameters,args)
    
    print(simulated_structures)

    steps = [[] for _ in simulated_structures]

    for x,step in enumerate(simulated_structures):
        for complex in step:
            steps[x].append(complex.kernel_string)


    set1 = {tuple(set(sublist)) for sublist in steps}
    set2 = {tuple(set(sublist)) for sublist in [['a b'],['a b a*','a( b )']]}
    
    assert set1 == set2

    pass 

def simulate_condensed_reactions():

    pass 





if __name__ == '__main__':

    pytest.main()
"""pytest moduel to test cocosim
"""
import pytest
import logging, os, re, subprocess
from unittest.mock import patch
from cocopaths.cocosim import simulate_system,map_transient_states,calc_macro_pop,enumerate_step,run_sim, enforce_cutoff_macrostate,apply_cutoff
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
def configure_logger(caplog):
    logger = logging.getLogger("cocosim")  # Use the same logger name as in your main code
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
    args, unknown = parser.parse_known_args()
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

    enum.k_slow = 0.001
    if bind21 in BI_REACTIONS:
        BI_REACTIONS.remove(bind21)
    enum.k_fast = 20 
    
    enum.enumerate()
    enum.condense()

    print(enum.to_pil(condensed=True,detailed = True))

    resting_complexes = [cplx for cplx in natsorted(enum.resting_complexes)] 
    
    transient_complexes = [cplx for cplx in natsorted(enum.transient_complexes)]
    
    


    parameters = {"k_slow": enum.k_slow,'k_fast': enum.k_fast, "cutoff": 0.0001,"d_length":None,"d_seq":None}


    return enum,resting_complexes,transient_complexes,all_complexes,complex_count,args,parameters

def test_map_transient_states_1(configure_logger,input_path):
    #2 transient to 1 resting state

    
    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args,parameters = create_input("map_transient_1.txt",[0.5,0.5],input_path)

    print(transient_complexes)
    print(resting_complexes)
    old_all_complexes = all_complexes.copy()
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)
    
    for id,complex in all_complexes.items():
        print(id,complex,complex[0].occupancy)

    
    new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
    
    assert len(new_entries) == 1 

    for complex in new_entries.values(): 
        assert round(complex[0].occupancy,6) == 1 

     

def test_map_transient_states_2(configure_logger,input_path):
    #2 transient to 1 resting state
    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args,parameters = create_input("map_transient_2.txt",[1],input_path)
    print(transient_complexes)
    print(resting_complexes)
    old_all_complexes = all_complexes.copy()
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)
    
    for id,complex in all_complexes.items():
        print(id,complex,complex[0].occupancy)

    
    new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
    
    assert len(new_entries) == 1 

    for complex in new_entries.values(): 
        assert complex[0].occupancy == 1 


def test_map_transient_states_3(configure_logger,input_path):
    #2 transient to 1 resting state
    enum,resting_complexes,transient_complexes,all_complexes, complex_count ,args,parameters = create_input("map_transient_3.txt",[0.3,0.3,0.3],input_path)

    print(transient_complexes)
    print(resting_complexes)
    old_all_complexes = all_complexes.copy()
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)
    
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
    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args,parameters = create_input("map_transient_4.txt",[0.2,0.2,0.2,0.2,0.2],input_path)

    print(transient_complexes)
    print(resting_complexes)
    old_all_complexes = all_complexes.copy()
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)
    
    for id,complex in all_complexes.items():
        print(id,complex,complex[0].occupancy)

    
    new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
    
    assert len(new_entries) == 1 
    solution_occs = [0.2,0.2,0.2,0.2,0.4]
    for complex,occ in zip(new_entries.values(),solution_occs): 
        assert float(complex[0].occupancy) in solution_occs
        solution_occs.remove(float(complex[0].occupancy))

def test_calc_macrostate_oc_1(configure_logger,input_path):


    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args,parameters = create_input("calc_macrostate1.txt",[0.5,0.5],input_path)
    old_all_complexes = all_complexes.copy()

    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)


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

    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args,parameters = create_input("calc_macrostate2.txt",[1],input_path)
    old_all_complexes = all_complexes.copy()

    
    print('\n\n\n\n complex count ',complex_count)
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)

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

    
    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args,parameters = create_input("calc_macrostate3.txt",[0.25,0.25,0.25,0.25],input_path)
    old_all_complexes = all_complexes.copy()
    print('\n\n\n\n complex count ',complex_count)
    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)


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





def test_enforce_cutoff_transient_1(input_path):
 
    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v","--verbosity",action="count")

    args, unknown = parser.parse_known_args()    
    args.cutoff = float(args.cutoff)
    args.condensed = True

    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args,parameters = create_input("calc_macrostate3.txt",[np.float128(0.25),np.float128(0.3),np.float128(0.2),np.float128(0.25)],input_path)
    
    
    for complex1,complex2 in zip(resting_complexes,transient_complexes):
        
        try:
            complex1.occupancy = round(complex1.occupancy,15)
            complex2.occupancy = round(complex2.occupancy,15)
        except:
            pass

    parameters = {"k_slow": None,'k_fast': None, "cutoff": 0.22,"d_length":None,}

    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)
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
   
   

def test_enforce_cutoff_transient_2(input_path):
 
    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v","--verbosity",action="count")

    args, unknown = parser.parse_known_args()
    args.cutoff = float(args.cutoff)
    args.condensed = True
    d_seq = "a* b* L0* c* d* S0 c L0 b S1 L0* S2 e d c L0 b a"

    d_seq_split = d_seq.split()

    d_length = {}

    for domain in d_seq.split():
        if domain[0] == "L" or domain[0] == "S":
            d_length[domain] = 8
        else: 
            d_length[domain] = 3 

    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args,parameters = create_input("macro_cutoff_2.txt",[np.float128(1)],input_path)
    
    
    for complex1,complex2 in zip(resting_complexes,transient_complexes):
        
        try:
            complex1.occupancy = round(complex1.occupancy,15)
            complex2.occupancy = round(complex2.occupancy,15)
        except:
            pass
    
    parameters = {"k_slow": 0.0001,'k_fast': 20, "cutoff": 0.05,"d_length":d_length}

    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)
    macrostates = enum._resting_macrostates
    calc_macro_pop(enum,all_complexes,resting_complexes,args)

    resting_complexes = simulate_system(enum,parameters,resting_complexes,all_complexes,parameters["d_length"][d_seq_split[17]])


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
    
    assert abs(occ_sum_all - occ_sum_resting) < 1e-10  
    
    assert below_t_count_resting == 0

    assert abs(1 - occ_sum_all) <= 1e-10

def test_enforce_cutoff_transient_3(input_path):
 
    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v","--verbosity",action="count")

    args, unknown = parser.parse_known_args()
    args.cutoff = float(args.cutoff)
    args.condensed = True
    d_seq = "a* b* L0* c* d* S0 c L0 b S1 L0* S2 e d c L0 b a"

    d_seq_split = d_seq.split()

    d_length = {}

    for domain in d_seq.split():
        if domain[0] == "L" or domain[0] == "S":
            d_length[domain] = 8
        else: 
            d_length[domain] = 3 

    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args,parameters = create_input("macro_cutoff_3.txt",[np.float128(1/3),np.float128(1/3),np.float128(1/3)],input_path)
    
    
    for complex1,complex2 in zip(resting_complexes,transient_complexes):
        
        try:
            complex1.occupancy = round(complex1.occupancy,15)
            complex2.occupancy = round(complex2.occupancy,15)
        except:
            pass
    
    parameters = {"k_slow": 0.0001,'k_fast': 20, "cutoff": 0.31,"d_length":d_length}

    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)
    macrostates = enum._resting_macrostates
    calc_macro_pop(enum,all_complexes,resting_complexes,args)

    resting_complexes = simulate_system(enum,parameters,resting_complexes,all_complexes,parameters["d_length"][d_seq_split[17]])


    below_t_count_all = 0
    below_t_count_resting = 0 
   
        
    apply_cutoff(enum,parameters,all_complexes)
    #print(enum)
    #calc_macro_pop(enum,all_complexes,resting_complexes,args)


    

    resting_complexes = set(enum._resting_complexes)
    occ_sum_all = 0
    for complex in all_complexes.values():
        if complex[0] in resting_complexes:
            print(complex)
            occ_sum_all += complex[1]
            if complex[1] == 0: 
                below_t_count_all += 1
            
    cut_macros = 0
    for macro in enum._resting_macrostates:
        if macro.occupancy == 0:
            cut_macros += 1 

    occ_sum_resting = 0
    for complex in resting_complexes:
        print(complex,complex.occupancy)
        occ_sum_resting += complex.occupancy
        if complex.occupancy == 0: 
                below_t_count_resting += 1
                
    print("\nMacrostates at the end of the test:")
    for macro in enum._resting_macrostates:
        print(macro,macro.occupancy)

    assert below_t_count_all == below_t_count_resting, f"Below t count all {below_t_count_all} == {below_t_count_resting} below t count resting "
    
    assert abs(occ_sum_all - occ_sum_resting) < 1e-10 
    
    assert cut_macros == 0, f"Macrocomplexes under the threshold {cut_macros}"

    assert abs(1 - occ_sum_all) <= 1e-10

    #resulting occupancies

    resulting_occs = [np.float64(0.96688541022548714067),np.float64(0.033114589774512491545)]
    resting_occs = [x.occupancy for x in enum._resting_complexes]
    resulting_occs.sort()
    resting_occs.sort()
    assert np.allclose(resulting_occs,resting_occs,atol = 1e-10),f"{resulting_occs = },{resting_occs = }"
        
 

def test_enforce_cutoff_transient_4(input_path):
 
    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v","--verbosity",action="count")

    args, unknown = parser.parse_known_args()
    args.cutoff = float(args.cutoff)
    args.condensed = True

    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args,parameters = create_input("cutoff_transient_4.txt",[np.float128(1/37) for _ in range(37)],input_path)
    
    
    for complex1,complex2 in zip(resting_complexes,transient_complexes):
        
        try:
            complex1.occupancy = round(complex1.occupancy,15)
            complex2.occupancy = round(complex2.occupancy,15)
        except:
            pass

    parameters = {"k_slow": None,'k_fast': None, "cutoff": 0.22,"d_length":None,}

    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)
    macrostates = enum._resting_macrostates
    calc_macro_pop(enum,all_complexes,resting_complexes,args)



    below_t_count_all = 0
    below_t_count_resting = 0 
   
        
    apply_cutoff(enum,parameters,all_complexes)

    resting_complexes = list(set(enum._resting_complexes))

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
                
    for reaction in enum.condensation._condensed_reactions:
        print(reaction)
    
    assert below_t_count_all == below_t_count_resting, f"Below t count all {below_t_count_all} == {below_t_count_resting} below t count resting "
    
    assert abs(occ_sum_all - occ_sum_resting) < 1e10 

    assert abs(1 - occ_sum_all) <= 1e-10
    
    resulting_occs = [0.038382816257001896757, 0.038770099035167028878, 0.0006144555402063971569, 8.7879528744050681395e-05, 0.0052310375663540544844, 0.005617991432653917495, 0.3519356958555787332, 0.00019380584501519420874, 6.4610613557798348055e-05, 0.0005378926658468193442, 0.00053600789885841550687, 0.0334075844267117457, 0.03348812881877702489, 0.03324272610860402426, 0.033325155267657873342, 0.00021855457948455538756, 0.004431760052938281573, 0.0037285041001922549505, 0.13884299678243543597, 0.13781190896046257759, 0.0006873918813152298702, 0.13884299678243551357]
    resting_occs = [x.occupancy for x in enum._resting_complexes]

    resulting_occs.sort()
    resting_occs.sort()


    assert np.allclose(resulting_occs,resting_occs,atol = 1e-10)
        


def test_enforce_v_transient_from_dseq(input_path,configure_logger):


    configure_logger

    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v","--verbosity",action="count")
    parser.add_argument("-l", "--logic", action="store_true", default=False,help="Visualizes Logic domain pairings. Best used when analyzing cocopaths generated sequences. (default = False)")

    args, unknown = parser.parse_known_args()
    args.cutoff = float(args.cutoff)
    args.condensed = True
    d_seq = "L0*  S0 a b L0 c d S1 c* L0* b* S2  L0  S3 d* c* L0* b* a* S4  L0  S5"
    d_length = {}

    for domain in d_seq.split():
        if domain[0] == "L":
            d_length[domain] = 12

        elif domain[0] == 'S':
            d_length[domain] = round(int(domain[1]) * 1.5)

        else: 
            d_length[domain] = 3


    parameters = {"k_slow": 0.001,'k_fast': 20, "cutoff": 0.05,"d_length":d_length,"logic":True}

    
    simulated_structures = run_sim(d_seq, parameters)

def test_wierd_outcome(input_path):
 
    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v","--verbosity",action="count")

    args, unknown = parser.parse_known_args()
    args.cutoff = float(args.cutoff)
    args.condensed = True

    enum,resting_complexes,transient_complexes,all_complexes, complex_count,args,parameters = create_input("cutoff_transient_1.txt",[np.float128(1/6) for _ in range(6)],input_path)
    
    
    for complex1,complex2 in zip(resting_complexes,transient_complexes):
        
        try:
            complex1.occupancy = round(complex1.occupancy,15)
            complex2.occupancy = round(complex2.occupancy,15)
        except:
            pass

    parameters = {"k_slow": None,'k_fast': None, "cutoff": 0.22,"d_length":None,"logic":True}

    map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,parameters)
    macrostates = enum._resting_macrostates
    calc_macro_pop(enum,all_complexes,resting_complexes,args)



    below_t_count_all = 0
    below_t_count_resting = 0 
   
        
    apply_cutoff(enum,parameters,all_complexes)

    print(" \n Done with cutoff")
    resting_complexes = list(set(enum._resting_complexes))


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
                
    for reaction in enum.condensation._condensed_reactions:
        print(reaction)
    
    assert below_t_count_all == below_t_count_resting, f"Below t count all {below_t_count_all} == {below_t_count_resting} below t count resting "
    
    assert abs(occ_sum_all - occ_sum_resting) < 1e-10 

    assert abs(1 - occ_sum_all) <= 1e-10
   
   


def test_run_sim():
    


    parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
    parser.add_argument("-cutoff","--cutoff", action= "store",default=float('-inf'), help="Cutoff value at which structures won't get accepted (default: -inf)")
    parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=0.0001)
    parser.add_argument("--k-fast", type=float, help="Specify k-fast. Determines the cutoffpoint for fast reactions.", default=20)
    parser.add_argument("-v","--verbosity",action="count")
    parser.add_argument("-l", "--logic", action="store_true", default=False,help="Visualizes Logic domain pairings. Best used when analyzing cocopaths generated sequences. (default = False)")

    args, unknown = parser.parse_known_args()
    args.cutoff = float(args.cutoff)
    args.condensed = True

    # create dictionary for domain lengths 
    d_length = {}

    for domain in 'a b a*'.split():
        if domain[0] == "L" or domain[0] == "S":
            d_length[domain] = 8
        else: 
            d_length[domain] = 3 

    parameters = {"k_slow": args.k_slow, 'k_fast': args.k_fast,"condensed":args.condensed, "cutoff": args.cutoff,'complexes': {},"d_length":d_length,"logic":True}

    simulated_structures = run_sim('a b a*',parameters)
    
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
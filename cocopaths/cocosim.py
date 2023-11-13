"""Compiler for cotranscriptional folding simulator

Simulates cotranscriptional folding of a domain level sequence.

uses Function from Peppercornenumerator to find reaction types/ possible products of reactions
"""
from peppercornenumerator import peppercorn 
from dsdobjects.objectio import read_pil as dsd_read_pil
from peppercornenumerator.input import read_pil
from peppercornenumerator import Enumerator 
from peppercornenumerator.enumerator import BI_REACTIONS
from peppercornenumerator.reactions import bind21
from natsort import natsorted
from utils import cv_db2kernel


def run_sim(d_seq):
    

    print("hello")
    
    
    #M0 S0 a M0 b S1 b M0 a S2

    global used_structure_names
    used_structure_names = 0

    d_seq_split = d_seq.split()
    complexes, reactions = input_parsing(d_seq_split[0:2],[cv_db2kernel(d_seq_split[0:2],"..")])



    final_structures = [[] for x in range(len(d_seq_split))]



    current_complexes = []
    current_complexes.append(enumerate_step(complexes=complexes,reactions=reactions))

    print("\n\n Loopy\n\n\n")
    for step in range(2,len(d_seq_split)):
        print("\n\n\nNext Step\n\n")
        print(f"Step: {step} Current Complexes: {current_complexes}")
        next_complexes = []



        for structure in current_complexes:
        
            next_struct = structure[0] + " " +  d_seq_split[step]

            print(f"Next Sturcture {next_struct}")
            next_complexes.append(next_struct)

        complexes, reactions = input_parsing(d_seq_split[0:step + 1],next_complexes)
        
        resting_complexes = enumerate_step(complexes=complexes,reactions=reactions)

        print("resulting resting comples",resting_complexes)

        final_structures[step].append(resting_complexes)
        current_complexes = [resting_complexes]

    print("Final Resting Complexes",resting_complexes)

    print("\n\n\nEnumerations after each step")
    for x in final_structures:
        
        if x and x[0][0][-2] == "S":
            print(x[0])

    




    

    
    print("hello")

def enumerate_step(complexes, reactions):

    k_slow = 0.1
    
    if bind21 in BI_REACTIONS:
            BI_REACTIONS.remove(bind21)


    init_cplxs = [x for x in complexes.values() if x.concentration is None or x.concentration[1] != 0]
    name_cplxs = list(complexes.values())
    enum = Enumerator(init_cplxs, reactions, named_complexes = name_cplxs)

    if k_slow: 
        enum.k_slow = k_slow



    #Start to enumerate
    enum.enumerate()

    output = enum.to_pil()
    
    resting_complexes = [cplx.kernel_string for cplx in natsorted(enum.resting_complexes)] 
    

    print(f"Resting complexes: {resting_complexes}")


    return resting_complexes





def input_parsing(d_seq,structures):
    global used_structure_names
    print("\n\n\n\n")
    print("Start Input Parsing")
    print("Structures", structures,"d_seq:",d_seq)
    toehold_length = 3 

    logic_length = 8

    space_length = 6



    unique_domains = set([domain.replace('*', '') for domain in d_seq])

    system_input = f""
    print(f"Unqiue domains: {unique_domains}")
    for unique_domain in unique_domains:

        if unique_domain.islower():
            system_input += f"length {unique_domain} = {toehold_length} \n"

        elif unique_domain.startswith("L"):
            system_input += f"length {unique_domain} = {logic_length} \n"


        elif unique_domain.startswith("S"):
            system_input += f"length {unique_domain} = {space_length} \n"

    system_input += "\n"


    for x,structure in enumerate(structures):
        system_input += f"E{used_structure_names} = {structure}\n"
        used_structure_names += 1 



    print(system_input)
    complexes, reactions = read_pil(system_input)

    print(f"Complexes: {complexes}")

    print(f"reactions: {reactions}")
    print("End Input Parsing")
    return complexes, reactions



def main():
    
    #Input parsing


    #simulation run 


    print("Nothing here come back later")


if __name__ == "__main__":

    #main()
    run_sim("L0*  S0 a L0 b S1 c* b* L0* a* d* S2 e d a L0 b c f S3 f* c* b* L0* a* d* e* S4")
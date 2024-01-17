#Script to analysis capabilities of cocopahts and cocosim 
#Should analyse all possible legal folding paths and do some statistics of the outcomes




from path_generator import generate_path
from cocopaths.cocopaths import build_graph
from cocopaths.cocosim import run_sim, write_output
from peppercornenumerator.objects import clear_memory
from cocopaths.cocosim import kernel_to_dot_bracket,only_logic_domain_struct
import logging 


def analyze_cocosim_output(simulated_structures,afp,d_seq):
    print("\n\nBeginning with analyzing cocosim output")
    print("Afp",afp)
    print("D_Seq",d_seq)
    occupancy_sum = 0
    target_dominant = False


    i = 0
    for step in simulated_structures:
        print("\n\nStep",step)

        occupancies = [complex.occupancy for complex in step]
        for complex in step: 
            kernel_string = kernel_to_dot_bracket(complex.kernel_string)
            db_struct = (only_logic_domain_struct(d_seq.split(),kernel_string))
            print(db_struct)
            print(afp,i)
            #part for avg occupancy of target
            try:
                if db_struct == afp[i] and complex.kernel_string.split()[-1][0] == "S":
                    print(complex,complex.kernel_string,complex.occupancy)
                    occupancy_sum += complex.occupancy
                    

                    #target dominant in ensemble? 
                    if max(occupancies) == complex.occupancy:
                        target_dominant = True
                    else:
                        target_dominant = False
                    
                    i += 1   
            except:
                continue
        
    print(occupancy_sum)

    avg_occupancy = occupancy_sum/(i)

    print(avg_occupancy)
    
    return f"0\t{afp}\t{target_dominant}\t{avg_occupancy}\n"




    




def main():


    #Generating all folding paths 
    
    n = 6

    print(f"Generating all possible Folding Paths upto a length of {n}")

    folding_paths = generate_path(n)


    #Calculating all domain level sequences for the folding paths

    domain_sequences = []
    
    real_paths = []
    for i,path in enumerate(folding_paths[-1]):
        print("\n\nCurrent Path",path)

        try:
            afp_graph = build_graph(path)

            domain_sequences.append(" ".join(afp_graph.get_domain_seq()))
            real_paths.append(path)
            print("Worked for : ",afp_graph.get_domain_seq())
        except: 
            
            print("Didn't work for:",path)

    for d_seq in domain_sequences:
        print(d_seq)


    tsv_header = "ID\tAFP\td_seq\tdominating_struct\tavg_occupancy\n"
    for d_seq,fp in zip(domain_sequences,real_paths):
        print("Current domain sequence",d_seq)
        print("Current Folding path",fp)


    #Now simulate the whole thing and evalute the output
        d_length = {}

        for domain in d_seq.split():
            if domain[0] == "L" or domain[0] == "S":
                d_length[domain] = 8
            else: 
                d_length[domain] = 3 


        parameters = {"k_slow": 0.001 , 'k_fast': 20, "cutoff": float('-inf'),"d_length":d_length,"d_seq":d_seq}


        simulated_structures = run_sim(d_seq,parameters)
        print(simulated_structures)


        with open(str(n) + "_steps_out.tsv", "a") as file: 
            new_data = analyze_cocosim_output(simulated_structures,fp,d_seq) 
            if file.tell() == 0:
                file.write(tsv_header)

            file.write(new_data)

        output = ""
        if 'S' in d_seq:
            output += write_output(simulated_structures,d_seq)
        print(output)

        clear_memory()

    print("\n\nend of script\n\n")





if __name__ == "__main__":
    main()
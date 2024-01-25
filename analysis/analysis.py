#Script to analysis capabilities of cocopahts and cocosim 
#Should analyse all possible legal folding paths and do some statistics of the outcomes




from path_generator import generate_path
from cocopaths.cocopaths import build_graph
from cocopaths.cocosim import run_sim, write_output
from cocopaths.utils import is_balanced_structure
from peppercornenumerator.objects import clear_memory
from cocopaths.cocosim import kernel_to_dot_bracket,only_logic_domain_struct
import logging 
import pandas as pd
import matplotlib.pyplot as plt
import os

def analyze_cocosim_output(simulated_structures,afp,d_seq):
    occupancy_sum = 0
    target_dominant = True

    dominant_path = ["."]
    i = 0

    target_occupancies = []
    for step in simulated_structures:
        

        occupancies = [complex.occupancy for complex in step]
        t_occ_append = False
        for complex in step: 
            print(complex,complex.kernel_string)
            kernel_string = kernel_to_dot_bracket(complex.kernel_string)
            db_struct = (only_logic_domain_struct(d_seq.split(),kernel_string))
            
            #part for avg occupancy of target
            if complex.kernel_string.split()[-1][0] == "S":
                if max(occupancies) == complex.occupancy and dominant_path[-1] != db_struct:
                    dominant_path.append(db_struct)
                    if  not is_balanced_structure(db_struct):
                        print("huh:",db_struct)
                        exit()
                try:
                    print("\n\n",db_struct,",",afp[i],"\n\n")
                    if db_struct == afp[i]:
                        occupancy_sum += complex.occupancy
                        target_occupancies.append(round(complex.occupancy,6))
                        t_occ_append = True
                        #target dominant in ensemble? 
                        if max(occupancies) != complex.occupancy:
                            target_dominant = False
                        i += 1
                      
                except:
                    continue

        if not t_occ_append and step[0].kernel_string.split()[-1][0] == "S":
            target_occupancies.append(0)
            i += 1
            target_dominant = False

    final_occupancy = None
    for complex in simulated_structures[-1]:
        kernel_string = kernel_to_dot_bracket(complex.kernel_string)
        db_struct = (only_logic_domain_struct(d_seq.split(),kernel_string))
        if db_struct == afp[-1]:
            final_occupancy = round(complex.occupancy,4)
        
    print("occ sum /len",occupancy_sum,len(target_occupancies))
    avg_occupancy = occupancy_sum/len(target_occupancies)

    if target_dominant:
        print("Target dominant",target_dominant)
        print("final occupancy")
        assert final_occupancy != None 
    
    return f"0\t{afp}\t{target_dominant}\t{round(avg_occupancy,4):8}\t{final_occupancy}\t{dominant_path}\t{target_occupancies}\t{d_seq}\n"



def statistical_analysis(folder_path,tsv,filename):


    if not os.path.exists("results/results.tsv"):
        print("\n\nMaking direktor")
        os.makedirs("results")
        with open("results/results.tsv","a") as file: 
            header = "Folder_Name     \tn_steps\ttotal\ttrue\tfalse\tavg_occ\tavg_occ_last\t0.0-0.1\t0.1-0.2\t0.2-0.3\t0.3-0.4\t0.4-0.5\t0.5-0.6\t0.6-0.7\t0.8-0.9\t1-1.1\n"
            file.write(header)
            print("wrote")


    print("start with analysis")
    # Read the TSV file into a DataFrame
    df = pd.read_csv(tsv, sep='\t', comment = "#")
    n = tsv[0]
    # Display the DataFrame
    print("Original DataFrame:")
    print(df.head())
    # Basic analysis
    print("\nBasic Analysis:")
    print("Total entries:", len(df))
    print("\nColumn Names:")
    print(df.columns)

    avg_occ = df['avg_occupancy'].mean()
    print("Average occupancy:", avg_occ)
    avg_last_step_occ = df['last_step_occ'].mean()

    # Count the total occurrences of True and False in dominant_fp
    total_counts = df['dominating_struct'].value_counts()
    dom_true = total_counts.get(True,0)
    dom_false = total_counts.get(False,0)
    print(total_counts,type(total_counts))
    print(dom_false)
    print(dom_true)
    print(df.columns)
    
    print(tsv)

    print("Total counts of True and False in dominant_fp:")
    print(total_counts)




    bin_edges = [i / 10 for i in range(0, 12)]  # [0.1, 0.2, ..., 1.0]

    # Create bins using pd.cut
    df['bin'] = pd.cut(df['avg_occupancy'], bins=bin_edges, right=False)

    # Count the number of rows in each bin
    counts = df['bin'].value_counts().sort_index()

    bins_dict = dict(zip(counts.index.astype(str), counts.values))
    print(bins_dict)

    with open("results/results.tsv","a") as file: 
        file.write(f"{folder_path:15}\t{filename[0]}\t{dom_false + dom_true:8}\t{dom_true:8}\t{dom_false:8}\t{avg_occ:8.4}\t{avg_last_step_occ:8.4}\t")
        for key,value in bins_dict.items():
            file.write(f"{value}\t")
        file.write("\n")



    print(counts)
    # Plot the results
    counts.plot(kind='bar', width=0.8, align='center')
    plt.xlabel('Threshold Bins')
    plt.ylabel('Number of Rows')
    plt.title(f'Number of Rows in Threshold Bins for {n} Steps')
    plt.xticks(rotation=45)
    #plt.show()


def find_next_folder_number(base_folder):
    existing_folders = [name for name in os.listdir('.') if os.path.isdir(name) and name.endswith('_run') and name.split('_')[0].isdigit()]
    print(existing_folders)
    print(os.path.abspath('.'))
    if not existing_folders:
        return 1
    existing_numbers = [int(name.split('_')[0]) for name in existing_folders]
    return max(existing_numbers) + 1

def write_paths_to_file(file_path, paths):
    print("beginn writing")
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    with open(file_path, 'w') as file:
        for step in paths:
            file.write(str(step) + '\n')

def read_paths_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            paths = [eval(line.strip()) for line in file.readlines()]
        return paths
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.")
        return []



def get_data(n,current_folder):

    #Generating all folding paths 
    

    print(f"Generating all possible Folding Paths upto a length of {n}")

    

    #Calculating all domain level sequences for the folding paths
    fp_locs = "./folding_paths"

    
    fp_locs = "./folding_paths" + "/" + str(n) + "_FPs.txt"

    if os.path.exists(fp_locs):
        folding_paths = read_paths_from_file(fp_locs)
        print(len(folding_paths),"1")
        
    else:    
        folding_paths = generate_path(n)
        print(len(folding_paths),"2")
    
    print(folding_paths)
    domain_sequences = []
    
    real_paths = []
    for i,path in enumerate(folding_paths):
        

        #print("\n\nCurrent Path",path)

        try:
            afp_graph = build_graph(path)

            domain_sequences.append(" ".join(afp_graph.get_domain_seq()))
            real_paths.append(path)
            #print("Worked for : ",afp_graph.get_domain_seq())
        except: 
            
            print("Didn't work for:",path)
    print("real paths",real_paths)
    if not os.path.exists(fp_locs):
        write_paths_to_file(fp_locs, real_paths)
    

    tsv_header = "ID\tAFP\tdominating_struct\tavg_occupancy\tlast_step_occ\tstep_occs\tdominant_fp\td_seq\n"
    for d_seq,fp in zip(domain_sequences,real_paths):
        


    #Now simulate the whole thing and evalute the output
        d_length = {}

        for domain in d_seq.split():
            if domain[0] == "L":
                d_length[domain] = 12

            elif domain[0] == 'S':
                d_length[domain] = 3

            else: 
                d_length[domain] = 3


        parameters = {"k_slow": 0.001 , 'k_fast': 20, "cutoff": 0.05,"d_length":d_length,"d_seq":d_seq}

    
        simulated_structures = run_sim(d_seq,parameters)
    
        with open(os.path.join(current_folder, f"{n}_steps_out.tsv"), "a") as file:
            new_data = analyze_cocosim_output(simulated_structures,fp,d_seq) 
            if file.tell() == 0:
                file.write("#" + str(parameters) + "\n")
                file.write(tsv_header)

            file.write(new_data)
    
        output = ""
        if 'S' in d_seq:
            output += write_output(simulated_structures,d_seq)
        print(output)
    
        
        #Clear memory of Peppercorn objects
        clear_memory()
           

    print("\n\nend of script\n\n")

def fill_data(n,tsv_file): 

    #Generating all folding paths 
    

    print(f"Generating all possible Folding Paths upto a length of {n}")



    #Calculating all domain level sequences for the folding paths
    fp_locs = "./folding_paths"

    with open(tsv_file,"r") as file:
        lines = file.readlines()

        start_line = len(lines) - 2 

    fp_locs = "./folding_paths" + "/" + str(n) + "_FPs.txt"

    if os.path.exists(fp_locs):
        folding_paths = read_paths_from_file(fp_locs)
    else:
        print("\n\nNew paths must be generated this can lead to duplicate fps and missing fps in the analysis\n\n")    
        folding_paths = generate_path(n)

    


    domain_sequences = []
    
    real_paths = []
    for i,path in enumerate(folding_paths):
        print("\n\nCurrent Path",path)

        try:
            afp_graph = build_graph(path)

            domain_sequences.append(" ".join(afp_graph.get_domain_seq()))
            real_paths.append(path)
            print("Worked for : ",afp_graph.get_domain_seq())
        except: 
            
            print("Didn't work for:",path)
    print("Start line:",start_line)
    print(domain_sequences)
    domain_sequences = domain_sequences[start_line:]
    real_paths = real_paths[start_line:]

    print(domain_sequences)
    for d_seq in domain_sequences:
        print(d_seq)

    write_paths_to_file



    tsv_header = "ID\tAFP\tdominating_struct\tavg_occupancy\tlast_step_occ\tstep_occs\tdominant_fp\td_seq\n"
    for d_seq,fp in zip(domain_sequences,real_paths):
        print("Current domain sequence")
        print("Current Folding path")


    #Now simulate the whole thing and evalute the output
        d_length = {}

        for domain in d_seq.split():
            if domain[0] == "L":
                d_length[domain] = 12

            elif domain[0] == 'S':
                d_length[domain] = 0

            else: 
                d_length[domain] = 3


        parameters = {"k_slow": 0.001 , 'k_fast': 20, "cutoff": 0.05,"d_length":d_length,"d_seq":d_seq}

        simulated_structures = run_sim(d_seq,parameters)
        print(simulated_structures)
        print("file:",tsv_file)
        with open(tsv_file, "a") as _file:

            new_data = analyze_cocosim_output(simulated_structures,fp,d_seq) 
            print("\n\n\nnew data",new_data)
            if _file.tell() == 0:
                _file.write("#" + str(parameters) + "\n")
                _file.write(tsv_header)

            _file.write(new_data)
    
        output = ""
        if 'S' in d_seq:
            output += write_output(simulated_structures,d_seq)
        print(output)

       

        
        
        

        clear_memory()
    

    print("\n\nend of script\n\n")


def main():

    base_folder = "run"
    next_folder_number = find_next_folder_number(base_folder)

    current_folder = f'{next_folder_number}_{base_folder}'
    
    os.makedirs(current_folder, exist_ok=True)
    for i in range(2,7):
        get_data(i,current_folder)
        
    print("S_3")


    #uncomment to fill up file if segfault happended 
    #check if parameters match 
    #fill_data(6,"10_run/6_steps_out.tsv")

    print("data is in ",current_folder)

    



if __name__ == "__main__":
    
    
    main()

    #analyze_folder = "Sx3+2_run"
    #for filename in os.listdir(analyze_folder):
    #        print(filename)
    #        if filename.endswith(".tsv") and os.path.isfile(os.path.join(analyze_folder, filename)):
    #            tsv_filepath = os.path.join(analyze_folder, filename)
    #            statistical_analysis(analyze_folder,tsv_filepath,filename)
    #statistical_analysis(folder_path,"7_steps_out.tsv")
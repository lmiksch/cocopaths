#Script to analysis capabilities of cocopahts and cocosim 
#Should analyse all possible legal folding paths and do some statistics of the outcomes




#from path_generator import generate_path
#from cocopaths.cocopaths import build_graph
#from cocopaths.cocosim import run_sim, write_output
#from peppercornenumerator.objects import clear_memory
#from cocopaths.cocosim import kernel_to_dot_bracket,only_logic_domain_struct
import logging 
import pandas as pd
import matplotlib.pyplot as plt

def analyze_cocosim_output(simulated_structures,afp,d_seq):
    print("\n\nBeginning with analyzing cocosim output")
    print("Afp",afp)
    print("D_Seq",d_seq)
    occupancy_sum = 0
    target_dominant = False

    dominant_path = ["."]
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
            if max(occupancies) == complex.occupancy and dominant_path[-1] != db_struct:
                dominant_path.append(db_struct)
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
    
    return f"0\t{afp}\t{target_dominant}\t{round(avg_occupancy,4):8}\t{dominant_path}\t{d_seq}\n"



def statistical_analysis(tsv):

        # Read the TSV file into a DataFrame
    df = pd.read_csv(tsv, sep='\t')

    # Display the DataFrame
    print("Original DataFrame:")
    print(df.head())
    # Basic analysis
    print("\nBasic Analysis:")
    print("Total entries:", len(df))
    
    print("Average occupancy:", df['avg_occupancy'].mean())

    # Count the total occurrences of True and False in dominant_fp
    total_counts = df['dominating_struct'].value_counts()

    print("Total counts of True and False in dominant_fp:")
    print(total_counts)

    """ Cumulative binning 
    # Add additional analysis based on your specific requirements
    threshold_range = [i / 10 for i in range(1, 10)] 
    # Store the counts for each threshold
    counts = []

    # Iterate over threshold values
    for threshold in threshold_range:
        filtered_rows = df[df['avg_occupancy'] > threshold]
        counts.append(len(filtered_rows))

    # Plot the results
    plt.bar(threshold_range, counts, width=0.1, align='center')
    plt.xlabel('Threshold')
    plt.ylabel('Entries')
    plt.title('Number of Rows vs. Threshold')
    """

    bin_edges = [i / 10 for i in range(1, 12)]  # [0.1, 0.2, ..., 1.0]

    # Create bins using pd.cut
    df['bin'] = pd.cut(df['avg_occupancy'], bins=bin_edges, right=False)

    # Count the number of rows in each bin
    counts = df['bin'].value_counts().sort_index()

    # Plot the results
    counts.plot(kind='bar', width=0.8, align='center')
    plt.xlabel('Threshold Bins')
    plt.ylabel('Number of Rows')
    plt.title('Number of Rows in Threshold Bins')
    plt.xticks(rotation=45)
    plt.show()




    




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


    tsv_header = "ID\tAFP\tdominating_struct\tavg_occupancy\tdominant_fp\td_seq\n"
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

        try:
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
        except KeyboardInterrupt:
            raise SystemExit
        except:
            pass


        clear_memory()
    
    

    print("\n\nend of script\n\n")





if __name__ == "__main__":
    #main()
    statistical_analysis("6_steps_out.tsv")
import argparse
from cocopaths.cocopath import build_graph
from cocopaths.cocodesign import rna_design,domain_path_to_nt_path,afp_to_domainfp,constrained_efe,objective_function,call_findpath,extend_domain_seq,score_sequence
from cocopaths.utils import only_logic_domain_struct
import RNA
import numpy as np 
import subprocess,os,re
import matplotlib.pyplot as plt 
import seaborn as sns
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from statistics import mean

#DrTransformer Parser




def coco_suite(folding_path,output,steps):
    """
    Uses the functions from CocoPaths to create a nucleotide sequence based on an abstract folding path.

    Args: 
        folding_path(list): List of folding steps each entry corresponds to a step in the folding path. 

    Returns: 
        nt_sequence(string): nucleotide sequence
        score(float)
    """
    #CocoPath

    afp_graph = build_graph(folding_path)

    domain_seq = ' '.join(afp_graph.get_domain_seq())

    print(f"Resulting Domain Level sequence: {domain_seq}")


    #creation of parameters using the default parameters 
    
    d_length = {}
    
    for domain in domain_seq.split():
        if domain[0] == "L":
            d_length[domain] = 8
        elif domain[0] == 'S':
            d_length[domain] =  round(int(domain[1]) * 4)  
        else: 
            d_length[domain] = 3 
                
    
    parameters = {"k_slow": 0.00001,'k_fast': 20, "cutoff": 0.05,"d_length":d_length,"d_seq":domain_seq,"logic":True,'steps':steps}

    
    #CocoDesign
    #for now force each folding path 

    domain_fp = afp_to_domainfp(folding_path,domain_seq)

    print(f"\n\nBegin RNA design\n")

    ext_fp = domain_path_to_nt_path(domain_fp,domain_seq,parameters)

    nt_sequence,score = rna_design(domain_seq,ext_fp,parameters)


    detailed_scores = score_sequence(nt_sequence,domain_seq,parameters,folding_path)[1]


    print(detailed_scores)
    #write output 
    with open(output + ".tsv","a") as file:
        file.write(f"\n{folding_path = }\t{parameters['steps'] = }\n{nt_sequence =}\n")
        file.write(detailed_scores)
    

    return nt_sequence,domain_seq,domain_fp,parameters,score


def evaluate_nt_fp(nt_sequence,domain_fp,domain_seq,parameters,output = "test"):
    """Function to evaluate the folding path of a nt sequence given an abstract foldingpath via a simulation with DrTransformer. 

    Args: 

        nt_sequence(str): nucleotide sequence
        domain_fp(list): the domain level folding path of the nucleotide sequence
        domain_seq(str): domain sequence from which the nucleotide sequence was designed separated by a space.
        parameteres(dict): parameters from which the nt seq was designed
        output(txt): output file which then gets printed

    Returns:
        populations(list): populations for each step 
    """

    print(nt_sequence)

    #create folder 
    # Check if the folder 'drt_out' exists, if not create it
    folder_name = output + '_drt_out'
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"Folder '{folder_name}' created.")
    else:
        print(f"Folder '{folder_name}' already exists.")

    # Regular expression pattern to match files with numbers
    pattern = re.compile(r'(\d+)_test\.drf')

    # Initialize current index
    current_index = 1

    # Check for numbered files and determine the highest number
    existing_files = os.listdir(folder_name)
    highest_number = 0

    for file in existing_files:
        match = pattern.match(file)
        if match:
            number = int(match.group(1))
            if number > highest_number:
                highest_number = number

    # Set current index
    if highest_number > 0:
        current_index = highest_number + 1
    

    print(f"\n{current_index = }  {nt_sequence = }")

    command = "DrTransformer --name " + str(current_index) + "_test"

    result = subprocess.run(command, input=nt_sequence, text=True, shell=True, cwd=folder_name, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    
    
    with open(output + "_drt_out/"+  str(current_index) + '_test.drf') as f:
        dr_out= [line.rstrip().split() for line in f]

    extended_fp = domain_path_to_nt_path(domain_fp,domain_seq,parameters)

    print(f"{extended_fp = }")


    populations = [[] for _ in range(len(extended_fp))]
    seen_drt_out = [[] for _ in range(len(extended_fp))]
    i = 0
    for x in range(1,len(dr_out)):
        if len(dr_out[x][3]) > len(extended_fp[i]):
            populations[i].append(0)
            i += 1
        if drt_match(dr_out[x][3],extended_fp[i]) and dr_out[x][3] not in seen_drt_out[i]:
            populations[i].append(float(dr_out[x][2]))
            seen_drt_out[i].append(dr_out[x][3])
        if i == len(extended_fp):
            break
    if i != len(extended_fp): # if last structure is not matching add 0 to populations only occurs if optimization fails
            populations[i].append(0)
    print("Populations", "Structure")
    for x in range(len(extended_fp)):
        print(sum(populations[x]),extended_fp[x])
    populations = [sum(sub_list) for sub_list in populations]
    if output != "test":
        with open(output + ".tsv","a") as file:
            file.write(f'\nCurrent Index: {current_index}\n')
            for x in range(len(extended_fp)):
                file.write(f'{populations[x]}\t{extended_fp[x]}\n')

    
    return populations

def drt_match(drt_out_struct,fp_struct):
    if len(drt_out_struct) != len(fp_struct)    :
        return False
    for drt_char,fp_char in zip(drt_out_struct,fp_struct):
        if fp_char == ".":
            continue
        elif drt_char != fp_char:
            return False 
    return True



def score_sim_nt_seq(nt_seq,afp,parameters):

    afp_graph = build_graph(afp)
    
    d_seq = ' '.join(afp_graph.get_domain_seq())

    score,output = score_sequence(nt_seq,domain_seq,parameters,folding_path_2)


    domain_fp = afp_to_domainfp(afp,d_seq)
    print(output)

    evaluate_nt_fp(nt_seq,domain_fp,d_seq,parameters)

def get_default_parameters():
    domain_seq_fp2 = 'a* L0* L1 L2 L3 L4 L1* L2* b* S0 g h i j k l m n o p q r s t u v w c L0 d S1 e b L0 a f S2 f* a* L0* b* e* S3 g* d* L0* c* h* S4 h c L0 d g S5'

    domain_seq = domain_seq_fp2
    d_length = {}
    
    for domain in domain_seq.split():
        if domain[0] == "L":
            d_length[domain] = 8
        elif domain[0] == 'S':
            d_length[domain] =  round(int(domain[1]) * 4)  
        else: 
            d_length[domain] = 3 
    
    parameters = {"k_slow": 0.00001,'k_fast': 20, "cutoff": 0.05,"d_length":d_length,"d_seq":domain_seq,"logic":True,'steps':None}

    return parameters

def plot_dist(dist,output_name):

    plt.figure(figsize=(8, 6))
    sns.histplot(x=dist, bins=int(1/0.05))
    plt.title(" ")
    plt.xlim(0, 1)
    plt.xlabel('min(Occupancy)')
    plt.ylabel('Counts')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(output_name)

def correlation_analysis(all_pop,all_score,obj_fun,out_folder):

    data = pd.DataFrame({
    'Occupancies': all_pop,
    'Scores': all_score
    })  
    print(f'{all_pop = }')

    print(f'{all_score = }')
    pearson_corr, pearson_p_value = pearsonr(data['Occupancies'], data['Scores'])

    spearman_corr, spearman_p_value = spearmanr(data['Occupancies'], data['Scores'])

    with open(out_folder + ".tsv","a") as f:
        f.write(f"\nPearson correlation coefficient: {pearson_corr}\n")
        f.write(f"Pearson p-value: {pearson_p_value}\n")
        f.write(f"Spearman correlation coefficient: {spearman_corr}\n")
        f.write(f"Spearman p-value: {spearman_p_value}\n")

    
    if not os.path.exists('corr_analysis.tsv'):
        with open('corr_analysis.tsv','a') as file:
            file.write(f'Objective Function                           \tcorr\tp_value\tmean_score\tmean_pop\n')

    with open('corr_analysis.tsv','a') as file: 
        file.write(f'{obj_fun}\t{spearman_corr}\t{spearman_p_value}\t{mean(all_score)}\t{mean(all_pop)}\n')

    plt.figure(figsize=(8, 6))
    plt.scatter(data['Occupancies'], data['Scores'])
    plt.xlabel('Occupancies')
    plt.ylabel('Scores')
    plt.title('Scatter plot of Occupancies vs. Scores')
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(out_folder + "/corr_plot")




def main(folding_path):

    folder_name = "stat_10"
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"Folder '{folder_name}' created.")
    else:
        print(f"Folder '{folder_name}' already exists.")
        raise SystemExit('Folder already exists')

    output = folder_name + "/stat_test"


    steps_5k = [5000 for x in range(0,15)]
    steps_10k = [10000 for x in range(0,10)]
    steps_20k = [20000 for x in range(0,5)]
    min_pop_5k = []
    min_pop_10k = []
    min_pop_20k = []
    score_5k = []
    score_10k = []
    score_20k = []


    for steps in steps_5k:
        nt_sequence,domain_sequence,domain_fp,parameters,score = coco_suite(folding_path,output,steps) 

        #DrTransformer 
        populations = evaluate_nt_fp(nt_sequence,domain_fp,domain_sequence,parameters,output)
        min_pop_5k.append(min(populations))
        score_5k.append(score)
    
    for steps in steps_10k:
        nt_sequence,domain_sequence,domain_fp,parameters,score = coco_suite(folding_path,output,steps) 


        #DrTransformer 
        populations = evaluate_nt_fp(nt_sequence,domain_fp,domain_sequence,parameters,output)
        min_pop_10k.append(min(populations))
        score_10k.append(score)


    for steps in steps_20k:
        nt_sequence,domain_sequence,domain_fp,parameters,score = coco_suite(folding_path,output,steps) 


        #DrTransformer 

        populations = evaluate_nt_fp(nt_sequence,domain_fp,domain_sequence,parameters,output)
        min_pop_20k.append((min(populations)))
        score_20k.append(score)

    
    
    plot_dist(min_pop_5k,folder_name + "/5kplot")
    plot_dist(min_pop_10k,folder_name + "/10kplot")
    plot_dist(min_pop_20k,folder_name + "/20kplot")

    all_pop = min_pop_5k + min_pop_10k + min_pop_20k
    all_score = score_5k + score_10k + score_20k
    plot_dist(all_pop,folder_name + "/allplot")
    
    obj_fun = objective_function(0,0,0,0,0,0)[1]

    correlation_analysis(all_pop,all_score,obj_fun,folder_name)
    
   

  
if __name__ == "__main__":
    
    print("stat_3")
    folding_path_1 = [".","()",".()","()()",".()()","()()()"]

    folding_path_2 = ['.', '()', '(.)', '()()', '.(())', '()()()']

    folding_path_3 = ['.', '()', '.()', '()()', '()().', '()(..)']
    test_path = [".","()",".()"]
    main(folding_path_3)




    #Extends the domain seq to nt length to check in DrForna
    #domain_seq_fp2 = 'a* L0* b* S0 c L0 d S1 e b L0 a f S2 f* a* L0* b* e* S3 g* d* L0* c* h* S4 h c L0 d g S5'



    domain_seq_fp3 = 'L0*  S0 a L0 b S1 c* d* b* L0* a* e* f* S2 e a L0 b d S3  L1*  S4 f e a L0 b d c S5'

    parameters = get_default_parameters()

    nt_seq = 'CUCGGACUAAAACAUCGUUUAGUUUUCAUUUUCAUUGUUUUAGUUCGAGGUUUCUUUUCUGACCUCGAGCUAAAGCAAUGUUCUUCUACUACGCUUGAGAACUGGGCGAAACGCAAUCAUUUCUCAUAGUUUCGCUCAGUUUUCGGGCCAAACCCACUCAAACAUCGU'

    #score_sim_nt_seq(nt_seq,folding_path_2,parameters)



    print(extend_domain_seq(domain_seq_fp3,parameters))

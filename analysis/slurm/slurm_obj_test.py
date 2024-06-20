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
import pandas as pd
import ast
from filelock import FileLock


def coco_suite(folding_path,output,steps,parameters):
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
    
    parameters['steps'] = steps
    
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
    print(f'{d_seq = }')

    score,output = score_sequence(nt_seq,d_seq,parameters,afp)


    domain_fp = afp_to_domainfp(afp,d_seq)
    print(output)

    evaluate_nt_fp(nt_seq,domain_fp,d_seq,parameters)

def get_default_parameters():
    domain_seq_fp2 = 'a* L0* L1 L2 L3 L4 L5 L5 L6 L7 L8 L9 L10 L1* L2* L3* L4* L5* L6* L7* L8* L9* L11 b* S0 g h i j k l m n o p q r s t u v w c L0 d S1 e b L0 a f S2 f* a* L0* b* e* S3 g* d* L0* c* h* S4 h c L0 d g S5'

    domain_seq = domain_seq_fp2
    d_length = {}
    
    for domain in domain_seq.split():
        if domain[0] == "L":
            d_length[domain] = 8
        elif domain[0] == 'S':
            d_length[domain] = 3 + round(int(domain[1]) * 3)  
        else: 
            d_length[domain] = 4
    
    parameters = {"k_slow": 0.00001,'k_fast': 20, "cutoff": 0.05,"d_length":d_length,"d_seq":domain_seq,"logic":True}

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
    '''calculataes correlation coeffiecient
    Returns: 
        spearman_coor(float): Spearman correlation coeffiecient
        spearman_p_value(float): calculated p-value
    '''

    data = pd.DataFrame({
    'Occupancies': all_pop,
    'Scores': all_score
    })  
    print(f'{all_pop = }')

    print(f'{all_score = }')
    if len(all_pop) <= 1 or len(all_score) <= 1:
        return -1,-1

    pearson_corr, pearson_p_value = pearsonr(data['Occupancies'], data['Scores'])

    spearman_corr, spearman_p_value = spearmanr(data['Occupancies'], data['Scores'])

    if out_folder:
        with open(out_folder + ".tsv","a") as f:
            f.write(f"\nPearson correlation coefficient: {pearson_corr}\n")
            f.write(f"Pearson p-value: {pearson_p_value}\n")
            f.write(f"Spearman correlation coefficient: {spearman_corr}\n")
            f.write(f"Spearman p-value: {spearman_p_value}\n")
        

        if not os.path.exists('corr_analysis.tsv'):
            with open('corr_analysis.tsv','a') as file:
                file.write(f'Objective Function                           \tcorr\tp_value\tmean_score\tmean_pop\n')

        with open('corr_analysis.tsv','a') as file: 
            file.write(f'{obj_fun}\t{spearman_corr}\t{spearman_p_value}\t{mean(all_score)}\t{mean(all_pop)}\t{out_folder}\n')

        plt.figure(figsize=(8, 6))
        plt.scatter(data['Occupancies'], data['Scores'])
        plt.xlabel('Occupancies')
        plt.ylabel('Scores')
        plt.title('Scatter plot of Occupancies vs. Scores')
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(out_folder + "/corr_plot")

    return spearman_corr,spearman_p_value

def analyze_all_fps(file_path,start_index):
    '''Goes through all possible folding paths as calculate in analysis and tries to generate a nt sequence until it formes the input folding path. 
    
    Args():
        analysis_result(tsv): file from analysis.py 
        start_index(int): index where the analysis should start

    Output():
        results(tsv): tsv file with results(fp_index,sucess,tries,avg_occ)
        fp_X(folder): X stands for the index at which the folding path occurs
            drt_out(folder): output of the drtransformer for each fp 
            X_result.tsv(): output of the nucleotide generation

    '''

    input_tsv = pd.read_csv(file_path, sep='\t', comment='#',skiprows = 1 )

    print(input_tsv.head())
    print(input_tsv.info())

    #convert AFP entries into list
    input_tsv['AFP'] = input_tsv['AFP'].apply(ast.literal_eval)

    parameters = get_default_parameters()
    obj_fun = objective_function(0,0,0,0,0,0)[1]
    

    out_folder = '/home/mescalin/miksch/Documents/cocopaths/analysis/obj_test/test/'

    output_file = '/home/mescalin/miksch/Documents/cocopaths/analysis/obj_test/test/all_analysis_test_out.tsv'

    if not os.path.exists(output_file):
        with open(output_file,'a') as out_file:
            out_file.write(f'Id\tAFP                    \tnt_success\ttries\tavg_occ\tdomain_sucess\tcorr_coeffiecient\tp_value\n')


    for index,row in input_tsv.iloc[start_index:].iterrows():
        print(row['AFP'],row['dominating_struct'])

        curr_path = row['AFP']

        result_pop = 0

        #create folder for each run 
        afp_folder = out_folder + str(index) + '_afp'
        os.makedirs(afp_folder)

        total_score = []
        total_pop = []

        score_list = []
        occ_list = []

        if row['dominating_struct']:
            num_tries = 20
        else:
            num_tries = 10


        
        mean_pop = 0
        result_pop = []

        for tries in range(0,num_tries):   
                output = afp_folder + '/' + str(index)
                try:
                    nt_sequence,domain_sequence,domain_fp,parameters,score = coco_suite(curr_path,output,3000,parameters) 
                    populations = evaluate_nt_fp(nt_sequence,domain_fp,domain_sequence,parameters,output)
                    cur_result_pop = min(populations)
                    result_pop.append(min(populations))

                    score_list.append(score)
                    occ_list.append(result_pop)



                    if cur_result_pop > 0.5:#checks if nt-design was sucessfull
                        mean_pop = mean(occ_list)
                        break

                except KeyboardInterrupt:
                    print("KeyboardInterrupt caught. Exiting gracefully.")
                    result_pop.append(-1)

                except:
                    total_score += score_list
                    total_pop += occ_list
                    with FileLock(output_file + ".lock"):
                        with open(output_file,'a') as output_f:
                            output_f.write(f'-')
                    result_pop.append(-1)

                    break
        total_score += score_list
        total_pop += occ_list         
        sp_corr, p_value = correlation_analysis(total_pop,total_score,obj_fun,False)
        
        with FileLock(output_file + ".lock"):
            with open(output_file,'a') as output_f:
                output_f.write(f'{index}\t{curr_path}\t{max(result_pop):.6f}\t{tries}\t{mean(result_pop):.6f}\t{row["dominating_struct"]}\t{sp_corr:.6f}\t{p_value:.6f}\n')

    #correlation analysis for whole data set

    end_sp_corr, end_p_value = correlation_analysis(total_pop,total_score,obj_fun,False)
    with FileLock(output_file + ".lock"):
        with open(output_file,'a') as output_f:
            output_f.write(f'\nTotal values:\t{end_sp_corr = }\t{end_p_value = }')




def analyze_single_fp(file_path, index, output_file,out_folder):
    '''Analyzes a single folding path given by the index.'''
    
    input_tsv = pd.read_csv(file_path, sep='\t', comment='#', skiprows=1)
    input_tsv['AFP'] = input_tsv['AFP'].apply(ast.literal_eval)
    
    parameters = get_default_parameters()
    obj_fun = objective_function(0, 0, 0, 0, 0, 0)[1]
    
    with FileLock(output_file + ".lock"):
        if not os.path.exists(output_file):
            with open(output_file, 'a') as out_file:
                out_file.write(f'Id\tAFP\tnt_success\ttries\tavg_occ\tdomain_success\tcorr_coeffiecient\tp_value\n')

    row = input_tsv.iloc[index]
    curr_path = row['AFP']
    result_pop = 0

    # Create folder for each run 
    afp_folder = out_folder + str(index) + '_afp'
    os.makedirs(afp_folder, exist_ok=True)

    total_score = []
    total_pop = []

    score_list = []
    occ_list = []

    if row['dominating_struct']:
        num_tries = 20
    else:
        num_tries = 10


    for tries in range(0, num_tries):   
        output = afp_folder + '/' + str(index)
        try:
            nt_sequence, domain_sequence, domain_fp, parameters, score = coco_suite(curr_path, output, 3000, parameters) 
            populations = evaluate_nt_fp(nt_sequence, domain_fp, domain_sequence, parameters, output)
            result_pop = min(populations)

            score_list.append(score)
            occ_list.append(result_pop)

            if result_pop > 0.5:  # checks if nt-design was successful
                break

        except Exception as e:
            print(f"Error processing {index}: {e}")
            total_score += score_list
            total_pop += occ_list
            with FileLock(output_file + ".lock"):
                with open(output_file, 'a') as output_f:
                    output_f.write(f'-')
            break

    total_score += score_list
    total_pop += occ_list         
    sp_corr, p_value = correlation_analysis(total_pop, total_score, obj_fun, False)

    with FileLock(output_file + ".lock"):
        with open(output_file, 'a') as output_f:
            output_f.write(f'{index}\t{curr_path}\t{max(occ_list):.6f}\t{tries}\t{mean(occ_list):.6f}\t{row["dominating_struct"]}\t{sp_corr:.6f}\t{p_value:.6f}\n')



def main(folding_path):

    folder_name = "6.3"
    if not os.path.exists(folder_name):
        os.makedirs(folder_name)
        print(f"Folder '{folder_name}' created.")
    else:
        print(f"Folder '{folder_name}' already exists.")
        raise SystemExit('Folder already exists')

    output = folder_name + "/stat_test"


    steps_1k = [1000 for x in range(0,30)]
    steps_5k = [5000 for x in range(0,15)]
    steps_10k = [10000 for x in range(0,5)]
    min_pop_1k = []
    min_pop_5k = []
    min_pop_10k = []
    score_1k = []
    score_5k = []
    score_10k = []

    parameters = get_default_parameters()
    for steps in steps_1k:
        nt_sequence,domain_sequence,domain_fp,parameters,score = coco_suite(folding_path,output,steps,parameters) 

        #DrTransformer 
        populations = evaluate_nt_fp(nt_sequence,domain_fp,domain_sequence,parameters,output)
        min_pop_1k.append(min(populations))
        score_1k.append(score)
    
    for steps in steps_5k:
        nt_sequence,domain_sequence,domain_fp,parameters,score = coco_suite(folding_path,output,steps,parameters) 


        #DrTransformer 
        populations = evaluate_nt_fp(nt_sequence,domain_fp,domain_sequence,parameters,output)
        min_pop_5k.append(min(populations))
        score_5k.append(score)


    for steps in steps_10k:
        nt_sequence,domain_sequence,domain_fp,parameters,score = coco_suite(folding_path,output,steps,parameters) 


        #DrTransformer 

        populations = evaluate_nt_fp(nt_sequence,domain_fp,domain_sequence,parameters,output)
        min_pop_10k.append((min(populations)))
        score_10k.append(score)

    
    
    plot_dist(min_pop_1k,folder_name + "/1kplot")
    plot_dist(min_pop_5k,folder_name + "/5kplot")
    plot_dist(min_pop_10k,folder_name + "/10kplot")

    all_pop = min_pop_1k + min_pop_5k + min_pop_10k
    all_score = score_1k + score_5k + score_10k
    plot_dist(all_pop,folder_name + "/allplot")
    
    obj_fun = objective_function(0,0,0,0,0,0)[1]

    correlation_analysis(all_pop,all_score,obj_fun,folder_name)
    
   


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Analyze a single folding path.')
    parser.add_argument('--file_path', type=str, required=True, help='Path to the input TSV file.')
    parser.add_argument('--index', type=int, required=True, help='Index for analysis.')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output TSV file.')

    args = parser.parse_args()

    
    analyze_single_fp(args.file_path, args.index, args.output_file,"test/")
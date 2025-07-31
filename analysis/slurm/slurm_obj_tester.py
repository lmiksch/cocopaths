import argparse
from cocopaths.cocopath import translate_acfp
from cocopaths.cocodesign import rna_design,domain_path_to_nt_path,acfp_to_domainfp,constrained_efe,objective_function,call_findpath,extend_domain_seq,score_sequence,toehold_structures
from cocopaths.utils import only_logic_domain_struct
import RNA
import numpy as np 
import subprocess,os,re
import matplotlib.pyplot as plt 
#import seaborn as sns
import pandas as pd
from scipy.stats import pearsonr, spearmanr
from statistics import mean
import pandas as pd
import ast
import traceback
from filelock import FileLock
import string
from check_drt_output import drt_match

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

    
    domain_list = translate_acfp(folding_path)
    domain_seq = ' '.join(domain_list)

    print(f"Resulting Domain Level sequence: {domain_seq}")
    #creation of parameters using the default parameters 
    
    parameters['steps'] = steps
    parameters['hard_constraint'] = False
    
    #CocoDesign
    #for now force each folding path 

    domain_fp = acfp_to_domainfp(folding_path,domain_seq)

    print(f"\n\nBegin RNA design\n")

    ext_fp = domain_path_to_nt_path(domain_fp,domain_seq,parameters)

    
    #toehold_domain_fp = toehold_structures(domain_fp,domain_seq,folding_path)

    nt_sequence,score = rna_design(domain_seq,ext_fp,parameters,folding_path,domain_fp)


    detailed_scores = score_sequence(nt_sequence,domain_seq,parameters,folding_path,domain_fp)[1]

    ext_fp = [x.replace('x','.') for x in ext_fp]
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
    
    #revert back to see wether or not success 
    extended_fp = [_.replace('x','.') for _ in extended_fp]

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

def score_sim_nt_seq(nt_seq,acfp,parameters):

    domain_list = translate_acfp(path)
    d_seq = ' '.join(domain_list)

    print(f'{d_seq = }')

    score,output = score_sequence(nt_seq,d_seq,parameters,acfp)


    domain_fp = acfp_to_domainfp(acfp,d_seq)
    print(output)

    evaluate_nt_fp(nt_seq,domain_fp,d_seq,parameters)

def get_default_parameters():

    # Combine lowercase and uppercase letters
    letters = string.ascii_lowercase 

    # Extend each letter with a '*' without space between the letter and '*'
    extended_letters = ' '.join(f"{letter}*" for letter in letters)

    # Generate L0 to L100 and S0 to S100
    l_values = ' '.join(f"L{num}" for num in range(101))
    l_star_values = ' '.join(f"L{num}*" for num in range(101))

    s_values = ' '.join(f"Z{num}" for num in range(101))

    # Combine the original letters, extended letters, and L/S sequences
    result = ' '.join(letters) + ' ' + extended_letters + ' ' + l_values + ' ' + s_values + ' ' + l_star_values

    domain_seq_fp2 = 'a b c d e f g h i j k l m n o p q r s t u v w x y  a* b* c* d* e* f* g* h* i* j* k* l* m* n* o* p* q* r* s* t* u* v* w* x* y* z* Z0 Z1 Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9 Z10 Z11 L12 L13 L14 L15 L16 L17 L18 L19 L20 L21 L22 L23 L24 S0 S1 S2 S3 S4 S5 S6 S7 S8 S9 S10 S11 S12 S13 S14 S15 S16 S17 S18 S19 S20 S21 S22 S23 S24 S25'

    domain_seq = result
    d_length = {}
    
    for domain in domain_seq.split():
        if domain[0] == "L":
            d_length[domain] = 8
        elif domain[0] == 'Z':
            d_length[domain] = 3 + round(int(domain[1]) * 3)  
        else: 
            d_length[domain] = 3
    
    parameters = {"k_slow": 0.00001,'k_fast': 20, "cutoff": 0.05,"d_length":d_length,"d_seq":domain_seq,"logic":True, "gap_length" : 3}

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

    #convert acfp entries into list
    input_tsv['acfp'] = input_tsv['acfp'].apply(ast.literal_eval)

    parameters = get_default_parameters()
    obj_fun = objective_function(0,0,0,0,0,0)[1]
    

    out_folder = '/home/mescalin/miksch/Documents/cocopaths/analysis/obj_test/test/'

    output_file = '/home/mescalin/miksch/Documents/cocopaths/analysis/obj_test/test/all_analysis_test_out.tsv'

    if not os.path.exists(output_file):
        with open(output_file,'a') as out_file:
            out_file.write(f'Id\tacfp                    \tnt_success\ttries\tavg_occ\tdomain_sucess\tcorr_coeffiecient\tp_value\n')


    for index,row in input_tsv.iloc[start_index:].iterrows():
        print(row['acfp'],row['domain_success'])

        curr_path = row['acfp']

        result_pop = 0

        #create folder for each run 
        acfp_folder = out_folder + str(index) + '_acfp'
        os.makedirs(acfp_folder)

        total_score = []
        total_pop = []

        score_list = []
        occ_list = []

        if row['domain_success']:
            num_tries = 20
        else:
            num_tries = 10


        
        mean_pop = 0
        result_pop = []

        for tries in range(0,num_tries):   
                output = acfp_folder + '/' + str(index)
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
        if not result_pop: 
            result_pop.append(-1)

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
    
    input_tsv = pd.read_csv(file_path, sep='\t', comment='#', skiprows=0)
    print(input_tsv.head())
    input_tsv['acfp'] = input_tsv['acfp'].apply(ast.literal_eval)
    
    parameters = get_default_parameters()
    obj_fun = objective_function(0, 0, 0, 0, 0, 0)[1]
    
    with FileLock(output_file + ".lock"):
        if not os.path.exists(output_file):
            with open(output_file, 'a') as out_file:
                out_file.write(f'Id\tacfp\tnt_success\ttries\tavg_occ\tdomain_success\tcorr_coeffiecient\tp_value\n')

    row = input_tsv.iloc[index]
    curr_path = row['acfp']
    result_pop = 0


    
    # Create folder for each run 
    acfp_folder = out_folder + str(index) + '_acfp'
    os.makedirs(acfp_folder, exist_ok=True)

    total_score = []
    total_pop = []

    score_list = []
    occ_list = []

    if row['domain_success']:
        num_tries = 20
    else:
        num_tries = 20


    for tries in range(0, num_tries):   
        output = acfp_folder + '/' + str(index)
        try:
            nt_sequence, domain_sequence, domain_fp, parameters, score = coco_suite(curr_path, output, 5000, parameters) 
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
    

def eval_nt_seq(nt_sequence,acfp,domain_seq):

    parameters = get_default_parameters()

    domain_fp = acfp_to_domainfp(acfp,domain_seq)


    print(evaluate_nt_fp(nt_sequence,domain_fp,domain_seq,parameters,output = "eval_nt_seq"))


def design_from_acfp_and_domain_seq(acfp, domain_seq, output_file, out_folder, file_path, index=0, tries_per_fp=20):
    """
    Designs an nt-sequence based on a given aCFP and domain sequence (skipping simulation),
    evaluates it with DrTransformer, and writes the result to output_file.
    
    Args:
        acfp (list): abstract cotranscriptional folding path
        domain_seq (str): precomputed domain-level sequence
        output_file (str): path to output .tsv
        out_folder (str): folder to store intermediate outputs
        index (int): ID of the folding path
        tries_per_fp (int): number of design attempts
    """
    from statistics import mean

    parameters = get_default_parameters()
    obj_fun = objective_function(0, 0, 0, 0, 0, 0)[1]

    os.makedirs(out_folder, exist_ok=True)
    acfp_folder = os.path.join(out_folder, f"{index}_acfp")
    os.makedirs(acfp_folder, exist_ok=True)


    # Init logs
    with FileLock(output_file + ".lock"):
        if not os.path.exists(output_file):
            with open(output_file, 'a') as out_file:
                out_file.write('Id\tacfp\tnt_success\ttries\tavg_occ\tdomain_success\n')

    domain_fp = acfp_to_domainfp(acfp, domain_seq)
    total_score = []
    total_pop = []

    score_list = []
    occ_list = []
    
    parameters['steps'] = 2000
    parameters['hard_constraint'] = False

    for tries in range(tries_per_fp):
        output = os.path.join(acfp_folder, f"{index}")
        try:
            ext_fp = domain_path_to_nt_path(domain_fp, domain_seq, parameters)
            nt_sequence, score = rna_design(domain_seq, ext_fp, parameters, acfp, domain_fp)

            populations = evaluate_nt_fp(nt_sequence, domain_fp, domain_seq, parameters, output)
            min_pop = min(populations)

            score_list.append(score)
            occ_list.append(min_pop)

            if min_pop > 0.5:
                break

        except Exception as e:
            print(f"Error in design try {tries}: {e}")
            traceback.print_exc()
            break

    if not occ_list:
        occ_list.append(-1)

    total_score += score_list
    total_pop += occ_list

    input_tsv = pd.read_csv(file_path, sep='\t', comment='#', skiprows=0)
    row = input_tsv.iloc[index]

    with FileLock(output_file + ".lock"):
        with open(output_file, 'a') as output_f:
            output_f.write(f'{index}\t{acfp}\t{max(occ_list):.6f}\t{tries}\t{mean(occ_list):.6f}\t{row["domain_success"]}\n')



if __name__ == "__main__":


    #analyze_single_fp("/home/mescalin/miksch/Documents/cocopaths/analysis/folding_paths/test.txt",8,'example','test/')
    
    
    import argparse
    
    parser = argparse.ArgumentParser(description='Analyze a single folding path.')
    parser.add_argument('--file_path', type=str, required=True, help='Path to the input TSV file.')
    parser.add_argument('--index', type=int, required=True, help='Index for analysis.')
    parser.add_argument('--output_file', type=str, required=True, help='Path to the output TSV file.')
    
    args = parser.parse_args()
    
    
    #analyze_single_fp(args.file_path, args.index, args.output_file,"test_relax/")


    df = pd.read_csv("/home/mescalin/miksch/Documents/data/cocopaths/paper_cocopath/sim_results_len6.tsv", sep="\t", )

    acfp_str = df.loc[args.index, "AFP"]
    acfp = ast.literal_eval(acfp_str)

    out_folder = "test_paper"


    output_file = args.output_file
    index = args.index
    print(f'{args.index = }\n {acfp =}')

    file_path = args.file_path
    
    domain_seq = df.loc[args.index, "d_seq"]
    design_from_acfp_and_domain_seq(acfp, domain_seq, output_file, out_folder, file_path, index,  tries_per_fp=20)
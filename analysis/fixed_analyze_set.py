import re
import string
import os
import pandas as pd

from cocopaths.cocosim import run_sim, write_output
from cocopaths.utils import is_balanced_structure
from peppercornenumerator.objects import clear_memory
from cocopaths.cocosim import kernel_to_dot_bracket, only_logic_domain_struct


def summarize_results(tsv_file: str):
    try:
        df = pd.read_csv(tsv_file, sep='\t', quotechar='"', comment='#', engine='python')
        total = len(df)
        dom_success = df['dominating_struct'].sum()
        domain_success = df['domain_success'].sum()
        both_success = ((df['dominating_struct']) & (df['domain_success'])).sum()

        avg_occ = df['avg_occupancy'].mean()
        avg_last_occ = df['last_step_occ'].mean()

        print(f"\n Summary of {tsv_file}")
        print("-" * 40)
        print(f"Total entries           : {total}")
        #print(f"Target dominant (True)  : {dom_success}")
        print(f"Domain success (True)   : {domain_success}")
        #print(f"Both conditions True    : {both_success}")
        print(f"Success rate            : {domain_success / total:.2%}")
        print(f"Average occupancy       : {avg_occ:.4f}")
        print(f"Average last step occ.  : {avg_last_occ:.4f}")
        print("-" * 40)

    except FileNotFoundError:
        print(f"[ERROR] File '{tsv_file}' not found.")
    except Exception as e:
        print(f"[ERROR] Could not process file: {e}")



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
            if complex.kernel_string.split()[-1][0] == "Z":
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

        if not t_occ_append and step[0].kernel_string.split()[-1][0] == "Z":
            target_occupancies.append(0)
            i += 1
            target_dominant = False

    final_occupancy = None
    for complex in simulated_structures[-1]:
        kernel_string = kernel_to_dot_bracket(complex.kernel_string)
        db_struct = (only_logic_domain_struct(d_seq.split(),kernel_string))

        if db_struct == afp[-1]:
            final_occupancy = round(complex.occupancy,4)
        else:
            print(f'{db_struct = }    {afp[-1] = }')
        
    print("occ sum /len",occupancy_sum,len(target_occupancies))
    avg_occupancy = occupancy_sum/len(target_occupancies)

    if target_dominant:
        print("Target dominant",target_dominant)
        print("final occupancy")
        assert final_occupancy != None 

    domain_success = all(x > 0.5 for x in target_occupancies)
    
    return f"0\t{afp}\t{target_dominant}\t{domain_success}\t{round(avg_occupancy,4):8}\t{final_occupancy}\t{dominant_path}\t{target_occupancies}\t{d_seq}\n"



def transform_sequence(sequence: str) -> str:
    # Replace T<number> → Z<number>
    sequence = re.sub(r'\bT(\d+)', r'Z\1', sequence)

    # Replace sX and sX* with a, b, ..., z
    support_domains = sorted(set(re.findall(r's\d+\*?', sequence)), key=lambda x: int(x[1:].rstrip('*')))
    letters = list(string.ascii_lowercase)

    domain_map = {}
    for i, dom in enumerate(support_domains):
        base = dom.rstrip('*')
        comp = dom.endswith('*')
        if base not in domain_map:
            domain_map[base] = letters[i]
        if comp:
            domain_map[dom] = domain_map[base] + '*'
        else:
            domain_map[dom] = domain_map[base]

    for original, replacement in domain_map.items():
        sequence = sequence.replace(original, replacement)

    return sequence


def simulate_from_sequence_file(input_file: str, output_file: str):
    if not os.path.exists(input_file):
        raise FileNotFoundError(f"Input file '{input_file}' does not exist.")

    with open(input_file, 'r') as f:
        lines = [
            line.strip().split('\t')
            for line in f
            if not line.startswith('#') and line.strip()
        ]

    tsv_header = "ID\tAFP\tdominating_struct\tdomain_success\tavg_occupancy\tlast_step_occ\tstep_occs\tdominant_fp\td_seq\n"
    with open(output_file, 'w') as out:
        out.write(tsv_header)

        for i, (afp, raw_seq) in enumerate(lines):
            d_seq = transform_sequence(raw_seq)
            afp_structs = afp.split()

            d_length = {}
            for domain in d_seq.split():
                if domain[0] == "L":
                    d_length[domain] = 8
                elif domain[0] == "Z":
                    d_length[domain] = 1 + round(int(domain[1]) * 3)
                else:
                    d_length[domain] = 4

            parameters = {
                "k_slow": 0.001,
                "k_fast": 20,
                "cutoff": 0.05,
                "d_length": d_length,
                "d_seq": d_seq,
                "logic": True,
            }

            try:
                sim_result = run_sim(d_seq, parameters)
                result_line = analyze_cocosim_output(sim_result, afp_structs, d_seq)
                out.write(result_line)
                clear_memory()

                if 'Z' in d_seq:
                    _ = write_output(sim_result, d_seq)

            except Exception as e:
                print(f"[ERROR] Simulation failed for entry {i}: {afp}")
                print("  →", e)
                continue

    print(f"\n Simulation complete. Results written to: {output_file}")


if __name__ == "__main__":



    input_tsv = "designable_sequences_len6.txt"
    output_tsv = "sim_results_len6.tsv"
    #simulate_from_sequence_file(input_tsv, output_tsv)


    summarize_results("sim_results_len6.tsv")
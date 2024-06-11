import argparse
from cocopaths.cocodesign import domain_path_to_nt_path, afp_to_domainfp
from peppercornenumerator.input import read_pil


def read_input_pil(args):
    afp = None
    d_length = {}  # Define d_length as an empty dictionary
    pil_input = read_pil(args.pil, True)
    print(pil_input)
    if len(pil_input[0]) == 1:
        input_complex = next(iter(pil_input[0].values()))
        print("Value of the entry:", input_complex)
        print(f"{input_complex._sequence = }  ")
        print(f"{input_complex.kernel_string = }  ")
        
        d_seq = input_complex.kernel_string


        if all(char.isalpha() or char == '*' for char in d_seq):
            raise SystemExit("SystemExit: Only a domain level sequence is accepted. (No structural information)")

        for domain in input_complex._sequence:
            length = domain._length
            name = domain._name
            if name not in d_length:
                d_length[name] = length

        print(f"{d_length = }")
        print(f'\nInput domain level sequence: {d_seq}\n')

        return d_seq, d_length
    else:
        raise ValueError("Invalid pil input format")


def get_drtout_pops(extended_fp,drt_out):
    """
    Parses the output of DrTransformer and checks if certain structures are present and returns their populations. 

    Args: 
        ext_fp(list): list of extended folding path 
        drt_out(filepath): 

    """

    with open(drt_out) as f:
        dr_out= [line.rstrip().split() for line in f]


    populations = []
    i = 0
    for x in range(1,len(dr_out)):
        print('\nDrout\n',dr_out[x][3],'\n',extended_fp[i])
        if len(dr_out[x][3]) > len(extended_fp[i]):
            populations.append("0")
            print('zero',extended_fp[i])
            i += 1
        if dr_out[x][3] == extended_fp[i]:
            populations.append(dr_out[x][2])
            print(dr_out[x][2],extended_fp[i])
            i += 1
            if i == len(extended_fp):
                break

    print("Populations", "Structure")
    for x in range(len(extended_fp)):
        print(populations[x],extended_fp[x])


def main(folding_path):
    parser = argparse.ArgumentParser(prog='DrT_parser', description='Parses the output of DrTransformer and checks if the desired structures at each transcription step are present')
    parser.add_argument("-drt_out", type=str, help="DrTransformer output file you want to analyze", default="NoName.drf")
    parser.add_argument("-pil", type=str, help="domain_seq_generator output file you want to analyze", default=None)

    args = parser.parse_args()

    if args.pil:
        try:
            d_seq, d_length = read_input_pil(args)
            print(f'Read input successfully: {d_seq}, {d_length}')
        except Exception as e:
            print(f"An error occurred: {e}")
    else:
        print("No pil file provided.")
    
    print(d_seq)
    print(d_length)
    domain_fp = afp_to_domainfp(folding_path,d_seq)        
    print('\n\nDseq',d_seq)
    parameters = {"d_length": d_length}

    extended_fp = domain_path_to_nt_path(domain_fp, d_seq, parameters)
    print(f'Extended folding path: {extended_fp}')

    get_drtout_pops(extended_fp,args.drt_out)



if __name__ == "__main__":
    print("begin")
    # Hard folding path 
    # folding_path_1 = [".","()",".()","()()",".()()","()()()"]
    folding_path_2 = ['.', '()', '(.)', '()()', '.(())', '()()()']
    main(folding_path_2)

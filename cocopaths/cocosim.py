"""Compiler for cotranscriptional folding simulator

Simulates cotranscriptional folding of a domain level sequence.

uses Function from Peppercornenumerator to find reaction types/ possible products of reactions
"""
import sys
from peppercornenumerator import peppercorn 
from dsdobjects.objectio import read_pil as dsd_read_pil
from peppercornenumerator.input import read_pil
from peppercornenumerator import Enumerator 
from peppercornenumerator.enumerator import BI_REACTIONS
from peppercornenumerator.reactions import bind21
from natsort import natsorted
from .utils import cv_db2kernel, kernel_to_dot_bracket, only_logic_domain_struct
import argparse
import logging 
from drtransformer.linalg import mx_simulate
import numpy as np 



logger = logging.getLogger('cocosim')
console_handler = logging.StreamHandler()
formatter = logging.Formatter('# %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)




def run_sim(d_seq, parameters):
	"""
	Runs a cotranscriptional simulation of the given d_seq using peppercorn from the peppercornenumerator

	Args:
		d_seq(str): domain level sequence 
		parameters(dict)
	"""    
	used_structure_names = 1


	d_seq_split = d_seq.split()

	# before input parsing define name of complexes 
	
	parameters["complexes"]["I0"] = [cv_db2kernel(d_seq_split[0:2], ".."),1]
	logger.debug(parameters["complexes"])
	complexes, reactions = input_parsing(d_seq_split[0:2], parameters["complexes"])

	final_structures = [[] for x in range(len(d_seq_split))]
	final_populations = [[] for x in range(len(d_seq_split))]
	current_complexes = []

	all_complexes = []


	final_structures[0].append([d_seq_split[0]])
	final_complexes, pops = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters)
	
	#current_complexes.append(first_complexes[0])
	final_structures[1].append([struct[0] for name,struct in parameters["complexes"].items()])

	all_complexes.append(final_complexes)
	#print("final Structures: ", final_structures)


	for step in range(2, len(d_seq_split)):
		logger.debug("\n\n\nNext Step\n\n")
		logger.debug(f"Step: {step} Current Complexes: {parameters}")
		next_complexes = []
		
		#put here add one domain to struct and give new name
		new_complexes = {}
		for name, struct in parameters["complexes"].items():
			#extend structure
			next_struct = struct[0] + " " + d_seq_split[step] 
			#parameters["complexes"][name][0] = next_struct
			logger.debug(f"Next Structure {next_struct}")
			
			#create new name for structure 
			# maybe change name to E to see which complexes were "manmade"
			new_complexes["E" + str(used_structure_names)] = [next_struct, struct[1]] 
			used_structure_names += 1 

		parameters["complexes"] = new_complexes
		logger.debug(f"\n\n\n\n New Parameters\n{parameters}\n\n\n\n")
			

		complexes, reactions = input_parsing(d_seq_split[0:step + 1], parameters["complexes"])
		#logger.debug(f"Following Complexes and reactions are put into enumerator \n {complexes} {reactions}")
		resting_complexes,pops_dict = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters)
		final_populations.append(pops_dict)
		#logger.debug(f"resulting resting complex: {resting_complexes}")

		final_structures[step].append(resting_complexes)
		
		all_complexes.append(resting_complexes)
	logger.debug("\n\n\n________Done with simulation__________\n\n")
	for step in final_populations:
		logger.info(step)
	print(f"\nComplexes: \n")

	for step in all_complexes:
		print("\n")
		for complex in step:
			print(f"Name {complex._name}, {complex.kernel_string}, {complex._concentration[1]}")

	return final_structures

	
	
def write_output(final_structures,d_seq,parameters = None):

	data_output = ""
	

	
	for step in final_structures:
		print(step)

	ts = 1
	data_output += ("\nResting Complexes after each Spacer:\n\n")
	data_output += "Transcription Step | Structure \n"
	for x in final_structures:
		print(x)
		if x and x[0][-1][-2] == "S": 
			for seq in x[0]:
				data_output += f"{ts}	{seq} \n"
			ts += 1 

			data_output += "\n"
	
	
	ts = 1
	data_output += ("\n\nOnly Logic Domain pairings:\n\n")
	data_output += "Transcription Step | Structure \n"
	for x in final_structures:
		if x and x[0][-1][-2] == "S": 
			for seq in x[0]:
				db_struct = (only_logic_domain_struct(d_seq.split(),kernel_to_dot_bracket(seq)))
				data_output += f"{ts}	{db_struct} \n"
			ts += 1
			data_output += "\n"

	return data_output

	


def enumerate_step(complexes, reactions, parameter):
	k_slow = parameter["k_slow"]
	
	if bind21 in BI_REACTIONS:
		BI_REACTIONS.remove(bind21)

	logger.debug(f"\n\n\nBeginning of Enumerate step:\nComplexes:{complexes}\nReactions:{reactions}\nParameters:{parameter}\n\n")
	init_cplxs = [x for x in complexes.values() if x.concentration is None or x.concentration[1] != 0]
	name_cplxs = list(complexes.values())
	enum = Enumerator(init_cplxs, reactions, named_complexes=name_cplxs)

	if k_slow: 
		enum.k_slow = k_slow
	
	# Start to enumerate
	enum.enumerate()
	
	
	condensed = parameter["condensed"]
	if condensed:
		rate_matrices, complex_indices = enum.condense()

	

	
	
	#Calculating Populations and updating the pop ins parameters
	logger.debug("\n\n\n\n\n\n Beginn rate matrixes\n\n_________________________________________________\n\n\n")
	pops_dict = {}
	for matrix,cplx_indices in zip(rate_matrices,complex_indices):
		
		if matrix is not None and len(matrix) > 1:
			


			#check to adjust parameter complexes
			new_complexes = {}
			for complex, index in cplx_indices.items():
				if complex._name not in parameter["complexes"]:
					
					new_complexes[complex._name] = [complex.kernel_string,None]
			for complex in new_complexes:
				parameter["complexes"][complex] = new_complexes[complex]

			keys_to_delete = []

			
			keys_to_keep = [key._name for key in cplx_indices if key._name  in parameter["complexes"]]


			for key in parameter["complexes"]:
				if key not in keys_to_keep:
					keys_to_delete.append(key)

			for key in keys_to_delete:
   				del parameter["complexes"][key]

			for complex_key in parameter["complexes"]:
				parameter["complexes"][complex_key][1] = 1 / len(parameter["complexes"])
							

			logger.debug(f"\n\n\nAfter Addition of new Structures\n\n{[parameter['complexes']]}")
			t1 = 0.1 
			t8 = 100000

			lin_times = np.array([t1], dtype='float128')
			log_times = np.array(np.logspace(np.log10(t1), np.log10(t8), num=10, dtype='float128'))
			log_times = np.delete(log_times, 0)
			
			log_times = []
			#print("lin times", lin_times)
			#print("log times",log_times)
			times = np.concatenate([lin_times, log_times])


			#print("times",times)

			#for row in matrix:
				#print(row)



			#change pops to real pops 
			curr_pops = []
			no_conc_complexes = []
			for complex,index in cplx_indices.items():
				if complex._concentration == None:
					no_conc_complexes.append(complex)
				else:
					curr_pops.append(complex._concentration[2])

			for complex in no_conc_complexes:
				complex._concentration = [0,1/len(no_conc_complexes),"M"]
				curr_pops.append(complex._concentration[1])
			
				


			logger.debug(f"Current populations {curr_pops}")
			
			for t, pt in mx_simulate(matrix,curr_pops,times):
				#logger.debug(f"\n\n\nFinal Populations\n{t:8.6e} {' '.join([f'{x:8.6e}' for x in abs(pt)])}")
				pops = pt


			logger.debug(f"Final Populations: {pops}")
			"""
			if cplx_indices is not None:
				#print(f"Cplx indices {cplx_indices}")
				for complex,index in cplx_indices.items():
					#print(f"Complex {complex._name}, Pop: {pops[index]}")
					#pops_dict[complex._name] = pops[index]
					#parameter["complexes"][complex._name][1] = pops[index]"""
			for complex,index in cplx_indices.items():
					logger.info(f"Complex {complex._name}, Pop: {pops[index]}")
					parameter["complexes"][complex._name][1] = pops[index]
					complex._concentration = [0,pops[index],"M"]

	
	for complex in parameter["complexes"]:				
		pops_dict[complex] = parameter["complexes"][complex][1]

	

	output = enum.to_pil(condensed=condensed,detailed = not condensed) # should be not condensed for same output format like peppercorn
	resting_complexes = [cplx for cplx in natsorted(enum.resting_complexes)] 
	

	if len(resting_complexes) == 1:
		for complex in resting_complexes:
			complex._concentration = [0,1,"M"]

	print("\n\n\n\n_________________________Complex Concentrations____________________\n\n")
	for complex in resting_complexes:
		print(complex._name,complex._concentration[1])
	

	#update parameters["complexes"]
	if len(rate_matrices) > 0:
		logger.info(f"\n\n\n\nOutput: \n {output} \n\n\n")
		return resting_complexes, pops_dict

	else:
		logger.info(f"\n\n\n\nOutput: \n {output} \n\n\n")
		print("Pops changed")
		return resting_complexes, pops_dict

	


def input_parsing(d_seq, complexes):
	"""Takes domain level sequence and structures and formats those into an input format for peppercorn

	Args:
		d_seq(str) : domain level sequence
		complexes(dict)	: dict of the complexes in the form of {Name:[kernel structure,population]}
	"""
	toehold_length = 5
	logic_length = 8
	space_length = 8
	logger.debug(f"Input Parsing: Structure input \n {complexes}")
	unique_domains = set([domain.replace('*', '') for domain in d_seq])

	system_input = f""
	for unique_domain in unique_domains:
		if unique_domain.islower():
			system_input += f"length {unique_domain} = {toehold_length} \n"
		elif unique_domain.startswith("L"):
			system_input += f"length {unique_domain} = {logic_length} \n"
		elif unique_domain.startswith("S"):
			system_input += f"length {unique_domain} = {space_length} \n"

	system_input += "\n"

	for name, lst in complexes.items():
		system_input += f"{name} = {lst[0]}\n"
		


	logger.debug(f"\n\n\nSystem Input \n\n{system_input}\n\n")
	complexes, reactions = read_pil(system_input)
	logger.debug(f"\nResulting complexes:{complexes}\nReactions{reactions}")
	return complexes, reactions


def extract_domain_sequence(input_lines):
	try:
	
		for line in input_lines:
			if not line.startswith("#") and not line.startswith("[") and not line.startswith("I") and len(line) > 1 and not line.startswith("Resulting Domain Level sequence:"):
				result_sequence = line
				return result_sequence
			elif line.startswith("Resulting Domain Level sequence:"):
				result_sequence = line.split(":",1)[1].strip()
				return result_sequence
	except:
		raise ImportError("Check Input: can't find anything")
		


def extract_afp(input_lines):
	try: 
		if input_lines[0] == "\n":
			return eval(next((line for line in input_lines if line.startswith("[")),None))
			
	except:
		return None 


def set_verbosity(console_handler,verbosity):
	if verbosity == 1:
		console_handler.setLevel(logging.INFO)
	elif verbosity >= 2:
		console_handler.setLevel(logging.DEBUG)


def main():
	parser = argparse.ArgumentParser(description="cocosim is a cotranscriptional folding path simulator using peppercornenumerate to simulate a domain level sequence during transcription.")
	
	parser.add_argument("-i", "--input_file", nargs='?', type=argparse.FileType('r'), default=sys.stdin,
						help="Input file. If not provided, reads from stdin.")
	parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=10)
	parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity level. -v only shows peppercorn output")
	parser.add_argument("-c","--condensed", action= "store_true",default=False, help="Condense reactions into only resting complexexes. (default: False)")



	args = parser.parse_args()
	
	set_verbosity(logger,args.verbose)


	if args.input_file.isatty():
		print("Please enter a domain level sequence:")
		d_seq = input()
		if len(d_seq) == 0:
			raise SystemExit("No Input given")
	else:
		input_lines = args.input_file.readlines()
		d_seq = extract_domain_sequence(input_lines)
		afp = extract_afp(input_lines)
		args.input_file.close()



	parameters = {"k_slow": args.k_slow, "condensed":args.condensed, "complexes": {}}



	
	logger.info(parameters)


	print("Given Domain Sequence:", d_seq)

	if afp:
		print("Following AFP was given for Sequence design:")
		for step in afp:
			print(step)
		print("\n")




	#______Running_Simulation___________# 
	simulated_structures = run_sim(d_seq, parameters)
	

	#_____Writing_and_printing_output___# 

	output = ""

	if args.verbose > 0:
		ts = 1
		output += ("\nResting Complexes after each step:\n\n")
		output += "Transcription Step | Structure \n"
		for x in simulated_structures: 
			for seq in x[0]:
				output += f"{ts}	{seq} \n"
			ts += 1 
			output += "\n"
	

	output += write_output(simulated_structures,d_seq)
	
	output += f"\n\nFollowing sequence was simulated:\n{d_seq}"

	if afp: 
		output += f"\n\nFollowing AFP was given for Sequence design"
		for step in afp:
			
			output += f"{step}\n"

	print(output)




if __name__ == "__main__":
	main()
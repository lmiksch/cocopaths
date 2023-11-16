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
from .utils import cv_db2kernel, kernel_to_dot_bracket, only_logic_domain_struct
import argparse
import sys
import logging 


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
	global used_structure_names
	used_structure_names = 0

	d_seq_split = d_seq.split()
	complexes, reactions = input_parsing(d_seq_split[0:2], [cv_db2kernel(d_seq_split[0:2], "..")])

	final_structures = [[] for x in range(len(d_seq_split))]

	current_complexes = []


	final_structures[0].append([d_seq_split[0]])
	first_complexes = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters)
	
	current_complexes.append(first_complexes[0])
	final_structures[1].append(first_complexes)


	
	for step in range(2, len(d_seq_split)):
		logger.debug("\n\n\nNext Step\n\n")
		logger.debug(f"Step: {step} Current Complexes: {current_complexes}")
		next_complexes = []
		
		for structure in current_complexes:
			next_struct = structure + " " + d_seq_split[step]
			logger.debug(f"Next Structure {next_struct}")
			next_complexes.append(next_struct)

		complexes, reactions = input_parsing(d_seq_split[0:step + 1], next_complexes)
		logger.debug(f"Following Complexes and reactions are put into enumerator \n {complexes} {reactions}")
		resting_complexes = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters)

		logger.debug(f"resulting resting complex: {resting_complexes}")

		final_structures[step].append(resting_complexes)
		current_complexes = resting_complexes

	return final_structures

	
	
def write_output(final_structures,d_seq,parameters = None):

	data_output = ""
	


	ts = 1
	data_output += ("\nResting Complexes after each Spacer:\n\n")
	data_output += "Transcription Step | Structure \n"
	for x in final_structures:
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

	init_cplxs = [x for x in complexes.values() if x.concentration is None or x.concentration[1] != 0]
	name_cplxs = list(complexes.values())
	enum = Enumerator(init_cplxs, reactions, named_complexes=name_cplxs)

	if k_slow: 
		enum.k_slow = k_slow
	condensed = parameter["condensed"]


	enum.enumerate()
	
	
	if condensed:
		enum.condense()



	# Start to enumerate
	

	output = enum.to_pil(condensed=condensed,detailed=condensed)
	resting_complexes = [cplx.kernel_string for cplx in natsorted(enum.resting_complexes)] 

	logger.info(f"\n\n\n\n Output: \n {output} \n\n\n")
	return resting_complexes


def input_parsing(d_seq, structures):
	global used_structure_names
	toehold_length = 5
	logic_length = 8
	space_length = 8
	logger.debug(f"Input Parsing: Structure input \n {structures}")
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

	for x, structure in enumerate(structures):
		system_input += f"S{used_structure_names} = {structure}\n"
		used_structure_names += 1 


	logger.debug(f"\n\n\nSystem Input \n{system_input}\n\n")
	complexes, reactions = read_pil(system_input)
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



	parameters = {"k_slow": args.k_slow, "condensed":args.condensed}
	
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
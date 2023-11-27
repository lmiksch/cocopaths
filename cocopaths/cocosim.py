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

	# Initiate datastructures
	used_structure_names = 1
	all_complexes = {}
	folding_step_complexes = []
	d_seq_split = d_seq.split()
	complex_count = 1



	# before input parsing define name of complexes 
	
	parameters["complexes"]["E0"] = [cv_db2kernel(d_seq_split[0:2], ".."),1]
	logger.debug(parameters["complexes"])
	complexes, reactions = input_parsing(d_seq_split[0:2], parameters["complexes"])
	

	resulting_complexes,transient_complexes = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters, all_complexes=all_complexes)
	
	
	

	# all complex keeps track of complexes throughtout the sim
	

	
	
	for complex in resulting_complexes:
		all_complexes["Id_" + str(complex_count)] = [complex,1/len(resulting_complexes)]
		complex_count += 1
	

	folding_step_complexes.append(resulting_complexes)	
	


	for step in range(2, len(d_seq_split)):
		logger.debug("\n\n\nNext Step\n\n")
		logger.info(f"Step: {step} Current Complexes: {parameters}")

		next_complexes = []
		

		#put here add one domain to struct and give new name
	
		new_complexes = {}
		

		#print("\n\nextendo\n\n")
		old_names = []
		logger.debug(f"All Complexes before:")
		for key,complex in all_complexes.items():
			logger.debug(f"{key}:{complex}\n")
		for complex in folding_step_complexes[-1]:
			print(complex)
			c_name = complex._name
			next_struct = complex.kernel_string + " " + d_seq_split[step]

			new_complexes["E" + str(used_structure_names)] = [next_struct,complex._concentration[1]]
			
			used_structure_names += 1
			old_names.append(c_name)
		#print(new_complexes)
		#print("\n\nextendo endo\n\n")

		
		logger.debug(f'Resulting new complexes: {new_complexes}')

		
		
		complexes, reactions = input_parsing(d_seq_split[0:step + 1], new_complexes)


		logger.info(f"New Complexes {complexes}")
		for key,complex in complexes.items():
			logger.info(f"{complex}:{complex._concentration}\n")


		logger.info(f"All Complexes before:")
		for key,complex in all_complexes.items():
			logger.info(f"{key}:{complex} {complex[0]._concentration	}\n")


		new_complexes = complexes.copy() 
		new_complex_items = []

		assert len(all_complexes) == len(new_complexes),"Not equal cant assign"

		"""
		for id, new_complex in zip(all_complexes,new_complexes):
			all_complexes[id] = [new_complex,99]"""
		
		for x, c_name in enumerate(old_names):
			for cid, complex_obj in all_complexes.items():
				
				if c_name == complex_obj[0]._name:
					next_complex = complexes.popitem()[1]
					next_complex._concentration = complex_obj[0]._concentration
					new_complex_items.append((cid, next_complex))

				elif len(complexes) == 1:
					next_complex = complexes.popitem()[1]
					next_complex._concentration = complex_obj[0]._concentration
					new_complex_items.append((cid,next_complex))

		
		#Add remaining ones to all_complexes should incluce concentration aswell
		if len(complexes) > 1:
			for key,complex in complexes.items():
				
				all_complexes["Id_" + str(complex_count)] = [complex,99] #dont know yet how to incorperate the concentration
				complex_count += 1

		

				
		# Add the new items to the dictionary
		for cid, complex_item in new_complex_items:
			#print("Complex Item",complex_item,complex_item._concentration)
			all_complexes[cid][0] = complex_item

		#for key, complex in all_complexes.items():
			#print(f"{key}	{complex[0].kernel_string}")
			#print(f"{complex[0]._concentration[1]}")
					
		
				
		logger.info(f"All Complexes after: ")
		for key,complex in all_complexes.items():
			logger.info(f"{key}:{complex}")

		#print("\n\n\nAll Complexes After extension ", all_complexes)

		complexes = new_complexes


		resting_complexes,transient_complexes = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters,all_complexes=all_complexes)
		
		folding_step_complexes.append(resting_complexes)

		#checks if there is a transient state which needs to be mapped
		for key,input_complex in complexes.items():
			for new_complex in transient_complexes:				
				if input_complex.kernel_string == new_complex.kernel_string:
						map_transient_states(resting_complexes,transient_complexes,all_complexes)
		
		
		map_resting_states(resting_complexes,transient_complexes,all_complexes,complex_count)



		#Check if unique ids match resting complexes should always be larger or equal 
		print("Current unique ids\nID	Structure	Concentration")
		for key, complex in all_complexes.items():
			print(f"{key} {complex[0]._name}	{complex[0].kernel_string}	{complex[0]._concentration}")

		assert len(all_complexes) >= len(resting_complexes),"More Resting complexes than unique_ids something went wrong."





	
	

	return folding_step_complexes

	
	



def enumerate_step(complexes, reactions, parameter, all_complexes):
	"""Takes complexes in form of PepperComplexes and uses Peppercorn Enumerator to find possible structures and their respective occupancies. 
	
	Args:
		complexes(dict): {Name:PepperComplex}
		reactions(class): PepperReaction not really used in our program
	Returns:
		resulting_complexes(list): Resulting PepperComplexes
	"""
	print("\n\nComplexes ", complexes)
	k_slow = parameter["k_slow"]
	
	if bind21 in BI_REACTIONS:
		BI_REACTIONS.remove(bind21)

	logger.debug(f"\n\n\nBeginning of Enumerate step:\nComplexes:{complexes}\nReactions:{reactions}\nParameters:{parameter}\n\n")
	init_cplxs = [x for x in complexes.values() if x.concentration is None or x.concentration[1] != 0]
	name_cplxs = list(complexes.values())
	enum = Enumerator(init_cplxs, reactions, named_complexes=name_cplxs)

	if k_slow: 
		enum.k_slow = k_slow


	#need to find suitable k_fast
	enum.k_fast = 100
	
	# Start to enumerate
	enum.enumerate()
	
	
	condensed = parameter["condensed"]
	if condensed:
		rate_matrices, complex_indices = enum.condense()


	resulting_complexes = [cplx for cplx in natsorted(enum.resting_complexes)] 
	
	transient_complexes = [cplx for cplx in natsorted(enum.transient_complexes)]

	"""
	if len(resulting_complexes) > 1 and len(enum.condensation._condensed_reactions) > 0: 
		update_occupancies(resulting_complexes,enum,all_complexes)
	else:
		for complex in resulting_complexes:
			if complex._concentration == None:
				complex._concentration = [0,1,"nM"]
			else:
				complex._concentration[1] = 1
	

	
	
	#Calculating Populations and updating the pop ins parameters
	pops_dict = {}
	for matrix,cplx_indices in zip(rate_matrices,complex_indices):
		
		if matrix is not None and len(matrix) > 1:
			

			#logger.debug(f"\nAfter Addition of new Structures\n{[parameter['complexes']]}")
			t1 = 0.1 
			t8 = 0.02

			lin_times = np.array([t1], dtype='float128')
			log_times = np.array(np.logspace(np.log10(t1), np.log10(t8), num=10, dtype='float128'))
			log_times = np.delete(log_times, 0)
			
			log_times = []

			times = np.concatenate([lin_times, log_times])



			#change pops to real pops 
			curr_pops = []
			no_conc_complexes = []
			for complex,index in cplx_indices.items():
				if complex._concentration == None:
					no_conc_complexes.append(complex)
				else:
					curr_pops.append(complex._concentration[2])

			for complex in no_conc_complexes:
				complex._concentration = [0,1/len(no_conc_complexes),"nM"]
				curr_pops.append(complex._concentration[1])
			
				


			logger.debug(f"Current populations {curr_pops}")
			
			#change to pilsimulator maybe use mx_simulate later 
			for t, pt in mx_simulate(matrix,curr_pops,times):
				pops = pt


			logger.debug(f"Final Populations: {pops}")


			for complex,index in cplx_indices.items():
					logger.info(f"Complex {complex._name}, Pop: {pops[index]}")
					complex._concentration = [0,pops[index],"nM"]

	

	#for complex in parameter["complexes"]:				
		#pops_dict[complex] = parameter["complexes"][complex][1]

	

	
	
	
	if len(resulting_complexes) == 1:
		for complex in resulting_complexes:
			complex._concentration = [0,1,"nM"]"""
	
	
	output = enum.to_pil(condensed=condensed,detailed = condensed) # should be not condensed for same output format like peppercorn
	

	


	#print("\n\n\n\n_________________________Complex Concentrations____________________\n\n")
	for complex in resulting_complexes:
		if complex._concentration is  None:
			complex._concentration = [0,9,'nM']
	
	print("\nAll Complexes: ")
	for key, complex in all_complexes.items():
		print(f"{key}	 {complex[0]._name}		{complex[0].kernel_string}		{complex[0]._concentration}")

	#update parameters["complexes"]
	if len(rate_matrices) > 0:
		logger.info(f"\n\n\n\nOutput: \n {output} \n\n\n")
		return resulting_complexes, transient_complexes

	else:
		logger.info(f"\n\n\n\nOutput: \n {output} \n\n\n")
		return resulting_complexes, transient_complexes



def map_transient_states(resting_complexes,transient_complexes,all_complexes):


	print(f"Map Transient states\n\n{transient_complexes}\n\n{all_complexes}")
	if len(resting_complexes) > 1 and len(transient_complexes) > 0:
		print("this")


def map_resting_states(resting_complexes,transient_complexes,all_complexes,complex_count):



	


	print("all_complexes",all_complexes)
	print("resting complexes", resting_complexes)
	for r_complex in resting_complexes:	
		if r_complex not in list(all_complexes.values())[0]:
			all_complexes["Id_" + str(complex_count)] = [r_complex,99]
			complex_count += 1 

	assert len(all_complexes) >= len(resting_complexes),"More Resting complexes than unique_ids something went wrong."

	print('map resting states')
	for key,complex in all_complexes.items():
		print(f"{key}	 {complex[0]._name}		{complex[0].kernel_string}		{complex[0]._concentration}")




def update_occupancies(resulting_complexes,enum,all_complexes):
	"""Updates Occupancies based on the resulting condensed 			resulting_complexes(list): Updated list of the resulting complexes with their respective concentrations
rates and stationary distribution of macrostates. 
		
		Args:
			resulting_complexes(list): List of the resulting resting complexes from enumeration step 
			enum(object): enumerator object which was used to calculate the resulting complexes

		Returns:
			resulting_complexes(list): Updated list of the resulting complexes with their respective concentrations
	
	"""	
	print("\n\n_______________Begin update Occupancy___________________________\n\n")
	# First calculate the distribution in the Macrostates to calculate the reactant occupancy of the outgoing reaction 
	
	macrostates = []
	exit_probs = {}


	# modify it so that it uses the stationary distribution 
	# Iterate over the outer dictionary 
	for key, inner_dict in enum.condensation.exit_prob.items():
		# Iterate over the inner dictionary 
		for reaction, value in inner_dict.items():
			# Extracting the float value 
			exit_prob = value
			

			# Extracting the product of the reaction
			reaction_product = reaction[1]._products  
			

			exit_probs[reaction_product[0]] = exit_prob

	#print(f"\n\nExit Probs: {exit_probs}\n\n")

	condensed_rates = {}
	stationary_dist = enum.condensation.stationary_dist
	for reaction in enum.condensation.condensed_reactions:
		#print(reaction)
		macrostates.append(reaction._reactants[0]._complexes)
		#print(reaction._reactants,reaction._reactants)
		for reactant,product in zip(reaction._reactants,reaction._products):
			#print("reactant",reactant,product)
			condensed_rates[(reactant._complexes[0], product._complexes[0])] = reaction._const

	#print(f"Condensed Rates: ",condensed_rates)
	#print(f"\n{macrostates}")
	#print(f"\n {exit_probs}\n")
	#print(f"\n {stationary_dist}\n\n\n")

	#print(resulting_complexes)
	
	#print("All Complexes: ")
	for key, complex in all_complexes.items():
		print(f"{key} {complex[0]._name}	{complex[0].kernel_string}	{complex[0]._concentration}")



	#create rate matrix based on condensed reactions
	rate_matrix = [[0 for _ in resulting_complexes] for _ in resulting_complexes]

	for y,y_complex in enumerate(resulting_complexes):
		
		for x,x_complex in enumerate(resulting_complexes):
			if (y_complex,x_complex) in condensed_rates:
				rate_matrix[y][x] = condensed_rates[(y_complex,x_complex)]
				#rate_matrix[x][y] = - condensed_rates[(y_complex,x_complex)]
			else:
				rate_matrix[y][y] -= rate_matrix[y][x]
	
	for i in range(len(resulting_complexes)):
		rate_matrix[i][i] = -sum(rate_matrix[i][:i] + rate_matrix[i][i+1:])

	
	


	# adjust populations based on exit probs

	
	populations = [] 

	for complex in resulting_complexes: 
		if complex in exit_probs and complex._concentration != None:
			populations.append(complex._concentration[1] * exit_probs[complex])
		else:
			populations.append(0)

	#logger.debug(f"\nAfter Addition of new Structures\n{[parameter['complexes']]}")
	t1 = 0.1 
	t8 = 0.02

	lin_times = np.array([t1], dtype='float128')
	log_times = np.array(np.logspace(np.log10(t1), np.log10(t8), num=10, dtype='float128'))
	log_times = np.delete(log_times, 0)
	
	log_times = []

	times = np.concatenate([lin_times, log_times])

	rate_matrix_np = np.array(rate_matrix)

	zero_mask = rate_matrix_np == 0

	rate_matrix_np[zero_mask] = None
	"""
	print("\n\nRate Matrix")
	for line in rate_matrix:
		print(line)"""

	###
	### Here Pilsimulator Stuff
	###


	"""
	print("Resulting Populations: ",pops)
	
	for complex,pop in zip(resulting_complexes,pops):
		if complex._concentration != None:
			complex._concentration[1] = pop
		else:
			complex._concentration = [0,pop,"nM"]

	"""

def write_output(final_structures,d_seq,parameters = None):

	data_output = ""
	

	ts = 1
	data_output += ("\nResting Complexes after each Spacer:\n\n")
	data_output += "Transcription Step | Structure | Occupancy \n"
	for x in final_structures:
		if x and x[0].kernel_string[-2] == "S": 
			for complex in x:
				data_output += f"{ts} {complex._name}	{complex.kernel_string}	{round(complex._concentration[1],4)} \n"
			ts += 1 

			data_output += "\n"
	
	'''
	ts = 1
	data_output += ("\n\nOnly Logic Domain pairings:\n\n")
	data_output += "Transcription Step | Structure	| Concentration \n"

	for x in final_structures:
		if x and x[0].kernel_string[-2] == "S": 
			for complex in x:
				kernel_string = kernel_to_dot_bracket(complex.kernel_string)
				#db_struct = (only_logic_domain_struct(d_seq.split(),kernel_string))
				data_output += f"{ts}	{db_struct} 	{round(complex._concentration[1],4)}\n"
			ts += 1
			#data_output += "\n"'''

	return data_output

	

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
		


	logger.info(f"\n\n\nSystem Input \n\n{system_input}\n\n")
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
		output += "Transcription Step | Structure | Occupancies \n"
		for x in simulated_structures: 
			for complex in x:
				output += f"{ts}	{complex.kernel_string}	{round(complex._concentration[1],4)}\n"
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
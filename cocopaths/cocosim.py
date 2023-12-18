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
from peppercornenumerator.objects import PepperComplex
from crnsimulator import ReactionGraph, get_integrator
from crnsimulator.odelib_template import add_integrator_args

from io import StringIO
from natsort import natsorted
from .utils import cv_db2kernel, kernel_to_dot_bracket, only_logic_domain_struct
import argparse
import logging 
import numpy as np 



logger = logging.getLogger('cocosim')
console_handler = logging.StreamHandler()
formatter = logging.Formatter('# %(levelname)s - %(message)s')
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)




def run_sim(d_seq, parameters,args):
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
	

	resulting_complexes,transient_complexes,enum = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters, all_complexes=all_complexes)
	
	
	
	map_transient_states(resulting_complexes,transient_complexes,all_complexes,enum,complex_count)
	
	
	
	# all complex keeps track of complexes throughtout the sim
	

	
	#only in first step this should only result in 1 structure
	for complex in resulting_complexes:
		all_complexes["Id_" + str(complex_count)] = [complex,1/len(resulting_complexes)]
		complex.occupancy = 1
		complex_count += 1
	

	folding_step_complexes.append(resulting_complexes)	
	


	for step in range(2, len(d_seq_split)):
		logger.info("\n\n\n______________________________________________")
		logger.info(f"Step: {step} Current Complexes: {parameters}")

		next_complexes = []
		

		#put here add one domain to struct and give new name
	
		new_complexes = {}
		

		#print("\n\nextendo\n\n")



		#______Function_for_extending_and_updating_names__ 

		old_names = []
		logger.debug(f"All Complexes before:")
		for key,complex in all_complexes.items():
			logger.debug(f"{key}:{complex} {complex[0].occupancy}\n")
		for complex in folding_step_complexes[-1]:
			#print("complex: ",complex)
			c_name = complex._name
			next_struct = complex.kernel_string + " " + d_seq_split[step]

			new_complexes["E" + str(used_structure_names)] = [next_struct,complex.occupancy]
			
			used_structure_names += 1
			old_names.append(c_name)

		#print(new_complexes)
		#print("\n\nextendo endo\n\n")


		logger.debug(f'Resulting new complexes: {new_complexes}')

		
		
		complexes, reactions = input_parsing(d_seq_split[0:step + 1], new_complexes)

		#Update all_complexes with new complexes 

		logger.debug(f"New Complexes")
		for key,complex in complexes.items():
			logger.debug(f"{complex}:{complex.occupancy}\n")


		logger.debug(f"All Complexes before:")
		for key,complex in all_complexes.items():
			logger.debug(f"{key}:{complex} {complex[0].occupancy	}\n")


		new_complexes = complexes.copy() 
		new_complex_items = []

		#assert len(all_complexes) == len(new_complexes),"Not equal cant assign"

		"""
		for id, new_complex in zip(all_complexes,new_complexes):
			all_complexes[id] = [new_complex,99]"""
		

		print("\n\ncomplexes ",len(complexes), complexes)
		print("Old Names", len(old_names), old_names)

		for x, c_name in enumerate(old_names):
			for cid, complex_obj in all_complexes.items():
				
				if c_name == complex_obj[0]._name:
					#print("\n",c_name,complex_obj[0]._name)
					#print("\n ",complexes,"\n")
					first_key = next(iter(complexes))
					next_complex = complexes.pop(first_key)
					next_complex.occupancy = complex_obj[0].occupancy
					new_complex_items.append((cid, next_complex))
				
		# Add the new items to the dictionary
		for cid, complex_item in new_complex_items:
			print("Complex Item",complex_item,complex_item.occupancy)
			all_complexes[cid] = [complex_item,complex_item.occupancy]	


		#print("New_complex_items:",new_complex_items)
		#Add remaining ones to all_complexes should incluce concentration aswell
		print("Remaining complexes",complexes)
		if len(complexes) > 0:
			for key,complex in complexes.items():
				print("Remaining complex in loop",complex)
				all_complexes["Id_" + str(complex_count)] = [complex,complex.occupancy] 
				complex_count += 1

		

				
		
		
		logger.debug(f"All Complexes after: ")
		for key,complex in all_complexes.items():
			logger.debug(f"{key}:{complex}")
			if complex[0].occupancy == None or complex[1] == None:
				raise SystemExit("Unoccupied Complex in all complex after updating")


		
		#print("\n\n\nAll Complexes After extension ", all_complexes)

		complexes = new_complexes


		resting_complexes,transient_complexes,enum = enumerate_step(complexes=complexes, reactions=reactions, parameter=parameters,all_complexes=all_complexes)
		
		
		
		#checks if there is a transient state which needs to be mapped
		
		map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,complex_count)
						
		
		#Additionally maps resting states 
		calc_macro_pop(enum,all_complexes,resting_complexes,complex_count)


		
		#simulate the condensed reactions
		resting_complexes = simulate_system(enum,args,resting_complexes,all_complexes)
		print(resting_complexes)

		#need to fix the duplication in update_macrostates -> should be fixed 
		resting_complexes = list(set(resting_complexes))
		oc_sum = 0
		for complex in resting_complexes:
			oc_sum += complex.occupancy 
		oc_sum = round(oc_sum,4) #can be later adjusted

		logger.debug(f"\n________________________________________\nFollowing Resting complexes result\n______________\n ")
		for complex in resting_complexes:
			logger.debug(f"{complex} {complex.occupancy}")
		logger.debug("-----")
		for complex in all_complexes.items():
			logger.debug(f"{complex[0]} {complex}")

		if oc_sum == 1:	
			folding_step_complexes.append(resting_complexes)
			print("\n\nFinal appending complexes: \n")
			for complex in resting_complexes:
				print(complex,complex.occupancy)
				if not is_complex_in_all_complexes(complex,all_complexes):
					raise SystemExit("resulting resting complex not in all complexes")
			print("\n")
		else:
			raise SystemExit(f"SystemExit: Occupancies summ up to {oc_sum}")
		
		#assert len(all_complexes) >= len(resting_complexes),"More Resting complexes than unique_ids something went wrong."

	

	return folding_step_complexes


def simulate_system(enum,args,resting_complexes,all_complexes):
	
	
	condensed_reactions = enum.condensation._condensed_reactions
	
	reactions = []

	oc_vector = []
	seen = set()
	if condensed_reactions:
		#create reactions list and occupancy_vector needed for crn simulations
		for reaction in condensed_reactions:
			print(reaction._reactants[0]._complexes)
			if reaction._reactants[0].name not in seen: 
				seen.add(reaction._reactants[0].name)
				
				reactants = [reaction._reactants[0].name]
				products = [reaction._products[0].name]

				occupancy1,occupancy2 = 0,0 

				print("\n\n calc macro occupancies\n")
				for r_complex in reaction._reactants[0]._complexes: 
						
					#set conc to 0 for undefined complexes 
					if r_complex.occupancy == None:
						r_complex.occupancy = 0
					
					
					
					
					occupancy1 += r_complex.occupancy
					print("occ1 ",occupancy1,r_complex)
					
					
				for p_complex in reaction._products[0]._complexes: 
					occupancy2 += p_complex.occupancy
					print("occ2 ",occupancy2,p_complex)


				rate = reaction._const
				reactions.append([reactants,products,[rate]])
				oc_vector.append(float(occupancy1))
				if occupancy2:
					oc_vector.append(float(occupancy2))
				else:	
					oc_vector.append(0)


				

		print("reactions: ",reactions)
		print("\n Oc vector: ",oc_vector)
		resulting_occupancies =	sim_condensed_rates(reactions,oc_vector,args)
		print("\n\nResulting occupancies after simulation",resulting_occupancies)
		#Update occupancies

		resting_complexes = update_macrostates(resulting_occupancies,all_complexes = all_complexes,enum= enum,resting_complexes=resting_complexes)

	return resting_complexes

def update_macrostates(result_occ,all_complexes,enum,resting_complexes):
	"""
	
	Args: 
		result_occ(dict): Dictionary, key = complex._name, value float
		all_complexes(dict): Id: [Complex,conc]
		enum(object): enumerated object from peppercorn
		resting_complexes(list)
	"""

	print("\n\nBegin update macrostates")

	print("\n\n_______________________________befor alle complexes\n")
	for complex in all_complexes.values():
		print(complex, complex[0].occupancy)

	print(result_occ)


	
	"""for reaction in result_occ.values():
		print(reaction) 
		for key,value in reaction.items():
			for cid,complex in all_complexes.items():
					if complex[0]._name == key:
						complex[0].occupancy = value
						update_complex_in_all_complexes(complex[0],all_complexes)
			for macrostate in enum._resting_macrostates:
				print(macrostate)"""

	for reaction in result_occ.values():
		print(reaction) 
		for key,value in reaction.items():
			for macrostate in enum._resting_macrostates:
				print(macrostate)
				if any(complex_obj._name == key for complex_obj in macrostate.complexes):
					macrostate.occupancy = value 
					print("\n\n\nresulting macrostate:  ",macrostate, macrostate.occupancy)


	print("\n\n_______________________________after alle complexes\n")
	for complex in all_complexes.values():
		print(complex, complex[0].occupancy)
	
		
	#print("After all complexes",all_complexes)

	print("before calc macro pop\n")
	for complex in resting_complexes:
		print(complex,complex.occupancy)
	print("Summe over resting complexes:", sum([complex.occupancy for complex in resting_complexes]))

	calc_macro_pop(enum,all_complexes,resting_complexes,complex_count)


	print("after calc macro pop\n")
	for complex in resting_complexes:
		print(complex,complex.occupancy)
	print("Summe over resting complexes:", sum([complex.occupancy for complex in resting_complexes]))
	#print("Resulting complexes",resulting_complexes)
	return resting_complexes

def is_complex_in_all_complexes(complex,all_complexes):

	for id,value in all_complexes.items():
		if complex == value[0]:	
			return True
	
	return False


def sim_condensed_rates(reactants,concvect,args):
	"""Script is based on Pillsimulator script from Peppercorn and uses CRNsimulator to simulate a chemical reaction network. 
	In this case it calculates the occupancies after accounting for the condensed rates.
	"""
	

	#need to hardcode all args:


	args.pyplot_labels = None
	
	args.list_labels = None

	args.nxy = False

	args.t8 = 20

	args.t0 = 0

	args.t_log = False

	args.t_lin = 2

	args.atol = None
	args.rtol = None
	args.mxstep = 0 
	args.labels_strict = False
	args.header = True

	args.pyplot = False
	args.labels = []
	args.pyplot_xlim = None 
	args.pyplot_ylim = None

	args.cutoff = 0.0000001

	p0 = []
	for i,s in enumerate(concvect,start=1):
		p0.append(str(i) + "=" + str(s))


	
	args.p0 = p0
	
	filename = "filename"

	odename = "odename"


	V = []
	for reaction in reactants:
		V.append(reaction[0][0])
		V.append(reaction[1][0])	

	

	# Split CRN into irreversible reactions
	crn = reactants
	new = []
	for [r, p, k] in crn:
		assert not (None in k)
		if len(k) == 2:
			new.append([r, p, k[0]])
			new.append([p, r, k[1]])
		else:
			new.append([r, p, k[0]])
	crn = new

	#print("Crn,species",crn)
	RG = ReactionGraph(crn)
	
	C = []
	Vars = []
	seen = set()
	for s,c in zip(V,concvect):
		if s not in seen:
			seen.add(s)
			#print("\n\nS,c\n\n",s,c)
			Vars.append(s)
			C.append(c)


	# set p0 
	p0 = []
	for i,s in enumerate(C,start=1):
		p0.append(str(i) + "=" + str(s))


	#print("p0",p0)
	args.p0 = p0

	#print("________V,C______",Vars,C)
	filename, odename = RG.write_ODE_lib(sorted_vars = Vars, concvect = C,
									  		jacobian= True,
											  const = [False for x in Vars],
											 filename = filename,
											 odename = odename)
		

	integrate = get_integrator(filename)

	#original_stdout = sys.stdout
	#sys.stdout = StringIO()

	time,occupancies = integrate(args) 
	print(occupancies)
	
	


	

	end_conc = [oc for oc in occupancies[1:]]
	print("end concs :",end_conc)
	

	resulting_concentrations = {}
	for i in range(len(reactants)):
		resulting_concentrations[i] = {}
		for complex, conc in zip(V,end_conc):
			if float(conc) >= args.cutoff:
				resulting_concentrations[i][complex] = conc
			else: 
				resulting_concentrations[i][complex] = 0


	#print(f"\n\n\nResulting Concentrations: {resulting_concentrations}\n\n\n\n")
	return resulting_concentrations



def calc_macro_pop(enum,all_complexes,resulting_complexes,complex_count):
	#print("\n______________\nCalc Macropop\n\n")
	resting_macrostates = enum._resting_macrostates

	stat_dist = enum.condensation.stationary_dist

	
	for macrostate in resting_macrostates:
		stat_dist_copy = dict(stat_dist[macrostate])

		try:
			macro_pop = macrostate.occupancy

			print("___________________yaaaaaas____",macro_pop,macrostate)


		except:
			macro_pop = 0
			for complex in macrostate._complexes:
					#print("\n\nComplex occupancy",complex.occupancy)

					# add to macrostate population if no occupancy is defined -> complex new 
					try:
						macro_pop += complex.occupancy

					except:
						macro_pop += 0
						
			if macro_pop == 0: 
				macro_pop = 1
			#print("Macro POP",macro_pop)
		
		for stat_complex,pop in stat_dist_copy.items():
			print("Stat Complex",stat_complex)
			print("pop,macropop",pop,macro_pop)
			new_dist = pop * macro_pop
			#stat_dist[macrostate][stat_complex] = new_dist
			print("stat_complex, new dist",stat_complex,new_dist)
			if is_complex_in_all_complexes(stat_complex,all_complexes):
				
				for c_id, all_complex in all_complexes.items():
			
					if stat_complex == all_complex[0]:
						#print("1:",complex)
						all_complexes[c_id] = [stat_complex,new_dist]
				
			else:
				#print("complex count",complex_count,stat_complex)
				all_complexes["Id_" + str(complex_count)] = [stat_complex,new_dist]
				complex_count += 1

			for r_complex in resulting_complexes:
				#print(complex,stat_complex)
				if stat_complex == r_complex: 
					print(stat_complex,float(new_dist))
					r_complex.occupancy = float(new_dist)
					
					
	#print("Resulting all_complexes")
	#for key,complex in all_complexes.items():
		#print(complex, complex[0].occupancy)
	#print("\n\nEnd Macropops\n\n")
	#print(all_complexes)
	#for each macrostate get input kernel and mulitply all results with 

	#return resting_complexes



def enumerate_step(complexes, reactions, parameter, all_complexes):
	"""Takes complexes in form of PepperComplexes and uses Peppercorn Enumerator to find possible structures and their respective occupancies. 
	
	Args:
		complexes(dict): {Name:PepperComplex}
		reactions(class): PepperReaction not really used in our program
	Returns:
		resulting_complexes(list): Resulting PepperComplexes
	"""
	#print("\n\nComplexes ", complexes)
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
	enum.k_fast = 20
	
	# Start to enumerate
	enum.enumerate()
	
	
	condensed = parameter["condensed"]
	if condensed:
		enum.condense()


	resulting_complexes = [cplx for cplx in natsorted(enum.resting_complexes)] 
	
	transient_complexes = [cplx for cplx in natsorted(enum.transient_complexes)]
	
	
	output = enum.to_pil(condensed=condensed,detailed = condensed) # should be not condensed for same output format like peppercorn
	

	

	macro_states = enum._resting_macrostates
	
	logger.info(f"\n\n\n\nOutput: \n {output} \n\n\n")
	return resulting_complexes, transient_complexes, enum

	



def map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,complex_count):

	
	logger.debug(f"\n\nMap Transient states\n\n")
	for complex in transient_complexes:
		logger.debug(f"{complex} {complex} ")
	logger.debug("\nAll Complexes Before\n")
	for key,complex in all_complexes.items():
		logger.debug(f"{key} {complex[0]} {complex[0].occupancy}")
		
	
	new_complexes = {}

	for key, complex in all_complexes.items():
		for t_complex in transient_complexes:
			# Check if the current complex matches a transient complex
			if complex[0] == t_complex:
				#print("\n\n________",complex[0])
				for tup, va in enum.condensation.cplx_decay_prob.items():
					
					if tup[0] == t_complex:
						macrostate = tup[1][0]
						stat_dist = enum.condensation.stationary_dist[macrostate]
						#print(tup,va)

						
						for stat_complex, value in stat_dist.items():
							print("Stat complex",stat_complex)
							##new_conc = curr_conc + (va(exit prob of transient) * value(stationary of resting state in macrostate) * t_c(concentration of transient state)) 
							try:
								if stat_complex.occupancy:
									new_conc = stat_complex.occupancy + (va * value * t_complex.occupancy) 
								else:
									new_conc = va * value * t_complex.occupancy
							except: 
									new_conc = va * value * t_complex.occupancy
							stat_complex.occupancy = new_conc

						
							# Check if the new complex is not already in all_complexes
							if not is_complex_in_all_complexes(stat_complex,all_complexes) and not is_complex_in_all_complexes(stat_complex,new_complexes):
								new_complexes["Id_" + str(complex_count)] = [stat_complex, new_conc]
								print("newconc stat complex ",stat_complex,new_conc,"complex count ",complex_count)
								complex_count += 1
							elif is_complex_in_all_complexes(stat_complex,new_complexes):
								update_complex_in_all_complexes(stat_complex,new_complexes)	

								
	
	all_complexes.update(new_complexes)
	logger.debug("\nAll Complexes After\n")
	for key,complex in all_complexes.items():
		logger.debug(f"{key} {complex[0]} {complex[0].occupancy}")
	logger.debug("\n\nEnd Mapping transient states")


def update_complex_in_all_complexes(complex,all_complexes):
	
	for key,a_complex in all_complexes.items():
		if a_complex[0] == complex:
			all_complexes[key] = [complex,complex.occupancy]



def update_occupancies(resulting_complexes,enum,all_complexes):
	"""Updates Occupancies based on the resulting condensed stationary distribution. 
		
		Args:
			resulting_complexes(list): List of the resulting resting complexes from enumeration step 
			enum(object): enumerator object which was used to calculate the resulting complexes

		Returns:
			resulting_complexes(list): Updated list of the resulting complexes with their respective concentrations
	
	"""	
	#print("\n\n_______________Begin update Occupancy___________________________\n\n")
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
	#for key, complex in all_complexes.items():
		#print(f"{key} {complex[0]._name}	{complex[0].kernel_string}	{complex[0].occupancy}")



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
		if complex in exit_probs and complex.occupancy != None:
			populations.append(complex.occupancy * exit_probs[complex])
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
	

	###
	### Here Pilsimulator Stuff
	###


	

def write_output(final_structures,d_seq,parameters = None):

	data_output = ""
	

	ts = 1
	data_output += ("\nResting Complexes after each Spacer:\n\n")
	data_output += "Transcription Step | Structure  \n"
	for x in final_structures:
		if x and x[0].kernel_string[-2] == "S": 
			for complex in x:
				data_output += f"{ts} 	{complex.kernel_string} {complex.occupancy}\n"
			ts += 1 

			data_output += "\n"
	
	
	ts = 1
	data_output += ("\n\nOnly Logic Domain pairings:\n\n")
	data_output += "Transcription Step | Structure	| Concentration \n"

	for x in final_structures:
		if x and x[0].kernel_string[-2] == "S": 
			data_output += "\n"
			for complex in x:
				kernel_string = kernel_to_dot_bracket(complex.kernel_string)
				db_struct = (only_logic_domain_struct(d_seq.split(),kernel_string))
				data_output += f"{ts}	{db_struct} 	{round(complex.occupancy,4)}\n"
			ts += 1
			#data_output += "\n"

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
	
	for key,complex in complexes.items():
		complex.occupancy = None
		complex.id = None

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
	parser.add_argument("--k-slow", type=float, help="Specify k-slow. Determines the cutoffpoint for slow reactions.", default=1)
	parser.add_argument("-v", "--verbose", action="count", default=0, help="Increase verbosity level. -v only shows peppercorn output")
	parser.add_argument("-c","--condensed", action= "store_true",default=True, help="Condense reactions into only resting complexexes. (default: False)")



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
	simulated_structures = run_sim(d_seq, parameters,args)
	

	#_____Writing_and_printing_output___# 

	output = ""

	if args.verbose > 0:
		ts = 1
		output += ("\nResting Complexes after each step:\n\n")
		output += "Transcription Step | Structure \n"
		for x in simulated_structures: 
			for complex in x:
				output += f"{ts}	{complex.kernel_string} {complex.occupancy}\n"
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
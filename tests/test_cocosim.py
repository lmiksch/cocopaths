"""pytest moduel to test cocosim
"""
import pytest
import logging, os, re, subprocess
from unittest.mock import patch
from cocopaths.cocosim import map_transient_states,calc_macro_pop,enumerate_step

from peppercornenumerator.input import read_pil 
from peppercornenumerator.enumerator import Enumerator
from peppercornenumerator.reactions import bind21
from peppercornenumerator.enumerator import BI_REACTIONS

from natsort import natsorted
#setup logger                                      
@pytest.fixture(scope="function")
def configure_logger():
	logger = logging.getLogger("copaths")  # Use the same logger name as in your main code
	console_handler = logging.StreamHandler()
	formatter = logging.Formatter("# %(levelname)s - %(message)s")
	console_handler.setFormatter(formatter)
	logger.addHandler(console_handler)
	logger.setLevel(logging.DEBUG)  # Set the desired logging level (e.g., DEBUG, INFO, etc.)
	yield
	logger.handlers = []

@pytest.fixture
def input_path():
	return os.path.join(os.path.dirname(__file__), 'test_cocosim_inputs')

def test_enumerate_step():

	pass 


def create_input(filename,occupancies,input_path):

	"""Creates complexes and structures for further testing
	

		Args():
			filename(str): filename of input file (.pil format)
			occupancies(list): list of occupancies of the corresponding complexes as in the input file 
	"""
	
	file_path = os.path.join(input_path,filename)
	with open(file_path, 'r') as file:
		file_content = file.read()
	complexes, reactions = read_pil(file_content)

	assert len(complexes) == len(occupancies) 
	complex_count = 1 

	all_complexes = {}

	

	for complex,occ in zip(complexes.values(),occupancies):
		
		complex.occupancy = occ
		all_complexes["Id_" + str(complex_count)] = [complex,occ]
		complex_count += 1

	init_cplxs = [x for x in complexes.values() if x.concentration is None or x.concentration[1] != 0]
	name_cplxs = list(complexes.values())
	enum = Enumerator(init_cplxs, reactions, named_complexes=name_cplxs)

	enum.k_slow = 1
	if bind21 in BI_REACTIONS:
		BI_REACTIONS.remove(bind21)
	enum.k_fast = 20 
	
	enum.enumerate()
	enum.condense()

	print(enum.to_pil(condensed=True,detailed = True))

	resting_complexes = [cplx for cplx in natsorted(enum.resting_complexes)] 
	
	transient_complexes = [cplx for cplx in natsorted(enum.transient_complexes)]


	return enum,resting_complexes,transient_complexes,all_complexes,complex_count

def test_map_transient_states_1(configure_logger,input_path):
	#2 transient to 1 resting state

	
	
	enum,resting_complexes,transient_complexes,all_complexes, complex_count = create_input("map_transient_1.txt",[0.5,0.5],input_path)

	print(transient_complexes)
	print(resting_complexes)
	old_all_complexes = all_complexes.copy()
	map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,complex_count)
	
	for id,complex in all_complexes.items():
		print(id,complex,complex[0].occupancy)

	
	new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
	
	assert len(new_entries) == 1 

	for complex in new_entries.values(): 
		assert complex[0].occupancy == 1 

	 

def test_map_transient_states_2(configure_logger,input_path):
	#2 transient to 1 resting state
	enum,resting_complexes,transient_complexes,all_complexes, complex_count = create_input("map_transient_2.txt",[1],input_path)

	print(transient_complexes)
	print(resting_complexes)
	old_all_complexes = all_complexes.copy()
	map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,complex_count)
	
	for id,complex in all_complexes.items():
		print(id,complex,complex[0].occupancy)

	
	new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
	
	assert len(new_entries) == 1 

	for complex in new_entries.values(): 
		assert complex[0].occupancy == 1 


def test_map_transient_states_3(configure_logger,input_path):
	#2 transient to 1 resting state
	enum,resting_complexes,transient_complexes,all_complexes, complex_count = create_input("map_transient_3.txt",[0.3,0.3,0.3],input_path)

	print(transient_complexes)
	print(resting_complexes)
	old_all_complexes = all_complexes.copy()
	map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,complex_count)
	
	for id,complex in all_complexes.items():
		print(id,complex,complex[0].occupancy)

	
	new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
	
	assert len(new_entries) == 5 
	solution_occs = [0.12803205580048818,0.0006771507479934317,0.12803205580048818,0.6331525165237,0.010106221127329963]
	for complex,occ in zip(new_entries.values(),solution_occs): 
		assert complex[0].occupancy in solution_occs
		solution_occs.remove(complex[0].occupancy)

def test_map_transient_states_4(configure_logger,input_path):
	#2 transient to 1 resting state
	enum,resting_complexes,transient_complexes,all_complexes, complex_count = create_input("map_transient_4.txt",[0.2,0.2,0.2,0.2,0.2],input_path)

	print(transient_complexes)
	print(resting_complexes)
	old_all_complexes = all_complexes.copy()
	map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,complex_count)
	
	for id,complex in all_complexes.items():
		print(id,complex,complex[0].occupancy)

	
	new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
	
	assert len(new_entries) == 1 
	solution_occs = [0.2,0.2,0.2,0.2,0.4]
	for complex,occ in zip(new_entries.values(),solution_occs): 
		assert complex[0].occupancy in solution_occs
		solution_occs.remove(complex[0].occupancy)

def test_calc_macrostate_oc_1(configure_logger,input_path):


	enum,resting_complexes,transient_complexes,all_complexes, complex_count = create_input("calc_macrostate1.txt",[0.5,0.5],input_path)
	old_all_complexes = all_complexes.copy()

	map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,complex_count)


	calc_macro_pop(enum,all_complexes,resting_complexes,complex_count)

	print(all_complexes)
	for key,complex in all_complexes.items():
		print(complex)
	
	new_entries = {key: all_complexes[key] for key in all_complexes if key not in old_all_complexes}
	
	assert len(new_entries) == 1 

	for complex in new_entries.values(): 
		assert complex[0].occupancy == 1

	pass 



def test_calc_macrostate_oc_2(configure_logger,input_path):

	enum,resting_complexes,transient_complexes,all_complexes, complex_count = create_input("calc_macrostate2.txt",[1],input_path)
	old_all_complexes = all_complexes.copy()

	
	print('\n\n\n\n complex count ',complex_count)
	map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,complex_count)

	print("\n\nResting Complexes before\n")
	for complex in resting_complexes:
		print(complex)
	calc_macro_pop(enum,all_complexes,resting_complexes,complex_count)
	print("\n\nResting Complexes after\n")
	for complex in resting_complexes:
		print(complex)


	print(all_complexes)
	for key,complex in all_complexes.items():
		print(complex)
	
	
	assert len(resting_complexes) == 1 

	for complex in resting_complexes: 
		assert complex.occupancy == 1 

	

def test_calc_macrostate_oc_3(configure_logger,input_path):

	
	enum,resting_complexes,transient_complexes,all_complexes, complex_count = create_input("calc_macrostate3.txt",[0.25,0.25,0.25,0.25],input_path)
	old_all_complexes = all_complexes.copy()
	print('\n\n\n\n complex count ',complex_count)
	map_transient_states(resting_complexes,transient_complexes,all_complexes,enum,complex_count)


	print("\n\nResting Complexes before\n")
	for complex in resting_complexes:
		print(complex,complex.occupancy)
	print("\n\nMAcrocalc")
	calc_macro_pop(enum,all_complexes,resting_complexes,complex_count)
	print("\n\nResting Complexes after\n")
	for complex in resting_complexes:
		print(complex,complex.occupancy)

	
	for key,complex in all_complexes.items():
		print(complex,complex)
	
	
	assert len(resting_complexes) == 4 

	solution_occs = [0.3740109438265814,0.0019781123468372314,0.3740109438265814,0.25]
	for complex,occ in zip(resting_complexes,solution_occs): 
		assert complex.occupancy in solution_occs
		solution_occs.remove(complex.occupancy)


def test_simulate_system():
	pass 

def simulate_condensed_reactions():

	pass 





if __name__ == '__main__':

	pytest.main()
"""pytest moduel to test cocosim
"""
import pytest
import logging, os, re, subprocess
from unittest.mock import patch
from cocopaths.cocosim import map_transient_states,map_resting_states,calc_macro_pop,enumerate_step
from peppercornenumerator.input import read_pil 
from peppercornenumerator.enumerator import Enumerator
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

	enum.k_fast = 20 
	
	enum.enumerate()
	enum.condense()

	print(enum.to_pil(condensed=True,detailed = True))

	resting_complexes = [cplx for cplx in natsorted(enum.resting_complexes)] 
	
	transient_complexes = [cplx for cplx in natsorted(enum.transient_complexes)]


	return enum,resting_complexes,transient_complexes,all_complexes,complex_count

def test_map_transient_states1(configure_logger,input_path):
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

	 

def test_map_transient_states2(configure_logger,input_path):
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


def test_map_transient_states3(configure_logger,input_path):
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






def test_calc_macrostate_oc():

	pass 


def simulate_condensed_reactions():

	pass 





if __name__ == '__main__':

	pytest.main()
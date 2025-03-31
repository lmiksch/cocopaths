import pytest
import logging, os, re
from unittest.mock import patch
from cocopaths.cocopath import (graph,                                     node,
                                     translate_acfp,
                                     input_parser,
                                     main)
from cocopaths.utils import find_connected_modules, path_to_pairtablepath, is_balanced_structure, pairtablepath_to_dotbracket

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
def sample_input(tmpdir):
    sample_content = """
        .
        ()
        .()
    """
    file_path = tmpdir.join("sample_input.txt")
    file_path.write_text(sample_content, encoding="utf-8")
    return str(file_path)

def test_input_parser_existing_file(sample_input):
    result = input_parser(sample_input)
    print(result)
    assert len(result) == 3
    assert result == ['.', '()', '.()']

def test_input_parser_empty_file(tmpdir):
    empty_file_path = tmpdir.join("empty_file.txt")
    empty_file_path.write_text("# Comment\n", encoding="utf-8")
    with pytest.raises(SystemExit):
        result = input_parser(str(empty_file_path))
    
    # Add more assertions if needed

def test_input_parser_nonexistent_file():
    non_existent_file = "nonexistent_file.txt"
    with pytest.raises(SystemExit):
        input_parser(non_existent_file)

def test_input_parser_invalid_content(tmpdir):
    invalid_content = """
        .
        abc
        ()
        .()
    """
    file_path = tmpdir.join("invalid_content.txt")
    file_path.write_text(invalid_content, encoding="utf-8")
    result = input_parser(str(file_path))
    assert result == ['.', '()', '.()']

def test_main_with_file_input(tmp_path,capsys):
    input_file = tmp_path / "input.txt"
    input_content = ".\n()\n.()"
    input_file.write_text(input_content)

    with patch("sys.argv", ["program_name", "--input", str(input_file)]):
        main()

    captured = capsys.readouterr()
    print(captured.out)

    expected_output = """\n\nInput folding path:\n['.', '()', '.()']\n\n\nResulting Domain seq: a* L0* b* Z0 c* b L0 a d* Z1 d a* L0* b* c Z2\n"""

    assert expected_output in captured.out

def test_is_balanced_structure():
    structure_1 = "()()"

    assert is_balanced_structure(structure_1)


    structure_2 = "(("

    assert not is_balanced_structure(structure_2)

    structure_3 = "()))"
    
    assert not is_balanced_structure(structure_3)




def test_get_edges(configure_logger):
    
    acfp = [[1, 0], [2, 2, 1], [3, 0, 3, 2], [4, 4, 3, 2, 1]]
    right_edges = {(1,2),(2,3),(1,4)}
    g = graph()
    g.create_nodes_from_pairtable(acfp)

    nodes = list(g.graph.keys())
    nodes.insert(0,0)
    g.create_edges(acfp,nodes)
    g.get_edges()
    

    for edge in g.edges:
        assert edge in right_edges
    pass

def test_bipartite_check(configure_logger):

    acfp = [[1, 0], [2, 2, 1], [3, 0, 3, 2], [4, 4, 3, 2, 1]]
    
    g = graph()
    g.create_nodes_from_pairtable(acfp)

    nodes = list(g.graph.keys())
    nodes.insert(0,0)
    g.create_edges(acfp,nodes)
    g.get_edges()
    g.print_nodes()
    connected_components = find_connected_modules(acfp)
    
    assert g.bipartite_check(connected_components=connected_components), f"Bipartite Check 1 failed result is {g.bipartite_check(connected_components)}"
    print("\nNodes after\n")
    g.print_nodes()
    print("\n")
    acfp = [[1,0],[2,2,1],[3,0,3,2],[4,3,4,1,2]]
    right_edges = {(1,2),(2,3),(1,4)}
    g = graph()
    g.create_nodes_from_pairtable(acfp)

    nodes = list(g.graph.keys())
    nodes.insert(0,0)
    g.create_edges(acfp,nodes)
    g.get_edges()
    
    connected_components = find_connected_modules(acfp)
    
    
    g.print_nodes()
    bipartite_result = g.bipartite_check(connected_components)
    assert bipartite_result == False,f"Bipartite Check 2 failed result is {g.bipartite_check(connected_components)}"


    acfp = [[1, 0], [2, 2, 1], [3, 0, 3, 2], [4, 4, 3, 2, 1], [5, 0, 3, 2, 5, 4], [6, 6, 3, 2, 5, 4, 1]]
    right_edges = {(1,2),(2,3),(1,4)}
    g = graph()
    g.create_nodes_from_pairtable(acfp)

    nodes = list(g.graph.keys())
    nodes.insert(0,0)
    g.create_edges(acfp,nodes)
    g.get_edges()
    
    connected_components = find_connected_modules(acfp)
    
    
    g.print_nodes()
    bipartite_result = g.bipartite_check(connected_components)
    assert bipartite_result == True,f"Bipartite Check 3 failed result is {g.bipartite_check(connected_components)}"

    acfp = [[1,0],[2,0,0],[3,0,3,2],[4,0,0,4,3],[5,0,5,4,3,2]]
    
    g = graph()
    g.create_nodes_from_pairtable(acfp)

    nodes = list(g.graph.keys())
    nodes.insert(0,0)
    g.create_edges(acfp,nodes)
    g.get_edges()
    
    connected_components = find_connected_modules(acfp)
    
    
    g.print_nodes()
    bipartite_result = g.bipartite_check(connected_components)
    assert bipartite_result == True,f"Bipartite Check 4 failed result is {g.bipartite_check(connected_components)}"


    pass


def test_pseudoknot_attack(configure_logger):

    with pytest.raises(SystemExit):
        acfp = [".","()",".()","(())","(()).","()(())"]

        domain_list = translate_acfp(acfp)
        domain_seq = ' '.join(domain_list)

        print(domain_seq)

    acfp = [".",'()',"().","(..)","().()"]

    domain_list = translate_acfp(acfp)
    domain_seq = ' '.join(domain_list)

    print(domain_seq)

    pass 

def test_detect_cycle(configure_logger):
    with pytest.raises(SystemExit):
        acfp = [".","()",".()","(())","(.())"]
        
        domain_list = translate_acfp(acfp)
        domain_seq = ' '.join(domain_list)

def test_graph(configure_logger):
    acfp = [[1, 0], [2, 2, 1], [3, 0, 3, 2], [4, 4, 3, 2, 1]]
    g = graph() 
    g.create_nodes_from_pairtable(acfp)
    assert len(g.graph) == len(acfp)

def test_1_different_substructures(configure_logger):
    with pytest.raises(SystemExit):
        acfp = [".","()","().",".()."]

        domain_list = translate_acfp(acfp)
        domain_seq = ' '.join(domain_list)

def test_2_different_substructures(configure_logger):

    with pytest.raises(SystemExit):
        acfp = [".","()","().","()()","()().","(())()"]

        domain_list = translate_acfp(acfp)
        domain_seq = ' '.join(domain_list)
        
def test_3_different_substructures(configure_logger):
     
    with pytest.raises(SystemExit):
        acfp = [".","()","().","(())","()()."]

        domain_list = translate_acfp(acfp)
        domain_seq = ' '.join(domain_list)

def test_4_different_substructures(configure_logger):
    
    with pytest.raises(SystemExit):
        acfp = [".","()",".()","()()","()().",".()()."]
        domain_list = translate_acfp(acfp)
        domain_seq = ' '.join(domain_list)


def test_5_different_substructures(configure_logger):
    
    with pytest.raises(SystemExit):
        acfp = [".","()",".()","()()","().()"]
        domain_list = translate_acfp(acfp)
        domain_seq = ' '.join(domain_list)

# List of favorable aCFPs (each is a list of dot-bracket strings)
favorable_acfps = [
    # Simple favorable path
    [".", "()", ".()"],
    # A more extended favorable path (example taken from one of your tests)
    [".", "()", ".()", "()()", ".()()"],
    

]

# List of unfavorable aCFPs
unfavorable_acfps = [
    # Example designed to trigger pseudoknot attack detection
    [".", "()", ".()", "(())", "()().", "()(())"],
    [".", "()", ".()", "().."],
    [".", "()", ".()", "..()"],
    [".", "()", "(.)", ".(.)"],
    [".", "()", ".()", "(())","(.())"],
    [".", "()", "(.)", ".(.)"],
    [".", "()", "(.)", "(.).", ".(())"],
    [".", "()", "(.)", "()()", "((.))"]
]

@pytest.mark.parametrize("acfp", favorable_acfps)
def test_acfp_favorable(acfp):
    """
    For each favorable aCFP the graph should be built successfully and a valid domain level sequence returned.
    """
    try:
        domain_list = translate_acfp(acfp)
        domain_seq = ' '.join(domain_list)
    except SystemExit as e:
        pytest.fail(f"Favorable aCFP {acfp} failed with SystemExit: {e}")
    
    # Check that a domain level sequence is produced.
    assert isinstance(domain_seq, str) and len(domain_seq) > 0, \
        f"Domain sequence for {acfp} is not valid: {domain_seq}"
    
    # Optionally, verify that each domain string contains expected markers (e.g., "L" and "S")
    #for ds in domain_seq.split():
    #    assert "L" in ds and "Z" in ds, f"Domain sequence{domain_seq} entry '{ds}' missing expected markers."

@pytest.mark.parametrize("acfp", unfavorable_acfps)
def test_acfp_unfavorable(acfp):
    """
    For each unfavorable aCFP the build_graph process should fail (raise a SystemExit).
    """
    with pytest.raises(SystemExit):
        domain_list = translate_acfp(acfp)
        domain_seq = ' '.join(domain_list)

def test_1_create_domain_seq(configure_logger):
       
    print("\n\n\nTest 1\n\n\n")
    
    acfp_1 = pairtablepath_to_dotbracket([[1,0],[2,2,1],[3,0,3,2],[4,4,3,2,1],[5,0,3,2,5,4],[6,6,3,2,5,4,1]])
    
    domain_list = translate_acfp(acfp_1)
    domain_seq_1 = ' '.join(domain_list)

    assert domain_seq_1 == 'i* e* a* L0* b* f* j* Z0 c* b L0 a d* Z1 d a* L0* b* c Z2 g* f b L0 a e h* Z3 h e* a* L0* b* f* g Z4 j f b L0 a e i Z5', f"Create Domain_seq 1 failed: Result: {domain_seq_1} Solution:  i* e* a* L0* b* f* j* Z0 c* b L0 a d* Z1 d a* L0* b* c Z2 g* f b L0 a e h* Z3 h e* a* L0* b* f* g Z4 j f b L0 a e i Z5"


    pass

def test_2_create_domain_seq(configure_logger):        
    
    print("\n\n\nTest 2\n\n\n")
    
    acfp_2 = pairtablepath_to_dotbracket([[1, 0], [2, 2, 1], [3, 2, 1, 0], [4, 4, 3, 2, 1]])

    domain_list = translate_acfp(acfp_2)
    domain_seq_2 = ' '.join(domain_list)



    assert domain_seq_2 == "c* a* L0* b* d* Z0 b L0 a Z1 L0* Z2 d b L0 a c Z3", f"Create domain_seq 2 failed Result: {domain_seq_2}, Solution c* a* L0* b* d* Z0 b L0 a Z1 L0* Z2 d b L0 a c Z3"

    pass

def test_3_create_domain_seq(configure_logger):

    print("\n\n\nTest 3\n\n\n")

    acfp_3 = [".","()","().","()()","().()","()(())"]

    domain_list = translate_acfp(acfp_3)
    domain_seq_3 = ' '.join(domain_list)

    assert domain_seq_3 == 'a* L0* b* Z0 b L0 a Z1 g* c* L1* d* h* Z2 e* d L1 c f* Z3 f c* L1* d* e Z4 h d L1 c g Z5', f"Create domain_seq 3 failed Result: {domain_seq_3}, Solution: a* L0* b* Z0 b L0 a Z1 g* c* L1* d* h* Z2 e* d L1 c f* Z3 f c* L1* d* e Z4 h d L1 c g Z5"

    pass

def test_4_create_domain_seq(configure_logger):

    print("\n\n\nTest 4\n\n\n")

    acfp_4 = [".","()","().","()()","()().","()()()","()()()."]

    domain_list = translate_acfp(acfp_4)
    domain_seq_4 = ' '.join(domain_list)
    

    assert domain_seq_4 == "a* L0* b* Z0 b L0 a Z1 c* L1* d* Z2 d L1 c Z3 e* L2* f* Z4 f L2 e Z5 L3* Z6", f"Create domain_seq 4 failed Result: {domain_seq_4},Solution: a* L0* b* Z0 b L0 a Z1 c* L1* d* Z2 d L1 c Z3 e* L2* f* Z4 f L2 e Z5 L3* Z6"

    pass 

def test_5_create_domain_seq(configure_logger):

    print("\n\n\nTest 5\n\n\n")
    
    path = [".","()",".()","(())","(()).","(()())"]
    
    domain_list = translate_acfp(path)
    domain_seq_5 = ' '.join(domain_list)

    assert domain_seq_5 == "g* e* a* L0* b* f* h* Z0 c* b L0 a d* Z1 d a* L0* b* c Z2 f b L0 a e Z3 L0* Z4 h f b L0 a e g Z5", f"Create domain_seq 5 failed Result: {domain_seq_5},Solution: ['g* e* a* L0* b* f* h* Z0 c* b L0 a d* Z1 d a* L0* b* c Z2 f b L0 a e Z3 L0* Z4 h f b L0 a e g Z5]"

    

    pass 


def test_6_create_domain_seq(configure_logger):
     

    print("\n\n\nTest 6\n\n\n")
        
    path = [".","()",".()","()()",".()()"]
    
    domain_list = translate_acfp(path)
    domain_seq_6 = ' '.join(domain_list)

    assert  domain_seq_6  == "a* L0* b* Z0 c* b L0 a d* Z1 e* d a* L0* b* c f* Z2 g* f c* b L0 a d* e h* Z3 h e* d a* L0* b* c f* g Z4", f"Create domain_seq_7 failed Result: {domain_seq_6},Solution: a* L0* b* Z0 c* b L0 a d* Z1 e* d a* L0* b* c f* Z2 g* f c* b L0 a d* e h* Z3 h e* d a* L0* b* c f* g Z4"

    pass

def test_7_create_domain_seq(configure_logger):
     

    print("\n\n\nTest 7\n\n\n")
        
    path = [".","()",".()","(())","(())."]

    domain_list = translate_acfp(path)
    domain_seq_7 = ' '.join(domain_list)

    assert  domain_seq_7  == 'e* a* L0* b* f* Z0 c* b L0 a d* Z1 d a* L0* b* c Z2 f b L0 a e Z3 L1* Z4'

    pass



def test_8_create_domain_seq(configure_logger):
    

    print("\n\n\nTest 8\n\n\n")
        
    path = [".","()","()."]
    
    domain_list = translate_acfp(path)
    domain_seq_8 = ' '.join(domain_list)

    assert  domain_seq_8  == 'a* L0* b* Z0 b L0 a Z1 L1* Z2'
    pass



if __name__ == '__main__':

    pytest.main()

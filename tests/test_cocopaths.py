import pytest
import logging, os, re
from unittest.mock import patch
from cocopaths.cocopath import (graph,
                                     node,
                                     build_graph,
                                     input_parser,
                                     main)
from cocopaths.utils import find_connected_modules, path_to_pairtablepath, is_balanced_structure

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
    assert "\n\nInput folding path:\n['.', '()', '.()']\n\n\nResulting Domain Level sequence:  L0*  S0 a L0 b S1 b* L0* a* S2\nLength of domain seq =  10\n" in captured.out 

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

        acfp_graph = build_graph(acfp)

        

        domain_seq = acfp_graph.get_domain_seq()

        print(domain_seq)

    acfp = [".",'()',"().","(..)","().()"]

    acfp_graph = build_graph(acfp)

    

    domain_seq = acfp_graph.get_domain_seq()

    print(domain_seq)

    pass 

def test_detect_cycle(configure_logger):
    with pytest.raises(SystemExit):
        acfp = [".","()",".()","(())","(.())"]
        
        acfp_graph = build_graph(acfp)

def test_graph(configure_logger):
    acfp = [[1, 0], [2, 2, 1], [3, 0, 3, 2], [4, 4, 3, 2, 1]]
    g = graph() 
    g.create_nodes_from_pairtable(acfp)
    assert len(g.graph) == len(acfp)

def test_1_different_substructures(configure_logger):
    with pytest.raises(SystemExit):
        path = [".","()","().",".()."]

        acfp_graph = build_graph(path)

def test_2_different_substructures(configure_logger):

    with pytest.raises(SystemExit):
        path = [".","()","().","()()","()().","(())()"]

        acfp_graph = build_graph(path)
        
def test_3_different_substructures(configure_logger):
     
    with pytest.raises(SystemExit):
        path = [".","()","().","(())","()()."]

        acfp_graph = build_graph(path)

def test_4_different_substructures(configure_logger):
    
    with pytest.raises(SystemExit):
        path = [".","()",".()","()()","()().",".()()."]

        acfp_graph = build_graph(path)
        acfp_graph.verify_weights(path)


def test_5_different_substructures(configure_logger):
    
    with pytest.raises(SystemExit):
        path = [".","()",".()","()()","().()"]

        acfp_graph = build_graph(path)
        acfp_graph.verify_weights(path)

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
        graph_obj = build_graph(acfp)
    except SystemExit as e:
        pytest.fail(f"Favorable aCFP {acfp} failed with SystemExit: {e}")
    
    # Check that a domain level sequence is produced.
    domain_seq = graph_obj.get_domain_seq()
    assert isinstance(domain_seq, list) and len(domain_seq) > 0, \
        f"Domain sequence for {acfp} is not valid: {domain_seq}"
    
    # Optionally, verify that each domain string contains expected markers (e.g., "L" and "S")
    for ds in domain_seq:
        assert "L" in ds and "S" in ds, f"Domain sequence entry '{ds}' missing expected markers."

@pytest.mark.parametrize("acfp", unfavorable_acfps)
def test_acfp_unfavorable(acfp):
    """
    For each unfavorable aCFP the build_graph process should fail (raise a SystemExit).
    """
    with pytest.raises(SystemExit):
        build_graph(acfp)

def test_1_create_domain_seq(configure_logger):
       
    print("\n\n\nTest 1\n\n\n")
    
    acfp_1 = [[1,0],[2,2,1],[3,0,3,2],[4,4,3,2,1],[5,0,3,2,5,4],[6,6,3,2,5,4,1]]

    acfp_graph = build_graph(acfp_1)

    
    domain_seq_1 = acfp_graph.get_domain_seq()

    assert domain_seq_1 == [' L0*  S0', 'a L0 b S1', 'b* L0* a* S2', 'c L0 d S3', 'd* L0* c* S4', ' L0  S5'], f"Create Domain_seq 1 failed: Result: {domain_seq_1} Solution:  [' m0* ', 'a m0 b', 'c* b* m0* a* e*', 'g e a m0 b c j', 'j* c* b* m0* a* e* g*']"


    pass

def test_2_create_domain_seq(configure_logger):        
    
    print("\n\n\nTest 2\n\n\n")
    
    acfp_2 = [[1, 0], [2, 2, 1], [3, 2, 1, 0], [4, 4, 3, 2, 1]]

    acfp_graph = build_graph(acfp_2)

    
    domain_seq_2 = acfp_graph.get_domain_seq()



    assert domain_seq_2 == ['a* b* L0* c* d* S0', 'c L0 b S1', ' L0*  S2', 'd c L0 b a S3'], f"Create domain_seq 2 failed Result: {domain_seq_2}, Solution ['a* b* m0* c* d*', 'c m0 a', ' m0* ', 'd c m0 b a']"

    pass

def test_3_create_domain_seq(configure_logger):

    print("\n\n\nTest 3\n\n\n")

    acfp_3 = [".","()","().","()()","().()","()(())"]

    acfp_graph = build_graph(acfp_3)

    
    domain_seq_3 = acfp_graph.get_domain_seq()


    assert domain_seq_3 == [' L0*  S0', ' L0  S1', ' L1*  S2', 'a L1 b S3', 'b* L1* a* S4', ' L1  S5'], f"Create domain_seq 3 failed Result: {domain_seq_3}, Solution: [' m0* ', ' m0 ', ' m1* ', 'a m1 b', 'b* m1* a*', ' m1 ']"

    pass

def test_4_create_domain_seq(configure_logger):

    print("\n\n\nTest 4\n\n\n")

    acfp_4 = [".","()","().","()()","()().","()()()","()()()."]


    acfp_graph = build_graph(acfp_4)

    
    domain_seq_4 = acfp_graph.get_domain_seq()
    

    assert domain_seq_4 == [' L0*  S0', ' L0  S1', ' L1*  S2', ' L1  S3', ' L2*  S4', ' L2  S5', ' L3*  S6'], f"Create domain_seq 4 failed Result: {domain_seq_4},Solution: [' m0* ', ' m0 ', ' m1* ', ' m1 ', ' m2* ', ' m2 ', ' m3* ']"

    pass 

def test_5_create_domain_seq(configure_logger):

    print("\n\n\nTest 5\n\n\n")
    
    path = [".","()",".()","(())","(()).","(()())"]
    
    acfp_graph = build_graph(path)

    
    domain_seq_5 = acfp_graph.get_domain_seq()
    

    assert domain_seq_5 == ['a* b* L0* c* d* S0', 'e L0 f S1', 'f* L0* e* S2', 'c L0 b S3', ' L0*  S4', 'd c L0 b a S5'], f"Create domain_seq 5 failed Result: {domain_seq_5},Solution: ['a* b* m0* c* d* l0', 'e m0 f l1', 'f* m0* e* l2', 'c m0 b l3', ' m0*  l4', 'd c m0 b a l5']"

    

    pass 


def test_6_create_domain_seq(configure_logger):
     

    print("\n\n\nTest 6\n\n\n")
        
    path = [".","()",".()","()()",".()()"]
    
    acfp_graph = build_graph(path)

    
    domain_seq_6 = acfp_graph.get_domain_seq()
    

    assert  domain_seq_6  == [' L0*  S0', 'a L0 b S1', 'c* b* L0* a* d* S2', 'e d a L0 b c f S3', 'f* c* b* L0* a* d* e* S4'], f"Create domain_seq_7 failed Result: {domain_seq_7},Solution: [' m0*  l0', 'a m0 b l1', 'c* b* m0* a* d* l2', 'e d a m0 b c f l3', 'f* c* b* m0* a* d* e* l4']"

    pass

def test_7_create_domain_seq(configure_logger):
     

    print("\n\n\nTest 7\n\n\n")
        
    path = [".","()",".()","(())","(())."]
    
    acfp_graph = build_graph(path)

    
    domain_seq_7 = acfp_graph.get_domain_seq()
    

    assert  domain_seq_7  == [' L0*  S0', 'a L0 b S1', 'b* L0* a* S2', ' L0  S3', ' L1*  S4']

    pass



def test_8_create_domain_seq(configure_logger):
     

    print("\n\n\nTest 8\n\n\n")
        
    path = [".","()","()."]
    
    acfp_graph = build_graph(path)

    
    domain_seq_8 = acfp_graph.get_domain_seq()
    

    assert  domain_seq_8  == [' L0*  S0', ' L0  S1', ' L1*  S2']

    pass



if __name__ == '__main__':

    pytest.main()

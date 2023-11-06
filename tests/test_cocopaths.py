import pytest
import logging
from cocopaths.cocopaths import (graph,
                                     node,
                                     build_graph)
from cocopaths.utils import find_connected_modules, path_to_pairtablepath

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


def test_get_edges(configure_logger):
    
    afp = [[1, 0], [2, 2, 1], [3, 0, 3, 2], [4, 4, 3, 2, 1]]
    right_edges = {(1,2),(2,3),(1,4)}
    g = graph()
    g.create_nodes_from_pairtable(afp)

    nodes = list(g.graph.keys())
    nodes.insert(0,0)
    g.create_edges(afp,nodes)
    g.get_edges()
    

    for edge in g.edges:
        assert edge in right_edges
    pass

def test_bipartite_check(configure_logger):

    afp = [[1, 0], [2, 2, 1], [3, 0, 3, 2], [4, 4, 3, 2, 1]]
    
    g = graph()
    g.create_nodes_from_pairtable(afp)

    nodes = list(g.graph.keys())
    nodes.insert(0,0)
    g.create_edges(afp,nodes)
    g.get_edges()
    g.print_nodes()
    connected_components = find_connected_modules(afp)
    
    assert g.bipartite_check(connected_components=connected_components), f"Bipartite Check 1 failed result is {g.bipartite_check(connected_components)}"
    print("\nNodes after\n")
    g.print_nodes()
    print("\n")
    afp = [[1,0],[2,2,1],[3,0,3,2],[4,3,4,1,2]]
    right_edges = {(1,2),(2,3),(1,4)}
    g = graph()
    g.create_nodes_from_pairtable(afp)

    nodes = list(g.graph.keys())
    nodes.insert(0,0)
    g.create_edges(afp,nodes)
    g.get_edges()
    
    connected_components = find_connected_modules(afp)
    
    
    g.print_nodes()
    bipartite_result = g.bipartite_check(connected_components)
    assert bipartite_result == False,f"Bipartite Check 2 failed result is {g.bipartite_check(connected_components)}"


    afp = [[1, 0], [2, 2, 1], [3, 0, 3, 2], [4, 4, 3, 2, 1], [5, 0, 3, 2, 5, 4], [6, 6, 3, 2, 5, 4, 1]]
    right_edges = {(1,2),(2,3),(1,4)}
    g = graph()
    g.create_nodes_from_pairtable(afp)

    nodes = list(g.graph.keys())
    nodes.insert(0,0)
    g.create_edges(afp,nodes)
    g.get_edges()
    
    connected_components = find_connected_modules(afp)
    
    
    g.print_nodes()
    bipartite_result = g.bipartite_check(connected_components)
    assert bipartite_result == True,f"Bipartite Check 3 failed result is {g.bipartite_check(connected_components)}"

    afp = [[1,0],[2,0,0],[3,0,3,2],[4,0,0,4,3],[5,0,5,4,3,2]]
    
    g = graph()
    g.create_nodes_from_pairtable(afp)

    nodes = list(g.graph.keys())
    nodes.insert(0,0)
    g.create_edges(afp,nodes)
    g.get_edges()
    
    connected_components = find_connected_modules(afp)
    
    
    g.print_nodes()
    bipartite_result = g.bipartite_check(connected_components)
    assert bipartite_result == True,f"Bipartite Check 4 failed result is {g.bipartite_check(connected_components)}"





    pass

def test_pseudoknot_attack(configure_logger):

    with pytest.raises(SystemExit):
        afp = [".","()",".()","(())","(()).","()(())"]

        afp_graph = build_graph(afp)

        

        domain_seq = afp_graph.get_domain_seq()

        print(domain_seq)

    pass 

def test_detect_cycle(configure_logger):
    with pytest.raises(SystemExit):
        afp = [".","()",".()","(())","(.())"]
        
        afp_graph = build_graph(afp)

def test_graph(configure_logger):
    afp = [[1, 0], [2, 2, 1], [3, 0, 3, 2], [4, 4, 3, 2, 1]]
    g = graph() 
    g.create_nodes_from_pairtable(afp)
    assert len(g.graph) == len(afp)

def test_1_different_substructures(configure_logger):
    with pytest.raises(SystemExit):
        path = [".","()","().",".()."]

        afp_graph = build_graph(path)

def test_2_different_substructures(configure_logger):

    with pytest.raises(SystemExit):
        path = [".","()","().","()()","()().","(())()"]

        afp_graph = build_graph(path)
        
def test_3_different_substructures(configure_logger):
     
    with pytest.raises(SystemExit):
        path = [".","()","().","(())","()()."]

        afp_graph = build_graph(path)

def test_1_create_domain_seq(configure_logger):
       
    print("\n\n\nTest 1\n\n\n")
    
    afp_1 = [[1,0],[2,2,1],[3,0,3,2],[4,4,3,2,1],[5,0,3,2,5,4],[6,6,3,2,5,4,1]]

    afp_graph = build_graph(afp_1)

    
    domain_seq_1 = afp_graph.get_domain_seq()

    assert domain_seq_1 == [' m0*  l0', 'a m0 b l1', 'b* m0* a* l2', 'c m0 d l3', 'd* m0* c* l4', ' m0  l5'], f"Create Domain_seq 1 failed: Result: {domain_seq_1} Solution:  [' m0* ', 'a m0 b', 'c* b* m0* a* e*', 'g e a m0 b c j', 'j* c* b* m0* a* e* g*']"


    pass

def test_2_create_domain_seq(configure_logger):        
    
    print("\n\n\nTest 2\n\n\n")
    
    afp_2 = [[1, 0], [2, 2, 1], [3, 2, 1, 0], [4, 4, 3, 2, 1]]

    afp_graph = build_graph(afp_2)

    
    domain_seq_2 = afp_graph.get_domain_seq()



    assert domain_seq_2 == ['a* b* m0* c* d* l0', 'c m0 b l1', ' m0*  l2', 'd c m0 b a l3'], f"Create domain_seq 2 failed Result: {domain_seq_2}, Solution ['a* b* m0* c* d*', 'c m0 a', ' m0* ', 'd c m0 b a']"

    pass

def test_3_create_domain_seq(configure_logger):

    print("\n\n\nTest 3\n\n\n")

    afp_3 = [".","()","().","()()","().()","()(())"]

    afp_graph = build_graph(afp_3)

    
    domain_seq_3 = afp_graph.get_domain_seq()


    assert domain_seq_3 == [' m0*  l0', ' m0  l1', ' m1*  l2', 'a m1 b l3', 'b* m1* a* l4', ' m1  l5'], f"Create domain_seq 3 failed Result: {domain_seq_3}, Solution: [' m0* ', ' m0 ', ' m1* ', 'a m1 b', 'b* m1* a*', ' m1 ']"

    pass

def test_4_create_domain_seq(configure_logger):

    print("\n\n\nTest 4\n\n\n")

    afp_4 = [".","()","().","()()","()().","()()()","()()()."]


    afp_graph = build_graph(afp_4)

    
    domain_seq_4 = afp_graph.get_domain_seq()
    

    assert domain_seq_4 == [' m0*  l0', ' m0  l1', ' m1*  l2', ' m1  l3', ' m2*  l4', ' m2  l5', ' m3*  l6'], f"Create domain_seq 4 failed Result: {domain_seq_4},Solution: [' m0* ', ' m0 ', ' m1* ', ' m1 ', ' m2* ', ' m2 ', ' m3* ']"

    pass 

def test_5_create_domain_seq(configure_logger):

    print("\n\n\nTest 5\n\n\n")
    
    path = [".","()",".()","(())","(()).","(()())"]
    
    afp_graph = build_graph(path)

    
    domain_seq_5 = afp_graph.get_domain_seq()
    

    assert domain_seq_5 == ['a* b* m0* c* d* l0', 'e m0 f l1', 'f* m0* e* l2', 'c m0 b l3', ' m0*  l4', 'd c m0 b a l5'], f"Create domain_seq 5 failed Result: {domain_seq_5},Solution: ['a* b* m0* c* d* l0', 'e m0 f l1', 'f* m0* e* l2', 'c m0 b l3', ' m0*  l4', 'd c m0 b a l5']"

    

    pass 

def test_6_create_domain_seq(configure_logger):
     

    print("\n\n\nTest 6\n\n\n")
        
    path = [".","..",".()","(())","(()).","(()())"]
    
    afp_graph = build_graph(path)

    
    domain_seq_6 = afp_graph.get_domain_seq()
    

    assert  domain_seq_6  == ['a* b* m0* c* d* l0', ' m1*  l1', ' m1  l2', 'c m0 b l3', ' m0*  l4', 'd c m0 b a l5'], f"Create domain_seq_6 failed Result: {domain_seq_6},Solution: ['a* b* m0* c* d* l0', ' m1*  l1', ' m1  l2', 'c m0 b l3', ' m0*  l4', 'd c m0 b a l5']"

    pass

def test_7_create_domain_seq(configure_logger):
     

    print("\n\n\nTest 6\n\n\n")
        
    path = [".","()",".()","()()",".()()"]
    
    afp_graph = build_graph(path)

    
    domain_seq_7 = afp_graph.get_domain_seq()
    

    assert  domain_seq_7  == [' m0*  l0', 'a m0 b l1', 'c* b* m0* a* d* l2', 'e d a m0 b c f l3', 'f* c* b* m0* a* d* e* l4'], f"Create domain_seq_7 failed Result: {domain_seq_7},Solution: [' m0*  l0', 'a m0 b l1', 'c* b* m0* a* d* l2', 'e d a m0 b c f l3', 'f* c* b* m0* a* d* e* l4']"

    pass
 


if __name__ == '__main__':

    pytest.main()
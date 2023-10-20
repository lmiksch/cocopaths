import pulp
import ast


def inequality_solver(inequalities, edges):
    """Takes a list of inequalities and searches for a numerical integer solution
    Args:
        inequalities(list): entries in list correspond to inequalities in the form of a string. Inequalities must have the form X > Y otherwise it won't work
        edges(set): names of the edges
    Return:
        inequality_solutions(dict): Dictionary of the numerical solutions
    """
    problem = pulp.LpProblem("Inequalities_Solver", pulp.LpMinimize)

    edge_vars = {}

    print("\nLP Solver \n")
    #pulp.LpSolverDefault.msg = 0 #silences message output
    print("\n Inequalities edges")
    print(inequalities)
    print(edges)
    print("\n\n")
    # Define Variables
    for edge in edges:
        edge_vars[edge] = pulp.LpVariable(str(edge), lowBound=1, cat="Integer")

    for inequality in inequalities:
        parts = inequality.split(">")
        lhs = parts[0].strip()
        rhs = parts[1].strip()
        
        lhs_split = lhs.split("+")
        rhs_split = rhs.split("+")
        
        lhs_split = [ast.literal_eval(edge) for edge in lhs_split]
        rhs_split = [ast.literal_eval(edge) for edge in rhs_split]
        # Create the left-hand side and right-hand side expressions
        lhs_expr = pulp.lpSum(edge_vars[edge] for edge in lhs_split)
        rhs_expr = pulp.lpSum(edge_vars[edge] for edge in rhs_split)

        # Add the constraint to the problem
        problem += lhs_expr - rhs_expr >= 1

    problem.solve()
    
    inequalities_solutions = {}

    if problem.status == 1:
        print("Solutions")
        # Display the integer solutions
        for edge in edges:
            if edge_vars[edge].varValue != None:
                inequalities_solutions[edge] = int(edge_vars[edge].varValue)
                print(edge, int(edge_vars[edge].varValue))
            
            else:
                inequalities_solutions[edge] = 1
                print("else",edge,"1")

    else:
        print("No feasible integer solution found.")


    return inequalities_solutions

if __name__ == "__main__":

    
    #edges = {"12", "23", "34"}
    #ineq = ["23 > 12", "34 > 23", "34 + 12 > 23"]

    #edges = {"AB","BC","CD","AD"}
    #ineq = ["BC > AB", "AD > AB","AD + BC > AB + CD"]

    #edges = {"AB","BC","AD","DE","CF","EF"}
    #ineq = ["BC > AB","AD > DE","AD + BC > AB","CF > BC", "AB > AD","AB + CF + DE > AD + BC + EF"]

    #edges = {"AB","BC","DE","CD","AF"}
    #ineq = ["AB > BC","CD > DE","AB + CD > BC + DE","AF > AB","BC > CD"]

    #edges = {"AB","BC","BE","EF","DE","CD"}
    #ineq = ["BC >AB","CD > BC","BE > BC","BE > AB","EF > DE"]

    
    #edges = {"12","23","34","25","16","14"}
    #ineq = ["23 > 12","14 > 34","16 > 12","16 > 14"]


    edges = {"12","23","34","28","15","14","47","56"}
    ineq = ["23 > 12","34 > 23","14 > 12","14 > 34", "15 > 12","56 > 15","47 > 14","47 > 34","28 > 23","28 > 12","15 > 14"]

    inequality_solver(inequalities=ineq, edges=edges)
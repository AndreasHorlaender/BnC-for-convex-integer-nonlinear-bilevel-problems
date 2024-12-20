from gurobipy import GRB

def check_callback(check, where):
    if where == GRB.Callback.MIPSOL:
        # maximization problem -> best_obj is largest lower bound and obj_bound smallest upper bound
        obj_lower_bound = check.cbGet(GRB.Callback.MIPSOL_OBJBST) # incumbent
        obj_upper_bound = check.cbGet(GRB.Callback.MIPSOL_OBJBND) 
        
        x = check.cbGetSolution(check._x)
        y = check.cbGetSolution(check._y)
        solve_subproblems_to_global_optimality = check._solve_subproblems_to_global_optimality
        
        #print("x = ", x, " y = ", y)
        
        #print("STATUS = ", check.status)
        #print("OBJ_UPPER_BOUND = ", obj_upper_bound)
        #print("OBJ_LOWER_BOUND = ", obj_lower_bound)
        
        if obj_upper_bound <= 0: # then there is no point with positive objective value
            #print("Subproblem gets terminated due to non-positive upper bound.")
            check.terminate()
            
        if not solve_subproblems_to_global_optimality:    
            if obj_lower_bound >= 1e-4: # then our incumbent solution already violates the hyperplane
                #print("Subproblem gets terminated due to positive lower bound. The incumbent solution is taken.")
                check.terminate()
                

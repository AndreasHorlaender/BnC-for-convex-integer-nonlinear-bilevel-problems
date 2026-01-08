import numpy as np
import time
from docplex.mp.model import Model
from refinement_procedure import refinement_procedure
import global_vars
from docplex.mp.model_reader import read_model


def modify_LL_constrs(self):  
    self.ll.parameters.timelimit = max(0, self.get_end_time() - self.get_time()) 
    #print("LL")
    #### DOcplexException: CPLEX cannot modify quadratic constraint
    #### -> build the quadratic constraint from scratch
    if global_vars.first_cnstr_quad:
        for i in range(self.b_len):
            rhs_i = self.b[i] - np.dot(self.C[i], self.x_j)
            if i >= 1:
                # modify linear cnstrs
                self.ll.get_constraint_by_name("ll_lin_cnstr"+str(i)).rhs=rhs_i
            else:
                if global_vars.quad_ll_cnstr:
                    # add quad cnstr
                    self.ll.add_constraint(self.ll.scal_prod([self.y_ll[j] for j in range(self.y_len)], [self.D[0][j] for j in range(self.y_len)]) 
                                           +  0.5 * (self.y_ll @ self.Q_cnstr @ self.y_ll) <= rhs_i, ctname="ll_quad_cnstr") 
                else:
                    # if we have no quad cnstr, the var first_cnstr_quad is set to True so that we only have to modify 1 part of the code
                    self.ll.add_constraint(self.ll.scal_prod([self.y_ll[j] for j in range(self.y_len)], [self.D[0][j] for j in range(self.y_len)]) 
                                           <= rhs_i, ctname="ll_quad_cnstr") 
    else:
        for i in range(self.b_len):
            rhs_i = self.b[i] - np.dot(self.C[i], self.x_j)
            if i <= self.b_len-2:
                # modify linear cnstrs
                self.ll.get_constraint_by_name("ll_lin_cnstr"+str(i)).rhs=rhs_i
            else:
                # add quad cnstr
                self.ll.add_constraint(self.ll.scal_prod([self.y_ll[j] for j in range(self.y_len)], [self.D[self.b_len-1][j] for j in range(self.y_len)]) 
                                       +  0.5 * (self.y_ll @ self.Q_cnstr @ self.y_ll) <= rhs_i, ctname="ll_quad_cnstr")             
             
       
def store_bilevel_feasible_point(self):
    found_bilevel_feas_point = False
    ll_obj_at_y_j = self.d_y @ self.y_j + 0.5 * (self.y_j @ self.Q_obj @ self.y_j)
    self.x_bf = self.x_j
    self.y_bf = []
    
    if abs(ll_obj_at_y_j - self.ll_obj_value) <= 1e-4:
        found_bilevel_feas_point = True
        self.node_sol_is_bilevel_feas = True
        self.y_bf = self.y_j
        
    elif (self.A @ self.x_j + self.B @ self.y_hat - self.a <= 1e-4).all():
        # (x_j, y_hat) satisfies UL cnstrs and, hence, is bilevel-feasible
        found_bilevel_feas_point = True
        self.y_bf = self.y_hat    
    
    # in any case we can use refinement procedure to improve y
    if global_vars.use_refinement_procedure:
        start_time_ref = time.time()
        refinement_has_solution, self.y_bf = refinement_procedure(self)
        global_vars.time_refinement += time.time() - start_time_ref
        if refinement_has_solution:
            found_bilevel_feas_point = True
            
    if found_bilevel_feas_point:
        # we store the list [x, y, x_bin, y_bin]
        var_names = []
        bil_feas_point = []
        
        # get x_j and y_bf var_names and values
        for i in range(self.x_len):
            var_names.append("x_" + str(i))
            bil_feas_point.append(int(self.x_bf[i]))
            
        for i in range(self.y_len):
            var_names.append("y_" + str(i))
            bil_feas_point.append(int(self.y_bf[i]))
        
        # get x_bin and y_bin var_names and values    
        for i in range(self.x_len):
            # get binary extension of x[i]
            aux_list_x_bin = [0]*self.number_of_x_binaries[i]
            if self.x_bf[i] <= -1e-3: # if the number is negative, then the last binary takes the value of 1
                aux_list_x_bin[-1] = 1
            k = 1
            while not bin(int(self.x_bf[i]))[-k] == 'b': # binary extension is something like 0b110 = 6 and reads from right to left
                if bin(int(self.x_bf[i]))[-k] == '1':
                    aux_list_x_bin[k-1] = 1
                elif bin(int(self.x_bf[i]))[-k] == '0':
                    aux_list_x_bin[k-1] = 0
                else:
                    raise Exception("Binary extension went wrong!") 
                k += 1  
                
            for j in range(self.number_of_x_binaries[i]):
                var_names.append("x_bin_" + str(i) + "_" + str(j))
                bil_feas_point.append(aux_list_x_bin[j])
        
        if not global_vars.use_INGC_only_on_x:        
            for i in range(self.y_len):
                # get binary extension of y[i]
                aux_list_y_bin = [0]*self.number_of_y_binaries[i]
                if self.y_bf[i] <= -1e-3: # if the number is negative, then the last binary takes the value of 1
                    aux_list_y_bin[-1] = 1
                k = 1
                while not bin(int(self.y_bf[i]))[-k] == 'b': # binary extension is something like 0b110 = 6 and reads from right to left
                    if bin(int(self.y_bf[i]))[-k] == '1':
                        aux_list_y_bin[k-1] = 1
                    elif bin(int(self.y_bf[i]))[-k] == '0':
                        aux_list_y_bin[k-1] = 0
                    else:
                        raise Exception("Binary extension went wrong!") 
                    k += 1  
                    
                for j in range(self.number_of_y_binaries[i]):
                    var_names.append("y_bin_" + str(i) + "_" + str(j))      
                    bil_feas_point.append(aux_list_y_bin[j])
    
        ##print("BILEVEL_FEASIBLE_POINT = ", var_names, bil_feas_point)
        # get UL objective of the bilevel-feasible point
        UL_obj = self.c_x @ self.x_bf + self.c_y @ self.y_bf
        
        global_vars.found_new_bil_feas_point = True
        global_vars.var_names = var_names
        global_vars.bil_feas_point = bil_feas_point
        global_vars.bil_feas_UL_obj_val = UL_obj
        
    else:
        global_vars.found_new_bil_feas_point = False
    
def solve_LL(self):  
    if self.get_time() >= self.get_end_time():
        self.abort()
        
    self.node_sol_is_bilevel_feas = False
    ll_modify_start_time = time.time()
    modify_LL_constrs(self)
    global_vars.ll_modify_time += time.time() - ll_modify_start_time
    ll_solve_start_time = time.time()
    #print("solve LL")
    self.ll.solve(clean_before_solve=True) 
    #print("LL_solved")
    global_vars.ll_solve_time += time.time() - ll_solve_start_time
    #self.ll.export_as_lp("/home/horlaender/Schreibtisch/CPLEX_model/src_HPC/ll.lp")
    #print(self.ll.solve_details.status)
    
    if self.ll.solve_details.status == "integer optimal with unscaled infeasibilities":
        #print(self.ll.solve_details.status)
        self.ll.parameters.read.scale = 1
        self.ll.solve(clean_before_solve=True) 
        #print(self.ll.solve_details.status)
        self.ll.parameters.read.scale = 0
    
    self.ll_obj_value = self.ll.objective_value
    
    self.y_hat = np.zeros(self.y_len)
    for i in range(self.y_len):
        self.y_hat[i] = np.round(self.y_ll[i].solution_value)
    
    # remove quad cnstr
    self.ll.remove_constraint("ll_quad_cnstr")
        
    store_bilevel_feasible_point(self)
    
    return(self.node_sol_is_bilevel_feas, self.y_hat)

import numpy as np
import global_vars
import time

    
def modify_scoop_constrs(self): 
    self.scoop.parameters.timelimit = max(0, self.get_end_time() - self.get_time()) 
    #print("SCOOP")
    #'''
    ### BOUNDS FOR VARS!
    for i in range(self.y_len):
        self.scoop.get_var_by_name("delta_y_"+str(i)).lb = np.array(self.y_bounds)[i,0] - np.array(self.y_j)[i]
        self.scoop.get_var_by_name("delta_y_"+str(i)).ub = np.array(self.y_bounds)[i,1] - np.array(self.y_j)[i]
    
    
    # OLD
    for i in range(self.y_len):
        self.scoop.get_constraint_by_name("aux"+str(i)).rhs=self.y_j[i]
        
    if global_vars.first_cnstr_quad:    
        # modify linear scoop cnstrs
        for i in range(1, self.b_len):
            rhs_i = self.b[i] - np.dot(self.C[i], self.x_j) - np.dot(self.D[i], self.y_j)
            self.scoop.get_constraint_by_name("scoop_lin_cnstr"+str(i)).rhs=rhs_i    
           
        ### last cnstr has quad term
        if self.b_len >= 1:
            if global_vars.quad_ll_cnstr:
                self.scoop.add_constraint(self.scoop.sum(self.C[0][j] * self.x_j[j] for j in range(self.x_len)) 
                                          + self.scoop.sum(self.D[0][j] * self.aux_var[j] for j in range(self.y_len)) 
                                          - self.b[0]  +  0.5 * self.scoop.sum(self.scoop.sum(self.aux_var[i] * self.Q_cnstr[i][j] 
                                          * self.aux_var[j] for j in range(self.y_len)) for i in range(self.y_len)) 
                                          + self.s[1] <= 0, ctname="scoop_quad_cnstr")    
            else:
                self.scoop.add_constraint(self.scoop.sum(self.C[0][j] * self.x_j[j] for j in range(self.x_len)) 
                                          + self.scoop.sum(self.D[0][j] * self.aux_var[j] for j in range(self.y_len)) 
                                          - self.b[0] + self.s[1] <= 0, ctname="scoop_quad_cnstr") 
            
    else:
        # modify linear scoop cnstrs
        for i in range(self.b_len-1):
            rhs_i = self.b[i] - np.dot(self.C[i], self.x_j) - np.dot(self.D[i], self.y_j)
            self.scoop.get_constraint_by_name("scoop_lin_cnstr"+str(i)).rhs=rhs_i    
        
            
        ### last cnstr has quad term
        if self.b_len >= 1:
            self.scoop.add_constraint(self.scoop.sum(self.C[self.b_len-1][j] * self.x_j[j] for j in range(self.x_len)) 
                                      + self.scoop.sum(self.D[self.b_len-1][j] * self.aux_var[j] for j in range(self.y_len)) 
                                      - self.b[self.b_len-1]  +  0.5 * self.scoop.sum(self.scoop.sum(self.aux_var[i] * self.Q_cnstr[i][j] 
                                      * self.aux_var[j] for j in range(self.y_len)) for i in range(self.y_len)) 
                                      + self.s[self.b_len] <= 0, ctname="scoop_quad_cnstr")            
       
        
    ### obj cnstr 
    self.scoop.add_constraint(self.scoop.scal_prod([self.Delta_y[j] for j in range(self.y_len)], [self.d_y[j] for j in range(self.y_len)]) 
                              + 0.5 * (self.aux_var @ self.Q_obj @ self.aux_var)
                              - 0.5 * (self.y_j @ self.Q_obj @ self.y_j)
                              + self.s[0] <= 0, ctname="scoop_quad_obj_cnstr")   
    
    
   
def solve_scoop(self):  
    if self.get_time() >= self.get_end_time():
        self.abort()
    #print("solve scoop")
    self.use_no_good_cut = False
    global_vars.node_sol_on_boundary_of_bf_set = False
    
    scoop_modify_start_time = time.time()
    modify_scoop_constrs(self)
    global_vars.scoop_modify_time += time.time() - scoop_modify_start_time
    scoop_solve_start_time = time.time()
    #self.scoop.export_as_lp("../src/scoop.lp")    
    self.scoop.solve(clean_before_solve=True)
    global_vars.scoop_solve_time += time.time() - scoop_solve_start_time
    
    
    self.active_cnstrs_of_bf_set = []
    
    if self.scoop.solve_details.status == 'integer optimal solution' or self.scoop.solve_details.status == 'integer optimal, tolerance':
        #print(self.scoop.solve_details.status)
        self.Delta_y_opt = np.zeros(self.y_len)
        for i in range(self.y_len):
            self.Delta_y_opt[i] = np.round(self.Delta_y[i].solution_value)
            
        # get the constraints of the bf-set that are satisfied with equality
        self.active_cnstrs_of_bf_set = [i for i, var in enumerate(self.s) if var.solution_value <= 1e-4]
        if self.active_cnstrs_of_bf_set:
            self.use_no_good_cut = True
            global_vars.node_sol_on_boundary_of_bf_set = True
        
       
    else:    
        # feasible with unscaled feasibilities
        self.Delta_y_opt = np.zeros(self.y_len)
        self.use_no_good_cut = True
        
    # remove quad cnstrs
    self.scoop.remove_constraint("scoop_quad_cnstr")
    self.scoop.remove_constraint("scoop_quad_obj_cnstr")
    
    
    
    return(self.Delta_y_opt, self.use_no_good_cut, self.active_cnstrs_of_bf_set)

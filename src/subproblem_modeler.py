import numpy as np
import time
import cplex
from docplex.mp.model import Model
import gurobipy as gp
from gurobipy import GRB, LinExpr
import global_vars

def get_data(self):
    self.x_len = self.packed_data[0]
    self.y_len = self.packed_data[1]
    self.x_bounds = self.packed_data[2]
    self.y_bounds = np.array(self.packed_data[3])
    self.number_of_x_binaries = self.packed_data[4]
    self.number_of_y_binaries = self.packed_data[5]
    
    self.a_len =self.packed_data[6]
    self.b_len = self.packed_data[7]
    
    self.c_x = np.array(self.packed_data[8])
    self.c_y = np.array(self.packed_data[9])
    self.d_y = np.array(self.packed_data[10])
    
    self.Q_obj = np.array(self.packed_data[11])
    self.Q_cnstr = np.array(self.packed_data[12])
    self.A = np.array(self.packed_data[13])
    self.B = np.array(self.packed_data[14])
    self.C = np.array(self.packed_data[15])
    self.D = np.array(self.packed_data[16])
    
    self.a = np.array(self.packed_data[17])
    self.b = np.array(self.packed_data[18])
    self.a_sense = np.array(self.packed_data[19])
    self.b_sense = np.array(self.packed_data[20])
    
    self.use_smart_dc_verification = self.packed_data[21]
    self.solve_subproblems_to_global_optimality = self.packed_data[22]
    
    
################ LOWER LEVEL #######################

def set_LL_parameters(self):
    self.ll = Model()
    
    self.ll.parameters.parallel = 1
    self.ll.parameters.threads = 1
    
def add_LL_vars(self):
    # x_j is fixed but we need to treat it as a variable to avoid constraints like 0 <= 0 which cplex does not
    # recognize as constraint but only as bool and hence creates an error
    self.y_ll = self.ll.integer_var_list(keys = self.y_len, lb=np.array(self.y_bounds)[:,0], ub=np.array(self.y_bounds)[:,1], name = "y_ll")
    
    
def add_LL_lin_constrs(self):   
    # add linear constraints
    if global_vars.first_cnstr_quad:
        self.ll.add_constraints([self.ll.scal_prod([self.y_ll[j] for j in range(self.y_len)], [self.D[i][j] for j in range(self.y_len)]) <= 0 for i in range(1, self.b_len)], ["ll_lin_cnstr"+str(i) for i in range(1, self.b_len)])
    else:
        self.ll.add_constraints([self.ll.scal_prod([self.y_ll[j] for j in range(self.y_len)], [self.D[i][j] for j in range(self.y_len)]) <= 0 for i in range(self.b_len-1)], ["ll_lin_cnstr"+str(i) for i in range(self.b_len-1)])
        
def add_LL_obj(self):   
    self.ll.set_objective('min', self.ll.scal_prod([self.y_ll[i] for i in range(self.y_len)], [self.d_y[i] for i in range(self.y_len)]) +  0.5 * (self.y_ll @ self.Q_obj @ self.y_ll))  
    
def build_lower_level(self):
    ll_model_start_time = time.time()
    set_LL_parameters(self)
    add_LL_vars(self)
    add_LL_lin_constrs(self)
    add_LL_obj(self)
    global_vars.ll_model_time = time.time() - ll_model_start_time
 
    
################ SCOOP PROBLEM #############################################    

def set_scoop_parameters(self):
    self.scoop = Model()
    
    self.scoop.parameters.parallel = 1
    self.scoop.parameters.threads = 1
    
def add_scoop_vars(self):
    self.Delta_y = self.scoop.integer_var_list(keys = self.y_len, lb=-cplex.infinity, name='delta_y')
    
    # the aux_vars need to have the same boundaries as the y vars, otherwise the problem may be anbounded
    self.aux_var = self.scoop.integer_var_list(keys = self.y_len)
    
    
    self.t = self.scoop.continuous_var(lb=-cplex.infinity, name='t')
    self.s = self.scoop.continuous_var_list(keys = self.b_len+1, lb=0, name='s')
    
def add_scoop_lin_constrs(self):  
    self.scoop.add_constraints([self.t <= self.s[i] for i in range(self.b_len+1)])
    if global_vars.first_cnstr_quad:
        self.scoop.add_constraints([self.scoop.scal_prod([self.Delta_y[j] for j in range(self.y_len)], [self.D[i][j] for j in range(self.y_len)]) + self.s[i+1] <= 0 for i in range(1, self.b_len)], ["scoop_lin_cnstr"+str(i) for i in range(1, self.b_len)])
    else:
        self.scoop.add_constraints([self.scoop.scal_prod([self.Delta_y[j] for j in range(self.y_len)], [self.D[i][j] for j in range(self.y_len)]) + self.s[i+1] <= 0 for i in range(self.b_len-1)], ["scoop_lin_cnstr"+str(i) for i in range(self.b_len-1)])
        
    self.scoop.add_constraints([self.aux_var[i] - self.Delta_y[i] == 0 for i in range(self.y_len)], ["aux"+str(i) for i in range(self.y_len)])
        
def add_scoop_obj(self):   
    self.scoop.set_objective('max', self.t)
    
def build_scoop_problem(self):
    scoop_model_start_time = time.time()
    set_scoop_parameters(self)
    add_scoop_vars(self)
    add_scoop_lin_constrs(self)
    add_scoop_obj(self)
    global_vars.scoop_model_time = time.time() - scoop_model_start_time
    

############## CVPs #######################################################
def set_cvp_parameters(self):
    self.cvp = gp.Model()
    self.cvp.Params.LogToConsole = 0  #no Gurobi output
    self.cvp.Params.Threads = 1

def add_cvp_vars(self):
    # change var bounds
    self.x_cvp = self.cvp.addMVar(self.x_len, vtype=GRB.INTEGER, name='x_cvp')
    self.y_cvp = self.cvp.addMVar(self.y_len, vtype=GRB.INTEGER, name='y_cvp')
   
    # binary variables to represent x and y
    self.x_bin_var_list = []
    self.y_bin_var_list = []
    for i in range(self.x_len):
        self.x_bin_var_list.append(self.cvp.addMVar(self.number_of_x_binaries[i], vtype=GRB.BINARY, name='x_bin_cvp_' + str(i)))
    for i in range(self.y_len):    
        self.y_bin_var_list.append(self.cvp.addMVar(self.number_of_y_binaries[i], vtype=GRB.BINARY, name='y_bin_cvp_' + str(i)))
    
    self.cvp.update()
    
def add_cvp_constrs(self):  
    
    quad_expr = gp.QuadExpr()
    for i in range(self.y_len):
        for j in range(self.y_len):
            if self.Q_cnstr[i, j] != 0:  # Only add non-zero terms for efficiency
                quad_expr.addTerms(self.Q_cnstr[i, j], self.y_cvp.tolist()[i], self.y_cvp.tolist()[j])
    
    #### HPR constraints
    self.cvp.addConstrs((LinExpr([self.A[i][j] for j in range(self.x_len)], self.x_cvp.tolist()) + 
                           LinExpr([self.B[i][j] for j in range(self.y_len)], self.y_cvp.tolist()) - 
                           self.a[i] <= 0 for i in range(self.a_len))) 
    
    if global_vars.first_cnstr_quad:
        self.cvp.addConstrs((LinExpr([self.C[i][j] for j in range(self.x_len)], self.x_cvp.tolist()) + 
                               LinExpr([self.D[i][j] for j in range(self.y_len)], self.y_cvp.tolist()) - 
                               self.b[i] <= 0 for i in range(1, self.b_len)))
        
        if self.b_len >= 1:
            if global_vars.quad_ll_cnstr:
                self.cvp.addQConstr(0.5 * quad_expr +
                                      LinExpr([self.C[0][j] for j in range(self.x_len)], self.x_cvp.tolist()) + 
                                      LinExpr([self.D[0][j] for j in range(self.y_len)], self.y_cvp.tolist()) - 
                                      self.b[0] <= 0) 
            else:
                self.cvp.addConstr(LinExpr([self.C[0][j] for j in range(self.x_len)], self.x_cvp.tolist()) + 
                                     LinExpr([self.D[0][j] for j in range(self.y_len)], self.y_cvp.tolist()) - 
                                     self.b[0] <= 0) 
    else:
        self.cvp.addConstrs((LinExpr([self.C[i][j] for j in range(self.x_len)], self.x_cvp.tolist()) + 
                           LinExpr([self.D[i][j] for j in range(self.y_len)], self.y_cvp.tolist()) - 
                           self.b[i] <= 0 for i in range(self.b_len-1)))
    
        if self.b_len >= 1:
            self.cvp.addQConstr(0.5 * quad_expr+
                                  LinExpr([self.C[self.b_len-1][j] for j in range(self.x_len)], self.x_cvp.tolist()) + 
                                  LinExpr([self.D[self.b_len-1][j] for j in range(self.y_len)], self.y_cvp.tolist()) - 
                                  self.b[self.b_len-1] <= 0) 
        
    self.cvp.update()
    
    #### binary expansion
    self.cvp.addConstrs((self.x_cvp[i] == - 2**(self.number_of_x_binaries[i] -1) * self.x_bin_var_list[i][self.number_of_x_binaries[i]-1] 
                           + LinExpr([2**k for k in range(self.number_of_x_binaries[i] - 1)], [self.x_bin_var_list[i].tolist()[k] for k in range(self.number_of_x_binaries[i] - 1)]) for i in range(self.x_len)))
    
    self.cvp.addConstrs((self.y_cvp[i] == - 2**(self.number_of_y_binaries[i] -1) * self.y_bin_var_list[i][self.number_of_y_binaries[i]-1] 
                           + LinExpr([2**k for k in range(self.number_of_y_binaries[i] - 1)], [self.y_bin_var_list[i].tolist()[k] for k in range(self.number_of_y_binaries[i] - 1)]) for i in range(self.y_len)))
    
    self.cvp.update()
    
def build_CVP_model(self):
    cvp_model_start_time = time.time()
    set_cvp_parameters(self)
    add_cvp_vars(self)
    add_cvp_constrs(self)    
    global_vars.cvp_model_time = time.time() - cvp_model_start_time



############################ snp #########################################################################################################
def set_snp_parameters(self):
    self.snp = gp.Model()
    self.snp.Params.LogToConsole = 0  #no Gurobi output 
    self.snp.Params.Threads = 1

def add_snp_vars(self):
    # change var bounds
    self.x_snp = self.snp.addMVar(self.x_len, vtype=GRB.INTEGER)
    self.y_snp = self.snp.addMVar(self.y_len, vtype=GRB.INTEGER)
   
    # binary variables to represent x and y
    self.x_bin_var_list_snp = []
    self.y_bin_var_list_snp = []
    for i in range(self.x_len):
        self.x_bin_var_list_snp.append(self.snp.addMVar(self.number_of_x_binaries[i], vtype=GRB.BINARY))
    for i in range(self.y_len):    
        self.y_bin_var_list_snp.append(self.snp.addMVar(self.number_of_y_binaries[i], vtype=GRB.BINARY))
    
    self.snp.update()
    
def add_snp_constrs(self):  
    #this model is only used when enabling sibling node pruning
    #### HPR constraints
    self.snp.addConstrs((LinExpr([self.A[i][j] for j in range(self.x_len)], [self.x_snp[j] for j in range(self.x_len)]) + 
                           LinExpr([self.B[i][j] for j in range(self.y_len)], [self.y_snp[j] for j in range(self.y_len)]) - 
                           self.a[i] <= 0 for i in range(self.a_len))) 
    
    if global_vars.first_cnstr_quad:
        self.snp.addConstrs((LinExpr([self.C[i][j] for j in range(self.x_len)], [self.x_snp[j] for j in range(self.x_len)]) + 
                               LinExpr([self.D[i][j] for j in range(self.y_len)], [self.y_snp[j] for j in range(self.y_len)]) - 
                               self.b[i] <= 0 for i in range(1, self.b_len)))
        
        if self.b_len >= 1:
            if global_vars.quad_ll_cnstr:
                self.snp.addQConstr(0.5 * (self.y_snp @ self.Q_cnstr @ self.y_snp) +
                                      LinExpr([self.C[0][j] for j in range(self.x_len)], [self.x_snp[j] for j in range(self.x_len)]) + 
                                      LinExpr([self.D[0][j] for j in range(self.y_len)], [self.y_snp[j] for j in range(self.y_len)]) - 
                                      self.b[0] <= 0) 
            else:
                self.snp.addConstr(LinExpr([self.C[0][j] for j in range(self.x_len)], [self.x_snp[j] for j in range(self.x_len)]) + 
                                   LinExpr([self.D[0][j] for j in range(self.y_len)], [self.y_snp[j] for j in range(self.y_len)]) - 
                                   self.b[0] <= 0) 
    else:
        self.snp.addConstrs((LinExpr([self.C[i][j] for j in range(self.x_len)], [self.x_snp[j] for j in range(self.x_len)]) + 
                           LinExpr([self.D[i][j] for j in range(self.y_len)], [self.y_snp[j] for j in range(self.y_len)]) - 
                           self.b[i] <= 0 for i in range(self.b_len-1)))
    
        if self.b_len >= 1:
            self.snp.addQConstr(0.5 * (self.y_snp @ self.Q_cnstr @ self.y_snp) +
                                  LinExpr([self.C[self.b_len-1][j] for j in range(self.x_len)], [self.x_snp[j] for j in range(self.x_len)]) + 
                                  LinExpr([self.D[self.b_len-1][j] for j in range(self.y_len)], [self.y_snp[j] for j in range(self.y_len)]) - 
                                  self.b[self.b_len-1] <= 0) 
                    
    self.snp.update()
    
    #### binary constraints 
    self.snp.addConstrs((self.x_snp[i] == - 2**(self.number_of_x_binaries[i] -1) * self.x_bin_var_list_snp[i][self.number_of_x_binaries[i]-1] 
                           + LinExpr([2**k for k in range(self.number_of_x_binaries[i] - 1)], [self.x_bin_var_list_snp[i][k] for k in range(self.number_of_x_binaries[i] - 1)]) for i in range(self.x_len)))
    
    self.snp.addConstrs((self.y_snp[i] == - 2**(self.number_of_y_binaries[i] -1) * self.y_bin_var_list_snp[i][self.number_of_y_binaries[i]-1] 
                           + LinExpr([2**k for k in range(self.number_of_y_binaries[i] - 1)], [self.y_bin_var_list_snp[i][k] for k in range(self.number_of_y_binaries[i] - 1)]) for i in range(self.y_len)))
    
    self.snp.update()

def add_snp_obj(self):
    # obj is zero because we only check feasibility
    self.snp.setObjective(0, gp.GRB.MAXIMIZE)

def build_snp_model(self):
    cvp_model_start_time = time.time()
    set_snp_parameters(self)
    add_snp_vars(self)
    add_snp_constrs(self) 
    add_snp_obj(self)
    global_vars.cvp_model_time = time.time() - cvp_model_start_time
    
    
###################################### Refinement Model ##################################################
def set_refinement_parameters(self):
    self.ref = Model()
    
    self.ref.parameters.parallel = 1
    self.ref.parameters.threads = 1
    
def add_refinement_vars(self):  
    self.y_ref = self.ref.integer_var_list(keys = self.y_len, lb=np.array(self.y_bounds)[:,0], ub=np.array(self.y_bounds)[:,1], name = "y")
    
def add_refinement_constrs(self):
    # lin UL cnstrs
    self.ref.add_constraints([self.ref.scal_prod([self.y_ref[j] for j in range(self.y_len)], [self.B[i][j] for j in range(self.y_len)]) <= 0 for i in range(self.a_len)], ["ul_lin_cnstr"+str(i) for i in range(self.a_len)])
    # lin LL cnstrs
    if global_vars.first_cnstr_quad:
        self.ref.add_constraints([self.ref.scal_prod([self.y_ref[j] for j in range(self.y_len)], [self.D[i][j] for j in range(self.y_len)]) <= 0 for i in range(1, self.b_len)], ["ll_lin_cnstr"+str(i) for i in range(1, self.b_len)])
    else:
        self.ref.add_constraints([self.ref.scal_prod([self.y_ref[j] for j in range(self.y_len)], [self.D[i][j] for j in range(self.y_len)]) <= 0 for i in range(self.b_len-1)], ["ll_lin_cnstr"+str(i) for i in range(self.b_len-1)])
        
def add_refinement_obj(self):
    # we min c_y @ y because x is fixed and, hence, the term c_x @ x would be constant
    self.ref.set_objective('min', self.ref.scal_prod([self.y_ref[i] for i in range(self.y_len)], [self.c_y[i] for i in range(self.y_len)])) 
    
def build_refinement_model(self):
    set_refinement_parameters(self)
    add_refinement_vars(self)
    add_refinement_constrs(self)
    add_refinement_obj(self)
import numpy as np
import time
# we have to take gurobi here because cplex cannot handle nonconvex quadratically-constrained problems
import gurobipy as gp
from gurobipy import GRB, LinExpr
import global_vars
from add_user_cuts_to_CVP_model import add_relevant_user_cuts_snp

def check_sibling_node(self):
    # subproblems
    self.subproblem_feasibility_snp = [True]*(self.b_len + 1)
    self.prune_sibling_node = False
    
    def modify_snp_model(self):
        # take the bounds of the parent node
        self.x_snp.lb = np.array(global_vars.parent_node_bounds)[0][0:self.x_len]
        self.x_snp.ub = np.array(global_vars.parent_node_bounds)[1][0:self.x_len]
        self.y_snp.lb = np.array(global_vars.parent_node_bounds)[0][self.x_len:self.x_len + self.y_len]
        self.y_snp.ub = np.array(global_vars.parent_node_bounds)[1][self.x_len:self.x_len + self.y_len]
        
        for i in range(self.x_len):
            idx = self.x_len + self.y_len + sum(self.number_of_x_binaries[k] for k in range(i))
            
            self.x_bin_var_list[i].lb = np.array(global_vars.parent_node_bounds)[0][idx:idx + self.number_of_x_binaries[i]]
            self.x_bin_var_list[i].ub = np.array(global_vars.parent_node_bounds)[1][idx:idx + self.number_of_x_binaries[i]]
        
        if not global_vars.use_INGC_only_on_x:
            for i in range(self.y_len):
                idx = self.x_len + self.y_len + sum(self.number_of_x_binaries[k] for k in range(self.x_len)) + sum(self.number_of_y_binaries[k] for k in range(i))
                
                self.y_bin_var_list[i].lb = np.array(global_vars.parent_node_bounds)[0][idx:idx + self.number_of_y_binaries[i]]
                self.y_bin_var_list[i].ub = np.array(global_vars.parent_node_bounds)[1][idx:idx + self.number_of_y_binaries[i]]
                
        self.snp.update()    
        
        constraints_to_remove = [constr for constr in self.snp.getConstrs() if constr.ConstrName.startswith("Cut")]
        for constr in constraints_to_remove:
            self.snp.remove(constr)
            
        self.snp.update()  
         
        add_relevant_user_cuts_snp(self)    
        
        self.snp.update()
       
        ##### discard previous solution ifo
        self.snp.reset(1)
        self.snp.update()

        
    def solve_snp(self):  
        y_Q_obj_y = gp.QuadExpr()
        for i in range(self.y_len):
            for j in range(self.y_len):
                if self.Q_obj[i, j] != 0:  # Only add non-zero terms for efficiency
                    y_Q_obj_y.addTerms(self.Q_obj[i, j], self.y_snp.tolist()[i], self.y_snp.tolist()[j])
        # solve subproblems
        #!!!!!!!!!!!
        for i in range(self.b_len + 1):
            self.snp.Params.TimeLimit = max(0, global_vars.end_time - time.time())
            if time.time() < global_vars.end_time:
                if global_vars.first_cnstr_quad:
                    if i == 0:
                        self.new_cnstr_snp = self.snp.addQConstr(0.5 * (self.y_j + self.Delta_y_opt) @ self.Q_obj @ (self.y_j + self.Delta_y_opt) -
                                                               0.5 * y_Q_obj_y + 
                                                               gp.quicksum([self.d_y[i] * self.y_j[i] for i in range(self.y_len)]) + 
                                                               gp.quicksum([self.d_y[i] * self.Delta_y_opt[i] for i in range(self.y_len)]) -
                                                               LinExpr([self.d_y[i] for i in range(self.y_len)], [self.y_snp.tolist()[i] for i in range(self.y_len)])>= 0)
                    elif i == 1:
                        # first cnstr is quad cnstr
                        if global_vars.quad_ll_cnstr:
                            self.new_cnstr_snp = self.snp.addConstr(0.5 * (self.y_j + self.Delta_y_opt) @ self.Q_cnstr @ (self.y_j + self.Delta_y_opt) +
                                                                  LinExpr([self.C[0][j] for j in range(self.x_len)], [self.x_snp.tolist()[j] for j in range(self.x_len)]) + 
                                                                  gp.quicksum([self.D[0][j] * (self.y_j[j] + self.Delta_y_opt[j]) for j in range(self.y_len)]) - 
                                                                  self.b[0]  >= 0)
                        else:
                            self.new_cnstr_snp = self.snp.addConstr(LinExpr([self.C[0][j] for j in range(self.x_len)], [self.x_snp.tolist()[j] for j in range(self.x_len)]) + 
                                                                  gp.quicksum([self.D[0][j] * (self.y_j[j] + self.Delta_y_opt[j]) for j in range(self.y_len)]) - 
                                                                  self.b[0]  >= 0)
                    else:
                        self.new_cnstr_snp = self.snp.addConstr(LinExpr([self.C[i-1][j] for j in range(self.x_len)], [self.x_snp.tolist()[j] for j in range(self.x_len)]) +
                                                              gp.quicksum(self.D[i-1][j] * (self.y_j[j] + self.Delta_y_opt[j]) for j in range(self.y_len)) - 
                                                              self.b[i-1] >= 0)
                            
                else:
                    if i == 0:
                        self.new_cnstr_snp = self.snp.addQConstr(0.5 * (self.y_j + self.Delta_y_opt) @ self.Q_obj @ (self.y_j + self.Delta_y_opt) -
                                                               0.5 * y_Q_obj_y + 
                                                               gp.quicksum([self.d_y[i] * self.y_j[i] for i in range(self.y_len)]) + 
                                                               gp.quicksum([self.d_y[i] * self.Delta_y_opt[i] for i in range(self.y_len)]) -
                                                               LinExpr([self.d_y[i] for i in range(self.y_len)], [self.y_snp.tolist()[i] for i in range(self.y_len)])>= 0)
                    elif i == self.b_len:
                        # last cnstr is quad cnstr
                        self.new_cnstr_snp = self.snp.addConstr(0.5 * (self.y_j + self.Delta_y_opt) @ self.Q_cnstr @ (self.y_j + self.Delta_y_opt) +
                                                              LinExpr([self.C[self.b_len-1][i] for i in range(self.x_len)], [self.x_snp.tolist()[i] for i in range(self.x_len)]) + 
                                                              gp.quicksum([self.D[self.b_len-1][i] * (self.y_j[i] + self.Delta_y_opt[i]) for i in range(self.y_len)]) - 
                                                              self.b[self.b_len-1]  >= 0)
                    else:
                        self.new_cnstr_snp = self.snp.addConstr(LinExpr([self.C[i-1][j] for j in range(self.x_len)], [self.x_snp.tolist()[j] for j in range(self.x_len)]) +
                                                              gp.quicksum(self.D[i-1][j] * (self.y_j[j] + self.Delta_y_opt[j]) for j in range(self.y_len)) - 
                                                              self.b[i-1] >= 0)
                    
                self.snp.update()
                
                    
                self.snp.optimize()
                
                    
                    
                if self.snp.status == GRB.INFEASIBLE or self.snp.status == GRB.INF_OR_UNBD:
                    self.subproblem_feasibility_snp[i] = False
                self.snp.remove(self.new_cnstr_snp)
                self.snp.update()   
                #print(self.snp.NumConstrs)
                
                    
        
        #print(self.subproblem_feasibility_2)
        if all(self.subproblem_feasibility_snp[i] == False for i in range(self.b_len+1)):
            self.prune_sibling_node = True   
            
            global_vars.x_node_bounds_prune = []
            global_vars.y_node_bounds_prune = []
            
            for i in range(self.x_len):
                global_vars.x_node_bounds_prune.append([self.x_snp[i].lb, self.x_snp[i].ub])
                
            for i in range(self.y_len):
                global_vars.y_node_bounds_prune.append([self.y_snp[i].lb, self.y_snp[i].ub])    
                
            
    
    start_time_snp = time.time()
    modify_snp_model(self)    
    solve_snp(self)    
    time_snp = time.time() - start_time_snp  
    global_vars.total_time_snp += time_snp  
    
    return self.prune_sibling_node

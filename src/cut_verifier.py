import time
import numpy as np
# we have to take gurobi here because cplex cannot handle nonconvex quadratically-constrained problems
import gurobipy as gp
from gurobipy import GRB, LinExpr

from gurobi_callback_for_CVPs import check_callback
from add_user_cuts_to_CVP_model import add_relevant_user_cuts
from sibling_node_pruning import check_sibling_node
import global_vars

        
def check_plane(self, hyperplane_was_valid_for_subproblem, alpha, beta, tau, Z):
    self.plane_feas_tol = 1e-5
    
    self.subproblem_feasibility = [True]*(self.b_len + 1)
    self.every_subproblem_is_infeasible = False
    self.hyperplane_was_valid_for_subproblem = hyperplane_was_valid_for_subproblem ### list of violated subproblems
    
    self.alpha_opt = alpha
    self.beta_opt = beta
    self.tau_opt = tau
    
    self.Z = Z
    
    def modify_cvp_model(self):
        self.x_cvp.lb = np.array(self.x_node_bounds)[:,0]
        self.x_cvp.ub = np.array(self.x_node_bounds)[:,1]
        self.y_cvp.lb = np.array(self.y_node_bounds)[:,0]
        self.y_cvp.ub = np.array(self.y_node_bounds)[:,1]
        
        
        for i in range(self.x_len):
            self.x_bin_var_list[i].lb = np.array(self.x_bin_node_bounds)[i,:,0]
            self.x_bin_var_list[i].ub = np.array(self.x_bin_node_bounds)[i,:,1]
        
        if not global_vars.use_INGC_only_on_x:
            for i in range(self.y_len):
                self.y_bin_var_list[i].lb = np.array(self.y_bin_node_bounds)[i,:,0]
                self.y_bin_var_list[i].ub = np.array(self.y_bin_node_bounds)[i,:,1]
                
        self.cvp.update()    
        
        constraints_to_remove = [constr for constr in self.cvp.getConstrs() if constr.ConstrName.startswith("Cut")]
        for constr in constraints_to_remove:
            self.cvp.remove(constr)
            
        self.cvp.update()     
        add_relevant_user_cuts(self)    
        self.cvp.update()
        
        
        
        ##### if we use the scoop problem and apply a DC on a node solution that is on the boundary of the bilevel-free set
        
        self.use_INGC_on_node_sol = False
        
        if not global_vars.use_max_t_scoop:
            if not global_vars.use_scoop or (global_vars.use_scoop and global_vars.node_sol_on_boundary_of_bf_set):
                self.use_INGC_on_node_sol = True
                self.x_j_bin = [list(map(int, bin(abs(num))[2:].zfill(self.number_of_x_binaries[i]-1)))[::-1] + [1 if num < 0 else 0]
                                for i, num in enumerate(self.x_j.astype(int))]
                self.y_j_bin = [list(map(int, bin(abs(num))[2:].zfill(self.number_of_y_binaries[i]-1)))[::-1] + [1 if num < 0 else 0]
                                for i, num in enumerate(self.y_j.astype(int))]
                
        ######################################################################################################################        
                
               
        
        self.cvp.setObjective(LinExpr([self.alpha_opt[i] for i in range(self.x_len)], [self.x_cvp.tolist()[i] for i in range(self.x_len)]) + 
                                LinExpr([self.beta_opt[i] for i in range(self.y_len)], [self.y_cvp.tolist()[i] for i in range(self.y_len)]) - 
                                self.tau_opt, gp.GRB.MAXIMIZE)
         
        self.cvp.update()
        
        
        ##### discard previous solution ifo
        self.cvp.reset(1)
        self.cvp.update()
        
        if self.use_INGC_on_node_sol:
            self.LHS_INGC = gp.quicksum( (1 - self.x_bin_var_list[i][j]) if self.x_j_bin[i][j] == 1 else self.x_bin_var_list[i][j] for i in range(self.x_len) for j in range(self.number_of_x_binaries[i])) + gp.quicksum( (1 - self.y_bin_var_list[i][j]) if self.y_j_bin[i][j] == 1 else self.y_bin_var_list[i][j] for i in range(self.y_len) for j in range(self.number_of_y_binaries[i]))
        
        
        
    def prepare_for_callback(self):
        self.cvp._x = self.x_cvp
        self.cvp._y = self.y_cvp
        self.cvp._solve_subproblems_to_global_optimality = self.solve_subproblems_to_global_optimality
        
    def solve_cvp(self):  
        y_Q_obj_y = gp.QuadExpr()
        for i in range(self.y_len):
            for j in range(self.y_len):
                if self.Q_obj[i, j] != 0:  # Only add non-zero terms for efficiency
                    y_Q_obj_y.addTerms(self.Q_obj[i, j], self.y_cvp.tolist()[i], self.y_cvp.tolist()[j])
                  
        self.latest_valid_subproblems = []
        
        # solve subproblems
        start_time_solve_all_CVPs = time.time()
        times_for_each_CVP = []
        
        for i in range(self.b_len + 1):
            self.cvp.Params.TimeLimit = max(0, global_vars.end_time - time.time())
            start_time_current_CVP = time.time()
            if time.time() < global_vars.end_time:
                if self.subproblem_feasibility[i]:
                    # only solve feasible subproblems
                    if not self.hyperplane_was_valid_for_subproblem[i]: # is initially always false
                        if self.use_INGC_on_node_sol and i in self.active_cnstrs_of_bf_set:
                            self.INGC_on_node_sol = self.cvp.addConstr(self.LHS_INGC >= 1)
                        if global_vars.first_cnstr_quad:
                            if i == 0:
                                self.new_cnstr = self.cvp.addQConstr(0.5 * (self.y_j + self.Delta_y_opt) @ self.Q_obj @ (self.y_j + self.Delta_y_opt) -
                                                                       0.5 * y_Q_obj_y + 
                                                                       gp.quicksum([self.d_y[i] * self.y_j[i] for i in range(self.y_len)]) + 
                                                                       gp.quicksum([self.d_y[i] * self.Delta_y_opt[i] for i in range(self.y_len)]) -
                                                                       LinExpr([self.d_y[i] for i in range(self.y_len)], [self.y_cvp.tolist()[i] for i in range(self.y_len)])>= 0)
                            elif i == 1:
                                # first cnstr is quad cnstr
                                if global_vars.quad_ll_cnstr:
                                    self.new_cnstr = self.cvp.addConstr(0.5 * (self.y_j + self.Delta_y_opt) @ self.Q_cnstr @ (self.y_j + self.Delta_y_opt) +
                                                                          LinExpr([self.C[0][j] for j in range(self.x_len)], [self.x_cvp.tolist()[j] for j in range(self.x_len)]) + 
                                                                          gp.quicksum([self.D[0][j] * (self.y_j[j] + self.Delta_y_opt[j]) for j in range(self.y_len)]) - 
                                                                          self.b[0]  >= 0)
                                else:
                                    self.new_cnstr = self.cvp.addConstr(LinExpr([self.C[0][j] for j in range(self.x_len)], [self.x_cvp.tolist()[j] for j in range(self.x_len)]) + 
                                                                          gp.quicksum([self.D[0][j] * (self.y_j[j] + self.Delta_y_opt[j]) for j in range(self.y_len)]) - 
                                                                          self.b[0]  >= 0)
                            else:
                                self.new_cnstr = self.cvp.addConstr(LinExpr([self.C[i-1][j] for j in range(self.x_len)], [self.x_cvp.tolist()[j] for j in range(self.x_len)]) +
                                                                      gp.quicksum(self.D[i-1][j] * (self.y_j[j] + self.Delta_y_opt[j]) for j in range(self.y_len)) - 
                                                                      self.b[i-1] >= 0)
                                    
                        else:
                            if i == 0:
                                self.new_cnstr = self.cvp.addQConstr(0.5 * (self.y_j + self.Delta_y_opt) @ self.Q_obj @ (self.y_j + self.Delta_y_opt) -
                                                                       0.5 * y_Q_obj_y + 
                                                                       gp.quicksum([self.d_y[i] * self.y_j[i] for i in range(self.y_len)]) + 
                                                                       gp.quicksum([self.d_y[i] * self.Delta_y_opt[i] for i in range(self.y_len)]) -
                                                                       LinExpr([self.d_y[i] for i in range(self.y_len)], [self.y_cvp.tolist()[i] for i in range(self.y_len)])>= 0)
                            elif i == self.b_len:
                                # last cnstr is quad cnstr
                                self.new_cnstr = self.cvp.addConstr(0.5 * (self.y_j + self.Delta_y_opt) @ self.Q_cnstr @ (self.y_j + self.Delta_y_opt) +
                                                                      LinExpr([self.C[self.b_len-1][i] for i in range(self.x_len)], [self.x_cvp.tolist()[i] for i in range(self.x_len)]) + 
                                                                      gp.quicksum([self.D[self.b_len-1][i] * (self.y_j[i] + self.Delta_y_opt[i]) for i in range(self.y_len)]) - 
                                                                      self.b[self.b_len-1]  >= 0)
                            else:
                                self.new_cnstr = self.cvp.addConstr(LinExpr([self.C[i-1][j] for j in range(self.x_len)], [self.x_cvp.tolist()[j] for j in range(self.x_len)]) +
                                                                      gp.quicksum(self.D[i-1][j] * (self.y_j[j] + self.Delta_y_opt[j]) for j in range(self.y_len)) - 
                                                                      self.b[i-1] >= 0)
                            
                        self.cvp.NumObj=0  #remove obj before defining it new, otherwise it wont set the obj correct
                        self.cvp.update()
                        
                        self.cvp.setObjective(self.alpha_opt @ self.x_cvp + self.beta_opt @ self.y_cvp - self.tau_opt, gp.GRB.MAXIMIZE)
                        
                        prepare_for_callback(self)
                        self.cvp.optimize(check_callback)
                        global_vars.total_nr_of_subproblems += 1
                        
                        # gather solution
                        if self.solve_subproblems_to_global_optimality:
                            if self.cvp.status == GRB.OPTIMAL: 
                                if self.cvp.getObjective().getValue() >= self.plane_feas_tol:
                                    # solution has positive objective value
                                    if not (np.array([*self.x_cvp.x, *self.y_cvp.x]).tolist() in [self.Z[i].tolist() for i in range(len(self.Z))]):
                                        if len(self.Z) == 0:
                                            self.Z = np.array([[*self.x_cvp.x, *self.y_cvp.x]], dtype=object) ### * transforms [1] to 1
                                        else:    
                                            self.Z = np.vstack([self.Z, np.array([*self.x_cvp.x, *self.y_cvp.x])])
                                else:
                                    self.hyperplane_was_valid_for_subproblem[i] = True
                                    self.latest_valid_subproblems.append(i)
                            else: # the subproblem is either infeasible or it got interrupted because of a non-positive global upper bound 
                                  # -> in both casees the hyperplane is valid
                                if self.cvp.status == GRB.INFEASIBLE or self.cvp.status == GRB.INF_OR_UNBD:
                                    # can never be unbounded, hence, infeasible
                                    self.subproblem_feasibility[i] = False
                                self.hyperplane_was_valid_for_subproblem[i] = True   
                                self.latest_valid_subproblems.append(i)
                        else: # stop early
                            if self.cvp.status == GRB.OPTIMAL or self.cvp.status == GRB.INTERRUPTED: # we got a feasible point (a new incumbent)
                                if self.cvp.getObjective().getValue() >= self.plane_feas_tol: # test if the incumbent value is positive
                                    if not (np.array([*self.x_cvp.x, *self.y_cvp.x]).tolist() in [self.Z[i].tolist() for i in range(len(self.Z))]):
                                        if len(self.Z) == 0:
                                            self.Z = np.array([[*self.x_cvp.x, *self.y_cvp.x]], dtype=object) ### * transforms [1] to 1
                                        else:    
                                            self.Z = np.vstack([self.Z, np.array([*self.x_cvp.x, *self.y_cvp.x])])
                                else:
                                    # either the optimal solution yields a non-positive objective value or we did terminate with a non-positive
                                    # upper bound for the optimal objective value
                                    self.hyperplane_was_valid_for_subproblem[i] = True 
                                    self.latest_valid_subproblems.append(i)
                                
                            else: # the subproblem is infeasible; unboundedness cannot occur because the HPR is bounded
                                self.subproblem_feasibility[i] = False
                                self.hyperplane_was_valid_for_subproblem[i] = True 
                                self.latest_valid_subproblems.append(i)
                        
                        ### remove constraints that only appear in subproblem i        
                        self.cvp.remove(self.new_cnstr)
                        if self.use_INGC_on_node_sol and i in self.active_cnstrs_of_bf_set:
                            self.cvp.remove(self.INGC_on_node_sol)
                        self.cvp.update()
            
                else:
                    # the subproblem is infeasible, set hyperplane_valid parameter to True
                    self.hyperplane_was_valid_for_subproblem[i] = True
                    self.latest_valid_subproblems.append(i)
            
            else:
                #### set everything to True because then the loop in CB_CP_Use_DC terminates
                self.hyperplane_was_valid_for_subproblem[i] = True
                self.every_subproblem_is_infeasible = True
                
            time_to_solve_current_CVP = time.time() - start_time_current_CVP
            times_for_each_CVP.append(time_to_solve_current_CVP)
        
        
        # store times to solve the CVPs
            # time for all together:
        time_to_solve_CVPs = time.time() - start_time_solve_all_CVPs   
        global_vars.total_time_to_solve_CVPs += time_to_solve_CVPs   
        
            # maximum time for a CVP in the current iteration:
        max_time_CVP = max(times_for_each_CVP)
        global_vars.idealized_time_to_solve_CVPs += max_time_CVP
        
        if all(self.subproblem_feasibility[i] == False for i in range(self.b_len+1)):
            self.every_subproblem_is_infeasible = True
            
    def check_sibling_node_feasibility(self):
        if all(self.subproblem_feasibility[i] == False for i in range(self.b_len+1)):
            self.every_subproblem_is_infeasible = True
    
            global_vars.prune_sibling_node = check_sibling_node(self) 
            
    check_modify_start_time = time.time()        
    modify_cvp_model(self)
    global_vars.cvp_modify_time += time.time() - check_modify_start_time
    cvp_solve_start_time = time.time()
    solve_cvp(self)    
    global_vars.cvp_solve_time += time.time() - cvp_solve_start_time

    if global_vars.use_sibling_node_pruning:
        if time.time() < global_vars.end_time:
            if global_vars.current_node_depth == global_vars.previous_node_depth +1:
                # then we have an unexplored sibling node
                check_sibling_node_feasibility(self) 

    return(self.hyperplane_was_valid_for_subproblem, self.Z, self.every_subproblem_is_infeasible, self.latest_valid_subproblems)

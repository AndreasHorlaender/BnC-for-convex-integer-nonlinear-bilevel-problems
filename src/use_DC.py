import time
from cut_generating_problem import generate_DC
from cut_verifier import check_plane

def use_disjunctive_cut(self):
    Z = []
    
    hyperplane_was_valid_for_subproblem = [False] * (self.b_len +1)
    hyperplane_final_check = False 
    
    while not all(item == True for item in hyperplane_was_valid_for_subproblem) or not hyperplane_final_check:
        # generate cut
        hyperplane_parameters, Z = generate_DC(self, Z)   
        alpha = hyperplane_parameters[0]
        beta = hyperplane_parameters[1]
        tau = hyperplane_parameters[2]
        
        # check cut
        hyperplane_was_valid_for_subproblem, Z, every_subproblem_is_infeasible, latest_valid_subproblems = check_plane(self, hyperplane_was_valid_for_subproblem, alpha, beta, tau, Z)
            
        if every_subproblem_is_infeasible:
            #we terminate the while loop
            hyperplane_final_check = True 
        else:
            # differ between verification methods
            if self.use_smart_dc_verification:
                # only solve the subproblems again which violate the hyperplane, needs a final check
                # if the hyerplane gets final checked and violates a subproblem we have to set this parameter to False again
                if not all(item == True for item in hyperplane_was_valid_for_subproblem):
                    hyperplane_final_check = False
                    
                else:
                    # final check
                    if not hyperplane_final_check: 
                        hyperplane_final_check = True
                        hyperplane_was_valid_for_subproblem = [False] * (self.b_len +1)
                        # the subproblems for which the current hyperplane is valid don't need to be checked again
                       
                        for i in latest_valid_subproblems:
                            hyperplane_was_valid_for_subproblem[i] = True
            else:
                # every subproblem get solved again and there is no final check
                if not all(item == True for item in hyperplane_was_valid_for_subproblem):
                    # set everything to False
                    hyperplane_was_valid_for_subproblem = [False] * (self.b_len +1)
                    #hyperplane_final_check = False # in this case this boolean never turns True
                else:
                    hyperplane_final_check = True # to terminate while-loop
    
    if time.time() >= self.get_end_time():
        self.abort()
    
    
    
    # store model var names in a list
    model_var_list = []
    for i in range(self.x_len):
        model_var_list.append("x_" + str(i))
    for i in range(self.y_len):
        model_var_list.append("y_" + str(i))
    
    # store hyperplane coefficients in a list
    coeff_list = []
    for i in range(self.x_len):
        coeff_list.append(alpha[i])
    for i in range(self.y_len):
        coeff_list.append(beta[i])    
    
    if every_subproblem_is_infeasible:
        # prune subtree
        self.add_local(constraint=[['x_0'],[0]], sense='E', rhs=1.0) 
        self.cuts.append([model_var_list,coeff_list, 'L', tau]) # to get a correct cut_list length
        self.nodes_pruned += 1
    else:    
        #print("RHS", tau)
        self.add_local(constraint=[model_var_list,coeff_list], sense='L', rhs=tau) 
        self.number_of_DC += 1
        self.cuts.append([model_var_list,coeff_list, 'L', tau])
        
        
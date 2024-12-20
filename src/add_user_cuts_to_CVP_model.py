from gurobipy import LinExpr

### Here we add the DCs and INGCs of the current subtree in the gurobi-model used for verifiying the new cut. The entry of self.cuts[i]
### is the i-th cut in the list and looks as follows.
###
### If we have a DC cut then 
### self.cuts[i][0] = ['x_0', ..., 'x_n_1', 'y_0', ..., 'y_n_2'],
### self.cuts[i][1] = [alpha[0], ..., alpha[n_1], beta[0], ..., beta[n_2]],
### self.cuts[i][2] = 'L'
### self.cuts[i][3] = tau
###
### If we have an INGC then
### self.cuts[i][0] = ['x_bin_0[0]', ..., 'x_bin_0[nr_of_x_bin[0]]', 'x_bin_n_1[0]', ..., 'x_bin_n_1[nr_of_x_bin[n_1]]', 
###                    'y_bin_0[0]', ..., 'y_bin_0[nr_of_y_bin[0]]', 'y_bin_n_2[0]', ..., 'y_bin_n_2[nr_of_x_bin[n_2]]'] 
### is a list which contains the name of all binary variables. I.e. self.cuts[i][0][j][k] = 'x_bin_j[k] is the k-th binary variable
### of the j-th x-variable (if j <= self.x_len, otherwise it is 'y_bin_(j-x_len)[k]', so the k-th binary variable of the (j-x_len)-th y-var.)
### self.cuts[i][1] contains the corresponding coefficients (0 or 1).
### self.cuts[i][2] = 'G'
### self.cuts[i][3] = 1 - (number of binary_vars with coefficient 1)

def add_relevant_user_cuts(self):
    for i in range(len(self.cuts)):
        if self.cuts[i][2] == 'L':   ## then we have a DC
            self.cvp.addConstr(LinExpr([self.cuts[i][1][j] for j in range(self.x_len)], [self.x_cvp.tolist()[j] for j in range(self.x_len)]) + 
                                 LinExpr([self.cuts[i][1][self.x_len + j] for j in range(self.y_len)], [self.y_cvp.tolist()[j] for j in range(self.y_len)]) <= self.cuts[i][3], name='Cut'+str(i))
        else:    
            self.cvp.addConstr(LinExpr([self.cuts[i][1][j][k] for j in range(self.x_len) 
                                          for k in range(self.number_of_x_binaries[j]) ],
                                         [self.x_bin_var_list[j].tolist()[k] for j in range(self.x_len) 
                                          for k in range(self.number_of_x_binaries[j])] ) 
                                 + LinExpr([self.cuts[i][1][self.x_len + j][k] for j in range(self.y_len)
                                            for k in range(self.number_of_y_binaries[j])],
                                           [self.y_bin_var_list[j].tolist()[k] for j in range(self.y_len)
                                            for k in range(self.number_of_y_binaries[j])] ) 
                                 >= self.cuts[i][3], name='Cut'+str(i))
            
def add_relevant_user_cuts_snp(self):
    for i in range(len(self.cuts)):
        if self.cuts[i][2] == 'L':   ## then we have a DC
            self.nnp.addConstr(LinExpr([self.cuts[i][1][j] for j in range(self.x_len)], [self.x_nnp[j] for j in range(self.x_len)]) + 
                                 LinExpr([self.cuts[i][1][self.x_len + j] for j in range(self.y_len)], [self.y_nnp[j] for j in range(self.y_len)]) <= self.cuts[i][3], name='Cut'+str(i))
        else:    
            self.nnp.addConstr(LinExpr([self.cuts[i][1][j][k] for j in range(self.x_len) 
                                          for k in range(self.number_of_x_binaries[j]) ],
                                         [self.x_bin_var_list_nnp[j][k] for j in range(self.x_len) 
                                          for k in range(self.number_of_x_binaries[j])] ) 
                                 + LinExpr([self.cuts[i][1][self.x_len + j][k] for j in range(self.y_len)
                                            for k in range(self.number_of_y_binaries[j])],
                                           [self.y_bin_var_list_nnp[j][k] for j in range(self.y_len)
                                            for k in range(self.number_of_y_binaries[j])] ) 
                                 >= self.cuts[i][3], name='Cut'+str(i))            
import global_vars 

def use_integer_no_good_cut(self):
    #print("USE NO-GOOD CUT")
    # store model var names and values in lists
    bin_var_names = []
    self.x_bin = [] 
    self.y_bin = [] 
    for i in range(self.x_len):
        for j in range(self.number_of_x_binaries[i]):
            bin_var_names.append("x_bin_" + str(i) + "_" + str(j))
            self.x_bin.append(self.get_values("x_bin_" + str(i) + "_" + str(j)))
    
    # store hyperplane coefficients in a list
    coeff_list = []
    lhs_constant = 0
    for i in range(len(self.x_bin)):
        if self.x_bin[i] <= 1e-4:
            coeff_list.append(1)
        else:
            coeff_list.append(-1)
            lhs_constant += 1
            
    if global_vars.use_INGC_only_on_x:
        # INGC only on x variables
        self.add_local(constraint=[bin_var_names, coeff_list], sense='G', rhs=1 - lhs_constant)
    else:     
        # INGC on both x and y
        for i in range(self.y_len):
            for j in range(self.number_of_y_binaries[i]):
                bin_var_names.append("y_bin_" + str(i) + "_" + str(j))
                self.y_bin.append(self.get_values("y_bin_" + str(i) + "_" + str(j)))
                
                
        for i in range(len(self.y_bin)):
            if self.y_bin[i] <= 1e-4:
                coeff_list.append(1)
            else:
                coeff_list.append(-1)   
                lhs_constant += 1   
                
        self.add_local(constraint=[bin_var_names, coeff_list], sense='G', rhs=1 - lhs_constant)
    
    
    self.number_of_INGC += 1
    
    # store bin coeffs of every variable in a separate list rather than all in one
    coeff_list_new = []
    idx = 0
    for i in range(self.x_len):
        coeff_list_new.append(coeff_list[idx:idx + self.number_of_x_binaries[i]])
        idx += self.number_of_x_binaries[i]
        
    idy = 0
    for i in range(self.y_len):
        coeff_list_new.append(coeff_list[idx + idy:idx + idy + self.number_of_y_binaries[i]])
        idy += self.number_of_y_binaries[i]    
        
    self.cuts.append([bin_var_names, coeff_list_new, 'G', 1-lhs_constant])

import numpy as np
from docplex.mp.model import Model
import global_vars

class Node_Builder:
    
    def __init__(self, packed_data):
        self.m = Model()
        
        self.m.parameters.mip.strategy.nodeselect = 0 ### depht-first search is needed to avoid "jumping" between the levels of the BnB tree
        self.m.parameters.preprocessing.presolve = False
        self.m.parameters.mip.strategy.heuristiceffort = 0
        
        self.m.parameters.parallel = 1
        self.m.parameters.threads = 1
        
        #self.m.parameters.mip.strategy.miqcpstrat = 2
       
        self.m.parameters.timelimit = 7200
        
        self.packed_data = packed_data
        # lenghts
        self.x_len = self.packed_data[0]
        self.y_len = self.packed_data[1]
        self.a_len = self.packed_data[6]
        self.b_len = self.packed_data[7]
        
        # bounds
        self.x_bounds = np.array(self.packed_data[2])
        self.y_bounds = np.array(self.packed_data[3])
        
        self.number_of_x_binaries = self.packed_data[4]
        self.number_of_y_binaries = self.packed_data[5]
        
        # objective info
        self.c_x = np.array(self.packed_data[8])
        self.c_y = np.array(self.packed_data[9])
        self.d_y = np.array(self.packed_data[10])
        self.Q_obj = self.packed_data[11]
        self.Q_cnstr = self.packed_data[12]
        
        # LHS constraints
        self.A = self.packed_data[13]
        self.B = self.packed_data[14]
        self.C = self.packed_data[15]
        self.D = self.packed_data[16]
        
        # RHS constraints
        self.a = np.array(self.packed_data[17])
        self.b = np.array(self.packed_data[18])
        self.a_sense = np.array(self.packed_data[19])
        self.b_sense = np.array(self.packed_data[20])
        

    def build(self):
    
        def add_node_vars(self):  
            self.x = self.m.integer_var_list(keys = self.x_len, lb=np.array(self.x_bounds)[:,0], ub=np.array(self.x_bounds)[:,1], name = "x")
            self.y = self.m.integer_var_list(keys = self.y_len, lb=np.array(self.y_bounds)[:,0], ub=np.array(self.y_bounds)[:,1], name = "y")
            
            # binary variables to represent x and y
            self.x_bin_var_list = []
            for i in range(self.x_len):
                self.x_bin_var_list.append(self.m.binary_var_list(keys = self.number_of_x_binaries[i], name = "x_bin_" + str(i)))
            
            if not global_vars.use_INGC_only_on_x:    
                self.y_bin_var_list = []    
                for i in range(self.y_len):    
                    self.y_bin_var_list.append(self.m.binary_var_list(keys = self.number_of_y_binaries[i], name = "y_bin_" + str(i)))
                
            ### modify the branching order of the variables cpx.order.set(var_name, priority, up/down/default)
            cpx = self.m.get_cplex()
            name_list = []
            order_list = []
            
            for i in range(self.x_len):
                name_list.append('x_' + str(i))
                order_list.append(1)
            for i in range(self.y_len):
                name_list.append('y_' + str(i))   
                order_list.append(1)
            ## no specifying for binary variables gives them a priority of 0 automatically
            cpx.order.set([(a, b, c) for a, b, c in zip(name_list, order_list, [0]*len(name_list))])
            
        def add_node_constrs(self):
            
            # lin UL cnstrs
            self.m.add_constraints((self.m.scal_prod([self.x[j] for j in range(self.x_len)], [self.A[i][j] for j in range(self.x_len)]) + 
                                    self.m.scal_prod([self.y[j] for j in range(self.y_len)], [self.B[i][j] for j in range(self.y_len)]) - 
                                    self.a[i] <= 0 for i in range(self.a_len)))
            
            if global_vars.first_cnstr_quad:
                # lin LL cnstrs
                self.m.add_constraints((self.m.scal_prod([self.x[j] for j in range(self.x_len)], [self.C[i][j] for j in range(self.x_len)]) + 
                                        self.m.scal_prod([self.y[j] for j in range(self.y_len)], [self.D[i][j] for j in range(self.y_len)]) - 
                                        self.b[i] <= 0 for i in range(1, self.b_len)))
            
                # quad LL cnstr
                if self.b_len >= 1:
                    if global_vars.quad_ll_cnstr:
                        self.m.add_quadratic_constraints([0.5 * (self.y @ self.Q_cnstr @ self.y) + 
                                                          self.m.scal_prod([self.x[j] for j in range(self.x_len)], [self.C[0][j] for j in range(self.x_len)]) + 
                                                          self.m.scal_prod([self.y[j] for j in range(self.y_len)], [self.D[0][j] for j in range(self.y_len)]) - 
                                                          self.b[0] <= 0])
                    else:
                        self.m.add_constraints([self.m.scal_prod([self.x[j] for j in range(self.x_len)], [self.C[0][j] for j in range(self.x_len)]) + 
                                                self.m.scal_prod([self.y[j] for j in range(self.y_len)], [self.D[0][j] for j in range(self.y_len)]) - 
                                                self.b[0] <= 0])
                    
            else: # last cnstr is quad
                # lin LL cnstrs
                self.m.add_constraints((self.m.scal_prod([self.x[j] for j in range(self.x_len)], [self.C[i][j] for j in range(self.x_len)]) + 
                                        self.m.scal_prod([self.y[j] for j in range(self.y_len)], [self.D[i][j] for j in range(self.y_len)]) - 
                                        self.b[i] <= 0 for i in range(self.b_len-1)))
                
                # quad LL cnstr
                if self.b_len >= 1:
                    self.m.add_quadratic_constraints([0.5 * (self.y @ self.Q_cnstr @ self.y) + 
                                                      self.m.scal_prod([self.x[j] for j in range(self.x_len)], [self.C[self.b_len-1][j] for j in range(self.x_len)]) + 
                                                      self.m.scal_prod([self.y[j] for j in range(self.y_len)], [self.D[self.b_len-1][j] for j in range(self.y_len)]) - 
                                                      self.b[self.b_len-1] <= 0])
        
            # binary expansion
            self.m.add_constraints((self.x[i] == - 2**(self.number_of_x_binaries[i] -1) * self.x_bin_var_list[i][self.number_of_x_binaries[i]-1]
                                               + self.m.scal_prod([self.x_bin_var_list[i][k] for k in range(self.number_of_x_binaries[i] - 1)], 
                                               [2**k for k in range(self.number_of_x_binaries[i] - 1)]) for i in range(self.x_len)))
            
            if not global_vars.use_INGC_only_on_x:    
                self.m.add_constraints((self.y[i] == - 2**(self.number_of_y_binaries[i] -1) * self.y_bin_var_list[i][self.number_of_y_binaries[i]-1]
                                                   + self.m.scal_prod([self.y_bin_var_list[i][k] for k in range(self.number_of_y_binaries[i] - 1)], 
                                                   [2**k for k in range(self.number_of_y_binaries[i] - 1)]) for i in range(self.y_len)))    
                    
        def add_node_obj(self):
            self.m.set_objective('min', self.m.scal_prod([self.x[i] for i in range(self.x_len)], [self.c_x[i] for i in range(self.x_len)]) + 
                                 self.m.scal_prod([self.y[i] for i in range(self.y_len)], [self.c_y[i] for i in range(self.y_len)])) 
            
        add_node_vars(self)
        add_node_constrs(self)
        add_node_obj(self)
        return(self.m, self.x)
    
   



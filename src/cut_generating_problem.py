import numpy as np
from docplex.mp.model import Model
import time
import global_vars
    
def generate_DC(self, Z):
    self.Z = Z
    
    self.plane = Model()
    
    self.plane.parameters.timelimit = max(0, self.get_end_time() - self.get_time()) 
    self.plane.parameters.parallel = 1
    self.plane.parameters.threads = 1
    
    def add_plane_vars(self):
        self.alpha = self.plane.continuous_var_list(keys = self.x_len, lb =-np.infty)
        self.beta = self.plane.continuous_var_list(keys = self.y_len, lb = -np.infty)
        self.tau = self.plane.continuous_var(lb = -np.infty)
        
    def add_plane_constrs(self):   
        ### 1-NORM
        self.plane.add_constraint(self.plane.sum(self.plane.abs(self.alpha[i]) for i in range(self.x_len)) + self.plane.sum(self.plane.abs(self.beta[i]) for i in range(self.y_len)) + self.plane.abs(self.tau) <= 1)
       
        def hyperplane(x,y):
            # alpha x + beta y - tau
            return(self.plane.scal_prod([self.alpha[i] for i in range(self.x_len)], [x[i] for i in range(self.x_len)]) + self.plane.scal_prod([self.beta[i] for i in range(self.y_len)], [y[i] for i in range(self.y_len)]) - self.tau)
        self.hyperplane = hyperplane
        
        for i in range(len(self.Z)):
            self.plane.add_constraint(self.hyperplane(self.Z[i][:self.x_len], self.Z[i][self.x_len:]) <= 0)
       
    def add_plane_obj(self):       
       self.plane.set_objective('max', self.plane.scal_prod([self.alpha[i] for i in range(self.x_len)], [self.x_j[i] for i in range(self.x_len)]) + self.plane.scal_prod([self.beta[i] for i in range(self.y_len)], [self.y_j[i] for i in range(self.y_len)]) - self.tau)
       
    def solve_plane(self):  
        self.plane.solve()
        
        # store the optimal plane coefficients in lists
        self.alpha_opt = np.zeros(self.x_len)
        for i in range(self.x_len):
            self.alpha_opt[i] = self.alpha[i].solution_value
            
        self.beta_opt = np.zeros(self.y_len)
        for i in range(self.y_len):
            self.beta_opt[i] = self.beta[i].solution_value
            
        self.tau_opt = self.tau.solution_value    
    
    cgp_model_start_time = time.time() 
    if time.time() >= global_vars.end_time:
        self.abort()
        
    add_plane_vars(self)
    add_plane_constrs(self)
    add_plane_obj(self)
    global_vars.cgp_model_time += time.time() - cgp_model_start_time
    cgp_solve_start_time = time.time()
    solve_plane(self)  
    global_vars.cgp_solve_time += time.time() - cgp_solve_start_time
    
    return([self.alpha_opt, self.beta_opt, self.tau_opt], self.Z)
import numpy as np
import global_vars

def refinement_procedure(self):
    ########################################################### first cnstr quad vs last cnstr quad vs no cnstr quad !!!!!!!!!!!!!!
    def modify_ref_constrs(self):  
        # modify HPR constraints
        for i in range(self.a_len):
            rhs_i = self.a[i] - np.dot(self.A[i], self.x_j)
            self.ref.get_constraint_by_name("ul_lin_cnstr"+str(i)).rhs=rhs_i
                                   
        if self.b_len >= 1:
            for i in range(self.b_len):
                rhs_i = self.b[i] - np.dot(self.C[i], self.x_j)
                if global_vars.first_cnstr_quad:
                    if i >= 1:
                        # modify linear cnstrs
                        self.ref.get_constraint_by_name("ll_lin_cnstr"+str(i)).rhs=rhs_i
                    else:
                        # add quad cnstr
                        self.ref.add_constraint(self.ref.scal_prod([self.y_ref[j] for j in range(self.y_len)], [self.D[0][j] for j in range(self.y_len)]) 
                                               +  0.5 * (self.y_ref @ self.Q_cnstr @ self.y_ref) <= rhs_i, ctname="ll_quad_cnstr") 
                       
                else:
                    if i <= self.b_len -2:
                        # modify linear cnstrs
                        self.ref.get_constraint_by_name("ll_lin_cnstr"+str(i)).rhs=rhs_i
                    else:
                        # add quad cnstr
                        self.ref.add_constraint(self.ref.scal_prod([self.y_ref[j] for j in range(self.y_len)], [self.D[self.b_len-1][j] for j in range(self.y_len)]) 
                                               +  0.5 * (self.y_ref @ self.Q_cnstr @ self.y_ref) <= rhs_i, ctname="ll_quad_cnstr") 
                
        # add refinement constraint
        self.ref.add_constraint(self.ref.scal_prod([self.y_ref[i] for i in range(self.y_len)], 
                                                   [self.d_y[i] for i in range(self.y_len)]) 
                                    +  0.5 * (self.y_ref @ self.Q_obj @ self.y_ref) <= self.ll_obj_value + 1e-4, ctname="ref_cnstr")
            
    def solve_ref_model(self):
        self.refinement_has_solution = True
        self.ref.solve(log_output=False)
        if self.ref.solve_details.status == 'integer optimal solution' or self.ref.solve_details.status == 'integer optimal, tolerance':
            self.y_bf = np.zeros(self.y_len)
            for i in range(self.y_len):
                self.y_bf[i] = np.round(self.y_ref[i].solution_value)
        else:
            self.refinement_has_solution = False
        
        
        
    modify_ref_constrs(self)
    solve_ref_model(self)
    
    # cant modify quad cnstrs in cplex -> build them from scratch
    self.ref.remove_constraint("ll_quad_cnstr")
    self.ref.remove_constraint("ref_cnstr")
    
    if self.refinement_has_solution:
        if not (abs(self.y_bf - self.y_hat) <= 1e-4).all():
            global_vars.nr_of_refinements += 1
            
    return self.refinement_has_solution, self.y_bf
    
    
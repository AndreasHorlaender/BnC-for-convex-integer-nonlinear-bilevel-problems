from lower_level import solve_LL
from scoop_problem import solve_scoop
from use_DC import use_disjunctive_cut
from use_INGC import use_integer_no_good_cut
import global_vars
import time

def cut_procedure(self):
    cut_procedure_start_time = time.time()
    
    # check bilevel feasibility
    node_sol_is_bil_feas, self.y_hat = solve_LL(self)
    #print("TIME LL: ", time.time() - start_time_ll)
    if not node_sol_is_bil_feas:
        if global_vars.use_only_INGCs:
            use_integer_no_good_cut(self)
        else: # we use DCs
            if global_vars.use_scoop:
                # solve scoop problem
                self.Delta_y_opt, use_no_good_cut, self.active_cnstrs_of_bf_set = solve_scoop(self)
                if global_vars.use_max_t_scoop and use_no_good_cut: # we use INGC if t = 0
                    use_integer_no_good_cut(self)
                else: # t /= 0 or we always use DCs
                    use_disjunctive_cut(self)
            else: # always use DCs based on lower-level response
                self.Delta_y_opt = self.y_hat - self.y_j
                self.active_cnstrs_of_bf_set = [i for i in range(self.b_len +1)]
                use_disjunctive_cut(self)
                
           
        
    else:
        #print("The point is bilevel-feasible.")
        # since the node gets pruned, we remove it from our node lists
        self.old_x_node_bounds = self.old_x_node_bounds[0:len(self.old_x_node_bounds)-1]
        self.old_y_node_bounds = self.old_y_node_bounds[0:len(self.old_y_node_bounds)-1]
        self.old_x_bin_node_bounds = self.old_x_bin_node_bounds[0:len(self.old_x_bin_node_bounds)-1]
        if not global_vars.use_INGC_only_on_x:
            self.old_y_bin_node_bounds = self.old_y_bin_node_bounds[0:len(self.old_y_bin_node_bounds)-1]


    global_vars.cut_procedure_time += time.time() - cut_procedure_start_time

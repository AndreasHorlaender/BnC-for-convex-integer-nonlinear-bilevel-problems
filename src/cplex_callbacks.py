import time

from cplex.callbacks import LazyConstraintCallback
from cplex.callbacks import BranchCallback
from cplex.callbacks import HeuristicCallback

from subproblem_modeler import get_data
from subproblem_modeler import build_lower_level
from subproblem_modeler import build_scoop_problem
from subproblem_modeler import build_CVP_model
from subproblem_modeler import build_snp_model
from subproblem_modeler import build_refinement_model

from get_node_info import get_node_info
from modify_cut_list import modify_cut_list
from handle_node_solution import cut_procedure
import global_vars

class update_incumbent_callback(HeuristicCallback):
    def __init__(self, env):
        HeuristicCallback.__init__(self, env)
        
    def __call__(self):
        if global_vars.found_new_bil_feas_point:
            inc_obj_val = self.get_incumbent_objective_value()
            
            if inc_obj_val >= global_vars.bil_feas_UL_obj_val + 1e-5:
                self.set_solution([global_vars.var_names, global_vars.bil_feas_point])
                
        global_vars.found_new_bil_feas_point = False # prevents the update from appearing in every single node       

class snp_callback(BranchCallback):
    def __init__(self, env):
        BranchCallback.__init__(self, env)
        
    def __call__(self):
        #track node depth
        global_vars.previous_node_depth = global_vars.current_node_depth
        global_vars.current_node_depth = self.get_current_node_depth()
        self.depth = self.get_current_node_depth()
        if self.get_num_branches() >= 2:
            global_vars.parent_node_bounds = [self.get_lower_bounds(), self.get_upper_bounds()]
           
        if global_vars.prune_sibling_node:
            x_len = len(global_vars.x_node_bounds_prune)
            y_len = len(global_vars.y_node_bounds_prune)
            
            self.x_node_bounds2 = []
            for i in range(x_len):
                self.x_node_bounds2.append([self.get_lower_bounds("x_" + str(i)), self.get_upper_bounds("x_" + str(i))])
               
            self.y_node_bounds2 = []
            for i in range(y_len):
                self.y_node_bounds2.append([self.get_lower_bounds("y_" + str(i)), self.get_upper_bounds("y_" + str(i))])
                
                
            if all(self.x_node_bounds2[i][0] >= global_vars.x_node_bounds_prune[i][0] for i in range(x_len)) and all(self.x_node_bounds2[i][1] <= global_vars.x_node_bounds_prune[i][1] for i in range(x_len)) and all(self.y_node_bounds2[i][0] >= global_vars.y_node_bounds_prune[i][0] for i in range(y_len)) and all(self.y_node_bounds2[i][1] <= global_vars.y_node_bounds_prune[i][1] for i in range(y_len)):
                # then we are now in the sibling node of the previous node  
                self.prune()
                global_vars.number_of_sibling_prunes += 1
                    
        global_vars.prune_sibling_node = False

class my_callback(LazyConstraintCallback):
    def __init__(self, env):
        LazyConstraintCallback.__init__(self, env)
       
    def read_info(self, packed_data):
        self.packed_data = packed_data
        self.cuts = []
        self.old_x_node_bounds = []
        self.old_y_node_bounds = []
        self.old_x_bin_node_bounds = []
        self.old_y_bin_node_bounds = []
        
        self.number_of_DC = 0
        self.number_of_INGC = 0
        self.nodes_pruned = 0
        self.sibling_nodes_pruned = 0
        
        self.b_len = self.packed_data[3]
        
    def build_models(self):
        get_data(self)
        build_lower_level(self) 
        if not global_vars.use_only_INGCs:
            build_scoop_problem(self)
            build_CVP_model(self)
            if global_vars.use_sibling_node_pruning:
                build_snp_model(self)
        if global_vars.use_refinement_procedure:
            build_refinement_model(self)
            
    def return_info(self):
        return self.number_of_DC, self.number_of_INGC, self.nodes_pruned, self.sibling_nodes_pruned
        
    def __call__(self):
        if global_vars.end_time == None:
            global_vars.end_time = self.get_end_time()
            
        cb_start_time = time.time()
        if self.get_solution_source() == 111: # then the point is actually a node solution
        # see https://www.ibm.com/docs/en/icos/20.1.0?topic=classes-cplexcallbackssolutionsource
             
            #track node depth
            global_vars.previous_node_depth = global_vars.current_node_depth
            global_vars.current_node_depth = self.get_current_node_depth()
            node_info_start_time = time.time()
            get_node_info(self)
            global_vars.node_info_time += time.time() - node_info_start_time
            
            
            if global_vars.quad_ll_cnstr:
                if self.get_quadratic_slacks()[0] >= -1e-3:
                    # we have an integer feasible point which also satisfies all constraints of the model ?????
                    modify_cut_list(self)
                    cut_procedure(self) 
            else:
                modify_cut_list(self)
                cut_procedure(self) 
                    
        global_vars.cb_time += time.time() - cb_start_time    
#paramterizations
use_only_INGCs = None
use_sibling_node_pruning = None
use_refinement_procedure = None
quad_ll_cnstr = None
first_cnstr_quad = None
prune_sibling_node = False
solve_max_ul_problems = None
use_INGC_only_on_x = None
use_scoop = None
use_max_t_scoop = None
use_1_norm_scoop = None

node_sol_on_boundary_of_bf_set = None

#measurements for enhancements
number_of_sibling_prunes = 0
nr_of_refinements = 0
number_of_useless_cuts = 0

number_of_iters_to_generate_DC = 0



#variables to measure the time/effort of the different steps in the algorithm
node_info_time = 0
modify_cut_list_time = 0
cut_procedure_time = 0
ll_model_time = 0
ll_modify_time = 0
ll_solve_time = 0
scoop_model_time = 0
scoop_modify_time = 0
scoop_solve_time = 0
cgp_model_time = 0
cgp_solve_time = 0
cvp_model_time = 0
cvp_modify_time = 0
cvp_solve_time = 0
cb_time = 0
total_time_to_solve_CVPs = 0
idealized_time_to_solve_CVPs = 0
total_time_NNP = 0
time_refinement = 0
total_nr_of_subproblems = 0

end_time = None

#helping variables
previous_node_depth = 0
current_node_depth = 0
parent_node_bounds = []
found_new_bil_feas_point = False
var_names = []    
bil_feas_point = []
bil_feas_UL_obj_val = float('inf')
inc_obj_val = float("inf")
x_node_bounds_prune = None
y_node_bounds_prune = None

DC_depth_list = [0]*1002
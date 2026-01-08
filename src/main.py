import os
import sys
import time
from data_reader import Data_Reader
from node_builder import Node_Builder
from node_solver import Node_Solver
from result_analyzer import Result_Analyzer
import global_vars
import parameter_settings


result_file = parameter_settings.result_filename

########### check cmd input
instance_name = sys.argv[1]

global_vars.use_only_INGCs = parameter_settings.use_only_INGCs

# scoop
global_vars.use_scoop = parameter_settings.use_scoop
global_vars.use_max_t_scoop = parameter_settings.use_max_t_scoop
global_vars.use_1_norm_scoop = parameter_settings.use_1_norm_scoop

# further algorithmic techniques
global_vars.use_sibling_node_pruning = parameter_settings.use_sibling_node_pruning
global_vars.use_refinement_procedure = parameter_settings.use_refinement_procedure
use_smart_dc_verification = parameter_settings.use_smart_dc_verification
solve_subproblems_to_global_optimality = parameter_settings.solve_subproblems_to_global_optimality

# instances
global_vars.quad_ll_cnstr = parameter_settings.quad_ll_cnstr
global_vars.first_cnstr_quad = parameter_settings.first_cnstr_quad
global_vars.solve_max_ul_problems = parameter_settings.solve_max_ul_problems
use_qkp_with_50_percent_ll_cnstrs = parameter_settings.use_qkp_with_50_percent_ll_cnstrs




global_vars.use_INGC_only_on_x = False        
#######################################################
# disable incompatible settings
if global_vars.use_only_INGCs:
    global_vars.use_max_t_scoop = False
    global_vars.use_scoop = False
    
if global_vars.use_max_t_scoop:
    global_vars.use_scoop = True    
    global_vars.use_1_norm_scoop = False

if global_vars.use_INGC_only_on_x:
    global_vars.use_refinement_procedure = True
    
   
        
########### modify result_file name
if use_qkp_with_50_percent_ll_cnstrs:
    instance_name = instance_name + "_50_50"
    result_file = result_file + "_50_50"
    ## for this setting we have ot use the last constraint
    global_vars.first_cnstr_quad = False
    
if global_vars.solve_max_ul_problems:
    result_file = result_file + "_max"
else:
    result_file = result_file + "_min"

if global_vars.quad_ll_cnstr:
    result_file = result_file + "_quad"
else:
    result_file = result_file + "_no_quad"    
    #global_vars.first_cnstr_quad = False
    
# prepare data parsing
data_directory_name = "../data/"
data_directory = os.fsencode(data_directory_name)

list_of_files = filter( lambda x: os.path.isfile(os.path.join(data_directory, x)), os.listdir(data_directory) )
for file in list_of_files:
    data_filename = os.fsdecode(file)
    if data_filename.endswith(".aux"):
        if os.path.splitext(data_filename)[0] == instance_name:
            current_file_name = os.path.splitext(data_filename)[0]
            current_mps_file = current_file_name + ".mps" 
            current_aux_file = current_file_name + ".aux"
            if global_vars.first_cnstr_quad:
                current_txt_file = current_file_name + "_first_cnstr.txt" # randomly generated matrices for quadratic program     
            else:
                current_txt_file = current_file_name + "_last_cnstr.txt"
                
            # for sepcific aux-files the name of the other files has to be adjusted
            if data_filename.endswith("_0_100.aux") or data_filename.endswith("_50_50.aux") or data_filename.endswith("_100_0.aux"):
                current_file_name = os.path.splitext(data_filename)[0]
                current_mps_file = current_file_name[:len(current_file_name)-6] + ".mps" 
                current_aux_file = current_file_name + ".aux"
                if global_vars.first_cnstr_quad:
                    current_txt_file = current_file_name + "_first_cnstr.txt" # randomly generated matrices for quadratic program     
                else:
                    current_txt_file = current_file_name + "_last_cnstr.txt"   
                if use_qkp_with_50_percent_ll_cnstrs:
                    ### the textfile is named first_cnstr but is now used for last cnstr
                    current_txt_file = current_file_name[:len(current_file_name)-6] + "_first_cnstr.txt" 
                
            #print(current_aux_file)
            # read csv file
            data_reader = Data_Reader()
            packed_data = data_reader.read_from_files(current_mps_file, current_aux_file, current_txt_file, use_smart_dc_verification, solve_subproblems_to_global_optimality)
            
            # reset global vars before every new instance
            global_vars.end_time = None
            global_vars.prune_sibling_node = False
            global_vars.previous_node_depth = 0
            global_vars.current_node_depth = 0
            global_vars.number_of_sibling_prunes = 0
            global_vars.parent_node_bounds = []
            global_vars.total_nr_of_subproblems = 0
            global_vars.found_new_bil_feas_point = False
            global_vars.var_names = []    
            global_vars.bil_feas_point = []
            global_vars.bil_feas_UL_obj_val = float('inf')
            global_vars.inc_obj_val = float('inf')
            global_vars.total_time_to_solve_CVPs = 0
            global_vars.idealized_time_to_solve_CVPs = 0
            global_vars.total_time_snp = 0
            global_vars.number_of_useless_cuts = 0
            
            global_vars.node_info_time = 0
            global_vars.modify_cut_list_time = 0
            global_vars.cut_procedure_time = 0
            global_vars.ll_model_time = 0
            global_vars.ll_modify_time = 0
            global_vars.ll_solve_time = 0
            global_vars.scoop_model_time = 0
            global_vars.scoop_modify_time = 0
            global_vars.scoop_solve_time = 0
            global_vars.cgp_model_time = 0
            global_vars.cgp_solve_time = 0
            global_vars.cvp_model_time = 0
            global_vars.cvp_modify_time = 0
            global_vars.cvp_solve_time = 0
            global_vars.cb_time = 0
            global_vars.time_refinement = 0
            global_vars.number_of_iters_to_generate_DC = 0
            
            global_vars.x_node_bounds_prune = None
            global_vars.y_node_bounds_prune = None
            global_vars.node_sol_on_boundary_of_bf_set = None
            
            global_vars.nr_of_refinements = 0
            global_vars.list_of_bf_sets = []
            global_vars.score_of_bf_sets = []
            global_vars.nr_of_bf_sets_recycled = 0
            global_vars.DC_depth_list = [0]*1002
            
            # create model
            #print("Building the node ...")
            node_builder = Node_Builder(packed_data)
            node_model, x = node_builder.build()
            
            #print("Solving the node ...")
            start_time = time.time()
            node_solver = Node_Solver(node_model, x, packed_data, instance_name)
            solution, rel_mip_gap, best_upper_bound, number_of_nodes, number_of_DC, number_of_INGC, nodes_pruned, sibling_nodes_pruned = node_solver.solve()
            solving_time = time.time() - start_time
            
            # save results
            result_analyzer = Result_Analyzer(packed_data, result_file, current_file_name, use_smart_dc_verification, solve_subproblems_to_global_optimality, solution, solving_time, rel_mip_gap, best_upper_bound, number_of_nodes, number_of_DC, number_of_INGC, nodes_pruned, sibling_nodes_pruned)
            result_analyzer.save_results()

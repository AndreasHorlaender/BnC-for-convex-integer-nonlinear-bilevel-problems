import csv 
import os
import global_vars
import numpy as np
from scipy.linalg import eigh

class Result_Analyzer():
    
    def __init__(self, packed_data, result_file, current_file_name, use_smart_dc_verification, solve_subproblems_to_global_optimality, solution, solving_time, rel_mip_gap, best_upper_bound, number_of_nodes, number_of_DC, number_of_INGC, nodes_pruned, sibling_nodes_pruned):
        self.packed_data = packed_data
        self.x_len = self.packed_data[0]
        self.y_len = self.packed_data[1]
        self.a_len = self.packed_data[6]
        self.b_len = self.packed_data[7]
        self.c_x = np.array(self.packed_data[8])
        self.c_y = np.array(self.packed_data[9])
        self.d_y = np.array(self.packed_data[10])
        
        # compute min/max objective coeffs
        self.min_ulo = min(min(self.c_x), min(self.c_y))
        self.max_ulo = max(max(self.c_x), max(self.c_y))
        self.min_llo = min(self.d_y)
        self.max_llo = max(self.d_y)
        
        # compute min/max eigenvalues of Q and P
        self.Q_obj = self.packed_data[11]
        self.Q_cnstr = self.packed_data[12]
        self.Q_eig_min = min(eigh(self.Q_obj, eigvals_only=True))
        self.Q_eig_max = max(eigh(self.Q_obj, eigvals_only=True))
        self.P_eig_min = min(eigh(self.Q_cnstr, eigvals_only=True))
        self.P_eig_max = max(eigh(self.Q_cnstr, eigvals_only=True))
        
        self.result_file = result_file
        self.current_file_name = current_file_name
        
        self.use_sibling_node_pruning = global_vars.use_sibling_node_pruning
        self.use_smart_dc_verification = use_smart_dc_verification
        self.solve_subproblems_to_global_optimality = solve_subproblems_to_global_optimality
        self.use_scoop = global_vars.use_scoop
        self.use_max_t_scoop = global_vars.use_max_t_scoop
        
        self.solution = solution
        self.solving_time = solving_time
        self.rel_mip_gap = rel_mip_gap
        self.best_upper_bound = best_upper_bound
        self.number_of_nodes = number_of_nodes
        self.number_of_DC = number_of_DC
        self.number_of_INGC = number_of_INGC
        self.nodes_pruned = nodes_pruned
        self.sibling_nodes_pruned = global_vars.number_of_sibling_prunes
        self.total_nr_of_subproblems = global_vars.total_nr_of_subproblems
        self.total_time_to_solve_CVPs = global_vars.total_time_to_solve_CVPs
        self.idealized_time_to_solve_CVPs = global_vars.idealized_time_to_solve_CVPs
        self.number_of_useless_cuts = global_vars.number_of_useless_cuts
        
        self.node_info_time = global_vars.node_info_time
        self.modify_cut_list_time = global_vars.modify_cut_list_time
        self.cut_procedure_time = global_vars.cut_procedure_time
        self.ll_model_time = global_vars.ll_model_time 
        self.ll_modify_time = global_vars.ll_modify_time 
        self.ll_solve_time = global_vars.ll_solve_time 
        self.scoop_model_time = global_vars.scoop_model_time 
        self.scoop_modify_time = global_vars.scoop_modify_time 
        self.scoop_solve_time = global_vars.scoop_solve_time 
        self.cgp_model_time = global_vars.cgp_model_time 
        self.cgp_solve_time = global_vars.cgp_solve_time 
        self.cvp_model_time = global_vars.cvp_model_time 
        self.cvp_modify_time = global_vars.cvp_modify_time 
        self.cvp_solve_time = global_vars.cvp_solve_time 
        self.snp_solve_time = global_vars.total_time_snp
        self.cb_time = global_vars.cb_time
        self.time_refinement = global_vars.time_refinement
        self.number_of_iters_to_generate_DC = global_vars.number_of_iters_to_generate_DC
        
        self.time_without_ref = self.solving_time - self.time_refinement
        
        self.use_refinement_procedure = global_vars.use_refinement_procedure
        self.nr_of_refinements = global_vars.nr_of_refinements
        self.nr_of_bf_sets_recycled = global_vars.nr_of_bf_sets_recycled
        
        self.DC_depth_list = global_vars.DC_depth_list
        
        self.idealized_solving_time = self.solving_time - self.total_time_to_solve_CVPs + self.idealized_time_to_solve_CVPs
        
        self.subproblems_time = (self.ll_model_time + self.ll_modify_time + self.ll_solve_time + 
                                 self.scoop_model_time + self.scoop_modify_time + self.scoop_solve_time + 
                                 self.cgp_model_time + self.cgp_solve_time + 
                                 self.cvp_model_time + self.cvp_modify_time + self.cvp_solve_time)
        
        
        if self.solving_time >= 7200.001:
            self.solution = 'time limit exceeded'
            self.solving_time = 9999999999999
            self.idealized_solving_time = 9999999999999
            self.number_of_nodes = 9999999999999
            self.number_of_DC = 9999999999999
            self.number_of_useless_cuts = 9999999999999
            self.number_of_INGC = 9999999999999
            self.nodes_pruned = 9999999999999
            self.sibling_nodes_pruned = 9999999999999
            self.solving_time = 9999999999999
            self.idealized_solving_time = 9999999999999
            self.total_nr_of_subproblems = 9999999999999
            self.time_without_ref = 9999999999999
            if self.rel_mip_gap == float('inf') or self.rel_mip_gap == 'nan' or self.rel_mip_gap == 'NaN' or self.rel_mip_gap == None:
                self.rel_mip_gap = 9999999999999
        
        
        
        self.method = ''
        if global_vars.use_only_INGCs:
            self.method = 'INGC'
        else:
            if self.use_max_t_scoop:
                self.method = 'INGC+DC+Scoop_min'
            elif self.use_scoop:
                if global_vars.use_1_norm_scoop:
                    self.method = 'DC+Scoop_sum'
                else:
                    self.method = 'DC+Scoop_min'
            else:
                self.method = 'DC_FR'
        
        
        if not global_vars.use_only_INGCs:
            if self.use_sibling_node_pruning:
                self.method += '+snp'
            if self.use_smart_dc_verification:
                self.method += '+smart'
            if not self.solve_subproblems_to_global_optimality:
                self.method += '+early'
            if self.use_refinement_procedure:
                self.method += '_ref'
            
    def save_results(self):
        header = ['name', 'objective value', 'method', 'time', 'idealized_time', 'nodes', 'DCs', 'useless DCs', 'INGCs', 
                  'sibling_node_pruning', 'smart_dc_verification', 'solve_cvp_globally', 
                  'Nr_total_CVPs', 'Nodes_pruned', 'sibling_nodes_pruned', 
                  'cb_time', 'node_info_time', 'modify_cut_list_time', 
                  'subproblems_time', 'll_model_time', 'll_modify_time', 'll_solve_time', 'scoop_model_time', 'scoop_modify_time', 
                  'scoop_solve_time', 'cgp_model_time', 
                  'cgp_solve_time', 'cvp_model_time', 'cvp_modify_time', 'cvp_solve_time', 'ref_time', 'time_without_ref', 'snp_time', 'rel_mip_gap', 
                  'total_time_CVPs',
                  'idealized_time_CVPs', 'total_nr_iter_adversarial', 'nr_refinements', 'best_upper_bound', 'nr_UL_vars',
                  'nr_LL_vars', 'nr_UL_cnstrs', 'nr_LL_cnstrs', 'min_UL_obj_coeff', 'max_UL_obj_coeff', 
                  'min_LL_obj_coeff', 'max_LL_obj_coeff', 'Q_eig_min', 'Q_eig_max', 'P_eig_min', 'P_eig_max', 'nr_bf_sets_recycled']
        
        data = [self.current_file_name, self.solution, self.method, self.solving_time, self.idealized_solving_time,
                self.number_of_nodes, self.number_of_DC, self.number_of_useless_cuts, self.number_of_INGC, 
                self.use_sibling_node_pruning, self.use_smart_dc_verification, self.solve_subproblems_to_global_optimality, 
                self.total_nr_of_subproblems, self.nodes_pruned, self.sibling_nodes_pruned, 
                self.cb_time, self.node_info_time, self.modify_cut_list_time, self.subproblems_time, self.ll_model_time, 
                self.ll_modify_time, self.ll_solve_time, self.scoop_model_time, self.scoop_modify_time, self.scoop_solve_time, 
                self.cgp_model_time, self.cgp_solve_time,
                self.cvp_model_time, self.cvp_modify_time, self.cvp_solve_time, self.time_refinement, self.time_without_ref, self.snp_solve_time,
                self.rel_mip_gap, 
                self.total_time_to_solve_CVPs, self.idealized_time_to_solve_CVPs, self.number_of_iters_to_generate_DC, self.nr_of_refinements,
                self.best_upper_bound, self.x_len, self.y_len, self.a_len, self.b_len, 
                self.min_ulo, self.max_ulo, self.min_llo, self.max_llo, self.Q_eig_min, self.Q_eig_max, 
                self.P_eig_min, self.P_eig_max, self.nr_of_bf_sets_recycled]
            
        with open('../results/csv/' + self.result_file + '.csv','a', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=';')
            
            if os.stat('../results/csv/' + self.result_file + '.csv').st_size == 0: #wenn die Datei leer ist
                writer.writerow(header)
            writer.writerow(data)
            
            
        '''
        header_2 = ['name', 'Nr DCs', *[i for i in range(1001)], '> 1000']
        
        data_2 = [self.current_file_name, self.number_of_DC, *self.DC_depth_list]
            
        with open('../results/csv/' + self.result_file + '_DC_depth_list.csv','a', newline='') as csvfile:
            writer = csv.writer(csvfile, delimiter=';')
            
            
            if os.stat('../results/csv/' + self.result_file + '_DC_depth_list.csv').st_size == 0: #wenn die Datei leer ist
                writer.writerow(header_2)
            writer.writerow(data_2)    
        '''
from cplex_callbacks import my_callback
from cplex_callbacks import snp_callback
from cplex_callbacks import update_incumbent_callback
            
class Node_Solver:
    
    def __init__(self, node_model, x, packed_data, instance_name):
        self.m = node_model
        self.x = x
        self.packed_data = packed_data
        self.instance_name = instance_name
       
    def solve(self):
        self.cb = self.m.register_callback(my_callback)
        self.m.register_callback(snp_callback) 
        self.m.register_callback(update_incumbent_callback)
        self.cb.read_info(self.packed_data)
        self.cb.build_models()
        
        solution = self.m.solve(log_output=False)
        
        number_of_DC, number_of_INGC, nodes_pruned, sibling_nodes_pruned = self.cb.return_info()
        number_of_nodes = self.m.solve_details.nb_nodes_processed
        
        rel_mip_gap = self.m.solve_details.mip_relative_gap
        best_upper_bound = self.m.solve_details.best_bound
        
        print(self.m.solve_details.status)
        if self.m.solve_details.status == 'integer optimal solution' or self.m.solve_details.status == 'integer optimal, tolerance':
            print(solution)
            
            return round(self.m.objective_value, 2), rel_mip_gap, best_upper_bound, number_of_nodes, number_of_DC, number_of_INGC, nodes_pruned, sibling_nodes_pruned 
        else:
            number_of_nodes = self.m.get_cplex().solution.progress.get_num_nodes_processed()
            return self.m.solve_details.status, rel_mip_gap, best_upper_bound, number_of_nodes, number_of_DC, number_of_INGC, nodes_pruned, sibling_nodes_pruned 

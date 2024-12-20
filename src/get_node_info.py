import numpy as np

def get_node_info(self): 
    self.node_depth = self.get_current_node_depth()
    
    # store values of the node solution (self.x^j, self.y^j)
    self.x_j = np.zeros(self.x_len)
    for i in range(self.x_len):
        self.x_j[i] = np.round(self.get_values("x_" + str(i)))
    
    self.y_j = np.zeros(self.y_len)
    for i in range(self.y_len):
        self.y_j[i] = np.round(self.get_values("y_" + str(i)))
    
    # get node bounds
    self.x_node_bounds = []
    for i in range(self.x_len):
        self.x_node_bounds.append([self.get_lower_bounds("x_" + str(i)), self.get_upper_bounds("x_" + str(i))])
       
    self.y_node_bounds = []
    for i in range(self.y_len):
        self.y_node_bounds.append([self.get_lower_bounds("y_" + str(i)), self.get_upper_bounds("y_" + str(i))])
    
    self.x_bin_node_bounds = []
    for i in range(self.x_len):
        self.x_aux = []
        for j in range(self.number_of_x_binaries[i]):
            self.x_aux.append([self.get_lower_bounds("x_bin_" + str(i) + "_" + str(j)), self.get_upper_bounds("x_bin_" + str(i) + "_" + str(j))])
        self.x_bin_node_bounds.append(self.x_aux)
           
    self.y_bin_node_bounds = []
    for i in range(self.y_len):
        self.y_aux = [] 
        for j in range(self.number_of_y_binaries[i]):
            self.y_aux.append([self.get_lower_bounds("y_bin_" + str(i) + "_" + str(j)), self.get_upper_bounds("y_bin_" + str(i) + "_" + str(j))])
        self.y_bin_node_bounds.append(self.y_aux)
import time
import global_vars
### We use depth-first search in our implementation so we check if the current belongs to the subtree we were in before by comparing the
### variable bounds of the current node with those of the nodes saved in the list old_.._node_bounds. First we check if the current node is 
### a child node of the last visited node. If yes, we keep all cuts from our list "cuts" and add them in the subproblems of the cut_verifier.
### If no, we remove the last node and the last cut of our lists and check the second-last visited node. And so on.

### self.old_x_node_bounds[n-1] = the last list inside the list self.old_x_node_bounds
### self.old_x_node_bounds[n-1][i] = the bounds of variable x[i] in the last-seen node
### self.old_x_node_bounds[n-1][i][1] = the upper bound of variable x[i] in the last-seen node

### self.x_bin_node_bounds[i] = list of the bounds of all binary variables of the i-th variable
### self.x_bin_node_bounds[i][j] = bounds of the j-th binary variable of the i-th variable
### self.x_bin_node_bounds[i][j][0] = lower bound of the j-th binary variable of the i-th variable

def modify_cut_list(self):
    modify_cut_list_start_time = time.time()
    ### check which is the parent node of the current node
    #print("x_bounds = ", self.x_node_bounds)
    #print("y_bounds = ", self.y_bounds)
    #print("y_bin=", self.y_bin_node_bounds)
    #print("")
    #print("y_bin_old=", self.old_y_bin_node_bounds)
    stop = False
    n = len(self.cuts) 
    #print("LEN CUTS BEFORE = ", n)
    #print("x", self.x_node_bounds)
    #print("xold",self.old_x_node_bounds)
    #print("xbin", self.x_bin_node_bounds)
    #print("xbinold",self.old_x_bin_node_bounds)
    if self.old_x_node_bounds:
        while stop == False:
            if all(self.x_node_bounds[i][0] >= self.old_x_node_bounds[n-1][i][0] for i in range(self.x_len)) and all(self.x_node_bounds[i][1] <= self.old_x_node_bounds[n-1][i][1] for i in range(self.x_len)):
                #print("x")
                if all(self.y_node_bounds[i][0] >= self.old_y_node_bounds[n-1][i][0] for i in range(self.y_len)) and all(self.y_node_bounds[i][1] <= self.old_y_node_bounds[n-1][i][1] for i in range(self.y_len)):
                    #print("y")
                    # if for all binary variables for all variables the bounds are ok
                    if all(all(self.x_bin_node_bounds[i][j][0] >= self.old_x_bin_node_bounds[n-1][i][j][0] for j in range(self.number_of_x_binaries[i])) and all(self.x_bin_node_bounds[i][j][1] <= self.old_x_bin_node_bounds[n-1][i][j][1] for j in range(self.number_of_x_binaries[i])) for i in range(self.x_len)):
                        #print("x_bin")
                        if all(all(self.y_bin_node_bounds[i][j][0] >= self.old_y_bin_node_bounds[n-1][i][j][0] for j in range(self.number_of_y_binaries[i])) and all(self.y_bin_node_bounds[i][j][1] <= self.old_y_bin_node_bounds[n-1][i][j][1] for j in range(self.number_of_y_binaries[i])) for i in range(self.y_len)):
                            #print("y_bin")
                            stop = True
                        #print("Its a child node!")
                        #print("N = ", n)
                        else:
                            n += -1
                    else:
                        n += -1     
                else:
                    n += -1        
            else:
                n += -1
            if n == 0:
                stop = True
                
    self.cuts = self.cuts[0:n]
    self.old_x_node_bounds = self.old_x_node_bounds[0:n]
    self.old_y_node_bounds = self.old_y_node_bounds[0:n]
    self.old_x_bin_node_bounds = self.old_x_bin_node_bounds[0:n]
    self.old_y_bin_node_bounds = self.old_y_bin_node_bounds[0:n]
    
    self.old_x_node_bounds.append(self.x_node_bounds)
    self.old_y_node_bounds.append(self.y_node_bounds)
    self.old_x_bin_node_bounds.append(self.x_bin_node_bounds)
    self.old_y_bin_node_bounds.append(self.y_bin_node_bounds)
    #print("LEN CUTS", len(self.cuts))
    global_vars.modify_cut_list_time += time.time() - modify_cut_list_start_time
    
    #assert all(self.get_linear_slacks()[i] >= -1e-5 for i in range(len(self.get_linear_slacks())))
    #assert all(self.get_quadratic_slacks()[i] >= -1e-5 for i in range(len(self.get_quadratic_slacks())))
     
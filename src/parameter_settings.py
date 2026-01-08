# every parameter (except result_filename) has to take either the value True or False

result_filename = "my_resultfile"   # name of the csv file containing the performance results
# different strings are added to the result_filename based on which of the other parameters are used (see main.py)


# methods
use_only_INGCs = False   # True: use INGCs, False: use DCs


# scoop problem
use_scoop = True   # True: solve scoop problem to compute Delta y, False: use follower's response
use_max_t_scoop = True   # True: max t in scoop objective, False: max ||s||_1
use_1_norm_scoop = False   # True: max ||s||_1 in scoop objective, False: max t


# further algorithmic techniques (see Section 3.4.3 in https://ubt.opus.hbz-nrw.de/frontdoor/deliver/index/docId/2786/file/Dissertation_Andreas_Horlaender.pdf)
use_sibling_node_pruning = False
use_refinement_procedure = False   # see Fischetti et al. 2018
use_smart_dc_verification = False   # True: only re-optimize the CVPs that violate the current hyperplane and perform a final check, False: re-optimize every feasible CVP
solve_subproblems_to_global_optimality = True   # True: solve CVPs to global optimality, False: stop as soon as a strictly positive lower bound is found

    
# instances
quad_ll_cnstr = True   # True: LL has a quadratic constraint, False: LL has only linear constraints
first_cnstr_quad = True   # True: first LL constraint is quadratic, False: last LL constraint is quadratic
solve_max_ul_problems = False   # True: max UL objective, False: min UL objective
use_qkp_with_50_percent_ll_cnstrs = False   # True: 50% of the QKP constraints are used as LL constraints, False: single LL constraint (see Gaar et al. 2023)
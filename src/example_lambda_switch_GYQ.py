"""
The BCN example for the $\lambda$ switch genetic network.
Solved using the algebraic method developed by Yuqian Guo etc.

Guo, Yuqian, Pan Wang, Weihua Gui, and Chunhua Yang. "Set stability and set stabilization of Boolean control networks
based on invariant subsets." Automatica 61 (2015): 106-112.

Please refer to "example_lambda_switch.py" for the results obtained using our proposed method.
"""

from algorithm.utils import read_network
from algorithm.related_work import *


if __name__ == '__main__':
    n, m, L = read_network('./networks/lambda_switch.txt')
    M_set = [11, 2, 30, 32, 31, 5, 20, 7, 24, 13]
    print('M_set = ', M_set)
    solver = GYQSolver(m, n, L, M_set)
    LCIS = solver.compute_largest_control_invariant_subset()
    print('LCIS = ', LCIS)
    print('Is globally set stabilizable? ', solver.is_set_stabilizable())
    print('Shortest transient period (T_M): ', solver.compute_shortest_transient_period())
    # time optimal state feedback (any logical sub-matrix of bF is a solution)
    # We can check that the one generated by the graphical method is indeed a sub-matrix of bF
    print('The bold F in Proposition 6 is:\n', solver.compute_time_optimal_stabilizer().astype(np.int8))
    print('(It can be validated that time-optimal F produced by our graphical method is a logical sub-matrix of the bold F here)')


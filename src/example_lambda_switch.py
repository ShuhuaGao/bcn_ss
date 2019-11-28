"""
The BCN example for the $\lambda$ switch genetic network.
Solved using the proposed approach.

Please refer to "example_lambda_switch_GYQ.py" for the results obtained using an existing algebraic method.
"""
from algorithm.proposed import GraphicalViewSolver
from algorithm.utils import read_network


def g(i: int, k: int) -> float:
    """
    Stage cost function for $\delta_N^i$ and $\delta_M^k$

    :return: cost of a transition
    """
    B = [0, 1]
    A = [14, 8, 14, 2, 17, 8, 1, 13, 3, 14, 11, 14, 19, 2, 19, 11,
         15, 10, 16, 6, 4, 8, 19, 12, 3, 14, 17, 8, 16, 11, 12, 11]
    return A[i - 1] + B[k - 1]


if __name__ == '__main__':
    to_visualize = False
    n, m, L = read_network('./networks/lambda_switch.txt')
    M_set = [11, 2, 30, 32, 31, 5, 20, 7, 24, 13]
    print('M_set = ', M_set)
    solver = GraphicalViewSolver(m, n, L, M_set)
    # STG
    stg = solver.build_stg()
    if to_visualize:
        from algorithm.utils import visualize
        visualize(stg, 'lambda_switch_stg.pdf', layout='dot')
    # LCIS
    LCIS = solver.compute_largest_control_invariant_subset()
    print('LCIS = ', LCIS)
    # BFT
    BFT = solver.build_BFT_for_M_extended_stg()
    if to_visualize and BFT:
        visualize(BFT, 'lambda_switch_bft.pdf', layout='dot')
    # Stabilizability
    print('Is globally set stabilizable? ', solver.is_set_stabilizable())
    # Shortest transient period
    print('Shortest transient period (T_M): ', solver.compute_shortest_transient_period())

    # Time optimal
    F = solver.compute_time_optimal_stabilizer()
    print('Time optimal: F = ', F)
    # General optimal and SPT
    SPT = solver.build_SPT_for_M_extended_stg(g)
    if to_visualize and SPT:
        visualize(SPT, 'lambda_switch_spt.pdf', layout='dot')
    F_opt = solver.compute_custom_optimal_stabilizer(g)
    print('General optimal: F = ', F_opt)
    # Get the path for \delta_{32}^{14} in the SPT and control sequence
    x0 = 14
    assert x0 in SPT
    p, u = solver.get_trajectory(SPT, x0)
    print(f'General optimal trajectory for initial state {x0}: ', p)
    print(f'\tachieved with control sequence: ', u)
    print(f'\twith a minimum cost: ', SPT.nodes[x0]['d_star'])


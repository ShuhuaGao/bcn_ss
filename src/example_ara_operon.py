"""
The Ara operon network example.
"""
from algorithm.proposed import GraphicalViewSolver
from algorithm.utils import read_network
import random


def inverse_map(i: int, n: int):
    """
    Accumulative STP of logical variables is bijective.
    Given a result i (\delta_{2^n}^i), find the corresponding logical values.

    :return a list of 0/1
    """
    r = []
    while n > 0:
        if i % 2 == 0:
            r.append(0)
            i = i // 2
        else:
            r.append(1)
            i = (i + 1) // 2
        n = n - 1
    r.reverse()
    return r


def g(i, k):
    """
    Stage cost

    :param i: state
    :param k: control
    :return: the cost
    """
    n = 9
    m = 4
    X = inverse_map(i, n)
    U = inverse_map(k, m)
    A = [0, 16, 40, 44, 28, 5, 18, 48, 24]
    B = [0, 48, 38, 12]
    return sum(a * x for a, x in zip(A, X)) + sum(b * u for b, u in zip(B, U))

if __name__ == '__main__':
    random.seed(18)  # to reproduce the result
    n, m, L = read_network('./networks/ara_operon.txt')
    expected_M_set_size = 40
    self_loops = {9, 265, 320, 352, 384, 472, 480, 504, 512}
    M_set = set()
    while len(M_set) < expected_M_set_size:
        M_set = self_loops.copy()
        M_set.update(random.sample(range(1, 2 ** n + 1), k=expected_M_set_size - len(M_set)))
    # M_set = {9, 265, 320, 352, 384, 472, 480, 504, 512}
    print('M set is ', M_set)
    solver = GraphicalViewSolver(m, n, L, M_set)
    stg = solver.build_stg()
    lcis = solver.compute_largest_control_invariant_subset()
    print(f'LCIS is of size {len(lcis)}: ', lcis)
    print('Is stabilizable? ', solver.is_set_stabilizable())
    t = solver.compute_shortest_transient_period()
    print('Network shortest transient period is: ', t)
    Ft = solver.compute_time_optimal_stabilizer()
    x0 = 11
    bft = solver.build_BFT_for_M_extended_stg()
    p_bft, u_bft = solver.get_trajectory(bft, x0)
    print(f'Time optimal trajectory from {x0}: ', p_bft, '\n\tdriven by: ', u_bft, '\n\twith time t*: ', bft.nodes[x0]['t_star'])
    spt = solver.build_SPT_for_M_extended_stg(g)
    p_spt, u_spt = solver.get_trajectory(spt, x0)
    print(f'General optimal trajectory measured by g from {x0}: ', p_spt, '\n\tdriven by: ', u_spt, '\n\twith cost J*: ', spt.nodes[x0]['d_star'])

    # visualization: generate Fig. 5 in the paper
    to_visualize = False
    if to_visualize:
        from networkx.drawing.nx_agraph import to_agraph
        import pygraphviz
        large_node_size = 0.09
        large_edge_width = 0.7
        stg.nodes[x0]['color'] = 'purple'
        stg.nodes[x0]['width'] = large_node_size
        for i in M_set:
            stg.nodes[i]['color'] = 'black'
            stg.nodes[i]['width'] = large_node_size
        for i in lcis:
            stg.nodes[i]['color'] = 'blue'
            stg.nodes[i]['width'] = large_node_size
        for i in range(len(p_bft) - 1):
            stg.edges[p_bft[i], p_bft[i + 1]]['color'] = 'green3'
            stg.edges[p_bft[i], p_bft[i + 1]]['penwidth'] = large_edge_width
            if i > 0:
                stg.nodes[p_bft[i]]['color'] = 'green3'
                stg.nodes[p_bft[i]]['width'] = large_node_size
        for i in range(len(p_spt) - 1):
            stg.edges[p_spt[i], p_spt[i + 1]]['color'] = 'red'
            stg.edges[p_spt[i], p_spt[i + 1]]['penwidth'] = large_edge_width
            if i > 0:
                stg.nodes[p_spt[i]]['color'] = 'red'
                stg.nodes[p_spt[i]]['width'] = large_node_size

        stg.graph['edge'] = {'arrowsize': '0', 'splines': 'curved', 'arrowhead': 'none', 'penwidth': 0.1, 'color': 'gray30', 'weight': 1.2}
        stg.graph['node'] = {'shape': 'point', 'fontsize': '10', 'width': '0.05', 'margin': 0, 'style': 'filled', 'color': 'gray52'}
        stg.graph['graph'] = {'scale': 1, 'nodesep': 0, 'ranksep': 0, 'margin': 0, 'ratio': 0.6}
        ag: pygraphviz.agraph.AGraph = to_agraph(stg)
        print('ag obtained')
        ag.layout('fdp')
        print('ag layout')
        ag.draw('ara_operon.pdf')



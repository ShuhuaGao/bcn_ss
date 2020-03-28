"""
Benchmark the performance of our approach against others using the TLGL network.

The GYQ approach refers to
Guo, Yuqian, Pan Wang, Weihua Gui, and Chunhua Yang. "Set stability and set stabilization of Boolean control networks
    based on invariant subsets." Automatica 61 (2015): 106-112.

The LRJ approach refers to 
R. Liu, J. Lu, W. X. Zheng, and J. Kurths, “Output feedback control
for set stabilization of boolean control networks,” IEEE transactions on
neural networks and learning systems, 2019.
"""

from algorithm.proposed import GraphicalViewSolver
from algorithm.related_work import GYQSolver, LRJSolver
from algorithm.utils import read_network
import random
from typing import List, Set, Callable
import time
import numpy as np


def compute_LCIS(m: int, n: int, L: List[int], M_set: Set[int], Solver):
    ts = time.time()
    solver = Solver(m, n, L, M_set)
    lcis = solver.compute_largest_control_invariant_subset()
    t = time.time() - ts
    return lcis, t


def check_global_stabilizability(m: int, n: int, L: List[int], M_set: Set[int], Solver):
    ts = time.time()
    solver = Solver(m, n, L, M_set)
    cs = solver.is_set_stabilizable()
    t = time.time() - ts
    return cs, t


def compute_time_optimal_stabilizer(m: int, n: int, L: List[int], M_set: Set[int], Solver):
    ts = time.time()
    solver = Solver(m, n, L, M_set)
    K = solver.compute_time_optimal_stabilizer()
    t = time.time() - ts
    return K, t


def generate_M_set(expected_M_set_size: int):
    M_set = set()
    while len(M_set) < expected_M_set_size:
        M_set = {50130, 50642, 58322, 58834, 65279, 65535}  # fix points
        M_set.update(random.sample(range(1, 2 ** n + 1),
                                   k=expected_M_set_size - len(M_set)))
    return M_set


def compute_fixed_points(m, n, L):
    """
    Find the fixed points of the BCN.
    A state is called a fixed point if it can transit to itself by a certain input.
    """
    M = 2 ** m
    N = 2 ** n
    fps = []
    for i in range(1, N + 1):
        for k in range(1, M + 1):
            blk = L[(k - 1) * N: k * N]
            j = blk[i - 1]
            if j == i:
                fps.append(i)
                print('FP: ', i)
                break
    return fps


N = 2 ** 16
M = 2 ** 3
Q = np.random.uniform(0, 10, N)
R = np.random.uniform(0, 10, M)


def g(i, k):
    """
    Stage cost

    :param i: state
    :param k: control
    :return: the cost
    """
    return Q[i - 1] + R[k - 1]


def compute_general_optimal_stabilizer(m: int, n: int, L: List[int], M_set: Set[int], g: Callable[[int, int], float], Solver):
    ts = time.time()
    solver = Solver(m, n, L, M_set)
    K = solver.compute_custom_optimal_stabilizer(g)
    t = time.time() - ts
    return K, t


if __name__ == '__main__':
    seed = 1
    random.seed(seed)    # to reproduce results
    n_runs = 20
    print(
        f'Experiment with random seed {seed} and {n_runs} runs for each task.')
    n, m, L = read_network('./networks/T_LGL.assr')
    expected_M_set_size = 1000
    print('The size of M set is ', expected_M_set_size)
    print("""NOTE: the GYQ approach will take extremely long time (possibly days).
        To test our approach only, set 'enable_GYQ = False'""")
    enable_GYQ = False
    print('Task 1: compute LCIS')
    t1 = []
    t2 = []
    t3 = []
    for r in range(n_runs):
        print(f' -Run {r}')
        M_set = generate_M_set(expected_M_set_size)
        LCIS1, t = compute_LCIS(m, n, L, M_set, GraphicalViewSolver)
        print(f'\tOur approach (seconds): {t}')
        t1.append(t)
        LCIS3, t = compute_LCIS(m, n, L, M_set, LRJSolver)
        print(f'\tLRJ approach (seconds): {t}')
        assert len(set(LCIS1) ^ LCIS3) == 0
        t3.append(t)
        if enable_GYQ:
            _, t = compute_LCIS(m, n, L, M_set, GYQSolver)
            print(f'\tGYQ approach (seconds): {t}')
            t2.append(t)

    print(' -Average time')
    print(f'\tOur approach (seconds): {sum(t1) / len(t1): .4f}')
    print(f'\tLRJ approach (seconds): {sum(t3) / len(t3): .4f}')
    if enable_GYQ:
        print(f'\tGYQ approach (seconds): {sum(t2) / len(t2): .4f}')

    print('Task 2: check stabilizability')
    t1 = []
    t2 = []
    for r in range(n_runs):
        print(f' -Run {r}')
        M_set = generate_M_set(expected_M_set_size)
        _, t = check_global_stabilizability(
            m, n, L, M_set, GraphicalViewSolver)
        print(f'\tOur approach (seconds): {t}')
        t1.append(t)
        if enable_GYQ:
            __, t = check_global_stabilizability(m, n, L, M_set, GYQSolver)
            print(f'\tGYQ approach (seconds): {t}')
            t2.append(t)
    print(' -Average time')
    print(f'\tOur approach (seconds): {sum(t1) / len(t1): .4f}')
    if enable_GYQ:
        print(f'\tGYQ approach (seconds): {sum(t2) / len(t2): .4f}')

    print('Task 3: time-optimal stabilization')
    t1 = []
    t2 = []
    for r in range(n_runs):
        print(f' -Run {r}')
        M_set = generate_M_set(expected_M_set_size)
        _, t = compute_time_optimal_stabilizer(
            m, n, L, M_set, GraphicalViewSolver)
        print(f'\tOur approach (seconds): {t}')
        t1.append(t)
        if enable_GYQ:
            _, t = compute_time_optimal_stabilizer(m, n, L, M_set, GYQSolver)
            print(f'\tGYQ approach (seconds): {t}')
            t2.append(t)
    print(' -Average time')
    print(f'\tOur approach (seconds): {sum(t1) / len(t1): .4f}')
    if enable_GYQ:
        print(f'\tGYQ approach (seconds): {sum(t2) / len(t2): .4f}')

    print('Task 4: general optimal set stabilization')
    t1 = []
    for r in range(n_runs):
        print(f' -Run {r}')
        M_set = generate_M_set(expected_M_set_size)
        _, t = compute_general_optimal_stabilizer(
            m, n, L, M_set, g, GraphicalViewSolver)
        print(f'\tOur approach (seconds): {t}')
        t1.append(t)
    print(' -Average time')
    print(f'\tOur approach (seconds): {sum(t1) / len(t1): .4f}')

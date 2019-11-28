"""
Benchmark the performance of our approach against others using the Ara operon network.

The GYQ approach refers to
Guo, Yuqian, Pan Wang, Weihua Gui, and Chunhua Yang. "Set stability and set stabilization of Boolean control networks
    based on invariant subsets." Automatica 61 (2015): 106-112.
"""

from algorithm.proposed import GraphicalViewSolver
from algorithm.related_work import GYQSolver
from algorithm.utils import read_network
import random
from typing import List, Set
import time


def compute_LCIS(m: int, n: int, L: List[int], M_set: Set[int], Solver):
    ts = time.time()
    solver = Solver(m, n, L, M_set)
    lcis = solver.compute_largest_control_invariant_subset()
    t = time.time() -ts
    return lcis, t

def check_global_stabilizability(m: int, n: int, L: List[int], M_set: Set[int], Solver):
    ts = time.time()
    solver = Solver(m, n, L, M_set)
    cs = solver.is_set_stabilizable()
    t = time.time() -ts
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
        M_set = {9, 265, 320, 352, 384, 472, 480, 504, 512}
        M_set.update(random.sample(range(1, 2 ** n + 1), k=expected_M_set_size - len(M_set)))
    return M_set


if __name__ == '__main__':
    random.seed(123)    # to reproduce results
    n, m, L = read_network('./networks/ara_operon.txt')
    expected_M_set_size = 40
    print('The size of M set is ', expected_M_set_size)
    n_runs = 10

    print('Task 1: compute LCIS')
    t1 = []
    t2 = []
    for r in range(n_runs):
        print(f' -Run {r}')
        M_set = generate_M_set(expected_M_set_size)
        _, t = compute_LCIS(m, n, L, M_set, GraphicalViewSolver)
        print(f'\tOur approach (seconds): {t}')
        t1.append(t)
        _, t = compute_LCIS(m, n, L, M_set, GYQSolver)
        print(f'\tGYQ approach (seconds): {t}')
        t2.append(t)
    print(' -Average time')
    print(f'\tOur approach (seconds): {sum(t1) / len(t1): .3f}')
    print(f'\tGYQ approach (seconds): {sum(t2) / len(t2): .3f}')

    print('Task 2: check stabilizability')
    t1 = []
    t2 = []
    for r in range(n_runs):
        print(f' -Run {r}')
        M_set = generate_M_set(expected_M_set_size)
        _, t = check_global_stabilizability(m, n, L, M_set, GraphicalViewSolver)
        print(f'\tOur approach (seconds): {t}')
        t1.append(t)
        __, t = check_global_stabilizability(m, n, L, M_set, GYQSolver)
        print(f'\tGYQ approach (seconds): {t}')
        t2.append(t)
    print(' -Average time')
    print(f'\tOur approach (seconds): {sum(t1) / len(t1): .3f}')
    print(f'\tGYQ approach (seconds): {sum(t2) / len(t2): .3f}')

    print('Task 3: time-optimal stabilization')
    t1 = []
    t2 = []
    for r in range(n_runs):
        print(f' -Run {r}')
        M_set = generate_M_set(expected_M_set_size)
        _, t = compute_time_optimal_stabilizer(m, n, L, M_set, GraphicalViewSolver)
        print(f'\tOur approach (seconds): {t}')
        t1.append(t)
        _, t = compute_time_optimal_stabilizer(m, n, L, M_set, GYQSolver)
        print(f'\tGYQ approach (seconds): {t}')
        t2.append(t)
    print(' -Average time')
    print(f'\tOur approach (seconds): {sum(t1) / len(t1): .3f}')
    print(f'\tGYQ approach (seconds): {sum(t2) / len(t2): .3f}')


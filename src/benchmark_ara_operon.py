"""
Benchmark the performance of our approach against others using the Ara operon network.

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

# the following two functions are used to calculate the stage cost according to our definition


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


def compute_general_optimal_stabilizer(m: int, n: int, L: List[int], M_set: Set[int], g: Callable[[int, int], float], Solver):
    ts = time.time()
    solver = Solver(m, n, L, M_set)
    K = solver.compute_custom_optimal_stabilizer(g)
    t = time.time() - ts
    return K, t


def generate_M_set(expected_M_set_size: int):
    M_set = set()
    while len(M_set) < expected_M_set_size:
        M_set = {9, 265, 320, 352, 384, 472, 480, 504, 512}
        M_set.update(random.sample(range(1, 2 ** n + 1),
                                   k=expected_M_set_size - len(M_set)))
    return M_set


if __name__ == '__main__':
    seed = 123
    random.seed(seed)    # to reproduce results
    n_runs = 20
    print(
        f'Experiment with random seed {seed} and {n_runs} runs for each task.')
    n, m, L = read_network('./networks/ara_operon.txt')
    expected_M_set_size = 40
    print('The size of M set is ', expected_M_set_size)

    print('Task 1: compute LCIS')
    t1 = []
    t2 = []
    t3 = []
    for r in range(n_runs):
        print(f' -Run {r}')
        M_set = generate_M_set(expected_M_set_size)
        _, t = compute_LCIS(m, n, L, M_set, GraphicalViewSolver)
        print(f'\tOur approach (seconds): {t}')
        t1.append(t)
        _, t = compute_LCIS(m, n, L, M_set, GYQSolver)
        print(f'\tGYQ approach (seconds): {t}')
        t2.append(t)
        _, t = compute_LCIS(m, n, L, M_set, LRJSolver)
        print(f'\tLRJ approach (seconds): {t}')
        t3.append(t)
    print(' -Average time')
    print(f'\tOur approach (seconds): {sum(t1) / len(t1): .4f}')
    print(f'\tGYQ approach (seconds): {sum(t2) / len(t2): .4f}')
    print(f'\tLRJ approach (seconds): {sum(t3) / len(t3): .4f}')

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
        _, t = compute_time_optimal_stabilizer(
            m, n, L, M_set, GraphicalViewSolver)
        print(f'\tOur approach (seconds): {t}')
        t1.append(t)
        _, t = compute_time_optimal_stabilizer(m, n, L, M_set, GYQSolver)
        print(f'\tGYQ approach (seconds): {t}')
        t2.append(t)
    print(' -Average time')
    print(f'\tOur approach (seconds): {sum(t1) / len(t1): .3f}')
    print(f'\tGYQ approach (seconds): {sum(t2) / len(t2): .3f}')

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
    print(f'\tOur approach (seconds): {sum(t1) / len(t1): .3f}')

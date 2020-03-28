"""
Implementation of existing algorithms.

"""
from typing import List, Set, Union
from .stp import *
from . import bool_algebra as ba
import time


class GYQSolver:
    """
    Guo, Yuqian, Pan Wang, Weihua Gui, and Chunhua Yang. "Set stability and set stabilization of Boolean control networks
    based on invariant subsets." Automatica 61 (2015): 106-112.
    """

    def __init__(self, m: int, n: int, L: Iterable, M_set: Iterable):
        """
        Initialize the solver.

        :param m: number of control inputs
        :param n: number of state variables
        :param L: network transition matrix in condensed form (an integer for each column)
        :param M_set: a given subset of the state space
        """
        self.m = m
        self.M = 2 ** m
        self.n = n
        self.N = 2 ** n
        self.q = len(M_set)
        if isinstance(M_set, set):
            self.M_set = M_set
        else:
            self.M_set = set(M_set)
        self.L = L
        # the one-step controllability matrix (N-by-N), computed by O(MN)
        self.C1 = np.zeros((self.N, self.N), dtype=np.bool_)
        for k in range(self.M):
            blk = L[k * self.N: (k + 1) * self.N]
            for i in range(self.N):
                self.C1[blk[i] - 1][i] = 1
        self.M0 = np.zeros((self.N, self.N), dtype=np.bool_)
        for j in range(self.N):
            if (j + 1) in self.M_set:
                self.M0[j, j] = 1

    def _get_logical_sub_vectors(self, f: Iterable) -> List[int]:
        """
        Given a Boolean vector, get is logical sub-vectors.

        :param f: a Boolean vector
        :return: all logical sub-vectors, each represented by an integer, i.e., i -> \delta_N^i
        """
        r = []
        for i in range(len(f)):
            if f[i]:
                r.append(i + 1)
        return r

    def _compute_Mqc(self) -> np.ndarray:
        """
        Eq. (19)
        """
        Mc = self.M0
        for _ in range(self.q):
            ts = time.time()
            Mc = ba.sp(ba.sp(Mc, self.C1), self.M0)
        return Mc

    def _compute_controllability_matrix(self) -> np.ndarray:
        """
        The matrix C: N-by-N
        """
        N = self.N
        C = self.C1
        Ck = self.C1
        for k in range(2, N + 1):
            Ck = ba.sp(Ck, self.C1)
            C = ba.add(C, Ck)
        return C

    def _get_predecessors(self, i: int) -> List[int]:
        """
        Get the states that can transit to `i` in one step.
        This is the $R^{-1}$ in the paper.

        :param i: a state
        :return: a list of states, each represented by an integer
        """
        ri = self.C1[i - 1]
        return self._get_logical_sub_vectors(ri)

    def _get_predecessors_for_set(self, ss: Iterable) -> Set:
        """
        Get the states that can transit to any state in `ss` in one step.
        This is the $R^{-1}$ in the paper.

        :param i: a state
        :return: a list of states, each represented by an integer
        """
        ps = set()
        for i in ss:
            ps.update(self._get_predecessors(i))
        return ps

    def compute_largest_control_invariant_subset(self) -> Iterable[int]:
        """
        Compute the LCIS, where each state is represented by an integer, i.e., i -> \delta_N^i.

        :return: LCIS
        """
        Mc = self._compute_Mqc()
        rs = ba.col_sum(Mc)
        return self._get_logical_sub_vectors(rs)

    def compute_shortest_transient_period(self, x0: int = None) -> Union[int, None]:
        """
        Get the shortest transient period.
        Proposition 5 (2).

        :param x0: initial state; if 'None', then get the overall shortest transient period.
            Note that the BCN should be stabilizable from `x0` or globally stabilizable.
        :return: the shortest transient period or `None`
        """
        N = self.N
        Mqc = self._compute_Mqc()
        N0 = ba.meet(np.eye(N, dtype=np.bool_), ba.sp(Mqc.T, Mqc))
        Ck = np.eye(N, dtype=np.bool_)
        if x0 is not None:
            for k in range(0, N):
                temp = ba.sp(N0, Ck)
                if np.any(temp[:, x0 - 1]):
                    return k
                Ck = ba.sp(Ck, self.C1)
            else:
                return None
        else:
            # globally reachable, get the maximum one among all states
            states = set(range(1, N + 1))
            ts = {}
            for k in range(0, N):
                temp = ba.sp(N0, Ck)
                determined = []
                for x0 in states:
                    if np.any(temp[:, x0 - 1]):
                        ts[x0] = k
                        determined.append(x0)
                Ck = ba.sp(Ck, self.C1)
                # we don't need to check the determined ones again
                for i in determined:
                    states.remove(i)
            if states:  # some states remain undetermined
                raise Warning(
                    f'BCN is not globally stabilizable! Only valid initial states are considered.')
            return max(ts.values()) if ts else None

    def is_set_stabilizable(self, x0: int = None) -> bool:
        """
        Check set stabilizability.
        Proposition 5 (1).

        :param x0: an initial state; if `None`, then check global set stabilizability.
        :return: a Boolean result
        """
        N = self.N
        Mqc = self._compute_Mqc()
        N0 = ba.meet(np.eye(N, dtype=np.bool_), ba.sp(Mqc.T, Mqc))
        assert N0.shape == (N, N)
        N0C = ba.sp(N0, self._compute_controllability_matrix())
        rs = ba.col_sum(N0C)
        return np.all(rs)

    def _exchange_x_u(self) -> List:
        """
        An equivalent form of the BCN with positions of x and u exchanged: Lm'xu
        See Eq. (30). However, it is very expensive to compute this naively using a swap matrix.
        Complexity: O(MN)

        :return: Lm', condensed form of a N-by-MN matrix
        """
        L_new = [None] * (self.M * self.N)
        for i in range(self.N):
            for k in range(self.M):
                # i transits to j with control k
                Lk = self.L[k * self.N: (k + 1) * self.N]
                j = Lk[i]
                L_new[i * self.M + k] = j
        return L_new

    def compute_time_optimal_stabilizer(self) -> np.ndarray:
        """
        Solve the time optimal set stabilization problem via state feedback.
        A M-by-N Boolean matrix is obtained, whose logical sub-matrix is a feedback matrix.
        First make sure that the LCIS is not empty.
        Proposition 6.

        :return: the logical matrix $\textbf(F)$
        """
        N = self.N
        M = self.M
        N0 = set(self.compute_largest_control_invariant_subset())
        TM = self.compute_shortest_transient_period()
        # compute a list of N's
        Ns = [None] * (TM + 1)
        Ns[0] = N0
        exclude = N0
        for i in range(1, TM + 1):
            Ns[i] = self._get_predecessors_for_set(Ns[i-1]) - exclude
            exclude = exclude | Ns[i]
        # define Boolean matrices
        NB = [None] * (TM + 1)
        for i in range(TM + 1):
            nb = np.zeros((N, N), dtype=np.bool_)
            for j in range(1, N + 1):
                if j in Ns[i]:
                    nb[j - 1, j - 1] = 1
            NB[i] = nb
        # convert BCN to an equivalent and prepare f(x)
        # NOTE: Eq. (30) and (31) are optimized here to reduce computational complexity
        # We don't use the naive matrix product, but exploit the special structure.
        L_tilde = self._exchange_x_u()
        NL = [None] * TM  # NL[i] = N_i * L_tilde
        for i in range(TM):
            NL[i] = np.empty((N, M*N), dtype=np.bool_)
            for col, j in enumerate(L_tilde):
                NL[i][:, col] = NB[i][j - 1]
        # compute F
        F = np.zeros((M, N), dtype=np.bool_)
        for x in range(1, N + 1):
            if x in Ns[0]:
                F[:, x - 1] = ba.col_sum(NL[0][:, (x - 1) * M: x*M]).T
            else:
                # find the set Ni that x lies in
                for i in range(1, TM + 1):
                    if x in Ns[i]:
                        break
                F[:, x - 1] = ba.col_sum(NL[i - 1][:, (x - 1) * M: x*M]).T
        return F


class LRJSolver:
    """
    R. Liu, J. Lu, W. X. Zheng, and J. Kurths, “Output feedback control
    for set stabilization of boolean control networks,” IEEE transactions on
    neural networks and learning systems, 2019.
    """

    def __init__(self, m: int, n: int, L: Iterable, M_set: Iterable):
        """
        Initialize the solver.

        :param m: number of control inputs
        :param n: number of state variables
        :param L: network transition matrix in condensed form (an integer for each column)
        :param M_set: a given subset of the state space
        """
        self.m = m
        self.n = n
        self.M = 2 ** m
        self.N = 2 ** n
        self.L = L
        self.M_set = M_set

    def compute_largest_control_invariant_subset(self) -> Set[int]:
        """
        Compute the LCIS, where each state is represented by an integer, i.e., i -> \delta_N^i.
        This function implements Algorithm 1 of the paper "R. Liu et. al."

        :return: LCIS
        """
        M, N = self.M, self.N
        # get F for each state (i.e., one-step reachable set)
        Fs = {}
        for i in self.M_set:
            F = set()
            for k in range(M):
                blk = self.L[k * N: (k + 1) * N]
                F.add(blk[i - 1])
            Fs[i] = F
        S = self.M_set  # Set
        while True:
            Q = set()
            for i in S:
                # check whether intersection is empty
                is_empty = True
                for j in Fs[i]:
                    if j in S:
                        is_empty = False
                        break
                if is_empty:
                    Q.add(i)
            # remove the states in Q from S
            if Q:
                for i in Q:
                    S.remove(i)
            else:
                break
        return S

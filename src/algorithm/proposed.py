"""
A Graphical View of Set Invariance and Set Stabilization of Boolean Control Networks
"""
from typing import Iterable, Tuple, List, Union, Dict, Callable
import networkx as nx
from collections import deque
from operator import itemgetter
from .fheap import FibonacciHeap

import time


class GraphicalViewSolver:
    """
    Algorithms for set stability/stabilization based on SCC and graph-theoretical approaches.
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
        self.L = L
        self.M_set = M_set if isinstance(M_set, set) else set(M_set)
        self._stg: nx.DiGraph = None
        self._bft: nx.DiGraph = None
        self._ssf: Dict = None
        self._spt: nx.DiGraph = None
        self.M_scc_list: List = None
        self.M_condensation: nx.DiGraph = None
        self.M_stg: nx.DiGraph = None
        self.M = 2 ** m
        self.N = 2 ** n
        self._lcis = None

    def build_stg(self) -> nx.DiGraph:
        """
        Build the STG for a given Boolean control network.

        :param m: number of inputs
        :param n: number of states
        :param Lm: network transition matrix (linear representation)
        :return: a directed graph
        """
        ts = time.time()
        if self._stg is not None:
            return self._stg
        L = self.L
        g = nx.DiGraph()
        N = self.N
        M = self.M
        g.add_nodes_from(i for i in range(1, N + 1))
        for k in range(M):
            blk = L[k * N: (k + 1) * N]
            for i in range(1, N + 1):
                j = blk[i - 1]  # state i -> j by control k
                g.add_edge(i, j)
                if 'Uij' in g.edges[i, j]:
                    g.edges[i, j]['Uij'].append(k + 1)
                else:
                    g.edges[i, j]['Uij'] = [k + 1]
        self._stg = g
        return g

    def build_BFT_for_M_extended_stg(self) -> Union[nx.DiGraph, None]:
        """
        Get the breadth-first tree for M-extended STG.
        If the LCIS is empty, i.e., this BCN is unstabilizable from any initial state, `None` is returned.

        :return: a tree, but represented as a directed graph; or `None`.
        """
        if self._bft is None:
            stg = self.build_stg()
            LCIS = self.compute_largest_control_invariant_subset()
            if not LCIS:
                self._bft = None
                return None
            # we don't need to build the M-extended STG explicitly
            # since only its BFT is useful
            bft = nx.DiGraph()
            Q = deque()
            bft.add_node(0, t_star=0)
            for i in LCIS:
                bft.add_node(i, t_star=1)
                bft.add_edge(0, i)
                Q.append(i)
            g = stg.reverse(copy=False)  # transpose
            while Q:
                i = Q.popleft()
                for j in g.adj[i]:
                    if j not in bft:
                        bft.add_node(j, t_star=bft.nodes[i]['t_star'] + 1)
                        bft.add_edge(i, j)
                        Q.append(j)
            self._bft = bft
        return self._bft

    def is_set_stabilizable(self, x0: int = None) -> bool:
        """
        Check set stabilizability.

        :param x0: an initial state; if `None`, then check global set stabilizability.
        :return: a Boolean result
        """
        bft = self.build_BFT_for_M_extended_stg()
        if bft is None:
            return False
        if x0 is None:
            stg = self.build_stg()
            return len(bft) == len(stg) + 1
        return x0 in bft

    def get_stability_domain(self) -> List[int]:
        """
        Compute the stability domain.

        :return a list of states
        """
        bft = self.build_BFT_for_M_extended_stg()
        if bft is None:
            return []
        # get all nodes in BFT except 0
        sd = list(bft)
        sd.remove(0)
        return sd

    def condense(self) -> Tuple[List, nx.DiGraph]:
        """
        Get the SCC and the condensation of the `M` set in the STG.

        :param _stg: an STG
        :param M: a set of states (vertices)
        :return: the list of SCCs and the acyclic SCC graph (condensation of `_stg`)
        """
        if self._stg is not None:    # we have already built STG G, just extract the induced graph
            sub_graph = self._stg.subgraph(self.M_set)
        else:
            sub_graph = self._build_M_graph()
        scc = nx.strongly_connected_components(sub_graph)
        scc = list(scc)
        self.M_scc_list = scc
        self.M_condensation = nx.condensation(sub_graph, scc)
        self.M_stg = sub_graph
        return self.M_scc_list, self.M_condensation

    def _build_M_graph(self) -> nx.DiGraph:
        """
        Construct the induced graph G[M]
        :return:
        """
        g = nx.DiGraph()
        g.add_nodes_from(self.M_set)
        for k in range(self.M):
            blk = self.L[k * self.N: (k + 1) * self.N]
            for i in self.M_set:
                j = blk[i - 1]  # state i -> j
                if j in self.M_set:
                    g.add_edge(i, j)
        return g

    def compute_largest_control_invariant_subset(self) -> List[int]:
        """
        Calculate the largest control invariant subset (LCIS) by excluding trivial SCCs.

        :return: the LCIS
        """
        if self._lcis is not None:
            return self._lcis
        self.condense()
        is_invariant = [None] * len(self.M_condensation)

        def check_invariance(i: int):
            """
            Check whether the ith scc belongs to the invariant set
            This is the $\kappa$ function in Algorithm 2 of the paper.
            """
            if is_invariant[i] is None:
                scc = self.M_scc_list[i]  # scc: a set
                # whether scc itself is nontrivial
                if len(scc) > 1:
                    is_invariant[i] = True
                else:   # a singleton
                    # has self-loop?
                    v = next(iter(scc))
                    if self.M_stg.has_edge(v, v):
                        is_invariant[i] = True
                    else:
                        # is in a invariance path?
                        for j in self.M_condensation.adj[i]:
                            if check_invariance(j):
                                is_invariant[i] = True
                                break
                        else:
                            is_invariant[i] = False
            return is_invariant[i]

        LCIS = []
        for i in self.M_condensation.nodes():
            if check_invariance(i):
                LCIS.extend(self.M_scc_list[i])
        self._lcis = LCIS
        return LCIS

    def _compute_steady_state_feedback(self) -> Dict:
        """
        Maintain BCN in the LCIS.

        :return: columns of the feedback matrix F (in form of an array) for states in the LCIS, i.e.,
            i -> v, if state i is in the LCIS, then v is the i-th column of F.
        """
        if self._ssf is None:
            self._ssf = {}
            stg = self.build_stg()
            LCIS = self.compute_largest_control_invariant_subset()
            g = stg.subgraph(LCIS)
            for i in g:
                js = g.adj[i]
                # we just choose one j if there are multiple
                j = next(iter(js))
                uij = stg.edges[i, j]['Uij'][0]
                self._ssf[i] = uij
        return self._ssf

    def compute_shortest_transient_period(self, x0: int = None) -> Union[int, None]:
        """
        Get the shortest transient period for set stabilization.
        Theorem 2

        :param x0: an initial state; if `None`, then the shortest transient period of the whole network.
        :return: shortest transient period; or `None`, if the network is not stabilizable from `x0`
        """
        bft: nx.DiGraph = self.build_BFT_for_M_extended_stg()
        if bft is None:
            return None
        if x0 is not None and self.is_set_stabilizable(x0):
            return bft.nodes[x0]['t_star'] - 1
        if x0 is None:
            if not self.is_set_stabilizable():
                raise Warning(
                    f'BCN is not globally stabilizable! Only valid initial states are considered.')
            return max(bft.nodes[i]['t_star'] for i in bft) - 1
        return None

    def compute_time_optimal_stabilizer(self) -> List:
        """
        Solve the time optimal set stabilization problem via state feedback.
        A M-by-N logical matrix (in form of an array) is designed as the feedback matrix.

        :return: the feedback matrix F (in form of an array). If the BCN is not globally stabilizable,
        then some columns of F are undetermined (denoted by `None`).
        """
        BFT: nx.DiGraph = self.build_BFT_for_M_extended_stg()
        F = [None] * self.N
        if BFT is None:
            return F
        STG = self.build_stg()
        # iterate over all edges
        for j, i in BFT.edges:
            if j > 0:
                Uij = STG.edges[i, j]['Uij']
                F[i - 1] = Uij[0]
        # steady state ones
        ssf = self._compute_steady_state_feedback()
        for i, u in ssf.items():
            F[i - 1] = u
        return F

    def build_SPT_for_M_extended_stg(self, g: Callable[[int, int], float]) -> Union[nx.DiGraph, None]:
        """
        Get the shortest path tree for M-extended STG.
        If the LCIS is empty, `None` is returned.

        :param g: the stage cost function
        :return: a tree, but represented as a directed graph; or `None`.
        """
        if self._spt is not None:
            return self._spt
        stg = self.build_stg()
        LCIS = self.compute_largest_control_invariant_subset()
        if not LCIS:
            return None
        # compute the stage cost for each transition in STG, i.e., to build a weighted STG
        for i, j in stg.edges:
            Uij = stg.edges[i, j]['Uij']
            u_star, weight = min(((u, g(i, u))
                                  for u in Uij), key=itemgetter(1))
            stg.edges[i, j].update(u_star=u_star, w=weight)
        tstg = stg.reverse(copy=False)  # transpose
        # we don't need to build the M-extended STG explicitly
        # since only its SPT is useful
        spt = nx.DiGraph()
        spt.add_node(0, d_star=0)
        Q = FibonacciHeap()
        nodes = {}
        rho = {}
        nodes[0] = Q.insert((0, 0))  # (priority, vertex)
        INF = float('inf')
        for v in stg:
            nodes[v] = Q.insert((INF, v))
        while Q.total_nodes > 0:
            di, i = Q.extract_min().data
            if i != 0:
                spt.add_node(i, d_star=di)
                p = rho[i]
                if p == 0:  # the edge connected to v0 needs no control and has zero weight
                    spt.add_edge(p, i, u_star=None, w=0)
                else:
                    spt.add_edge(
                        p, i, u_star=tstg.edges[p, i]['u_star'], w=tstg.edges[p, i]['w'])
            successors = LCIS if i == 0 else tstg.adj[i]
            for j in successors:
                if i == 0:  # the edge connected to v0 needs no control and has zero weight
                    d = di
                else:
                    d = di + tstg.edges[i, j]['w']
                if d < nodes[j].data[0]:
                    Q.decrease_key(nodes[j], (d, j))
                    rho[j] = i
        self._spt = spt
        return spt

    def compute_custom_optimal_stabilizer(self, g: Callable[[int, int], float]) -> List[int]:
        """
        Solve the optimal set stabilization problem via state feedback for the given stage cost.
        A M-by-N logical matrix (in form of an array) is designed as the feedback matrix.

        :return: the feedback matrix F (in form of an array). If the BCN is not globally stabilizable,
        then some columns of F are undetermined (denoted by `None`).
        """
        spt: nx.DiGraph = self.build_SPT_for_M_extended_stg(g)
        F = [None] * self.N
        if spt is None:
            return F
        # transient ones
        for j, i in spt.edges:
            if j > 0:
                u_star = spt.edges[j, i]['u_star']
                F[i - 1] = u_star
        # steady state ones
        ssf = self._compute_steady_state_feedback()
        for i, u in ssf.items():
            F[i - 1] = u
        return F

    def get_trajectory(self, tree: nx.DiGraph, x0: int) -> Tuple[List[int], List[int]]:
        """
        Get a trajectory for state `x0` in the BFT/SPT for set stabilization.

        :param tree: a BFT or an SPT
        :param x0: initial state
        :return: a trajectory and the associated input sequence
        """
        if tree is None or x0 not in tree:
            raise RuntimeError(f'The BCN is not stabilizable from state {x0}')
        j = x0
        p = [x0]
        u = []
        while True:
            # edge (i, j) in the tree: j has only one parent
            i = next(iter(tree.predecessors(j)))
            if i == 0:
                break
            p.append(i)
            try:
                u.append(tree.edges[i, j]['u_star'])
            except:
                u.append(self._stg.edges[j, i]['Uij'][0])
            j = i
        return p, u

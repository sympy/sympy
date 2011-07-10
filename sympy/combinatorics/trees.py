from sympy.core import Basic
from sympy.matrices import Matrix
from sympy.abc import x

from numpy.matlib import matrix

import networkx as nx

class Graph(Basic):
    """
    A Graph is a mathematical construct consisting of
    nodes and edges that connect these nodes. They are
    represented as G = (V, E), where V is the set of
    vertices and E is the set of edges.
    """

    _g = None
    _graph_type = 'undefined'
    _mst = None

    @property
    def graph(self):
        return self._g

    @property
    def graph_type(self):
        return self.graph_type

    @property
    def vertices(self):
        return self.graph.nodes()

    @property
    def edges(self):
        return self.graph.edges()

    @property
    def adjacency_matrix(self):
        return Matrix(nx.adj_matrix(self.graph))

    @property
    def laplacian(self):
        return Matrix(nx.laplacian(self.graph))

    @property
    def normalized_laplacian(self):
        return Matrix(nx.normalized_laplacian(self.graph))

    @property
    def laplacian_spectrum(self):
        return Matrix(nx.laplacian_spectrum(self.graph))

    @property
    def adjacency_spectrum(self):
        return Matrix(nx.adjacency_spectrum(self.graph))

    @property
    def spectral_eigenvector(self):
        return self.adjacency_matrix.eigenvects()

    @property
    def spectral_eigenvalues(self):
        return self.adjacency_matrix.berkowitz_eigenvals()

    @property
    def spectral_polynomial(self):
        return self.adjacency_matrix.berkowitz_charpoly(x)

    @property
    def degree(self):
        return nx.degree(self.graph)

    @property
    def degree_histogram(self):
        return nx.degree_histogram(self.graph)

    @property
    def density(self):
        return nx.density(self.graph)

    @property
    def is_directed(self):
        return nx.is_directed(self.graph)

    @property
    def m(self):
        return nx.number_of_nodes(self.graph)

    @property
    def n(self):
        return nx.number_of_edges(self.graph)

    @property
    def node_iter(self):
        return nx.node_iter(self.graph)

    @property
    def edge_iter(self):
        return nx.edge_iter(self.graph)

    def from_sympy_matrix(self, mat):
        return nx.from_numpy_matrix(matrix(mat))

    def from_dict_of_dicts(self, d):
        return Graph(nx.from_dict_of_dicts(d))

    def from_edge_list(self, e):
        return Graph(nx.from_edgelist(e))

    def from_dict_of_lists(self, l):
        return Graph(nx.from_dict_of_lists(l))

    @property
    def to_dict_of_lists(self):
        return nx.to_dict_of_lists(self.graph)

    @property
    def to_edge_list(self):
        return nx.to_edgelist(self.graph)

    @property
    def to_dict_of_dicts(self):
        return nx.to_dict_of_dicts(self.graph)

    @property
    def average_neighbour_degree(self):
        return nx.average_neighbour_degree(self.graph)

    @property
    def average_neighbour_in_degree(self):
        return nx.average_neighbour_in_degree(self.graph)

    @property
    def average_neighbour_out_degree(self):
        return nx.average_neighbour_out_degree(self.graph)

    @property
    def average_degree_connectivity(self):
        return nx.average_degree_connectivity(self.graph)

    @property
    def average_in_degree_connectivity(self):
        return nx.average_in_degree_connectivity(self.graph)

    @property
    def averate_out_degree_connectivity(self):
        return nx.averate_out_degree_connectivity(self.graph)

    @property
    def k_nearest_neighbours(self):
        return nx.k_nearest_neighbours(self.graph)

    @property
    def dfs_edges(self):
        return nx.dfs_edges(self.graph)

    @property
    def dfs_tree(self):
        return nx.dfs_tree(self.graph)

    @property
    def dfs_predecessors(self):
        return nx.dfs_predecessors(self.graph)

    @property
    def dfs_successors(self):
        return nx.dfs_successors(self.graph)

    @property
    def dfs_preorder_nodes(self):
        return nx.dfs_preorder_nodes(self.graph)

    @property
    def dfs_postorder_nodes(self):
        return nx.dfs_postorder_nodes(self.graph)

    @property
    def dfs_labeled_edges(self):
        return nx.dfs_labeled_edges(self.graph)

    @property
    def bfs_edges(self):
        return nx.bfs_edges(self.graph)

    @property
    def bfs_tree(self):
        return nx.bfs_tree(self.graph)

    @property
    def bfs_successors(self):
        return nx.bfs_successors(self.graph)

    @property
    def bfs_predecessors(self):
        return nx.bfs_predecessors(self.graph)

    @property
    def shortest_path(self):
        return nx.shortest_path(self.graph)

    @property
    def shortest_path_length(self):
        return nx.shortest_path_length(self.graph)

    @property
    def average_shortest_path_length(self):
        return nx.average_shortest_path_length(self.graph)

    @property
    def floyd_warshall(self):
        return nx.floyd_warshall(self.graph)

    @property
    def floyd_warshall_predecessor_and_distance(self):
        return nx.floyd_warshall_predecessor_and_distance(self.graph)

    def astar_path(self, source, target):
        return nx.astar_path(self.graph, source, target)

    def astar_path_length(self, source, target):
        return nx.astar_path_length(self.graph, path, target)

    @property
    def mst(self):
        if self._mst is None:
            self._mst = nx.minimum_spanning_tree(self.graph)
        return self._mst

    @property
    def mst_edges(self):
        return self.mst.edges()

    def maximal_independent_set(self, nodes):
        return nx.maximal_independent_set(self.graph, nodes)

    @property
    def pagerank(self):
        return nx.pagerank(self.graph)

    @property
    def hub_matrix(self):
        return Matrix(nx.hub_matrix(self.graph))

    @property
    def authority_matrix(self):
        return Matrix(nx.authority_matrix(self.graph))

    def __mul__(self, other):
        return Graph(nx.cartesian_product(self.graph, other.graph))

    def __add__(self, other):
        return Graph(nx.union(self.graph, other.graph))

    def __inv__(self):
        return Graph(nx.complement(self.graph))

    def __sub__(self, other):
        return Graph(nx.difference(self.graph, other.graph))

    def __div__(self, other):
        return Graph(nx.intersection(self.graph, other.graph))

    def __mod__(self, other):
        return nx.is_isomorphic(self.graph, other.graph)

    def is_isomorphic(self, other):
        return self % other

    def could_be_isomorphic(self, other):
        """
        Isomorphism checks are expensive and sometimes we just need
        to know if two graphs are definitely NOT isomorphic.
        """
        return nx.faster_could_be_isomorphic(self.graph, other.graph)

    def is_isolate(self, node):
        return nx.is_isolate(self.graph, node)

    @property
    def isolates(self):
        return nx.isolates(self.graph)

    @property
    def center(self):
        return nx.center(self.graph)

    @property
    def diameter(self):
        return nx.diameter(self.graph)

    @property
    def eccentricity(self):
        return nx.eccentricity(self.graph)

    @property
    def periphery(self):
        return nx.periphery(self.graph)

    @property
    def radius(self):
        return nx.radius(self.graph)

    @property
    def is_eulerian(self):
        return nx.is_eulerian(self.graph)

    @property
    def eulerian_circuit(self):
        return nx.eulerian_circuit(self.graph)

    @property
    def is_distance_regular(self):
        return nx.is_distance_regular(self.graph)

    @property
    def intersection_array(self):
        return nx.intersection_array(self.graph)

    def cycle_basis(self):
        return nx.cycle_basis(self.graph)

    @property
    def simple_cycles(self):
        return nx.simple_cycles(self.graph)

    @property
    def is_connected(self):
        return nx.is_connected(self.graph)

    @property
    def number_connected_components(self):
        return nx.number_connected_components(self.graph)

    @property
    def connected_components(self):
        return nx.connected_components(self.graph)

    def node_connected_component(self, n):
        return nx.node_connected_component(self.graph, n)

    @property
    def is_strongly_connected(self):
        return nx.is_strongly_connected(self.graph)

    @property
    def number_of_strongly_connected_components(self):
        return nx.number_strongly_connected_components(self.graph)

    @property
    def strongly_connected_components(self):
        return nx.strongly_connected_components(self.graph)

    @property
    def strongly_connected_component_subgraphs(self):
        return nx.strongly_connected_component_subgraphs(self.graph)

    @property
    def strongly_connected_components_recursive(self):
        return nx.strongly_connected_components_recursive(self.graph)

    @property
    def kosaraju_strongly_connected_components(self):
        return nx.kosaraju_strongly_connected_components(self.graph)

    @property
    def condensation(self):
        return nx.condensation(self.graph)

    @property
    def is_weakly_connected(self):
        return nx.is_weakly_connected(self.graph)

    @property
    def number_weakly_connected_components(self):
        return nx.number_weakly_connected_components(self.graph)

    @property
    def weakly_connected_components(self):
        return nx.weakly_connected_components(self.graph)

    @property
    def weakly_connected_components_subgraphs(self):
        return nx.weakly_connected_components_subgraphs(self.graph)

    def __new__(cls, *args, **kw_args):
        """
        The arguments given are graph type and parameters.
        """
        ret_obj = Basic.__new__(cls, *args, **kw_args)
        if isinstance(args[0], nx.Graph):
            ret_obj._g = args[0]
            return ret_obj
        graph_type = kw_args['graph_type']

        if graph_type == 'star':
            ret_obj._g = nx.star_graph(args[0])
        elif graph_type == 'path':
            ret_obj._g = nx.path_graph(args[0])
        elif graph_type == 'lollipop':
            ret_obj._g = nx.lollipop_graph(args[0], args[1])
        elif graph_type == 'ladder':
            ret_obj._g = nx.ladder_graph(args[0])
        elif graph_type == 'hypercube':
            ret_obj._g = nx.hypercube_graph(args[0])
        elif graph_type == 'grid_2d':
            ret_obj._g = nx.grid_2d_graph(args[0], args[1])
        elif graph_type == 'grid':
            ret_obj._g = nx.grid_graph(args[0])
        elif graph_type == 'empty':
            ret_obj._g = nx.empty_graph(args[0])
        elif graph_type == 'cycle':
            ret_obj._g = nx.cycle_graph(args[0])
        elif graph_type == 'dorogovetsev_goltsev_mendes':
            ret_obj._g = nx.dorogovetsev_goltsev_mendes_graph(args[0])
        elif graph_type == 'circular_ladder':
            ret_obj._g = nx.circular_ladder_graph(args[0])
        elif graph_type == 'complete_bipartite':
            ret_obj._g = nx.complete_bipartite_graph(args[0], args[1])
        elif graph_type == 'barbell':
            ret_obj._g = nx.barbell_graph(args[0], args[1])
        elif graph_type == 'balanced':
            ret_obj._g = nx.balanced_graph(args[0], args[1])
        else:
            ret_obj._g = nx.empty_graph()

        return ret_obj

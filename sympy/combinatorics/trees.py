from sympy.core import Basic
from sympy.matrices import Matrix

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

    def __new__(cls, *args, **kw_args):
        """
        The arguments given are graph type and parameters.
        """
        ret_obj = Basic.__new__(cls, *args, **kw_args)
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

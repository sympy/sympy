"""
This module can be used to solve problems related
to 2D Trusses.
"""

from sympy.core.sympify import sympify



class Truss:
    """
    A Truss is an assembly of members such as beams,
    connected by nodes, that create a rigid structure.
    In engineering, a truss is a structure that
    consists of two-force members only.

    .. image:: truss_explanation.png

    Trusses are extremely important in engineering applications
    and can be seen in numerous real-world applications like bridges.

    Examples
    ========

    There is a Truss consisting of four nodes and five
    members connecting the nodes. A force P acts
    downward on the node D and there also exist pinned
    and roller joints on the nodes A and B respectively.

    .. image:: truss_example.png

    >>> from sympy.physics.continuum_mechanics.truss import Truss
    >>> from sympy import (Symbol, symbols)
    >>> A = Symbol('A')
    >>> B = Symbol('B')
    >>> C = Symbol('C')
    >>> AB, BC, AC = symbols("AB, BC, AC")
    >>> P = Symbol('P')
    >>> t = Truss()
    >>> t.add_node(A, 0, 0)
    >>> t.add_node(B, 2, 2)
    >>> t.add_node(C, 3, 0)
    >>> t.add_member(AB, A, B)
    >>> t.add_member(BC, B, C)
    >>> t.add_member(AC, A, C)
    >>> t.apply_load(A, magnitude=P, direction=90)
    >>> t.add_support(A, type="roller")
    >>> t.add_support(B, type="fixed")
    """

    def __init__(self):
        """
        Initializes the class
        """
        self._nodes = []
        self._members = []
        self._loads = {}
        self._supports = {}
        self._node_labels = [self._nodes[i][0] for i in range(len(self._nodes))]
        self._node_positions = [[self._nodes[i][1], self._nodes[i][2]] for i in range(len(self._nodes))]
        self._node_position_x = [self._node_positions[0] for i in range(len(self._nodes))]
        self._node_position_y = [self._node_positions[1] for i in range(len(self._nodes))]
        self._member_labels = [self._members[i][0] for i in range(len(self._members))]
        self._member_nodes = {}
        self._nodes_occupied = {}
        self._reaction_forces = {}
        self._internal_forces = {}

    @property
    def nodes(self):
        """
        Returns the nodes of the truss along with their positions.
        """
        return self._nodes

    @property
    def node_labels(self):
        """
        Returns the node labels of the truss.
        """
        return self._node_labels

    @property
    def node_positions(self):
        """
        Returns the positions of the nodes of the truss.
        """
        return self._node_positions

    @property
    def members(self):
        """
        Returns the members of the truss along with the start and end points.
        """
        return self._members

    @property
    def member_labels(self):
        """
        Returns the members of the truss along with the start and end points.
        """
        return self._member_labels

    @property
    def supports(self):
        """
        Returns the nodes with provided supports along with the kind of support provided i.e.
        pinned or roller.
        """
        return self._supports

    @property
    def loads(self):
        """
        Returns the loads acting on the truss.
        """
        return self._loads

    @property
    def reaction_forces(self):
        """
        Returns the reaction forces for all supports which are all initialized to 0.
        """
        return self._reaction_forces

    @property
    def internal_forces(self):
        """
        Returns the internal forces for all members which are all initialized to 0.
        """
        return self._internal_forces

    def add_node(self, label, x, y):
        """
        This method adds a node to the truss along with its name/label and its location.

        Parameters
        ==========
        label:  String or a Symbol
            The label for a node. It is the only way to identify a particular node.

        x: Sympifyable
            The x-coordinate of the position of the node.

        y: Sympifyable
            The y-coordinate of the position of the node.

        Examples
        ========

        >>> from sympy.physics.continuum_mechanics.truss import Truss
        >>> t = Truss()
        >>> t.add_node('A', 0, 0)
        >>> t.nodes
        [('A', 0, 0)]
        >>> t.add_node('B', 3, 0)
        >>> t.nodes
        [('A', 0, 0), ('B', 3, 0)]
        """
        x = sympify(x)
        y = sympify(y)

        if label in self._node_labels:
            raise ValueError("Node needs to have a unique label")

        elif x in self._node_position_x and y in self._node_position_y:
            raise ValueError("A node already exists at the given position")

        else :
            self._nodes.append((label, x, y))
            self._node_labels.append(label)
            self._node_positions.append((x, y))
            self._node_position_x.append(x)
            self._node_position_y.append(y)
            self._loads[label] = [[0, 90]]
            self._supports[label] = "none"

    def remove_node(self, label):
        """
        This method removes a node from the truss.

        Parameters
        ==========
        label:  String or Symbol
            The label of the node to be removed.

        Examples
        ========

        >>> from sympy.physics.continuum_mechanics.truss import Truss
        >>> t = Truss()
        >>> t.add_node('A', 0, 0)
        >>> t.nodes
        [('A', 0, 0)]
        >>> t.add_node('B', 3, 0)
        >>> t.nodes
        [('A', 0, 0), ('B', 3, 0)]
        >>> t.remove_node('A')
        >>> t.nodes
        [('B', 3, 0)]
        """
        for i in range(len(self.nodes)):
            if self._node_labels[i] == label:
                x = self._node_position_x[i]
                y = self._node_position_y[i]

        if label not in self._node_labels:
            raise ValueError("No such node exists in the truss")

        else:
            members_duplicate = self._members.copy()
            for member in members_duplicate:
                if label == self._member_nodes[member[0]][0] or label == self._member_nodes[member[0]][1]:
                    self.remove_member(member[0])
            self._nodes.remove((label, x, y))
            self._node_labels.remove(label)
            self._node_positions.remove((x, y))
            self._node_position_x.remove(x)
            self._node_position_y.remove(y)
            self._loads.pop(label)
            self._supports.pop(label)

    def add_member(self, label, start, end):
        """
        This method adds a member between any two nodes in the given truss.

        Parameters
        ==========
        label: String or Symbol
            The label for a member. It is the only way to identify a particular member.

        start: String or Symbol
            The label of the starting point/node of the member.

        end: String or Symbol
            The label of the ending point/node of the member.

        Examples
        ========

        >>> from sympy.physics.continuum_mechanics.truss import Truss
        >>> t = Truss()
        >>> t.add_node('A', 0, 0)
        >>> t.add_node('B', 3, 0)
        >>> t.add_node('C', 2, 2)
        >>> t.add_member('AB', 'A', 'B')
        >>> t.members
        [['AB', 'A', 'B']]
        """

        if start not in self._node_labels or end not in self._node_labels or start==end:
            msg = ("The start and end points of the member must be unique nodes")
            raise ValueError(msg)

        elif label in self._member_labels:
            raise ValueError("A member with the same label already exists for the truss")

        elif self._nodes_occupied.get(tuple([start, end])):
            raise ValueError("A member already exists between the two nodes")

        else:
            self._members.append([label, start, end])
            self._member_labels.append(label)
            self._member_nodes[label] = [start, end]
            self._nodes_occupied[start, end] = True
            self._nodes_occupied[end, start] = True
            self._internal_forces[label] = 0

    def remove_member(self, label):
        """
        This method removes a member from the given truss.

        Parameters
        ==========
        label: String or Symbol
            The label for the member to be removed.

        Examples
        ========

        >>> from sympy.physics.continuum_mechanics.truss import Truss
        >>> t = Truss()
        >>> t.add_node('A', 0, 0)
        >>> t.add_node('B', 3, 0)
        >>> t.add_node('C', 2, 2)
        >>> t.add_member('AB', 'A', 'B')
        >>> t.add_member('AC', 'A', 'C')
        >>> t.add_member('BC', 'B', 'C')
        >>> t.members
        [['AB', 'A', 'B'], ['AC', 'A', 'C'], ['BC', 'B', 'C']]
        >>> t.remove_member('AC')
        >>> t.members
        [['AB', 'A', 'B'], ['BC', 'B', 'C']]
        """
        if label not in self._member_labels:
            raise ValueError("No such member exists in the Truss")

        else:
            self._members.remove([label, self._member_nodes[label][0], self._member_nodes[label][1]])
            self._member_labels.remove(label)
            self._nodes_occupied.pop(tuple([self._member_nodes[label][0], self._member_nodes[label][1]]))
            self._nodes_occupied.pop(tuple([self._member_nodes[label][1], self._member_nodes[label][0]]))
            self._member_nodes.pop(label)
            self._internal_forces.pop(label)

    def change_label(self, label, new_label):
        """
        This method changes the label of a node or a member.

        Parameters
        ==========
        label: String or Symbol
            The label of the node/member for which the label has
            to be changed.

        new_label: String or Symbol
            The new label of the node/member.

        Examples
        ========

        >>> from sympy.physics.continuum_mechanics.truss import Truss
        >>> t = Truss()
        >>> t.add_node('A', 0, 0)
        >>> t.add_node('B', 3, 0)
        >>> t.nodes
        [('A', 0, 0), ('B', 3, 0)]
        >>> t.change_label('A', 'C')
        >>> t.nodes
        [('C', 0, 0), ('B', 3, 0)]
        >>> t.add_member('BC', 'B', 'C')
        >>> t.members
        [['BC', 'B', 'C']]
        >>> t.change_label('BC', 'BC_new')
        >>> t.members
        [['BC_new', 'B', 'C']]
        """
        if label not in self._node_labels and label not in self._member_labels:
            raise ValueError("No such node/member exists")
        if label in self._node_labels:
            for node in self._nodes:
                if node[0] == label:
                    self._node_labels[self._node_labels.index(node[0])] = new_label
                    self._loads[new_label] = self._loads[node[0]]
                    self._loads.pop(node[0])
                    self._supports[new_label] = self._supports[node[0]]
                    self._supports.pop(node[0])
                    if self._supports[new_label] == 'pinned':
                        self._reaction_forces['R_'+str(new_label)+'_x'] = self._reaction_forces['R_'+str(label)+'_x']
                        self._reaction_forces['R_'+str(new_label)+'_y'] = self._reaction_forces['R_'+str(label)+'_y']
                    elif self._supports[new_label] == 'roller':
                        self._reaction_forces['R_'+str(new_label)+'_y'] = self._reaction_forces['R_'+str(label)+'_y']
                    for member in self._members:
                        if member[1] == node[0]:
                            member[1] = new_label
                            self._member_nodes[member[0]] = [new_label, member[2]]
                            self._nodes_occupied[(new_label, member[2])] = True
                            self._nodes_occupied[(member[2], new_label)] = True
                            self._nodes_occupied.pop(tuple([label, member[2]]))
                            self._nodes_occupied.pop(tuple([member[2], label]))
                        elif member[2] == node[0]:
                            member[2] = new_label
                            self._member_nodes[member[0]] = [member[1], new_label]
                            self._nodes_occupied[(member[1], new_label)] = True
                            self._nodes_occupied[(new_label, member[1])] = True
                            self._nodes_occupied.pop(tuple([member[1], label]))
                            self._nodes_occupied.pop(tuple([label, member[1]]))
                    self._nodes[self._nodes.index((label, node[1], node[2]))] = (new_label, node[1], node[2])

        if label in self._member_labels:
            members_duplicate = self._members.copy()
            for member in members_duplicate:
                if member[0] == label:
                    self._member_labels[self.member_labels.index(member[0])] = new_label
                    self._member_nodes[new_label] = [self._member_nodes[label][0], self._member_nodes[label][1]]
                    self._member_nodes.pop(label)
                    self._internal_forces[new_label] = self._internal_forces[label]
                    self._internal_forces.pop(label)
                    self._members[self._members.index([label, member[1], member[2]])] = [new_label, member[1], member[2]]

    def apply_load(self, location, magnitude, direction):
        """
        This method applies an external load at a particular node

        Parameters
        ==========
        location: String or Symbol
            Label of the Node at which load is applied.

        magnitude: Sympifyable
            Magnitude of the load applied. It must always be positive and any changes in
            the direction of the load are not reflected here.

        direction: Sympifyable
            The angle, in degrees, that the load vector makes with the horizontal
            in the counter-clockwise direction . It takes the values 0 to 360,
            inclusive.

        Examples
        ========

        >>> from sympy.physics.continuum_mechanics.truss import Truss
        >>> from sympy import symbols
        >>> t = Truss()
        >>> t.add_node('A', 0, 0)
        >>> t.add_node('B', 3, 0)
        >>> P = symbols('P')
        >>> t.apply_load('A', P, 90)
        >>> t.apply_load('A', P/2, 45)
        >>> t.apply_load('A', P/4, 90)
        >>> t.loads
        {'A': [[5*P/4, 90], [P/2, 45]], 'B': [[0, 90]]}
        """
        magnitude = sympify(magnitude)
        direction = sympify(direction)

        if location not in self.node_labels:
            raise ValueError("Load must be applied at a known node")

        else:
            if direction in [self._loads[location][i][1] for i in range(len(self._loads[location]))]:
                for i in range(len(self._loads[location])):
                    if self._loads[location][i][1] == direction:
                        self._loads[location][i][0] += magnitude
                        break
            else:
                self._loads[location].append([magnitude, direction])

    def add_support(self, location, type):
        """
        This method adds a pinned or roller support at a particular node

        Parameters
        ==========

        location: String or Symbol
            Label of the Node at which support is added.

        type: String
            Type of the support being provided at the node.

        Examples
        ========

        >>> from sympy.physics.continuum_mechanics.truss import Truss
        >>> t = Truss()
        >>> t.add_node('A', 0, 0)
        >>> t.add_node('B', 3, 0)
        >>> t.add_support('A', 'pinned')
        >>> t.supports
        {'A': 'pinned', 'B': 'none'}
        """
        if location not in self._node_labels:
            raise ValueError("Support must be added on a known node")

        else:
            self._supports[location] = type

        if type == "pinned":
            self._reaction_forces['R_'+str(location)+'_x'] = 0
            self._reaction_forces['R_'+str(location)+'_y'] = 0

        elif type == "roller":
            self._reaction_forces['R_'+str(location)+'_y'] = 0
            if 'R_'+str(location)+'_x' in list(self._reaction_forces):
                self._reaction_forces.pop('R_'+str(location)+'_x')

    def remove_support(self, location):
        """
        This method removes support from a particular node

        Parameters
        ==========

        location: String or Symbol
            Label of the Node at which support is to be removed.

        Examples
        ========

        >>> from sympy.physics.continuum_mechanics.truss import Truss
        >>> t = Truss()
        >>> t.add_node('A', 0, 0)
        >>> t.add_node('B', 3, 0)
        >>> t.add_support('A', 'pinned')
        >>> t.supports
        {'A': 'pinned', 'B': 'none'}
        >>> t.remove_support('A')
        >>> t.supports
        {'A': 'none', 'B': 'none'}
        """
        if location not in self._node_labels:
            raise ValueError("No such node exists in the Truss")

        else:
            self._reaction_forces.pop('R_'+str(location)+'_y')
            if self._supports[location] == "pinned":
                self._reaction_forces.pop('R_'+str(location)+'_x')
            self._supports[location] = "none"

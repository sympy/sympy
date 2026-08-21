"""
This module can be used to solve 2D problems with
singularity functions in mechanics.
"""

from sympy.core import Basic, Symbol, symbols
from sympy.core.relational import Eq
from sympy.core.numbers import pi
from sympy.external import import_module
from sympy.functions import SingularityFunction, Piecewise
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import atan2, cos, sin, tan
from sympy.geometry.polygon import deg, rad
from sympy.physics.continuum_mechanics.beam import Beam
from sympy.physics.continuum_mechanics.column import Column
from sympy.plotting import plot
from sympy.simplify import nsimplify, simplify
from sympy.solvers import linsolve
from sympy.integrals import integrate
from sympy.series import limit
from sympy.utilities.lambdify import lambdify

plt = import_module(
    "matplotlib.pyplot",
    import_kwargs={
        "fromlist": [
            "pyplot",
        ]
    },
)

patches = import_module(
    "matplotlib.patches",
    import_kwargs={
        "fromlist": [
            "FancyArrow",
            "Circle",
            "Rectangle",
            "Polygon",
        ]
    },
)

np = import_module(
    "numpy",
    import_kwargs={
        "fromlist": [
            "numpy",
        ]
    },
)


class Member:
    def __init__(self, x1, y1, x2, y2, E, I, A, member_id):
        self.x1 = x1
        self.y1 = y1
        self.x2 = x2
        self.y2 = y2
        self.E = E
        self.I = I
        self.A = A
        self.member_id = member_id

        # properties that are auto computed
        self.length = self._compute_length()
        self.angle = self._compute_angle()
        self.angle_deg = deg(self.angle)
        self.member_loads = []
        self.member_eq = self._compute_eq()

    def _compute_eq(self):
        x1, y1, x2, y2 = self.x1, self.y1, self.x2, self.y2
        x_mem = Symbol("x")
        if x2 == x1:
            return "vertical"
        a = (y2 - y1) / (x2 - x1)
        b = y1 - a * x1
        return a * x_mem + b

    def _compute_length(self):
        """Compute the length of the member."""
        return sqrt((self.x2 - self.x1) ** 2 + (self.y2 - self.y1) ** 2)

    def _compute_angle(self):
        """Compute the angle of the member in radians."""
        return atan2(self.y2 - self.y1, self.x2 - self.x1)

    def __repr__(self):
        return f"Member(ID={self.member_id}, Length={self.length}, Global_Angle={float(self.angle_deg):.02f}deg)"


class Node:
    def __init__(self, x, y, node_type, node_id):
        self.x = x
        self.y = y
        self.node_type = node_type
        self.node_id = node_id

        # loads
        self.node_loads = []

    def __repr__(self):
        return f"Node(ID={self.node_id}, Type={self.node_type}, X={self.x}, Y={self.y})"


class Load:
    def __init__(
        self,
        start_x,
        start_y,
        value,
        global_angle,
        order,
        end_x=None,
        end_y=None,
        load_id=None,
    ):
        self.start_x = start_x
        self.start_y = start_y
        self.value = value  # simplify this breaks plotter
        # self.value = simplify(value)
        self.global_angle = rad(global_angle)
        self.order = order
        self.end_x = end_x
        self.end_y = end_y
        self.load_id = load_id

        self.x_component = self._compute_x_component()
        self.y_component = self._compute_y_component()
        self.applied_to = None
        self.local_start = None
        self.local_end = None

    def __repr__(self):
        return f"Load(applied_to={self.applied_to},X1={self.start_x}, Y1={self.start_y},X2={self.end_x},Y2={self.end_y} Load={self.value}, Global_Angle={self.global_angle}, Order={self.order}, X-Component={self.x_component}, Y-Component={self.y_component})"

    def _compute_x_component(self):
        return self.value * cos(self.global_angle)

    def _compute_y_component(self):
        return self.value * sin(self.global_angle) * -1

    def _compute_local_loc(self):
        pass


class Structure2d:
    """
    Represents a 2D structure composed of members, supports, and loads.

    .. note::
        - A consistent sign convention should be maintained for all forces and reactions.
            The positive direction for the x-axis is to the right, and the positive direction for the y-axis is downwards.
        - The angles for applying loads follow an anticlockwise convention.
            Load applied along the positive x-axis is 0 degrees, along the upward y-axis direction is 90 degrees,
            along the negative x-axis is 180 degrees, and along the downward y-axis direction is 270 degrees.


        Limitations:
            - Only non-branching or non-intersecting structures are supported. (All members must either be connected to ONE other member at each end or to nothing)
            - All E, A, I values must be the same for all members.
            - Members must be added in order of the unwrapping of the structure, from left to right or right to left.
            - Support must be added at the end of the a member.


    Examples
    ========
    There is a structure containing 2 members. A constant distributed load is applied over the length of the first member.
    The structure is supported by a pin support at one end and a roller support at the other end.

    .. plot::
        :context: close-figs
        :format: doctest
        :include-source: True

        >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
        >>> E = 3e4
        >>> I = 1
        >>> A = 1e4
        >>> F = 15
        >>> s = Structure2d()
        >>> s.add_member(x1=0, y1=0, x2=3, y2=4, E=E, I=I, A=A)
        >>> s.add_member(x1=3, y1=4, x2=7, y2=-1, E=E, I=I, A=A)
        >>> s.apply_load(start_x=1.5, start_y=2, value=F, global_angle=0, order=-1)
        >>> s.apply_load(
        ...     start_x=5,
        ...     start_y=1.5,
        ...     value=F / 2,
        ...     global_angle=s.members[1].angle_deg + 270,
        ...     order=0,
        ...     end_x=7,
        ...     end_y=-1,
        ... )
        >>> s.apply_load(
        ...     start_x=0,
        ...     start_y=0,
        ...     value=F * 0.8,
        ...     global_angle=270,
        ...     order=0,
        ...     end_x=3,
        ...     end_y=4,
        ... )
        >>> s.apply_support(x=7, y=-1, type="roller")
        >>> s.apply_support(x=0, y=0, type="pin")
        >>> s.solve_for_reaction_loads()
        {R_h (x=0,y=0): 15/4, R_v (x=0,y=0): -5115/112, R_v (x=7,y=-1): -3285/112}
        >>> s.draw(show_load_values=True) #doctest: +SKIP

    There is a structure containing 3 members. A point load is applied at 1/4 L of the first member.
    A constant distributed load is applied over the second half of the first member and a second distributed load is applied over first hald of the third member.
    The structure is supported by a pin support at one end and a roller support at the other end.

    .. plot::
        :context: close-figs
        :format: doctest
        :include-source: True

        >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
        >>> E = 3e4
        >>> I = 1
        >>> A = 1e4
        >>> F = 15
        >>> s = Structure2d()
        >>> s.add_member(x1=0, y1=0, x2=12, y2=5, E=E, I=I, A=A)
        >>> s.add_member(x1=12, y1=5, x2=12, y2=2, E=E, I=I, A=A)
        >>> s.add_member(x1=12, y1=2, x2=15, y2=2, E=E, I=I, A=A)
        >>> s.apply_load(
        ...     start_x=6,
        ...     start_y=2.5,
        ...     value=20,
        ...     global_angle=270,
        ...     order=0,
        ...     end_x=12,
        ...     end_y=5,
        ... )
        >>> s.apply_load(
        ...     start_x=12,
        ...     start_y=2,
        ...     value=50,
        ...     global_angle=270,
        ...     order=0,
        ...     end_x=13.5,
        ...     end_y=2,
        ... )
        >>> s.apply_load(
        ...     start_x=3,
        ...     start_y=1.25,
        ...     value=100,
        ...     global_angle=s.members[0].angle_deg + 270,
        ...     order=-1,
        ... )
        >>> s.apply_support(x=15, y=2, type="roller")
        >>> s.apply_support(x=0, y=0, type="pin")
        >>> s.solve_for_reaction_loads()
        {R_h (x=0,y=0): -500/13, R_v (x=0,y=0): -20887/156, R_v (x=15,y=2): -1961/12}
        >>> s.draw(show_load_values=True, forced_load_size=2) #doctest: +SKIP

    """

    def __init__(self):
        """
        Initializes the Structure2d class with empty lists for members, supports, nodes, and loads.

        This method sets up the internal structure and prepares it for adding members,
        supports, nodes, and loads.
        """
        self.members = []
        self.supports = []
        self.support_symbols = []
        self.rotation_jumps = []
        self.nodes = []
        self.loads = []
        self.unwrapped_bendpoints = []
        self.unwrapped_loadpoints = []
        self.beam = self._init_beam()
        self.column = self._init_column()
        self.reaction_loads = {}
        self.load_qz = 0
        self.load_qx = 0
        self.V   = None
        self.M   = None
        self.phi = None
        self.uz  = None
        self.N   = None
        self.ux  = None
        self._is_solved = False

    def __repr__(self):
        return f"Structure2d(Members={len(self.members)}, Nodes={len(self.nodes)}, Supports={len(self.supports)})"

    def _init_beam(
        self,
        length=1,
        elastic_modulus=1,
        second_moment=1,
        area=Symbol("A"),
        variable=Symbol("x"),
        base_char="C",
        ild_variable=Symbol("a"),
    ):
        return Beam(
            length=length,
            elastic_modulus=elastic_modulus,
            second_moment=second_moment,
        )

    def _init_column(
            self,
            length=1,
            elastic_modulus=1,
            area=1,
            variable=Symbol('x'),
            base_char='C',
    ):
        return Column(
        length=length,
        elastic_modulus=elastic_modulus,
        area=area

        )


    def add_member(self, x1, y1, x2, y2, E, I, A):
        """
        Adds a member to the structure on a x and y grid.

        This method creates a Member object with the specified coordinates and properties.
        It also adds the corresponding nodes to the structure and updates their types based on the member's location.

        Parameters
        ==========

        x1 : Sympifyable
            X coordinate of the start of the member.

        y1 : Sympifyable
            Y coordinate of the start of the member.

        x2 : Sympifyable
            X coordinate of the end of the member.

        y2 : Sympifyable
            Y coordinate of the end of the member.

        E : Sympifyable
            A SymPy expression representing the Beam's Modulus of Elasticity (Young's modulues).

        I : Sympifyable
            Describes the cross-section of the member via a SymPy expression
            representing the members's second moment of area.

        A : Sympifyable
            Describes the cross-section of the member via a SymPy expression
            representing the members's area.

        Returns:
            Member: The newly created Member object that has been added to the structure.

        Examples
        ========
        A simple horizontal member with coordinates x1=0, y1=0, x2=4, y2=0 is added to the structure.

        >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
        >>> E = 3e4
        >>> I = 1
        >>> A = 1e4
        >>> F = 15
        >>> s = Structure2d()
        >>> s.add_member(0, 0, 4, 0, E, I, A)
        """

        member_id = len(self.members)
        member = Member(x1, y1, x2, y2, E, I, A, member_id)
        self.members.append(member)
        if len(self.members) == 1:
            self.beam.elastic_modulus = E
            self.column.elastic_modulus = E
            self.beam.second_moment = I
            self.column.area = A

        self._add_or_update_node(x1, y1, "fixed", overwrite_type=False)
        self._add_or_update_node(x2, y2, "fixed", overwrite_type=False)


    def _add_or_update_node(self, x, y, new_node_type, overwrite_type=True):
        """Adds a node to the structure at a specified location."""

        for node in self.nodes:
            if node.x == x and node.y == y:
                # If you overwrite here, make sure the correct coordinates are being passed
                if overwrite_type:
                    node.node_type = new_node_type
                return node

        # If the node is new, create it with the correct coordinates
        node_id = len(self.nodes)
        new_node = Node(x, y, new_node_type, node_id)
        self.nodes.append(new_node)
        return new_node


    def _find_applied_to(self, x, y, x_end=None, y_end=None):
        """ Finds the member or node to which the load is applied. """
        # Calculate the midpoint if end coordinates are provided
        if x_end is not None and y_end is not None:
            x_start, y_start = x, y
            x = (x + x_end) / 2
            y = (y + y_end) / 2

        # Check if the point corresponds to a node
        for node in self.nodes:
            if simplify(node.x - x) == 0 and simplify(node.y - y) == 0:
                return f"n_{node.node_id}", 0, None

        # Check if the point corresponds to a member
        for member in self.members:
            # Handle non-vertical members
            if member.member_eq != "vertical":
                member_y = member.member_eq.subs({"x": x})
                if (
                    min(member.x1, member.x2) <= x <= max(member.x1, member.x2)
                ) and simplify(member_y - y) == 0:
                    local_x_start = sqrt((member.x1 - x) ** 2 + (member.y1 - y) ** 2)
                    if x_end is not None and y_end is not None:
                        local_x_start = sqrt(
                            (member.x1 - x_start) ** 2 + (member.y1 - y_start) ** 2
                        )
                        local_x_end = member.length - sqrt(
                            (member.x2 - x_end) ** 2 + (member.y2 - y_end) ** 2
                        )
                    else:
                        local_x_end = None
                    return f"m_{member.member_id}", local_x_start, local_x_end

            # Handle vertical members
            elif simplify(member.x1 - x) == 0 and min(member.y1, member.y2) <= y <= max(
                member.y1, member.y2
            ):
                local_x_start = Abs(member.y1 - y)
                if x_end is not None and y_end is not None:
                    local_x_start = Abs(member.y1 - y_start)
                    local_x_end = member.length - Abs(member.y2 - y_end)
                else:
                    local_x_end = None
                return f"m_{member.member_id}", local_x_start, local_x_end


    def apply_load(
        self, start_x, start_y, value, global_angle, order, end_x=None, end_y=None
    ):
        """
        This method adds up the force load to a particular beam object.

        .. note::
           The ``value`` parameter can be negative to represent a force in the opposite direction. However,
           note that the ``global_angle`` parameter automatically adjusts the sign. This means that when ``value = 10`` and ``global_angle = 0``, it is equivalent to ``value = -10`` and ``global_angle = 180``.

        Parameters
        ==========
        start_x : Sympifyable
            X coordinate of the start of the load.

        start_y : Sympifyable
            Y coordinate of the start of the load.

        value : Sympifyable
            The magnitude of an applied load.

        global_angle : Sympifyable
            The angle of the applied load with respect to the positive global x-axis.

        order : Integer
            The order of the applied load.
            - For point loads, order=-1
            - For constant distributed load, order=0
            - For ramp loads, order=1
            - For parabolic ramp loads, order=2
            - ... so on.

        end_x : Sympifyable, optional
            X coordinate of the end of the load. Defaults to None.

        end_y : Sympifyable, optional
            Y coordinate of the end of the load. Defaults to None.

        Returns:
            Load: The newly created Load object that has been added to the structure.


        Examples
        ========
        A point load of 15 units is applied at coordinates (2, 0) with an angle of 225 degrees from the positive x-axis.

        >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
        >>> E = 3e4
        >>> I = 1
        >>> A = 1e4
        >>> F = 15
        >>> s = Structure2d()
        >>> s.add_member(0, 0, 4, 0, E, I, A)
        >>> s.apply_load(2, 0, F, global_angle=270, order=-1)

        A constant distributed load of 15 units is applied from (2, 0) to (3, 0) with an angle of 225 degrees from the positive x-axis.

        >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
        >>> E = 3e4
        >>> I = 1
        >>> A = 1e4
        >>> F = 15
        >>> s = Structure2d()
        >>> s.add_member(0, 0, 4, 0, E, I, A)
        >>> s.apply_load(2, 0, F, global_angle=225, order=0, end_x=3, end_y=0)
        """


        # Create a new load object
        load_id = len(self.loads)
        load = Load(start_x, start_y, value, global_angle, order, end_x, end_y, load_id)
        self.loads.append(load)

        # Find the member or node to which the load is applied
        if end_x is not None and end_y is not None:
            load.applied_to, load.local_start, load.local_end = self._find_applied_to(
                start_x, start_y, end_x, end_y
            )

        else:
            load.applied_to, load.local_start, load.local_end = self._find_applied_to(
                start_x, start_y
            )

        if "n" in load.applied_to:
            self.nodes[int(load.applied_to.split("_")[1])].node_loads.append(load)
        elif "m" in load.applied_to:
            self.members[int(load.applied_to.split("_")[1])].member_loads.append(load)

        # build load equation
        ########################

        bend_points, member_angles, L = self._unwrap_structure()
        # self._unwrap_structure returns
        """
        load_types supported in the module are -2(moment loads), -1(point loads), 0(distributed loads)
        # Convert these load_types to alexes load_types
        load types in alex algorithm
                # Moment(-2)           -->      'moment_load'
                # Point load(-1)       -->      'vertical_point_load', 'horizontal_point_load'
                # Distributed load(0)  -->      'vertical_distributed_load', 'horizontal_distributed_load'
                # hinges(-3 only for internal usage) --> 'hinge'

        As the members are inclined with some angles with global axis
        every load is projected into members axial and vertical so every
        load order will be split and projected into both vertical, horizontal
        after resolving into their components
        moment is purely applied no resolutions
        point loads are resolved into their respective components and split to both horizontal point load
        and vertical point load
        similarly with distributed loads they are resolved into their respective components and split into vertical
        distributed load and horizontal distributed load
        """
        # for hinge
        if load.order == -3:
            load_p = load.value
            load_values = [nsimplify(load_p)]
            load_points = [nsimplify(self.unwrapped_loadpoints[-1]["locals"][0])]
            load_type = ['hinge']
        # for moment loads and for fixed support moment reaction from apply_support().
        elif load.order == -2:
            qz = load.value
            load_values = [nsimplify(qz)]
            load_points = [nsimplify(self.unwrapped_loadpoints[-1]["locals"][0])]
            load_type = ['moment_load']
        # for point loads and for reaction loads from apply_support()
        elif load.order == -1:
            Fv = load.y_component
            Fh = load.x_component
            load_values = [nsimplify(Fv), nsimplify(Fh)]
            load_points = [
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][0]),
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][0]),
            ]
            load_type = ['vertical_point_load', 'horizontal_point_load']
        # for distributed loads for order 0
        else:
            qv = load.y_component
            qh = load.x_component
            load_values = [nsimplify(qv), nsimplify(-qv), nsimplify(qh), nsimplify(-qh)]
            load_points = [
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][0]),
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][1]),
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][0]),
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][1]),
            ]
            load_type = ['vertical_distributed_load', 'vertical_distributed_load', 'horizontal_distributed_load', 'horizontal_distributed_load']

        # ALEX algo ############################################################################################
        # bendpoints
        for i in range(len(load_values)):
            for j in range(len(bend_points)):
                if load_points[i] == bend_points[-1]:
                    # moment load.
                    if load_type[i] == 'moment_load':
                        self.beam.apply_load(load_values[i], load_points[i], -2)

                    # vertical point load.
                    elif load_type[i] == 'vertical_point_load':
                        self.beam.apply_load(load_values[i] * cos(member_angles[-1]), load_points[i], -1)
                        self.column.apply_load(load_values[i] * -sin(member_angles[-1]), load_points[i], -1)

                    # horizontal point load
                    elif load_type[i] == 'horizontal_point_load':
                        self.beam.apply_load(load_values[i] * sin(member_angles[-1]), load_points[i], -1)
                        self.column.apply_load(load_values[i] * cos(member_angles[-1]), load_points[i], -1)

                    # vertical distributed load
                    elif load_type[i] == 'vertical_distributed_load':
                        self.beam.apply_load(load_values[i] * cos(member_angles[-1]), load_points[i], 0)
                        self.column.apply_load(load_values[i] * -sin(member_angles[-1]), load_points[i], 0)

                    # horizontal distributed load
                    elif load_type[i] == 'horizontal_distributed_load':
                        self.beam.apply_load(load_values[i] * sin(member_angles[-1]), load_points[i], 0)
                        self.column.apply_load(load_values[i] * cos(member_angles[-1]), load_points[i], 0)

                    break
                else:
                    if load_points[i] < bend_points[j]:
                        # moment load.
                        if load_type[i] == 'moment_load':
                            self.beam.apply_load(load_values[i], load_points[i], -2)

                        # vertical point load.
                        elif load_type[i] == 'vertical_point_load':
                            self.beam.apply_load(load_values[i] * cos(member_angles[j - 1]), load_points[i], -1)
                            self.column.apply_load(load_values[i] * -sin(member_angles[j - 1]), load_points[i], -1)

                        # horizontal point load.
                        elif load_type[i] == 'horizontal_point_load':
                            self.beam.apply_load(load_values[i] * sin(member_angles[j - 1]), load_points[i], -1)
                            self.column.apply_load(load_values[i] * cos(member_angles[j - 1]), load_points[i], -1)

                        # vertical distributed load.
                        elif load_type[i] == 'vertical_distributed_load':
                            self.beam.apply_load(load_values[i] * cos(member_angles[j - 1]), load_points[i], 0)
                            self.column.apply_load(load_values[i] * -sin(member_angles[j - 1]), load_points[i], 0)

                        # horizontal distributed load.
                        elif load_type[i] == 'horizontal_distributed_load':
                            self.beam.apply_load(load_values[i] * sin(member_angles[j - 1]), load_points[i], 0)
                            self.column.apply_load(load_values[i] * cos(member_angles[j - 1]), load_points[i], 0)

                        # for hinge inspired from beam module and alex algorithm.
                        elif load_type[i] == 'hinge':
                            E = self.beam.elastic_modulus
                            I = self.beam.second_moment
                            self.beam.apply_load(load_values[i]*E*I,load_points[i],-3)
                        break

        for i in range(len(load_values)):
            for j in range(len(bend_points) - 1):
                if load_points[i] < bend_points[j]:
                    # vertical point load.
                    if load_type[i] == 'vertical_point_load':
                        self.beam.apply_load(
                            load_values[i] * (cos(member_angles[j]) - cos(member_angles[j - 1])), bend_points[j], -1
                        )
                        self.column.apply_load(
                            load_values[i] * (-sin(member_angles[j]) + sin(member_angles[j - 1])), bend_points[j], -1
                        )

                    # horizontal point load.
                    elif load_type[i] == 'horizontal_point_load':
                        self.beam.apply_load(
                            load_values[i] * (sin(member_angles[j]) - sin(member_angles[j - 1])), bend_points[j], -1
                        )
                        self.column.apply_load(
                            load_values[i] * (cos(member_angles[j]) - cos(member_angles[j - 1])), bend_points[j], -1
                        )

                    # vertical distributed load.
                    elif load_type[i] == 'vertical_distributed_load':
                        self.beam.apply_load(
                            load_values[i] * (cos(member_angles[j]) - cos(member_angles[j - 1])), bend_points[j], 0
                        )
                        self.beam.apply_load(
                            load_values[i] * (bend_points[j] - load_points[i]) * (cos(member_angles[j]) - cos(member_angles[j - 1])),
                            bend_points[j],
                            -1,
                        )
                        self.column.apply_load(
                            load_values[i] * (-sin(member_angles[j]) + sin(member_angles[j - 1])), bend_points[j], 0
                        )
                        self.column.apply_load(
                            load_values[i] * (bend_points[j] - load_points[i]) * (-sin(member_angles[j]) + sin(member_angles[j - 1])),
                            bend_points[j],
                            -1,
                        )

                    # horizontal distributed load.
                    elif load_type[i] == 'horizontal_distributed_load':
                        self.beam.apply_load(
                            load_values[i] * (sin(member_angles[j]) - sin(member_angles[j - 1])), bend_points[j], 0
                        )
                        self.beam.apply_load(
                            load_values[i] * (bend_points[j] - load_points[i]) * (sin(member_angles[j]) - sin(member_angles[j - 1])),
                            bend_points[j],
                            -1,
                        )
                        self.column.apply_load(
                            load_values[i] * (cos(member_angles[j]) - cos(member_angles[j - 1])), bend_points[j], 0
                        )
                        self.column.apply_load(
                            load_values[i] * (bend_points[j] - load_points[i]) * (cos(member_angles[j]) - cos(member_angles[j - 1])),
                            bend_points[j],
                            -1,
                        )
        # all vertical projections on to the beam
        self.load_qz = self.beam._load
        # all horizontal projections on to the beam
        self.load_qx = self.column._load

    ###################################################################################################################

    def _unwrap_structure(self):
        """
        Unwraps the structure to a 1D beam for analysis and
          returns all the bend points or joints, angles of members
          with global horizontal axis, length of the un wrapped 1d structure.
        """

        unwrapped_len = 0
        unwrapped_bendpoints = []
        unwrapped_loadpoints = []

        # Process members
        for member in self.members:
            unwrapped_bendpoints.append(
                {
                    "bend_point": [unwrapped_len, unwrapped_len + member.length],
                    "angle": member.angle,
                    "2positions": [(member.x1, member.y1), (member.x2, member.y2)]
                }
            )

            # Process loads on the member
            for load in member.member_loads:
                if load.end_x is not None and load.end_y is not None:
                    unwrapped_loadpoints.append(
                        {
                            "l_id": load.load_id,
                            "locals": [
                                nsimplify(unwrapped_len + load.local_start),
                                nsimplify(unwrapped_len + load.local_end)],
                            "2positions": [(load.start_x, load.start_y), (load.end_x, load.end_y)]
                        }
                    )
                else:
                    unwrapped_loadpoints.append(
                        {
                            "l_id": load.load_id,
                            "locals": [nsimplify(unwrapped_len + load.local_start)],
                            "2positions": [(load.start_x, load.start_y)]
                        }
                    )

            # Process loads at the start node
            x_start, y_start = member.x1, member.y1
            for node in self.nodes:
                if simplify(node.x - x_start) == 0 and simplify(node.y - y_start) == 0:
                    for load in node.node_loads:
                        unwrapped_loadpoints.append(
                            {"l_id": load.load_id, "locals": [nsimplify(unwrapped_len)],
                            "2positions": [(load.start_x, load.start_y)]
                             }
                        )
                    break

            unwrapped_len += member.length

            # Process loads at the end node
            x_end, y_end = member.x2, member.y2
            for node in self.nodes:
                if simplify(node.x - x_end) == 0 and simplify(node.y - y_end) == 0:
                    for load in node.node_loads:
                        unwrapped_loadpoints.append(
                            {"l_id": load.load_id, "locals": [nsimplify(unwrapped_len)],
                             "2positions": [(load.start_x, load.start_y)]
                             }
                        )
                    break

        self.unwrapped_bendpoints = unwrapped_bendpoints
        self.unwrapped_loadpoints = sorted(
            unwrapped_loadpoints, key=lambda x: x["l_id"]
        )
        ### Emulate Alexes input for his algorithm
        bend_points = []
        member_angles = []

        for unwrapped_bendpoint in unwrapped_bendpoints:
            bend_points.append(unwrapped_bendpoint["bend_point"][0])
            member_angles.append(unwrapped_bendpoint["angle"])
        bend_points.append(unwrapped_bendpoints[-1]["bend_point"][1])

        L = unwrapped_len
        self.beam.length = L
        self.column.length = L

        # returns all bend points, member angles , length.
        return bend_points, member_angles, L
        ### Emulate Alexes input for his algorithm


    def _find_unwrapped_position(self, x, y):
        """Finds the unwrapped position based on the x and y from the grid."""
        self._unwrap_structure()

        unwrapped_position_x = 0
        current_length = 0

        for member in self.members:
            if (member.x1 <= x <= member.x2 or member.x2 <= x <= member.x1) and (
                member.y1 <= y <= member.y2 or member.y2 <= y <= member.y1
            ):
                local_x = sqrt((member.x1 - x) ** 2 + (member.y1 - y) ** 2)
                unwrapped_position_x = current_length + local_x
                break

            current_length += member.length

        if unwrapped_position_x == 0:
            unwrapped_position_x = current_length
        # print("Debug - unwrapped position:", unwrapped_position_x)
        return nsimplify(unwrapped_position_x)  # ALSO NEEDS SIMPLIFY


    def apply_support(self, x, y, type="pin"):
        """
        Applies a support at a specified location in the structure.

        This method adds a support node at the given coordinates and applies the corresponding loads
        based on the type of support (pin, roller, or fixed). It updates the beam's boundary conditions
        and marks the loads as support reactions.

        Parameters
        ==========
        x : Sympifyable
            X coordinate of the support location.

        y : Sympifyable
            Y coordinate of the support location.

        type : str
            The type of support to be applied. It can be one of the following (Default is "pin").
            - "pin": Pin support that restricts vertical and horizontal movement.
            - "roller": Roller support that restricts vertical movement.
            - "fixed": Fixed support that restricts vertical, horizontal, and rotational movement.

        Examples
        ========
        A pin support is applied at coordinates (0, 0) and a roller support is applied at coordinates (4, 0).

        >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
        >>> E = 3e4
        >>> I = 1
        >>> A = 1e4
        >>> F = 15
        >>> s = Structure2d()
        >>> s.add_member(0, 0, 4, 0, E, I, A)
        >>> s.apply_load(2, 0, F, global_angle=225, order=0, end_x=3, end_y=0)
        >>> s.apply_support(x=0, y=0, type="pin")
        >>> s.apply_support(x=4, y=0, type="roller")
        """

        support = self._add_or_update_node(x, y, type)
        self.supports.append(support)

        unwarap_x = self._find_unwrapped_position(x, y)

        Rh = Symbol(f"R_h (x={round(x,2)},y={round(y,2)})")
        Rv = Symbol(f"R_v (x={round(x,2)},y={round(y,2)})")
        T = Symbol(f"T (x={round(x,2)},y={round(y,2)})")

        if type == "pin" or type == "roller" or type == "fixed":
            pass
        else:
            raise ValueError(
                "Invalid support type. Choose from 'pin', 'roller', or 'fixed'."
            )

        # Apply support loads based on the type of support
        if type == "pin":
            load_v = -1 * Rv
            load_h = Rh
            # apply vertical support reaction load  as a point  load with 90 degrees which is upwards
            self.apply_load(
                start_x=x, start_y=y, value=load_v, global_angle=90, order=-1
            )
            # apply horizontal support reaction load as a point load with 0 degrees which is right wards
            self.apply_load(
                start_x=x, start_y=y, value=load_h, global_angle=0, order=-1
            )

            # Mark these loads as support reactions
            self.loads[-1].is_support_reaction = True
            self.loads[-2].is_support_reaction = True

            # apply boundary conditions for pin support.
            self.beam.bc_deflection.append((unwarap_x, 0))
            self.column._bc_extension.append(unwarap_x)

            self.support_symbols.append(Rv)
            self.support_symbols.append(Rh)

        elif type == "roller":
            load_v = -1 * Rv
            # apply vertical support reaction load  as a point  load with 90 degrees which is upwards
            self.apply_load(
                start_x=x, start_y=y, value=load_v, global_angle=90, order=-1
            )

            # Mark this load as a support reaction
            self.loads[-1].is_support_reaction = True

            # apply boundary conditions for roller support.
            self.beam.bc_deflection.append((unwarap_x, 0))
            self.support_symbols.append(Rv)

        elif type == "fixed":
            load_t = T
            load_v = -1 * Rv
            load_h = Rh
            # apply support reaction load as a moment load.
            self.apply_load(
                start_x=x, start_y=y, value=load_t, global_angle=0, order=-2
            )
            # apply vertical support reaction load  as a point  load with 90 degrees which is upwards
            self.apply_load(
                start_x=x, start_y=y, value=load_v, global_angle=90, order=-1
            )
            # apply horizontal support reaction load as a point load with 0 degrees which is right wards
            self.apply_load(
                start_x=x, start_y=y, value=load_h, global_angle=0, order=-1
            )

            # Mark these loads as support reactions
            self.loads[-1].is_support_reaction = True
            self.loads[-2].is_support_reaction = True
            self.loads[-3].is_support_reaction = True

            # apply boundary conditions for fixed support.
            self.beam.bc_deflection.append((unwarap_x, 0))
            # This i think sould be slope of the beam at this point beacause 0 is assuming supports and horizontal members
            # unwrap is already called so it should be extaracatble from the unwrapped position
            self.beam.bc_slope.append((unwarap_x, 0))
            self.column._bc_extension.append(unwarap_x)

            self.support_symbols.append(Rv)
            self.support_symbols.append(Rh)
            self.support_symbols.append(T)


    def apply_rotation_hinge(self, x, y):

        unwarap_x = self._find_unwrapped_position(x, y)
        P = Symbol(f"P (x={round(x,2)},y={round(y,2)})")
        load_p = P
        self.apply_load(
                start_x=x, start_y=y, value=load_p, global_angle=0, order=-3
            )
        self.rotation_jumps.append(P)
        self.beam.bc_bending_moment.append((unwarap_x,0))


    def _build_local_displacements(self, uz, ux):

        x = self.beam.variable
        bend_points, member_angles, _ = self._unwrap_structure()
        o0 = member_angles[0] if member_angles else 0

        # local v = projection of (uz, ux) on each segment's local v-axis
        uvz = uz.subs(x, bend_points[0]) * cos(o0)
        uvx = -ux.subs(x, bend_points[0]) * sin(o0)
        for i in range(len(member_angles)):
            a, b, th = bend_points[i], bend_points[i+1], member_angles[i]
            uvz += ((uz - uz.subs(x, a)) * SingularityFunction(x, a, 0) - (uz - uz.subs(x, b)) * SingularityFunction(x, b, 0)) * cos(th)
            uvx += -((ux - ux.subs(x, a)) * SingularityFunction(x, a, 0) - (ux - ux.subs(x, b)) * SingularityFunction(x, b, 0)) * sin(th)
        uv = uvz + uvx

        # local h = projection of (uz, ux) on each segment's local h-axis
        uhz = uz.subs(x, bend_points[0]) * sin(o0)
        uhx = ux.subs(x, bend_points[0]) * cos(o0)
        for i in range(len(member_angles)):
            a, b, th = bend_points[i], bend_points[i+1], member_angles[i]
            uhz += ((uz - uz.subs(x, a)) * SingularityFunction(x, a, 0) - (uz - uz.subs(x, b)) * SingularityFunction(x, b, 0)) * sin(th)
            uhx += ((ux - ux.subs(x, a)) * SingularityFunction(x, a, 0) - (ux - ux.subs(x, b)) * SingularityFunction(x, b, 0)) * cos(th)
        uh = uhz + uhx
        return uv, uh


    def solve_for_reaction_loads(self):
        """
        Solves for the reaction loads in the structure based on the applied loads and applied support reactions.

        Returns:
            dict: A dictionary containing the solved reaction loads for the structure

        Examples
        ========
        >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
        >>> E = 3e4
        >>> I = 1
        >>> A = 1e4
        >>> s = Structure2d()
        >>> s.add_member(x1=0, y1=0, x2=12, y2=5, E=E, I=I, A=A)
        >>> s.add_member(x1=12, y1=5, x2=12, y2=2, E=E, I=I, A=A)
        >>> s.add_member(x1=12, y1=2, x2=15, y2=2, E=E, I=I, A=A)
        >>> s.apply_load(
        ...     start_x=6,
        ...     start_y=2.5,
        ...     value=20,
        ...     global_angle=270,
        ...     order=0,
        ...     end_x=12,
        ...     end_y=5,
        ... )
        >>> s.apply_load(
        ...     start_x=12,
        ...     start_y=2,
        ...     value=50,
        ...     global_angle=270,
        ...     order=0,
        ...     end_x=13.5,
        ...     end_y=2,
        ... )
        >>> s.apply_load(
        ...     start_x=3,
        ...     start_y=1.25,
        ...     value=10,
        ...     global_angle=s.members[0].angle_deg + 270,
        ...     order=-1,
        ... )
        >>> s.apply_support(x=15, y=2, type="roller")
        >>> s.apply_support(x=0, y=0, type="pin")
        >>> s.solve_for_reaction_loads()
        {R_h (x=0,y=0): -3.84615384615385, R_v (x=0,y=0): -70.3141025641026, R_v (x=15,y=2): -143.916666666667}
        """

        C_V, C_M, C_phi, C_uz, C_N, C_ux = symbols('C_V C_M C_phi C_uz C_N C_ux')
        b, c = self.beam, self.column
        x = b.variable
        L = b.length


        E = b.elastic_modulus
        I = b.second_moment
        A = c.area

        qz = self.load_qz
        qx = self.load_qx

        self.V   = integrate(-qz, x) + C_V
        self.M   =  integrate(self.V,  x) + C_M
        self.phi =  integrate(self.M/(E*I), x) + C_phi
        self.uz  = integrate(-self.phi, x) + C_uz


        self.N   = integrate(-qx, x) + C_N
        self.ux  =  integrate(self.N/(E*A), x) + C_ux


        uv, uh = self._build_local_displacements(self.uz, self.ux)


        eqs = []
        eqs.append(limit(self.V, x, 0,  dir='-'))   # V(0-) = 0
        eqs.append(limit(self.V, x, L,  dir='+'))   # V(L+) = 0
        eqs.append(limit(self.N, x, 0,  dir='-'))   # N(0-) = 0
        eqs.append(limit(self.N, x, L,  dir='+'))   # N(L+) = 0

        support_pos = []
        for s in self.supports:
            s_pos = self._find_unwrapped_position(float(s.x), float(s.y))
            support_pos.append((s_pos, s.node_type))

        for s_pos, s_type in support_pos:
            if s_type == "roller":
                eqs.append(uv.subs(x, s_pos))
            elif s_type == "pin":
                eqs.append(uv.subs(x, s_pos))
                eqs.append(uh.subs(x, s_pos))
            elif s_type == "fixed":
                eqs.append(uv.subs(x, s_pos))
                eqs.append(uh.subs(x, s_pos))
                eqs.append(self.phi.subs(x, s_pos))

        left_is_fixed  = any((abs(float(s.x) - float(self.members[0].x1)) < 1e-12 and
                            abs(float(s.y) - float(self.members[0].y1)) < 1e-12 and
                            s.node_type == "fixed") for s in self.supports)

        right_is_fixed = any((abs(float(s.x) - float(self.members[-1].x2)) < 1e-12 and
                            abs(float(s.y) - float(self.members[-1].y2)) < 1e-12 and
                            s.node_type == "fixed") for s in self.supports)

        # left end
        if left_is_fixed:
            eqs.append(limit(self.M, x, 0, dir='-'))
        else:
            eqs.append(self.M.subs(x, 0))

        # right end
        if right_is_fixed:
            eqs.append(limit(self.M, x, L, dir='+'))
        else:
            eqs.append(self.M.subs(x, L))


        for pos, val in b.bc_bending_moment:

            if val != 0:
                eqs.append(self.M.subs(x, pos) - val)
                continue


            if pos > 0:
                eqs.append(limit(self.M, x, pos, dir='-'))
            if pos < L:
                eqs.append(limit(self.M, x, pos, dir='+'))

        support_positions = {p for p, _ in support_pos}
        for pos, val in getattr(b, "bc_slope", []):
            if pos not in support_positions:
                eqs.append(self.phi.subs(x, pos) - val)
        for pos, val in getattr(b, "bc_deflection", []):
            if pos not in support_positions:
                eqs.append(self.uz.subs(x, pos) - val)


        for pos in getattr(c, "_bc_extension", []):
            if pos not in support_positions:
                eqs.append(self.ux.subs(x, pos))


        reaction_syms  = tuple(self.support_symbols)
        rotation_jumps = tuple(self.rotation_jumps)
        unknowns = (C_V, C_M, C_phi, C_uz) + reaction_syms + rotation_jumps + (C_N, C_ux)


        flat_eqs = []
        for e in eqs:
            flat_eqs.append(e.lhs - e.rhs if isinstance(e, Eq) else e)


        sol_tuple = list((linsolve(flat_eqs, unknowns).args)[0])

        sol_map   = dict(zip(unknowns, sol_tuple))


        b._reaction_loads = {s: sol_map[s] for s in reaction_syms if s in sol_map}
        c._reaction_loads = {s: sol_map[s] for s in reaction_syms if s in sol_map}
        self.reaction_loads = dict(b._reaction_loads)



        b._load      = b._load.subs(sol_map)
        self.load_qz = self.load_qz.subs(sol_map)
        c._load      = c._load.subs(sol_map)
        self.load_qx = self.load_qx.subs(sol_map)

        self.V   = self.V.subs(sol_map)
        self.M   = self.M.subs(sol_map)
        self.phi = self.phi.subs(sol_map)
        self.uz  = self.uz.subs(sol_map)
        self.N   = self.N.subs(sol_map)
        self.ux  = self.ux.subs(sol_map)

        self._is_solved = True
        return self.reaction_loads


    def shear_force(self, x=None, y=None):
        """
        Calculates the shear force at a specified point on the beam.

        This method returns the shear force equation if no arguments are provided.
        When both x and y coordinates are provided, it calculates the shear force at the specified (x, y) coordinates.
        (If the x-coordinate is provided, it returns the shear force at the unwrapped location, this for debugging purposes)

        Parameters
        ==========
        x : Sympifyable, optional
            X coordinate of the point to calculate the shear force. Defaults to None.
        y : Sympifyable, optional
            Y coordinate of the point to calculate the shear force. Defaults to None.

        Returns:
            If no arguments are provided: The shear force equation.
            If both x and y coordinates are provided: The shear force at the specified (x, y) coordinates.
            (If only the x-coordinate is given: The shear force at the unwrapped location.)
        """

        if x is None:
            return self.V
        if y is None:
            x_symbol = Symbol("x")
            return self.V.subs(x_symbol, x)
        else:
            unwarap_x = self._find_unwrapped_position(x, y)
            return self.V.subs(Symbol("x"), unwarap_x)


    def plot_shear_force(self):
        """
        Plots the shear force diagram for the beam.

        Returns:
            A plot showing the shear force distribution along the structure.

        Examples
        ========
        Plot the shear force diagram of unwrapped structure.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
            >>> from sympy import symbols
            >>> s = Structure2d()
            >>> E, I, A = symbols('E I A')
            >>> E = 10**4
            >>> I = 10**4
            >>> A = 10**4
            >>> s.add_member(0, 0, 4, 0, E, I, A)
            >>> s.add_member(4, 0, 8, 3, E, I, A)
            >>> s.add_member(8, 3, 11, -1, E, I, A)
            >>> s.apply_load(0, 0, 15, 0, -1)
            >>> s.apply_load(2, 0, 16, 270, -1)
            >>> s.apply_load(0, 0, 6, 270, 0, 4, 0)
            >>> s.apply_load(4, 0, 6, 270, 0, 6, 1.5)
            >>> s.apply_support(11,-1,"fixed")
            >>> s.solve_for_reaction_loads()
            {R_h (x=11,y=-1): -15, R_v (x=11,y=-1): -55, T (x=11,y=-1): -435}
            >>> s.plot_shear_force()  # doctest: +SKIP
        """
        if not self._is_solved :
            raise RuntimeError("Call solve_for_reaction_loads() first.")
        x = self.column.variable
        L = self.beam.length
        return plot(self.V, (x, 0, L),
                    title='Shear Force',
                    xlabel=r'$\mathrm{x}$',
                    ylabel=r'$\mathrm{V}$',
                    line_color='g')


    def axial_force(self, x=None, y=None):
        """
        Calculates the axial force at a specified point on the structure.

        This method returns the axial force equation if no arguments are provided.
        When both x and y coordinates are provided, it calculates the axial force at the specified (x, y) coordinates.
        (If the x-coordinate is provided, it returns the axial force at the unwrapped location)

        Parameters
        ==========
        x : Sympifyable, optional
            X coordinate of the point to calculate the axial force. Defaults to None.
        y : Sympifyable, optional
            Y coordinate of the point to calculate the axial force. Defaults to None.

        Returns:
            If no arguments are provided: The axial force equation.
            If both x and y coordinates are provided: The axial force at the specified (x, y) coordinates.
            (If only the x-coordinate is given: The axial force at the unwrapped location.)
        """

        if x is None:
            return self.N
        if y is None:
            x_symbol = Symbol("x")
            return self.N.subs(x_symbol, x)
        else:
            unwarap_x = self._find_unwrapped_position(x, y)
            return self.N.subs(Symbol("x"), unwarap_x)


    def plot_axial_force(self):
        """
        Plots the axial force diagram for the structure.

        Returns:
            A plot showing the axial force distribution along the structure.

        Examples
        ========
        Plot the axial force diagram on the unwrapped 1d structure after solving reactions.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
            >>> from sympy import symbols
            >>> s = Structure2d()
            >>> E, I, A = symbols('E I A')
            >>> E = 10**4
            >>> I = 10**4
            >>> A = 10**4
            >>> s.add_member(0, 0, 4, 0, E, I, A)
            >>> s.add_member(4, 0, 8, 3, E, I, A)
            >>> s.add_member(8, 3, 11, -1, E, I, A)
            >>> s.apply_load(0, 0, 15, 0, -1)
            >>> s.apply_load(2, 0, 16, 270, -1)
            >>> s.apply_load(0, 0, 6, 270, 0, 4, 0)
            >>> s.apply_load(4, 0, 6, 270, 0, 6, 1.5)
            >>> s.apply_support(11,-1,"fixed")
            >>> s.solve_for_reaction_loads()
            {R_h (x=11,y=-1): -15, R_v (x=11,y=-1): -55, T (x=11,y=-1): -435}
            >>> s.plot_axial_force()  # doctest: +SKIP
        """

        if not self._is_solved :
            raise RuntimeError("Call solve_for_reaction_loads() first.")
        x = self.column.variable
        L = self.beam.length
        return plot(self.N, (x, 0, L),
                    title='Axial Force',
                    xlabel=r'$\mathrm{x}$',
                    ylabel=r'$\mathrm{N(x)}$',
                    line_color='c')


    def bending_moment(self, x=None, y=None):
        """
        Calculates the bending moment at a specified point on the beam.

        This method returns the bending moment equation if no arguments are provided.
        When both x and y coordinates are provided, it computes the bending moment at the specified (x, y) coordinates.
        (If only the x-coordinate is given, it calculates the bending moment at the unwrapped location.)

        Parameters
        ==========
        x : Sympifyable, optional
            X coordinate of the point to calculate the bending moment. Defaults to None.
        y : Sympifyable, optional
            Y coordinate of the point to calculate the bending moment. Defaults to None.

        Returns:
            If no arguments are provided: The bending moment equation.
            If both x and y coordinates are provided: The bending moment at the specified (x, y) coordinates.
            (If only the x-coordinate is given: The bending moment at the unwrapped location.)
        """

        if x is None:
            return self.M
        if y is None:
            x_symbol = Symbol("x")
            return self.M.subs(x_symbol, x)
        else:
            unwarap_x = self._find_unwrapped_position(x, y)
            return self.M.subs(Symbol("x"), unwarap_x)


    def plot_bending_moment(self):
        """
        Plots the bending moment diagram for the beam.

        Returns:
            A plot showing the bending moment distribution along the structure.

        Examples
        ========
        Plot the bending moment diagram on the unwrapped 1d structure after solving reactions.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
            >>> from sympy import symbols
            >>> s = Structure2d()
            >>> E, I, A = symbols('E I A')
            >>> E = 10**4
            >>> I = 10**4
            >>> A = 10**4
            >>> s.add_member(0, 0, 3, 4, E, I, A)
            >>> s.add_member(3, 4, 6, 0, E, I, A)
            >>> s.apply_load(0, 0, 60, 0, 0, 3, 4)
            >>> s.apply_load(3, 4, 60, 0, 0, 6, 0)
            >>> s.apply_support(0, 0,"pin")
            >>> s.apply_support(6, 0,"pin")
            >>> s.solve_for_reaction_loads()
            {R_h (x=0,y=0): -300, R_h (x=6,y=0): -300, R_v (x=0,y=0): 200, R_v (x=6,y=0): -200}
            >>> s.plot_bending_moment()  # doctest: +SKIP
        """
        if not self._is_solved :
            raise RuntimeError("Call solve_for_reaction_loads() first.")
        x = self.column.variable
        L = self.beam.length
        return plot(self.M, (x, 0, L),
                    title='Bending Moment',
                    xlabel=r'$\mathrm{x}$',
                    ylabel=r'$\mathrm{M}$',
                    line_color='b')



    def extension(self, x=None, y=None):
        """
        Calculates the extension at a specified point on the beam.

        This method returns the extension equation if no arguments are provided.
        When both x and y coordinates are provided, it computes the extension at the specified (x, y) coordinates.
        (If only the x-coordinate is given, it calculates the extension at the unwrapped location.)

        Parameters
        ==========
        x : Sympifyable, optional
            X coordinate of the point to calculate the extension. Defaults to None.
        y : Sympifyable, optional
            Y coordinate of the point to calculate the extension. Defaults to None.

        Returns:
            If no arguments are provided: The extension equation.
            If both x and y coordinates are provided: The extension at the specified (x, y) coordinates.
            (If only the x-coordinate is given: The extension at the unwrapped location.)
        """

        if x is None:
            return self.ux
        if y is None:
            x_symbol = Symbol("x")
            return self.ux.subs(x_symbol, x)
        else:
            unwarap_x = self._find_unwrapped_position(x, y)
            return self.ux.subs(Symbol("x"), unwarap_x)


    def plot_extension(self):
        """
        Plots the extension diagram for the beam.

        Returns:
            A plot showing the extension along the structure.


        Examples
        ========
        Plot the extension diagram on the unwrapped 1d structure after solving reactions.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
            >>> from sympy import symbols
            >>> s = Structure2d()
            >>> E, I, A = symbols('E I A')
            >>> E = 10**4
            >>> I = 10**4
            >>> A = 10**4
            >>> s.add_member(0, 0, 3, 4, E, I, A)
            >>> s.add_member(3, 4, 6, 0, E, I, A)
            >>> s.apply_load(0, 0, 60, 0, 0, 3, 4)
            >>> s.apply_load(3, 4, 60, 0, 0, 6, 0)
            >>> s.apply_support(0, 0,"pin")
            >>> s.apply_support(6, 0,"pin")
            >>> s.solve_for_reaction_loads()
            {R_h (x=0,y=0): -300, R_h (x=6,y=0): -300, R_v (x=0,y=0): 200, R_v (x=6,y=0): -200}
            >>> s.plot_extension()  # doctest: +SKIP
        """

        if not self._is_solved :
            raise RuntimeError("Call solve_for_reaction_loads() first.")

        x = self.column.variable
        L = self.beam.length
        return plot(self.ux, (x, 0, L),
                    title='Extension',
                    xlabel=r'$\mathrm{x}$',
                    ylabel=r'$\mathrm{u(x)}$',
                    line_color='m')


    def deflection(self, x=None, y=None):
        """
        Calculates the deflection at a specified point on the beam.

        This method returns the deflection equation if no arguments are provided.
        When both x and y coordinates are provided, it computes the deflection at the specified (x, y) coordinates.
        (If only the x-coordinate is given, it calculates the deflection at the unwrapped location.)

        Parameters
        ==========
        x : Sympifyable, optional
            X coordinate of the point to calculate the deflection. Defaults to None.
        y : Sympifyable, optional
            Y coordinate of the point to calculate the deflection. Defaults to None.

        Returns:
            If no arguments are provided: The deflection equation.
            If both x and y coordinates are provided: The deflection at the specified (x, y) coordinates.
            (If only the x-coordinate is given: The deflection at the unwrapped location.)
        """

        if x is None:
            return self.uz
        if y is None:
            x_symbol = Symbol("x")
            return self.uz.subs(x_symbol, x)
        else:
            unwarap_x = self._find_unwrapped_position(x, y)
            return self.uz.subs(Symbol("x"), unwarap_x)


    def plot_deflection(self):
        """
        Plots the deflection diagram for the beam.

        Returns:
            A plot showing the deflection along the structure.

        Examples
        ========
        Plot the deflection diagram on the unwrapped 1d structure after solving reactions.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
            >>> from sympy import symbols
            >>> s = Structure2d()
            >>> E, I, A = symbols('E I A')
            >>> E = 10**4
            >>> I = 10**4
            >>> A = 10**4
            >>> s.add_member(0, 0, 3, 4, E, I, A)
            >>> s.add_member(3, 4, 6, 0, E, I, A)
            >>> s.apply_load(0, 0, 60, 0, 0, 3, 4)
            >>> s.apply_load(3, 4, 60, 0, 0, 6, 0)
            >>> s.apply_support(0, 0,"pin")
            >>> s.apply_support(6, 0,"pin")
            >>> s.solve_for_reaction_loads()
            {R_h (x=0,y=0): -300, R_h (x=6,y=0): -300, R_v (x=0,y=0): 200, R_v (x=6,y=0): -200}
            >>> s.plot_deflection()  # doctest: +SKIP
        """
        if not self._is_solved :
            raise RuntimeError("Call solve_for_reaction_loads() first.")

        x = self.column.variable
        L = self.beam.length
        return plot(self.uz, (x, 0, L),
                    title='Deflection',
                    xlabel=r'$\mathrm{x}$',
                    ylabel=r'$\mathrm{u(z)}$',
                    line_color='r')


    def _plot_expression_on_structure(self, expr, *, factor=75.0, show_values=True, title=None, _color='tab:red', scale_text=None):
        # helper plotting function to reduce duplication of code

        SAMPLES_PER_MEMBER = 60
        KEEP_FRACTION      = 1/10
        LINE_WIDTH         = 1.0
        ZERO_TOL           = 1e-9
        LABEL_DECIMALS     = 2
        LABEL_FMT          = "[{v:.2f}]"
        LABEL_OFFSET       = 0.2

        label_mode = 'ends' if bool(show_values) else 'none'

        x = self.beam.variable
        self._unwrap_structure()

        fig, ax = plt.subplots()
        ax.set_aspect("equal")

        # draw structure
        for m in self.members:
            ax.plot([float(m.x1), float(m.x2)],
                    [float(m.y1), float(m.y2)],
                    color="black", linewidth=2, zorder=10)

        def _eval_sub(e, xi):
            return float(e.subs(x, xi).evalf())

        def _eval_limit(e, xi, side):
            return float(limit(e, x, xi, dir=side))

        stride = max(1, int(np.ceil(1.0 / float(KEEP_FRACTION))))

        cum = 0.0
        for m in self.members:
            Lm = float(m.length)
            if Lm <= 0:
                cum += Lm
                continue

            dx = float(m.x2 - m.x1); dy = float(m.y2 - m.y1)
            mag = (dx*dx + dy*dy)**0.5 or 1e-15
            tx, ty = dx/mag, dy/mag
            nx, ny = -ty, tx

            n   = max(2, int(SAMPLES_PER_MEMBER))
            xs  = np.linspace(cum, cum + Lm, n)
            bx  = float(m.x1) + (xs - cum) * tx
            by  = float(m.y1) + (xs - cum) * ty

            vals       = np.array([_eval_sub(expr, s) for s in xs], dtype=float)
            vals[0]    = _eval_limit(expr, cum,    '+')
            vals[-1]   = _eval_limit(expr, cum+Lm, '-')

            scaled     = vals / float(factor) if factor else vals
            scaled[np.abs(scaled) < float(ZERO_TOL)] = 0.0

            txs = bx + scaled * nx
            tys = by + scaled * ny

            # Plot individual line segments
            for j in range(0, n, stride):
                ax.plot([bx[j], txs[j]], [by[j], tys[j]], color=_color, linewidth=1.0, zorder=20)
            ax.plot([bx[0],  txs[0]],  [by[0],  tys[0]],  color=_color, linewidth=1.0, zorder=20)
            ax.plot([bx[-1], txs[-1]], [by[-1], tys[-1]], color=_color, linewidth=1.0, zorder=20)

            # Plot the expression line
            ax.plot(txs, tys, color=_color, linewidth=LINE_WIDTH, zorder=21)

            # Add filled polygon between member line (bx, by) and expression line (txs, tys)
            polygon_coords = (
                [(bx[i], by[i]) for i in range(n)] +  # Member line points
                [(txs[i], tys[i]) for i in range(n-1, -1, -1)]  # Expression line points in reverse
            )
            polygon = plt.Polygon(
                polygon_coords,
                closed=True,
                fill=True,
                color=_color,
                edgecolor=None,
                alpha=0.25,  # Semi-transparent fill, similar to distributed load in draw method
                zorder=19,   # Place below tip and expression line
            )
            ax.add_patch(polygon)

            if label_mode == 'ends':
                for j in (0, n-1):
                    v_raw = vals[j]
                    if abs(scaled[j]) < float(ZERO_TOL):
                        v_raw = 0.0
                    sign = 1 if (j % 2 == 0) else -1
                    lx = txs[j] + LABEL_OFFSET * nx
                    ly = tys[j] + sign*LABEL_OFFSET * ny
                    txt = LABEL_FMT.format(v=round(v_raw, int(LABEL_DECIMALS)))
                    ax.text(lx, ly, txt, fontsize=8, color="black",
                            zorder=30, ha='center', va='center', bbox={
                            "facecolor": "lightgrey",
                            "alpha": 0.6,
                            "edgecolor": "none",
                        },
                    )

            cum += Lm

        core = title
        ax.set_title(core)
        ax.set_xlabel("x")
        ax.set_ylabel("y")
        ax.grid(True, zorder=5)
        scale_text_ = scale_text
        if scale_text:
            ax.text(
                0.99, 0.99, scale_text_,
                transform=ax.transAxes,
                fontsize=7,
                ha="right", va="top",
            )

        return fig, ax


    def plot_shear_force_on_structure(self, *, factor=75.0, show_values=True):
        """
        Plots the shear force on the structure geometry (member-wise).

        This method draws short ticks representing magnitude of shear force/factor at the point perpendicular to each member
          the tick length is V/factor.


        Parameters
        ==========
        factor : float(optional)
            Visual scale for tick length (offset = V / factor). Default is 75.0.
            This is to adjust the magnitude lengths for better visual analysis.
        show_values : bool, optional
            If True, prints shear values at member ends. Default is True.


        Examples
        ========
        Plot the shear diagram on the structure after solving reactions.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
            >>> from sympy import symbols
            >>> s = Structure2d()
            >>> E, I, A = symbols('E I A')
            >>> E = 10**4
            >>> I = 10**4
            >>> A = 10**4
            >>> s.add_member(0, 0, 3, 4, E, 10**4, 10**4)
            >>> s.add_member(3, 4, 7, 1, E, I, A)
            >>> s.apply_load(3, 4, 60, 270, -1)
            >>> s.apply_load(3, 4, 18, 0, 0, 7, 1)
            >>> s.apply_support(0, 0,"fixed")
            >>> s.apply_support(7, 1,"roller")
            >>> s.solve_for_reaction_loads()
            {R_h (x=0,y=0): -90, R_v (x=0,y=0): -44421/1120, R_v (x=7,y=1): -22779/1120, T (x=0,y=0): 42021/160}
            >>> s.plot_shear_force_on_structure(factor=75.0, show_values=True)  # doctest: +SKIP
        """
        if not self._is_solved :
            raise RuntimeError("Call solve_for_reaction_loads() first.")
        return self._plot_expression_on_structure(self.shear_force(),
                                                factor=factor,
                                                show_values=show_values,
                                                title="Shear diagram",_color='tab:green',scale_text="(x, y in metre)\n (Shear in kN)"
)


    def plot_bending_moment_on_structure(self, *, factor=150.0, show_values=True):
        """
        Plots the bending moment on the structure geometry (member-wise).

        This method draws short ticks representing magnitude of bending moment/factor at the point perpendicular to each member
          the tick length is M/factor.


        Parameters
        ==========
        factor : float(optional)
            Visual scale for tick length (offset = M / factor). Default is 150.0.
            This is to adjust the magnitude lengths for better visual analysis.
        show_values : bool, optional
            If True, prints moment values at member ends. Default is True.


        Examples
        ========
        Plot the bending moment diagram on the structure after solving reactions.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
            >>> from sympy import symbols
            >>> s = Structure2d()
            >>> E, I, A = symbols('E I A')
            >>> E = 10**4
            >>> I = 10**4
            >>> A = 10**4
            >>> s.add_member(0, 0, 3, 4, E, 10**4, 10**4)
            >>> s.add_member(3, 4, 7, 1, E, I, A)
            >>> s.apply_load(3, 4, 60, 270, -1)
            >>> s.apply_load(3, 4, 18, 0, 0, 7, 1)
            >>> s.apply_support(0, 0,"fixed")
            >>> s.apply_support(7, 1,"roller")
            >>> s.solve_for_reaction_loads()
            {R_h (x=0,y=0): -90, R_v (x=0,y=0): -44421/1120, R_v (x=7,y=1): -22779/1120, T (x=0,y=0): 42021/160}
            >>> s.plot_bending_moment_on_structure(factor=150.0, show_values=True)  # doctest: +SKIP
        """
        if not self._is_solved :
            raise RuntimeError("Call solve_for_reaction_loads() first.")
        return self._plot_expression_on_structure(self.bending_moment(),
                                                factor=factor,
                                                show_values=show_values,
                                                title="Bending moment diagram", _color='tab:blue', scale_text="(x, y in metre)\n (Bending in kN/metre)")


    def plot_axial_force_on_structure(self, *, factor=75.0, show_values=True):
        """
        Plots the axial force on the structure geometry (member-wise).

        This method draws short ticks representing magnitude of axial force/factor at the point perpendicular to each member
          the tick length is N/factor.


        Parameters
        ==========
        factor : float(optional)
            Visual scale for tick length (offset = V / factor). Default is 75.0.
            This is to adjust the magnitude lengths for better visual analysis.
        show_values : bool, optional
            If True, prints axial force values at member ends. Default is True.


        Examples
        ========
        Plot the axial force diagram on the structure after solving reactions.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
            >>> from sympy import symbols
            >>> s = Structure2d()
            >>> E, I, A = symbols('E I A')
            >>> E = 10**4
            >>> I = 10**4
            >>> A = 10**4
            >>> s.add_member(0, 0, 3, 4, E, 10**4, 10**4)
            >>> s.add_member(3, 4, 7, 1, E, I, A)
            >>> s.apply_load(3, 4, 60, 270, -1)
            >>> s.apply_load(3, 4, 18, 0, 0, 7, 1)
            >>> s.apply_support(0, 0,"fixed")
            >>> s.apply_support(7, 1,"roller")
            >>> s.solve_for_reaction_loads()
            {R_h (x=0,y=0): -90, R_v (x=0,y=0): -44421/1120, R_v (x=7,y=1): -22779/1120, T (x=0,y=0): 42021/160}
            >>> s.plot_axial_force_on_structure(factor=75.0, show_values=True)  # doctest: +SKIP
        """
        if not self._is_solved :
            raise RuntimeError("Call solve_for_reaction_loads() first.")
        return self._plot_expression_on_structure(self.axial_force(),
                                                factor=factor,
                                                show_values=show_values,
                                                title="Axial force diagram", _color='tab:orange', scale_text="(x, y in metre)\n (Axial in kN)")


    def _build_geometry_functions(self):

        x = self.beam.variable
        aa, oo, _ = self._unwrap_structure()


        x0 = self.members[0].x1
        y0 = self.members[0].y1

        h = x0
        v = y0

        for i in range(len(oo)):
            a = aa[i]
            b = aa[i+1]
            th = oo[i]
            gate_len = (x - a) * SingularityFunction(x, a, 0) - (x - b) * SingularityFunction(x, b, 0)
            h += gate_len * cos(th)
            v += gate_len * sin(th)
        return h, v


    def plot_deformation_on_structure(self, factor=1/10000.0):
        """
        Plots the deformation on the structure geometry.

        Parameters
        ==========
        factor : float(optional)
            Visual scale for deformation.
            This is to make shift in the deformed structure for better visual analysis.

        Examples
        ========
        Plot the deformation on the structure after solving reactions.

        >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
        >>> E, I, A = 3e4, 1, 1e4
        >>> F = 15
        >>> s = Structure2d()
        >>> s.add_member(x1=0, y1=0, x2=4, y2=0, E=E, I=I, A=A)
        >>> s.apply_load(start_x=2, start_y=0, value=F, global_angle=270, order=-1)
        >>> s.apply_support(x=0, y=0, type="pin")
        >>> s.apply_support(x=4, y=0, type="roller")
        >>> s.solve_for_reaction_loads()
        {R_h (x=0,y=0): 0, R_v (x=0,y=0): -7.5, R_v (x=4,y=0): -7.5}
        >>> s.plot_deformation_on_structure(factor=75.0)  # doctest: +SKIP
        """
        if not self._is_solved :
            raise RuntimeError("Call solve_for_reaction_loads() first.")


        x = self.beam.variable
        L = float(self.beam.length)

        # symbolic geometry + local displacements
        h, v = self._build_geometry_functions()
        uv, uh = self._build_local_displacements(self.uz, self.ux)

        h_np  = lambdify(x, h.rewrite(Piecewise))
        v_np  = lambdify(x, v.rewrite(Piecewise))
        uv_np = lambdify(x, uv.rewrite(Piecewise))
        uh_np = lambdify(x, uh.rewrite(Piecewise))

        s = np.linspace(0.0, L, int(800))

        base_h = np.array(h_np(s), dtype=float).reshape(-1)
        base_v = np.array(v_np(s), dtype=float).reshape(-1)


        def_h  = base_h + uh_np(s) / float(factor if factor else 1.0)
        def_v  = base_v - uv_np(s) / float(factor if factor else 1.0)

        # plot
        fig, ax = plt.subplots()
        ax.set_aspect("equal")

        # draw members as solid black for reference
        for m in self.members:
            ax.plot([float(m.x1), float(m.x2)],
                    [float(m.y1), float(m.y2)],
                    linewidth=2)

        # baseline polyline & deformed
        ax.plot(base_h, base_v, linewidth=2, label='structure')
        ax.plot(def_h,  def_v,  linewidth=2, label='deformed')

        ax.set_xlabel('x'); ax.set_ylabel('y')

        ax.spines['right']
        ax.spines['top']
        ax.spines['bottom']
        ax.spines['left']
        ax.grid(True); ax.axis('scaled')

        ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), ncol=2)
        return fig, ax


    def summary(self, verbose=True, round_digits=None):
        """
        Provides a summary of the structure, including reaction loads and points of interest.

        This method prints a formatted summary of the structure's members, nodes, and supports,
        along with the reaction loads and points of interest related to bending moments shear forces and axial forces.
        It allows for a quick overview of the structural state and key values for analysis.

        Parameters
        ==========
        verbose : bool, optional
            If True, the summary will be printed. Defaults to True.
        round_digits : int, optional
            The number of decimal places to round the values in the summary. Defaults to None (exact solution).
        Returns:
            A formatted summary of the structure's key values.

        Examples
        ========
        >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
        >>> E = 3e4
        >>> I = 1
        >>> A = 1e4
        >>> F = 10
        >>> s = Structure2d()
        >>> s.add_member(0, 0, 4, 0, E, I, A)
        >>> s.apply_load(2, 0, F, global_angle=270, order=-1)
        >>> s.apply_support(x=0, y=0, type="pin")
        >>> s.apply_support(x=4, y=0, type="roller")
        >>> s.solve_for_reaction_loads()
        {R_h (x=0,y=0): 0, R_v (x=0,y=0): -5.0, R_v (x=4,y=0): -5.0}
        >>> s.summary(round_digits=2)
        ===================== Structure Summary =====================
        <BLANKLINE>
        Reaction Loads:
        R_v   (x=0.00,y=0.00)  (unwrapped x=0.00)          = -5.0
        R_h   (x=0.00,y=0.00)  (unwrapped x=0.00)          = 0.0
        R_v   (x=4.00,y=0.00)  (unwrapped x=4.00)          = -5.0
        <BLANKLINE>
        Points of Interest - Bending Moment:
        bending_moment at (x=0.00,y=0.00)  (unwrapped x=0.00) = 0.0
        bending_moment at (x=2.00,y=0.00)  (unwrapped x=2.00) = 10.0
        bending_moment at (x=4.00,y=0.00)  (unwrapped x=4.00-) = 0.0
        <BLANKLINE>
        Points of Interest - Shear Force:
        shear_force at (x=0.00,y=0.00)  (unwrapped x=0.00+) = 5.0
        shear_force at (x=2.00,y=0.00)  (unwrapped x=2.00-) = 5.0
        shear_force at (x=2.00,y=0.00)  (unwrapped x=2.00+) = -5.0
        shear_force at (x=4.00,y=0.00)  (unwrapped x=4.00-) = -5.0
        <BLANKLINE>
        Points of Interest - Axial Force:
        axial_force at (x=0.00,y=0.00)  (unwrapped x=0.00+) = 0.0
        axial_force at (x=2.00,y=0.00)  (unwrapped x=2.00-) = 0.0
        axial_force at (x=2.00,y=0.00)  (unwrapped x=2.00+) = 0.0
        axial_force at (x=4.00,y=0.00)  (unwrapped x=4.00-) = 0.0
        """

        title = "Structure Summary"
        line_length = 60

        print(
            f'{"=" * ((line_length - len(title)) // 2)} {title} {"=" * ((line_length - len(title)) // 2)}'
        )

        if not self._is_solved :
            raise RuntimeError("Call solve_for_reaction_loads() first.")

        if verbose:
            self._print_reaction_loads(round_digits)
            self._print_points_of_interest(round_digits)


    def _print_reaction_loads(self, round_digits):
        import re
        print("\nReaction Loads:")
        for key, value in self.reaction_loads.items():
            match = re.match(r"(\w+.*?)\s*\(x=([\d.-]+),y=([\d.-]+)\)", str(key))

            support_name = match.group(1).strip()
            support_x_loc = float(match.group(2))
            support_y_loc = float(match.group(3))
            unwrapped_xl_loc = self._find_unwrapped_position(support_x_loc, support_y_loc)
            if round_digits is not None and not value.has(Symbol):
                value = round(float(value), round_digits)
            support_str = f"{support_name:<5} (x={support_x_loc:.2f},y={support_y_loc:.2f})  (unwrapped x={unwrapped_xl_loc:.2f})"
            print(f"{support_str:<50} = {value}")


    def _print_points_of_interest(self, round_digits):
        dx = 1e-6
        print("\nPoints of Interest - Bending Moment:")
        bend_points = []
        for item in self.unwrapped_bendpoints:
            bend_points.extend(item["bend_point"])
        moment_points = sorted([
            nsimplify(p)
            for item in self.unwrapped_bendpoints
            for p in (item["bend_point"][0], item["bend_point"][1], (item["bend_point"][0] + item["bend_point"][1]) / 2)
        ])
        for point in moment_points:
            coords = None
            for item in self.unwrapped_bendpoints:
                if nsimplify(point) in [nsimplify(item["bend_point"][0]), nsimplify(item["bend_point"][1])]:
                    idx = 0 if nsimplify(point) == nsimplify(item["bend_point"][0]) else 1
                    coords = item["2positions"][idx]
                    break
                elif nsimplify(point) == nsimplify((item["bend_point"][0] + item["bend_point"][1]) / 2):
                    x_mid = (item["2positions"][0][0] + item["2positions"][1][0]) / 2
                    y_mid = (item["2positions"][0][1] + item["2positions"][1][1]) / 2
                    coords = (x_mid, y_mid)
                    break
            if coords is None:
                coords = (float(point), 0.0)
            x_loc, y_loc = coords
            unwrapped_x_loc = float(point)
            if point == moment_points[-1]:
                bending_moment_value = self.bending_moment(point - dx)
                if round_digits is not None and not bending_moment_value.has(Symbol):
                    bending_moment_value = round(float(bending_moment_value), round_digits)
                string_text = f"bending_moment at (x={x_loc:.2f},y={y_loc:.2f})  (unwrapped x={unwrapped_x_loc:.2f}-)"
                print(f"{string_text:<50} = {bending_moment_value}")
            else:
                bending_moment_value = self.bending_moment(nsimplify(point))
                if round_digits is not None and not bending_moment_value.has(Symbol):
                    bending_moment_value = round(float(bending_moment_value), round_digits)
                string_text = f"bending_moment at (x={x_loc:.2f},y={y_loc:.2f})  (unwrapped x={unwrapped_x_loc:.2f})"
                print(f"{string_text:<50} = {bending_moment_value}")

        print("\nPoints of Interest - Shear Force:")
        load_points = []
        for item in self.unwrapped_loadpoints:
            if len(item["locals"]) == 1:
                load_points.append((item["locals"][0], item["2positions"][0]))
        for item in self.unwrapped_bendpoints:
            load_points.extend([(p, c) for p, c in zip(item["bend_point"], item["2positions"])])
        load_points = sorted(set(load_points), key=lambda x: x[0])
        for point, coords in load_points:
            x_loc, y_loc = coords
            unwrapped_x_loc = float(point)
            shear_force_value_minus = self.shear_force(point - dx)
            shear_force_value_plus = self.shear_force(nsimplify(point))
            if round_digits is not None:
                if isinstance(shear_force_value_minus, Basic) and not shear_force_value_minus.has(Symbol):
                    shear_force_value_minus = round(float(shear_force_value_minus), round_digits)
                if isinstance(shear_force_value_plus, Basic) and not shear_force_value_plus.has(Symbol):
                    shear_force_value_plus = round(float(shear_force_value_plus), round_digits)
            if point == load_points[0][0]:
                string_text = f"shear_force at (x={x_loc:.2f},y={y_loc:.2f})  (unwrapped x={unwrapped_x_loc:.2f}+)"
                print(f"{string_text:<50} = {shear_force_value_plus}")
            elif point == load_points[-1][0]:
                string_text = f"shear_force at (x={x_loc:.2f},y={y_loc:.2f})  (unwrapped x={unwrapped_x_loc:.2f}-)"
                print(f"{string_text:<50} = {shear_force_value_minus}")
            else:
                string_text = f"shear_force at (x={x_loc:.2f},y={y_loc:.2f})  (unwrapped x={unwrapped_x_loc:.2f}-)"
                print(f"{string_text:<50} = {shear_force_value_minus}")
                string_text = f"shear_force at (x={x_loc:.2f},y={y_loc:.2f})  (unwrapped x={unwrapped_x_loc:.2f}+)"
                print(f"{string_text:<50} = {shear_force_value_plus}")

        print("\nPoints of Interest - Axial Force:")
        axial_points = load_points
        for point, coords in axial_points:
            x_loc, y_loc = coords
            unwrapped_x_loc = float(point)
            axial_value_minus = self.axial_force(point - dx)
            axial_value_plus = self.axial_force(nsimplify(point))
            if round_digits is not None:
                if isinstance(axial_value_minus, Basic) and not axial_value_minus.has(Symbol):
                    axial_value_minus = round(float(axial_value_minus), round_digits)
                if isinstance(axial_value_plus, Basic) and not axial_value_plus.has(Symbol):
                    axial_value_plus = round(float(axial_value_plus), round_digits)
            if point == axial_points[0][0]:
                string_text = f"axial_force at (x={x_loc:.2f},y={y_loc:.2f})  (unwrapped x={unwrapped_x_loc:.2f}+)"
                print(f"{string_text:<50} = {axial_value_plus}")
            elif point == axial_points[-1][0]:
                string_text = f"axial_force at (x={x_loc:.2f},y={y_loc:.2f})  (unwrapped x={unwrapped_x_loc:.2f}-)"
                print(f"{string_text:<50} = {axial_value_minus}")
            else:
                string_text = f"axial_force at (x={x_loc:.2f},y={y_loc:.2f})  (unwrapped x={unwrapped_x_loc:.2f}-)"
                print(f"{string_text:<50} = {axial_value_minus}")
                string_text = f"axial_force at (x={x_loc:.2f},y={y_loc:.2f})  (unwrapped x={unwrapped_x_loc:.2f}+)"
                print(f"{string_text:<50} = {axial_value_plus}")


    #################################################################################
    def draw(
        self, forced_load_size=None, show_load_values=False, draw_support_icons=False
    ):
        """
        Draws the structure, including members, loads, and supports.

        This method creates a visual representation of the structure using matplotlib, displaying
        the members, applied loads, and support types. Loads that are symbolic and remail symbolic after solving are displayed in a grey color.

        Parameters
        ==========
        forced_load_size : float, optional
            This caps the size of all loads to this value, uses the grid units not value units.
            Default is None which scales the loads on the plot based on their value.

        show_load_values : bool, optional
            If True, the load values will be displayed on the plot. Default is False.

        Returns:
            Plot of the structure including members, loads, and supports.

        Examples
        ========
        A structure with two members and multiple loads symbolic loads is created and visualized using the draw method.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
            >>> from sympy.core.symbol import symbols
            >>> E = 3e4
            >>> I = 1
            >>> A = 1e4
            >>> F = symbols("F")
            >>> s = Structure2d()
            >>> s.add_member(x1=0, y1=0, x2=3, y2=4, E=E, I=I, A=A)
            >>> s.add_member(x1=3, y1=4, x2=7, y2=-1, E=E, I=I, A=A)
            >>> s.apply_load(start_x=1.5, start_y=2, value=F, global_angle=0, order=-1)
            >>> s.apply_load(
            ...     start_x=5,
            ...     start_y=1.5,
            ...     value=F / 2,
            ...     global_angle=s.members[1].angle_deg + 270,
            ...     order=0,
            ...     end_x=7,
            ...     end_y=-1,
            ... )
            >>> s.apply_load(
            ...     start_x=0,
            ...     start_y=0,
            ...     value=F * 0.8,
            ...     global_angle=270,
            ...     order=0,
            ...     end_x=3,
            ...     end_y=4,
            ... )
            >>> s.apply_support(x=7, y=-1, type="roller")
            >>> s.apply_support(x=0, y=0, type="pin")
            >>> s.solve_for_reaction_loads()
            {R_h (x=0,y=0): F/4, R_v (x=0,y=0): -341*F/112, R_v (x=7,y=-1): -219*F/112}
            >>> s.draw(show_load_values=True) #doctest: +SKIP

        The same plot can be generated without symbols using nummeric values.

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
            >>> E = 3e4
            >>> I = 1
            >>> A = 1e4
            >>> F = 15
            >>> s = Structure2d()
            >>> s.add_member(x1=0, y1=0, x2=3, y2=4, E=E, I=I, A=A)
            >>> s.add_member(x1=3, y1=4, x2=7, y2=-1, E=E, I=I, A=A)
            >>> s.apply_load(start_x=1.5, start_y=2, value=F, global_angle=0, order=-1)
            >>> s.apply_load(
            ...     start_x=5,
            ...     start_y=1.5,
            ...     value=F / 2,
            ...     global_angle=s.members[1].angle_deg + 270,
            ...     order=0,
            ...     end_x=7,
            ...     end_y=-1,
            ... )
            >>> s.apply_load(
            ...     start_x=0,
            ...     start_y=0,
            ...     value=F * 0.8,
            ...     global_angle=270,
            ...     order=0,
            ...     end_x=3,
            ...     end_y=4,
            ... )
            >>> s.apply_support(x=7, y=-1, type="roller")
            >>> s.apply_support(x=0, y=0, type="pin")
            >>> s.solve_for_reaction_loads()
            {R_h (x=0,y=0): 15/4, R_v (x=0,y=0): -5115/112, R_v (x=7,y=-1): -3285/112}
            >>> s.draw(show_load_values=True) #doctest: +SKIP

        The optional parameter ``draw_support_icons`` will draw the following icons for the supports:

        .. plot::
            :context: close-figs
            :format: doctest
            :include-source: True

            >>> from sympy.physics.continuum_mechanics.structure2d import Structure2d
            >>> E = 3e4
            >>> I = 1
            >>> A = 1e4
            >>> s = Structure2d()
            >>> s.add_member(x1=0, y1=0, x2=7, y2=0, E=E, I=I, A=A)
            >>> s.apply_support(x=0, y=0, type="pin")
            >>> s.apply_support(x=7/2, y=0, type="roller")
            >>> s.apply_support(x=7, y=0, type="fixed")
            >>> s.draw(show_load_values=True, forced_load_size=2, draw_support_icons=True) #doctest: +SKIP
        """

        fig, ax = plt.subplots()
        ax.set_aspect("equal")

        # Define colors
        colors = {
            "point_load": "blue",
            "point_load_symbolic": "lightblue",
            "distributed_load": "#A30000",
            "distributed_load_symbolic": "lightred",
            "support_reaction": "#e6772d",
            "support_reaction_symbolic": "yellow",
            "uncomputed": "#696969",
        }  # symbolic colors not yet implemented uncomputed is used for now
        scale = 0.75

        # Draw members
        for member in self.members:
            x1, y1 = float(member.x1), float(member.y1)
            x2, y2 = float(member.x2), float(member.y2)
            ax.plot([x1, x2], [y1, y2], color="black", linewidth=3, zorder=10)

        # Draw loads
        for load in self.loads:
            x, y = float(load.start_x), float(load.start_y)
            is_support_reaction = getattr(load, "is_support_reaction", False)

            load_color = (
                colors["support_reaction"]
                if is_support_reaction
                else (
                    colors["point_load"]
                    if load.order < 0 or load.order == 1
                    else colors["distributed_load"]
                )
            )

            if isinstance(load.value, Basic) and load.value.has(Symbol):
                load_color = colors["uncomputed"]

            if load.order == -2:
                # Moment load
                radius = 0.5
                angle_start, angle_end, angle_moment_icon, angle2moment_icon = (
                    (180 + 30, 270 - 10, 180 + 30, 180 + 30 - 90)
                    if load.value.is_negative
                    else (180 + 10, 270 - 30, 270 - 30, 270 - 30 + 90)
                )

                arc = patches.Arc(
                    (x, y),
                    width=2.5 * radius,
                    height=2.5 * radius,
                    theta1=angle_start,
                    theta2=angle_end,
                    color=load_color,
                    alpha=1.0,
                    linewidth=2,
                )
                ax.add_patch(arc)

                arrowhead = patches.FancyArrow(
                    x + radius * 1.25 * np.cos(np.radians(angle_moment_icon)),
                    y + radius * 1.25 * np.sin(np.radians(angle_moment_icon)),
                    0.1 * np.cos(np.radians(angle2moment_icon)),
                    0.1 * np.sin(np.radians(angle2moment_icon)),
                    width=0.02,
                    head_width=0.1,
                    head_length=0.2,
                    color=load_color,
                    alpha=1.0,
                )
                ax.add_patch(arrowhead)

                if show_load_values:
                    color_text = (
                        "black" if load_color == colors["uncomputed"] else load_color
                    )
                    value_text = (
                        f"{float(abs(load.value)):0.2f}"
                        if not isinstance(load.value, Basic)
                        or not load.value.has(Symbol)
                        else f'{str((load.value)).split("__")[0]}'
                    )
                    plt.text(
                        x - radius * 2,
                        y - radius,
                        value_text,
                        fontsize=8,
                        color=color_text,
                        zorder=100,
                        ha="center",
                        va="center",
                        bbox={
                            "facecolor": "lightgrey",
                            "alpha": 0.75,
                            "edgecolor": "none",
                        },
                    )

            elif load.order < 0 or load.order == 1:
                # Point load
                load_length = (
                    forced_load_size
                    if forced_load_size is not None
                    else (
                        abs(float(load.value)) / 10.0
                        if isinstance(load.value, (int, float))
                        else 1.0
                    )
                )

                if is_support_reaction:
                    if draw_support_icons:
                        support_offset = scale * 0.65
                    else:
                        support_offset = 0

                    if load.y_component.has(Symbol):
                        y += (
                            support_offset
                            if load.y_component.is_negative
                            else -support_offset
                        )
                    elif load.x_component.has(Symbol):
                        x += (
                            support_offset
                            if load.x_component.is_negative
                            else -support_offset
                        )
                    elif not load.y_component.has(Symbol) and not load.x_component.has(
                        Symbol
                    ):
                        if float(load.y_component) == 0:
                            x += (
                                -support_offset
                                if float(load.x_component) > 0
                                else support_offset
                            )
                        if float(load.x_component) == 0:
                            y += (
                                support_offset
                                if float(load.y_component) > 0
                                else -support_offset
                            )
                if not isinstance(load.value, Basic) or not load.value.has(Symbol):
                    if float(load.value) < 0:
                        angle = float(load.global_angle + pi)
                    else:
                        angle = float(load.global_angle)
                else:
                    angle = float(load.global_angle)
                # angle = float(load.global_angle)
                dx = load_length * np.cos(angle)
                dy = load_length * np.sin(angle)

                arrow_patch = patches.FancyArrow(
                    x - dx,
                    y - dy,
                    dx,
                    dy,
                    length_includes_head=True,
                    width=0.02,
                    head_width=0.1,
                    head_length=0.2,
                    color=load_color,
                    alpha=1.0,
                    zorder=100,
                )
                ax.add_patch(arrow_patch)

                if show_load_values:
                    color_text = (
                        "black" if load_color == colors["uncomputed"] else load_color
                    )
                    value_text = (
                        f"{float(abs(load.value)):0.2f}"
                        if not isinstance(load.value, Basic)
                        or not load.value.has(Symbol)
                        else f'{str((load.value)).split("__")[0]}'
                    )
                    plt.text(
                        x - dx,
                        y - dy,
                        value_text,
                        fontsize=8,
                        color=color_text,
                        zorder=100,
                        ha="center",
                        va="center",
                        bbox={
                            "facecolor": "lightgrey",
                            "alpha": 0.75,
                            "edgecolor": "none",
                        },
                    )

            else:
                # Distributed load
                load_length = (
                    forced_load_size
                    if forced_load_size is not None
                    else (
                        abs(float(load.value)) / 10.0
                        if isinstance(load.value, (int, float))
                        else 1.0
                    )
                )

                x1, y1 = float(load.start_x), float(load.start_y)
                if load.end_x is not None and load.end_y is not None:
                    x2, y2 = float(load.end_x), float(load.end_y)
                else:
                    member = next(
                        (
                            m
                            for m in self.members
                            if (float(m.x1) == x1 and float(m.y1) == y1)
                            or (float(m.x2) == x1 and float(m.y2) == y1)
                        ),
                        None,
                    )
                    if member is None:
                        continue
                    x2, y2 = float(member.x2), float(member.y2)

                member_length = np.hypot(x2 - x1, y2 - y1)
                num_arrows = max(int(member_length / 0.5), 1)

                arrow_tails_x = []
                arrow_tails_y = []

                arrow_heads_x = []
                arrow_heads_y = []

                for i in range(num_arrows + 1):
                    t = i / num_arrows
                    x_point = x1 + t * (x2 - x1)
                    y_point = y1 + t * (y2 - y1)

                    if not isinstance(load.value, Basic) or not load.value.has(Symbol):
                        if float(load.value) < 0:
                            angle = float(load.global_angle + pi)
                        else:
                            angle = float(load.global_angle)

                    dx = load_length * np.cos(angle)
                    dy = load_length * np.sin(angle)

                    arrow_patch = patches.FancyArrow(
                        x_point - dx,
                        y_point - dy,
                        dx,
                        dy,
                        length_includes_head=True,
                        width=0.02,
                        head_width=0.1,
                        head_length=0.2,
                        color=load_color,
                        alpha=1.0,
                        zorder=90,
                    )
                    ax.add_patch(arrow_patch)

                    arrow_tails_x.append(x_point - dx)
                    arrow_tails_y.append(y_point - dy)

                    arrow_heads_x.append(x_point)
                    arrow_heads_y.append(y_point)
                    plt.plot(arrow_tails_x, arrow_tails_y, color=load_color, alpha=1)

                    arrow_coords = [
                        (arrow_tails_x[0], arrow_tails_y[0]),  # First tail point
                        (arrow_tails_x[-1], arrow_tails_y[-1]),  # Last tail point
                        (arrow_heads_x[-1], arrow_heads_y[-1]),  # Last head point
                        (arrow_heads_x[0], arrow_heads_y[0]),  # First head point
                    ]

                polygon = plt.Polygon(
                    arrow_coords,
                    closed=True,
                    fill=True,
                    color=load_color,
                    edgecolor=None,
                    alpha=0.25,
                    zorder=89,
                )
                ax.add_patch(polygon)

                if show_load_values:
                    dx = load_length * np.cos(angle)
                    dy = load_length * np.sin(angle)
                    avrage_x = (x1 - dx + x2 - dx) / 2
                    avrage_y = (y1 - dy + y2 - dy) / 2
                    color_text = (
                        "black" if load_color == colors["uncomputed"] else load_color
                    )
                    value_text = (
                        f"{float(abs(load.value)):0.2f}"
                        if not isinstance(load.value, Basic)
                        or not load.value.has(Symbol)
                        else f'{str(abs(load.value)).split("__")[0]}'
                    )
                    plt.text(
                        avrage_x,
                        avrage_y,
                        value_text,
                        fontsize=8,
                        color=color_text,
                        zorder=100,
                        ha="center",
                        va="center",
                        bbox={
                            "facecolor": "lightgrey",
                            "alpha": 0.75,
                            "edgecolor": "none",
                        },
                    )

        # Draw supports
        if draw_support_icons:
            for support in self.supports:
                x, y = float(support.x), float(support.y)

                if support.node_type == "pin":
                    triangle_height = 0.5 * scale
                    angle = 34 * (pi / 180)
                    half_width = triangle_height * float(tan(angle))

                    triangle_vertices = [
                        [x, y],
                        [x - half_width, y - triangle_height],
                        [x + half_width, y - triangle_height],
                    ]

                    triangle = patches.Polygon(
                        triangle_vertices,
                        edgecolor="k",
                        facecolor="none",
                        linewidth=1.5,
                        zorder=20,
                    )
                    ax.add_patch(triangle)

                    pivot_circle = patches.Circle(
                        (x, y),
                        0.1 * scale,
                        edgecolor="k",
                        facecolor="white",
                        linewidth=1.5,
                        zorder=20,
                    )
                    ax.add_patch(pivot_circle)

                elif support.node_type == "roller":
                    triangle_height = 0.5 * scale
                    angle = 34 * (pi / 180)
                    half_width = triangle_height * float(tan(angle))

                    triangle_vertices = [
                        [x, y],
                        [x - half_width, y - triangle_height],
                        [x + half_width, y - triangle_height],
                    ]

                    triangle = patches.Polygon(
                        triangle_vertices,
                        edgecolor="k",
                        facecolor="none",
                        linewidth=1.5,
                        zorder=20,
                    )
                    ax.add_patch(triangle)

                    pivot_circle = patches.Circle(
                        (x, y),
                        0.1 * scale,
                        edgecolor="k",
                        facecolor="white",
                        linewidth=1.5,
                        zorder=20,
                    )
                    ax.add_patch(pivot_circle)

                    plt.plot(
                        [x - half_width, x + half_width],
                        [
                            y - triangle_height - (scale * 0.1),
                            y - triangle_height - (scale * 0.1),
                        ],
                        color="k",
                        zorder=20,
                    )

                elif support.node_type == "fixed":
                    box_size = scale * 0.5
                    x_box = x - box_size / 2
                    y_box = y - box_size / 2

                    scaled_vertices = [
                        [x_box, y_box],
                        [(x_box + box_size), y_box],
                        [(x_box + box_size), (y_box + box_size)],
                        [x_box, (y_box + box_size)],
                    ]

                    polygon = patches.Polygon(
                        scaled_vertices, edgecolor="k", facecolor="none", linewidth=1.5
                    )
                    ax.add_patch(polygon)

                    num_lines = 4
                    for i in range(num_lines):
                        x_start = x_box + i * (box_size / num_lines)
                        y_start = y_box
                        x_end = x_box + (i + 1) * (box_size / num_lines)
                        y_end = y_box + (box_size / num_lines)
                        ax.plot(
                            [x_start, x_end],
                            [y_start, y_end],
                            color="black",
                            linewidth=1,
                        )

        x_ticks = plt.xticks()[0]
        y_ticks = plt.yticks()[0]

        if len(x_ticks) > 1 and len(y_ticks) > 1:
            x_step = x_ticks[1] - x_ticks[0]
            y_step = y_ticks[1] - y_ticks[0]

            step = min(x_step, y_step)

            if step < 1:
                step = 1

            x_min, x_max = plt.xlim()
            y_min, y_max = plt.ylim()

            # Generate new ticks using the updated step
            new_x_ticks = [
                i * step for i in range(int(x_min // step), int(x_max // step) + 1)
            ]
            new_y_ticks = [
                i * step for i in range(int(y_min // step), int(y_max // step) + 1)
            ]

            plt.xticks(new_x_ticks)
            plt.yticks(new_y_ticks)

        ax.grid(True, zorder=10)

"""
This module can be used to solve 2D problems with
singularity functions in mechanics.
"""

from sympy.core import Basic, Symbol
from sympy.core.numbers import pi
from sympy.external import import_module
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.trigonometric import atan2, cos, sin, tan
from sympy.geometry.polygon import deg, rad
from sympy.physics.continuum_mechanics.beam import Beam
from sympy.physics.continuum_mechanics.column import Column
from sympy.simplify import nsimplify, simplify

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
        A consistent sign convention should be maintained for all forces and reactions.
        The positive direction for the x-axis is to the right, and the positive direction for the y-axis is downwards.

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
        >>> Rv1 = s.apply_support(x=7, y=-1, type="roller")
        >>> Rv2, Rh2 = s.apply_support(x=0, y=0, type="pin")
        >>> s.solve_for_reaction_loads(Rv1, Rv2, Rh2)
        {R_h__0,__0: 3.75, R_v__0,__0: -45.6696428571429, R_v__7,__-1: -29.3303571428571}
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
        >>> Rv1 = s.apply_support(x=15, y=2, type="roller")
        >>> Rv2, Rh2 = s.apply_support(x=0, y=0, type="pin")
        >>> s.solve_for_reaction_loads(Rv1, Rv2, Rh2)
        {R_h__0,__0: -38.4615384615385, R_v__0,__0: -133.891025641026, R_v__15,__2: -163.416666666667}
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
        self.nodes = []
        self.loads = []
        self.unwrapped_bendpoints = []
        self.unwrapped_loadpoints = []
        self.beam = self._init_beam()
        self.column = self._init_column()
        self.reaction_loads = {}
        self.load_qz = 0

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
        variable=Symbol("x"),
        base_char="C",
        ):
        return Column(
            length=length,
            elastic_modulus=elastic_modulus,
            area=area,
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

        aa, oo, L = self._unwrap_structure()
        # Convert to alexes input
        # T = 1, Fv = 2, Fh = 3, qv = 4, qh = 5
        if load.order == -2:
            qz = load.value
            B = [nsimplify(qz)]
            bb = [nsimplify(self.unwrapped_loadpoints[-1]["locals"][0])]
            nn = [1]
        elif end_x is None and end_y is None:
            Fv = load.y_component
            Fh = load.x_component
            B = [nsimplify(Fv), nsimplify(Fh)]
            bb = [
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][0]),
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][0]),
            ]
            nn = [2, 3]
        else:
            qv = load.y_component
            qh = load.x_component
            B = [nsimplify(qv), nsimplify(-qv), nsimplify(qh), nsimplify(-qh)]
            bb = [
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][0]),
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][1]),
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][0]),
                nsimplify(self.unwrapped_loadpoints[-1]["locals"][1]),
            ]
            nn = [4, 4, 5, 5]

        # ALEX algo ############################################################################################
        # bendpoints
        for i in range(len(B)):
            for j in range(len(aa)):
                if bb[i] == aa[-1]:
                    if nn[i] == 1:
                        self.beam.apply_load(B[i], bb[i], -2)

                    if nn[i] == 2:
                        self.beam.apply_load(B[i] * cos(oo[-1]), bb[i], -1)

                    if nn[i] == 3:
                        self.column.apply_load(B[i] * cos(oo[-1]), bb[i], -1)

                    if nn[i] == 4:
                        self.beam.apply_load(B[i] * cos(oo[-1]), bb[i], 0)

                    if nn[i] == 5:
                        self.column.apply_load(B[i] * cos(oo[-1]), bb[i], 0)

                    break
                else:
                    if bb[i] < aa[j]:
                        if nn[i] == 1:
                            self.beam.apply_load(B[i], bb[i], -2)

                        if nn[i] == 2:
                            self.beam.apply_load(B[i] * cos(oo[j - 1]), bb[i], -1)

                        if nn[i] == 3:
                            self.column.apply_load(B[i] * cos(oo[j - 1]), bb[i], -1)

                        if nn[i] == 4:
                            self.beam.apply_load(B[i] * cos(oo[j - 1]), bb[i], 0)

                        if nn[i] == 5:
                            self.column.apply_load(B[i] * cos(oo[j - 1]), bb[i], 0)
                        break

        for i in range(len(B)):
            for j in range(len(aa) - 1):
                if bb[i] < aa[j]:
                    if nn[i] == 2:
                        self.beam.apply_load(
                            B[i] * (cos(oo[j]) - cos(oo[j - 1])), aa[j], -1
                        )

                    if nn[i] == 3:
                        self.column.apply_load(
                            B[i] * (cos(oo[j]) - cos(oo[j - 1])), aa[j], -1
                        )

                    if nn[i] == 4:
                        self.beam.apply_load(
                            B[i] * (cos(oo[j]) - cos(oo[j - 1])), aa[j], 0
                        )
                        self.beam.apply_load(
                            B[i] * (aa[j] - bb[i]) * (cos(oo[j]) - cos(oo[j - 1])),
                            aa[j],
                            -1,
                        )

                    if nn[i] == 5:
                        self.column.apply_load(
                            B[i] * (cos(oo[j]) - cos(oo[j - 1])), aa[j], 0
                        )
                        self.column.apply_load(
                            B[i] * (aa[j] - bb[i]) * (cos(oo[j]) - cos(oo[j - 1])),
                            aa[j],
                            -1,
                        )

        self.load_qz = self.beam._load

    ###################################################################################################################

    def _unwrap_structure(self):
        """Unwraps the structure to a 1D beam for analysis."""
        unwrapped_len = 0
        unwrapped_bendpoints = []
        unwrapped_loadpoints = []

        # Process members
        for member in self.members:
            unwrapped_bendpoints.append(
                {
                    "bend_point": [unwrapped_len, unwrapped_len + member.length],
                    "angle": member.angle,
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
                                nsimplify(unwrapped_len + load.local_end),
                            ],
                        }
                    )
                else:
                    unwrapped_loadpoints.append(
                        {
                            "l_id": load.load_id,
                            "locals": [nsimplify(unwrapped_len + load.local_start)],
                        }
                    )

            # Process loads at the start node
            x_start, y_start = member.x1, member.y1
            for node in self.nodes:
                if simplify(node.x - x_start) == 0 and simplify(node.y - y_start) == 0:
                    for load in node.node_loads:
                        unwrapped_loadpoints.append(
                            {"l_id": load.load_id, "locals": [nsimplify(unwrapped_len)]}
                        )
                    break

            unwrapped_len += member.length

            # Process loads at the end node
            x_end, y_end = member.x2, member.y2
            for node in self.nodes:
                if simplify(node.x - x_end) == 0 and simplify(node.y - y_end) == 0:
                    for load in node.node_loads:
                        unwrapped_loadpoints.append(
                            {"l_id": load.load_id, "locals": [nsimplify(unwrapped_len)]}
                        )
                    break

        self.unwrapped_bendpoints = unwrapped_bendpoints
        self.unwrapped_loadpoints = sorted(
            unwrapped_loadpoints, key=lambda x: x["l_id"]
        )
        ### Emulate Alexes input for his algorithm
        aa = []
        oo = []

        for unwrapped_bendpoint in unwrapped_bendpoints:
            aa.append(unwrapped_bendpoint["bend_point"][0])
            oo.append(unwrapped_bendpoint["angle"])
        aa.append(unwrapped_bendpoints[-1]["bend_point"][1])

        L = unwrapped_len
        self.beam.length = L

        return aa, oo, L
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

        Returns:
            SymPy Symbol(s): The reaction loads at the support location. The return value depends on the type of support applied.

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
        >>> Rv1, Rh1 = s.apply_support(x=0, y=0, type="pin")
        >>> Rv2 = s.apply_support(x=4, y=0, type="roller")
        """

        support = self._add_or_update_node(x, y, type)
        self.supports.append(support)

        unwarap_x = self._find_unwrapped_position(x, y)

        Rh = Symbol(f"R_h__{round(x,2)},__{round(y,2)}")
        Rv = Symbol(f"R_v__{round(x,2)},__{round(y,2)}")
        T = Symbol(f"T__{round(x,2)},__{round(y,2)}")

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
            self.apply_load(
                start_x=x, start_y=y, value=load_v, global_angle=90, order=-1
            )
            self.apply_load(
                start_x=x, start_y=y, value=load_h, global_angle=0, order=-1
            )

            # Mark these loads as support reactions
            self.loads[-1].is_support_reaction = True
            self.loads[-2].is_support_reaction = True

            self.beam.bc_deflection.append((unwarap_x, 0))
            self.column._bc_deflection.append(unwarap_x)
            return Rv, Rh

        elif type == "roller":
            load_v = -1 * Rv
            self.apply_load(
                start_x=x, start_y=y, value=load_v, global_angle=90, order=-1
            )

            # Mark this load as a support reaction
            self.loads[-1].is_support_reaction = True

            self.beam.bc_deflection.append((unwarap_x, 0))
            return Rv

        elif type == "fixed":
            load_t = T
            load_v = -1 * Rv
            load_h = Rh
            self.apply_load(
                start_x=x, start_y=y, value=load_t, global_angle=0, order=-2
            )
            self.apply_load(
                start_x=x, start_y=y, value=load_v, global_angle=90, order=-1
            )
            self.apply_load(
                start_x=x, start_y=y, value=load_h, global_angle=0, order=-1
            )

            # Mark these loads as support reactions
            self.loads[-1].is_support_reaction = True
            self.loads[-2].is_support_reaction = True
            self.loads[-3].is_support_reaction = True

            self.beam.bc_deflection.append((unwarap_x, 0))
            self.column._bc_deflection.append(unwarap_x)
            # This i think sould be slope of the beam at this point beacause 0 is assuming supports and horizontal members
            # unwrap is already called so it should be extaracatble from the unwrapped position
            self.beam.bc_slope.append((unwarap_x, 0))
            return T, Rv, Rh

    def solve_for_reaction_loads(self, *args):
        """
        Solves for the reaction loads in the structure based on the applied loads and applied support reactions.

        Parameters
        ==========
        *args : dict
            Variables representing the reaction loads.

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
        >>> Rv1 = s.apply_support(x=15, y=2, type="roller")
        >>> Rv2, Rh2 = s.apply_support(x=0, y=0, type="pin")
        >>> s.solve_for_reaction_loads(Rv1, Rv2, Rh2)
        {R_h__0,__0: -3.84615384615385, R_v__0,__0: -70.3141025641026, R_v__15,__2: -143.916666666667}
        """


        # Split arguments into vertical and horizontal reaction loads
        reaction_loads_vertical = [arg for arg in args if "R_v" in str(arg)]
        reaction_loads_horizontal = [arg for arg in args if "R_h" in str(arg)]
        reaction_moments = [arg for arg in args if "T_" in str(arg)]

        args_for_beam_solver = tuple(
            reaction_loads_vertical + reaction_moments + reaction_loads_horizontal
        )

        # display(self.beam.load)
        # display(self.beam.shear_force(),'--')
        # Solve for moment and vertical reaction loads using the beam solver
        # print(args_for_beam_solver)
        self.beam.solve_for_reaction_loads(*args_for_beam_solver)

        # Compute the horizontal reaction load by summing up all horizontal forces
        sum_horizontal = 0
        for load in self.loads:
            if isinstance(load.value, (int, float)):
                if load.order == -1:
                    sum_horizontal += load.x_component
                elif load.order == 0:
                    length = sqrt(
                        (load.end_y - load.start_y) ** 2
                        + (load.end_x - load.start_x) ** 2
                    )
                    sum_horizontal += length * load.x_component

        # Substitute horizontal reactions with their solved values
        horizontal_key = {
            arg: float(-1 * sum_horizontal) for arg in args if "R_h" in str(arg)
        }

        # Update solved reaction loads dictionary with horizontal reaction values
        for key in self.beam._reaction_loads.items():
            key = key[0]
            self.beam._reaction_loads[key] = self.beam._reaction_loads[key].subs(
                horizontal_key
            )

        # Store reaction loads in the structure's state
        self.reaction_loads = self.beam.reaction_loads

        # Substitute solved reactions into the beam's load equation
        self.beam._load = self.beam._load.subs(self.beam.reaction_loads)

        # Check for symbolic or numerical reaction load solution
        list_of_symbols = []
        for key in self.reaction_loads.items():
            list_of_symbols.append(str(key[0]))

        list_of_symbols_reactions = []
        for load in self.loads:
            if isinstance(load.value, Basic):
                list_of_symbols_reactions.append(str(load.value.as_coeff_Mul()[1]))

        # If all symbols are resolved
        to_ignore_list = set(list_of_symbols).union(list_of_symbols_reactions)
        if len(to_ignore_list) == len(list_of_symbols):
            # print('Debug - Reaction loads are numerical')
            vertical_key = {
                arg: self.beam._reaction_loads[arg] for arg in args if "R_v" in str(arg)
            }
            bending_moment_key = {
                arg: self.beam._reaction_loads[arg] for arg in args if "T_" in str(arg)
            }

            for load in self.loads:
                if isinstance(load.value, Basic) and load.order != -2:
                    # Substitute horizontal and vertical reaction loads
                    load.value = load.value.subs(vertical_key).subs(horizontal_key)
                    load.y_component = load.y_component.subs(vertical_key).subs(
                        horizontal_key
                    )
                    load.x_component = load.x_component.subs(vertical_key).subs(
                        horizontal_key
                    )
                    if isinstance(load.value, Basic):
                        pass
                    else:
                        # Adjust load direction
                        if load.y_component > 0:
                            load.global_angle = load.global_angle + pi
                        if load.x_component < 0:
                            load.global_angle = load.global_angle + pi

                elif isinstance(load.value, Basic) and load.order == -2:
                    load.value = load.value.subs(bending_moment_key)

        else:
            # print('Debug - Reaction loads are symbolic')
            pass
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
            return self.beam.shear_force()
        if y is None:
            x_symbol = Symbol("x")
            return self.beam.shear_force().subs(x_symbol, x)
        else:
            unwarap_x = self._find_unwrapped_position(x, y)
            return self.beam.shear_force().subs(Symbol("x"), unwarap_x)

    def plot_shear_force(self):
        """
        Plots the shear force diagram for the beam.

        Returns:
            Matplotlib plot: A plot showing the shear force distribution along the structure.
        """

        self.beam.plot_shear_force()

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
            return self.beam.bending_moment()
        if y is None:
            x_symbol = Symbol("x")
            return self.beam.bending_moment().subs(x_symbol, x)
        else:
            unwarap_x = self._find_unwrapped_position(x, y)
            return self.beam.bending_moment().subs(Symbol("x"), unwarap_x)

    def plot_bending_moment(self):
        """
        Plots the bending moment diagram for the beam.

        Returns:
            Matplotlib plot: A plot showing the bending moment distribution along the structure.
        """

        self.beam.plot_bending_moment()

    def summary(self, verbose=True, round_digits=None):
        """
        Provides a summary of the structure, including reaction loads and points of interest.

        This method prints a formatted summary of the structure's members, nodes, and supports,
        along with the reaction loads and points of interest related to bending moments and shear forces.
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
        >>> Rv1, Rh1 = s.apply_support(x=0, y=0, type="pin")
        >>> Rv2 = s.apply_support(x=4, y=0, type="roller")
        >>> s.solve_for_reaction_loads(Rh1, Rv1, Rv2)
        {R_h__0,__0: 0.0, R_v__0,__0: -5, R_v__4,__0: -5}
        >>> s.summary(round_digits=2)
        ===================== Structure Summary =====================
        <BLANKLINE>
        Reaction Loads:
        R_v   [0.00,0.00]  (0.00)                = -5.0
        R_v   [4.00,0.00]  (4.00)                = -5.0
        R_h   [0.00,0.00]  (0.00)                = 0.0
        <BLANKLINE>
        Points of Interest - Bending Moment:
        bending_moment at [x.xx,y.yy]  (0.00)    = 0.0
        bending_moment at [x.xx,y.yy]  (4.00)-   = 0.0
        <BLANKLINE>
        Points of Interest - Shear Force:
        shear_force at [x.xx,y.yy]  (0.00+)      = 5.0
        shear_force at [x.xx,y.yy]  (2.00-)      = 5.0
        shear_force at [x.xx,y.yy]  (2.00+)      = -5.0
        shear_force at [x.xx,y.yy]  (4.00-)      = -5.0
        """

        title = "Structure Summary"
        line_length = 60

        print(
            f'{"=" * ((line_length - len(title)) // 2)} {title} {"=" * ((line_length - len(title)) // 2)}'
        )

        if not self.reaction_loads:
            print("\nPlease solve for reaction loads first")
            return

        if verbose:
            self._print_reaction_loads(round_digits)
            self._print_points_of_interest(round_digits)

    def _print_reaction_loads(self, round_digits):
        """Prints the reaction loads in the structure."""
        print("\nReaction Loads:")
        for key, value in self.reaction_loads.items():
            support_name = str(key).split("__")[0]
            support_x_loc = float(str(key).split("__")[1].strip(","))
            support_y_loc = float(str(key).split("__")[2].strip(","))
            unwrapped_xl_loc = self._find_unwrapped_position(
                support_x_loc, support_y_loc
            )

            if round_digits is not None and not value.has(Symbol):
                value = round(float(value), round_digits)

            support_str = f"{support_name:<5} [{support_x_loc:.2f},{support_y_loc:.2f}]  ({unwrapped_xl_loc:.2f})"
            print(f"{support_str:<40} = {(value)}")

    def _print_points_of_interest(self, round_digits):
        """Prints the points of interest for shear force and bending moment."""
        dx = 1e-6
        # dx = 0

        print("\nPoints of Interest - Bending Moment:")
        bend_points = sorted(
            {
                nsimplify(point)
                for item in self.unwrapped_bendpoints
                for point in item["bend_point"]
            }
        )
        for point in bend_points:
            if point == bend_points[-1]:
                bending_moment_value = self.bending_moment(point - dx)
                if round_digits is not None and not bending_moment_value.has(Symbol):
                    bending_moment_value = round(
                        float(bending_moment_value), round_digits
                    )
                string_text = f"bending_moment at [x.xx,y.yy]  ({point:.02f})-"
                print(f"{string_text:<40} = {bending_moment_value}")

            else:
                bending_moment_value = self.bending_moment(nsimplify(point))

                if round_digits is not None and not bending_moment_value.has(Symbol):
                    bending_moment_value = round(
                        float(bending_moment_value), round_digits
                    )
                string_text = f"bending_moment at [x.xx,y.yy]  ({point:.02f})"
                print(f"{string_text:<40} = {bending_moment_value}")

        print("\nPoints of Interest - Shear Force:")
        load_points = sorted(
            {
                nsimplify(local)
                for point in self.unwrapped_loadpoints
                for local in point["locals"]
            }
        )
        load_points = []
        for item in self.unwrapped_loadpoints:
            # print(item["locals"])
            if len(item["locals"]) == 1:
                load_points.append(item["locals"][0])
        load_points = sorted(set(load_points + bend_points))
        # print(load_points)

        for point in sorted(load_points):
            shear_force_value_minus = self.shear_force(point - dx)
            shear_force_value_plus = self.shear_force(nsimplify(point))

            if round_digits is not None:
                if isinstance(
                    shear_force_value_minus, Basic
                ) and not shear_force_value_minus.has(Symbol):
                    shear_force_value_minus = round(
                        float(shear_force_value_minus), round_digits
                    )
                if isinstance(
                    shear_force_value_plus, Basic
                ) and not shear_force_value_plus.has(Symbol):
                    shear_force_value_plus = round(
                        float(shear_force_value_plus), round_digits
                    )

            if point == load_points[0]:
                string_text = f"shear_force at [x.xx,y.yy]  ({point:.02f}+)"
                print(f"{string_text:<40} = {shear_force_value_plus}")
            elif point == load_points[-1]:
                string_text = f"shear_force at [x.xx,y.yy]  ({point:.02f}-)"
                print(f"{string_text:<40} = {shear_force_value_minus}")
            else:
                string_text = f"shear_force at [x.xx,y.yy]  ({point:.02f}-)"
                print(f"{string_text:<40} = {shear_force_value_minus}")
                string_text = f"shear_force at [x.xx,y.yy]  ({point:.02f}+)"
                print(f"{string_text:<40} = {shear_force_value_plus}")

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
            >>> Rv1 = s.apply_support(x=7, y=-1, type="roller")
            >>> Rv2, Rh2 = s.apply_support(x=0, y=0, type="pin")
            >>> s.solve_for_reaction_loads(Rv1, Rv2, Rh2)
            {R_h__0,__0: 0.0, R_v__0,__0: -345*F/112, R_v__7,__-1: -125*F/56}
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
            >>> Rv1 = s.apply_support(x=7, y=-1, type="roller")
            >>> Rv2, Rh2 = s.apply_support(x=0, y=0, type="pin")
            >>> s.solve_for_reaction_loads(Rv1, Rv2, Rh2)
            {R_h__0,__0: 3.75, R_v__0,__0: -45.6696428571429, R_v__7,__-1: -29.3303571428571}
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
            >>> Rv1, Rh1 = s.apply_support(x=0, y=0, type="pin")
            >>> Rv2 = s.apply_support(x=7/2, y=0, type="roller")
            >>> Rv3, Rh3, T1 = s.apply_support(x=7, y=0, type="fixed")
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

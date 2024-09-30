from sympy.physics.continuum_mechanics.beam import Beam
from sympy.core import symbols,Symbol

from sympy.core.numbers import pi
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.trigonometric import sin, cos, atan2
from sympy.simplify import nsimplify
from sympy.simplify.simplify import simplify
from sympy.geometry.polygon import deg, rad
from sympy.core import Expr
from sympy.external import import_module

plt = import_module(
    "matplotlib.pyplot",
    import_kwargs={
        "fromlist": [
            "pyplot",
        ]
    },
)

class Member:
    def __init__(self,x1, y1, x2, y2, E, I, A, member_id):

        self.x1 = simplify(x1)
        self.y1 = simplify(y1)
        self.x2 = simplify(x2)
        self.y2 = simplify(y2)
        self.E = simplify(E)
        self.I = simplify(I)
        self.A = simplify(A)
        self.member_id = member_id

        # properties that are auto computed
        self.length = self._compute_length()
        self.angle = self._compute_angle()
        self.angle_deg = deg(self.angle)
        self.member_loads = []
        self.member_eq = self._compute_eq()

    def _compute_eq(self):
        x1, y1, x2, y2 = self.x1, self.y1, self.x2, self.y2
        x_mem = Symbol('x')

        a = (y2 - y1) / (x2 - x1)
        b = y1 - a * x1
        eq = a * x_mem + b
        return eq

        # loads
        # self.loads = []

    def _compute_length(self):
        """Compute the length of the member."""
        return sqrt((self.x2 - self.x1) ** 2 + (self.y2 - self.y1) ** 2)

    def _compute_angle(self):
        """Compute the angle of the member in radians."""
        return atan2(self.y2 - self.y1, self.x2 - self.x1)

    def __repr__(self):
        return f"Member(ID={self.member_id}, Length={self.length}, Global_Angle={float(self.angle_deg):.02f}°)"

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
    def __init__(self, start_x, start_y, value, global_angle, order,end_x=None,end_y=None,load_id=None):
        self.start_x = start_x
        self.start_y = start_y
        self.value = value
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

        # self.relative_loc = None

    def __repr__(self):
        return f"Load(applied_to={self.applied_to},X1={self.start_x}, Y1={self.start_y},X2={self.end_x},Y2={self.end_y} Load={self.value}, Global_Angle={self.global_angle}, Order={self.order}, X-Component={self.x_component}, Y-Component={self.y_component})"

    def _compute_x_component(self):
        return self.value * cos(self.global_angle) * -1

    def _compute_y_component(self):
        return self.value * sin(self.global_angle)

    def _compute_local_loc(self):
        pass



class Structure2d:
    def __init__(self):
        self.members = []
        self.supports = []
        self.nodes = []
        self.loads = []
        # self.unwrapped_len = 0
        self.unwrapped_bendpoints = []
        self.unwrapped_loadpoints = []
        self.beam = self._init_beam()

    def __repr__(self):
        return f"Structure2d(Members={len(self.members)}, Nodes={len(self.nodes)}, Supports={len(self.supports)})"

    def _init_beam(self, length=1, elastic_modulus=1, second_moment=1, area=Symbol('A'), variable=Symbol('x'), base_char='C', ild_variable=Symbol('a')):
        beam = Beam(length=length, elastic_modulus=elastic_modulus, second_moment=second_moment)
        return beam

    def add_member(self, x1, y1, x2, y2, E, I, A):
        member_id = len(self.members)
        member = Member(x1, y1, x2, y2, E, I, A, member_id)
        self.members.append(member)
        self.add_or_update_node(x1, y1, "fixed", overwrite_type=False)
        self.add_or_update_node(x2, y2, "fixed", overwrite_type=False)
        return member

    def add_or_update_node(self, x, y, new_node_type, overwrite_type=True):
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


    def _find_applied_to(self, x, y,x_end=None,y_end=None):

        if x_end is not None and y_end is not None:
            x_start = x
            y_start = y
            x = (x + x_end) / 2
            y = (y + y_end) / 2

        for node in self.nodes:
            if node.x == x and node.y == y:
                return f'n_{node.node_id}',0,None


        for member in self.members:
            member_y = member.member_eq.subs({'x':x})
            if member.x1 <= x <= member.x2 or member.x2 <= x <= member.x1:
                if simplify(member_y - y)== 0:
                    local_x_start = sqrt((member.x1 - x)**2 + (member.y1 - y)**2)
                    if x_end is not None and y_end is not None:
                        local_x_start = sqrt((member.x1 - x_start)**2 + (member.y1 - y_start)**2)
                        local_x_end = member.length - sqrt((member.x2 - x_end)**2 + (member.y2 - y_end)**2)
                    else:
                        local_x_end = None
                    return f'm_{member.member_id}',local_x_start,local_x_end



    def apply_load(self, start_x, start_y, value, global_angle, order, end_x=None, end_y=None):
        load_id = len(self.loads)
        load = Load(start_x, start_y, value, global_angle, order, end_x, end_y, load_id)
        self.loads.append(load)

        if end_x is not None and end_y is not None:
            load.applied_to,load.local_start,load.local_end = self._find_applied_to(start_x, start_y,end_x,end_y)

        else:
            load.applied_to,load.local_start,load.local_end = self._find_applied_to(start_x, start_y)

        if 'n' in load.applied_to:
            self.nodes[int(load.applied_to.split('_')[1])].node_loads.append(load)
        elif 'm' in load.applied_to:
            self.members[int(load.applied_to.split('_')[1])].member_loads.append(load)

        # self._unwrap_structure()
        ########################

    def _unwrap_structure(self):
        unwrapped_len = 0
        unwrapped_bendpoints = []
        unwrapped_loadpoints = []

        for member in self.members:
            unwrapped_bendpoints.append({'bend_point':unwrapped_len,'angle':member.angle})

            for load in member.member_loads:
                if load.end_x is not None and load.end_y is not None:
                    unwrapped_loadpoints.append({'l_id':load.load_id,'locals':[unwrapped_len+load.local_start,unwrapped_len+load.local_end]})
                else:
                    unwrapped_loadpoints.append({'l_id':load.load_id,'locals':[unwrapped_len+load.local_start]})

            x = member.x1
            y = member.y1
            for node in self.nodes:
                if node.x == x and node.y == y:
                    for load in node.node_loads:
                        unwrapped_loadpoints.append({'l_id':load.load_id,'locals':[unwrapped_len]})
                    break
            unwrapped_len += member.length

        ### UGLY CODE
        #manualy do the last node
        unwrapped_bendpoints.append({'bend_point':unwrapped_len,'angle':0})

        x = member.x1
        y = member.y1
        for node in self.nodes:
            if node.x == x and node.y == y:
                for load in node.node_loads:
                    unwrapped_loadpoints.append({'l_id':load.load_id,'locals':[unwrapped_len]})
                break

        self.unwrapped_bendpoints = unwrapped_bendpoints
        #sort based on load id
        self.unwrapped_loadpoints = sorted(unwrapped_loadpoints,key=lambda x: x['l_id'])
        # self.unwrapped_loadpoints = unwrapped_loadpoints



    def apply_support(self, x, y, type="pin"):
        support = self.add_or_update_node(x, y, type)
        self.supports.append(support)
        pass







#################################################################################
    def draw(self):
        fig, ax = plt.subplots()
        # Set aspect ratio to 'equal' so that the members are plotted proportionally
        ax.set_aspect('equal')

        # Draw each member
        for member in self.members:
            x1, y1, x2, y2 = float(member.x1), float(member.y1), float(member.x2), float(member.y2)
            ax.plot([x1, x2], [y1, y2], color='black')

        # Draw each node
        for node in self.nodes:
            x, y = float(node.x), float(node.y)
            ax.plot(x, y, 'o', color='k')

        # Draw each load (point and distributed loads)
        for load in self.loads:
            if load.order < 0:
                # Point load
                x, y = float(load.start_x), float(load.start_y)
                load_value = float(load.value)
                angle = float(load.global_angle)

                load_value = load_value / 10
                # Compute the start point for the arrow (away from the load's application point)
                arrow_x = x + load_value * cos(angle)
                arrow_y = y + load_value * sin(angle)

                # Add arrow annotations for the load
                ax.annotate(
                    "",
                    xy=(x, y),
                    xytext=(arrow_x, arrow_y),
                    arrowprops=dict(arrowstyle='->', color='red',shrinkA=0, shrinkB=0)  # noqa: C408
                )
            else:
                # Distributed load
                x1, y1 = float(load.start_x), float(load.start_y)
                if load.end_x is not None and load.end_y is not None:
                    x2, y2 = float(load.end_x), float(load.end_y)
                else:
                    member = next(m for m in self.members if m.x1 == load.start_x and m.y1 == load.start_y)
                    x2, y2 = float(member.x2), float(member.y2)

                member_length = sqrt((x2 - x1)**2 + (y2 - y1)**2)
                arrow_spacing = 0.5  # Define a consistent arrow spacing
                num_arrows = int(member_length / arrow_spacing)

                dx = (x2 - x1) / num_arrows
                dy = (y2 - y1) / num_arrows

                # Plot the load line that the arrows will originate from

                arrow_coords = []
                for i in range(num_arrows + 1):
                    x = x1 + i * dx
                    y = y1 + i * dy

                    arrow_length = load.value * (x) ** load.order
                    # print(arrow_length)

                    arrow_length = arrow_length / 10  # Scale the arrow length for better visualization
                    # Arrows originate from the load line, pointing in the direction of the load
                    arrow_x = x + arrow_length * cos(float(load.global_angle))
                    arrow_y = y + arrow_length * sin(float(load.global_angle))

                    arrow_coords.append((arrow_x, arrow_y))


                    ax.annotate(
                        "",
                        xy=(x, y),
                        xytext=(arrow_x, arrow_y),
                        arrowprops=dict(arrowstyle='->', color='green',shrinkA=0, shrinkB=0)  # noqa: C408
                    )
                for i in range(len(arrow_coords)-1):
                    ax.plot([arrow_coords[i][0], arrow_coords[i+1][0]], [arrow_coords[i][1], arrow_coords[i+1][1]], color='green')

        # Draw each support
        for support in self.supports:
            x, y = float(support.x), float(support.y)
            ax.plot(x, y, 's', color='blue')

    def draw_unwrapped(self):
        fig, ax = plt.subplots()
        # Set aspect ratio to 'equal' for correct proportions
        ax.set_aspect('equal', adjustable='datalim')

        x_start = 0

        # Draw each member in the unwrapped configuration
        for member in self.members:
            plt.plot([x_start, x_start + member.length], [0, 0], 'o-', color='black')
            x_start += member.length


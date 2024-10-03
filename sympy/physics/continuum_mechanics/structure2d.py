from sympy.physics.continuum_mechanics.beam import Beam
from sympy.core import Symbol

# from sympy.core.numbers import pi
from sympy.functions.elementary.miscellaneous import sqrt
# from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.trigonometric import sin, cos, atan2
# from sympy.simplify import nsimplify
from sympy.simplify.simplify import simplify
from sympy.geometry.polygon import deg, rad
# from sympy.core import Expr
from sympy.external import import_module

from sympy.functions import SingularityFunction

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
        if x2 == x1:
            return 'vertical'
        else:
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

        self.load_qz = 0

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


    def _find_applied_to(self, x, y, x_end=None, y_end=None):
        if x_end is not None and y_end is not None:
            x_start = x
            y_start = y
            x = (x + x_end) / 2
            y = (y + y_end) / 2


        for node in self.nodes:
            if node.x == x and node.y == y:
                return f'n_{node.node_id}', 0, None

        # Check if point corresponds to a member
        for member in self.members:
            if member.member_eq != 'vertical':
                # Handle non-vertical members
                member_y = member.member_eq.subs({'x': x})
                if member.x1 <= x <= member.x2 or member.x2 <= x <= member.x1:
                    if simplify(member_y - y) == 0:
                        local_x_start = sqrt((member.x1 - x) ** 2 + (member.y1 - y) ** 2)
                        if x_end is not None and y_end is not None:
                            local_x_start = sqrt((member.x1 - x_start) ** 2 + (member.y1 - y_start) ** 2)
                            local_x_end = member.length - sqrt((member.x2 - x_end) ** 2 + (member.y2 - y_end) ** 2)
                        else:
                            local_x_end = None
                        return f'm_{member.member_id}', local_x_start, local_x_end
            else:
                # Handle vertical members
                if simplify(member.x1 - x) ==0:
                    if member.y1 <= y <= member.y2 or member.y2 <= y <= member.y1:
                        local_x_start = abs(member.y1 - y)
                        if x_end is not None and y_end is not None:
                            local_x_start = abs(member.y1 - y_start)
                            local_x_end = member.length - abs(member.y2 - y_end)
                        else:
                            local_x_end = None
                        return f'm_{member.member_id}', local_x_start, local_x_end






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

        # build load equation
        ########################

        aa,oo,L = self._unwrap_structure()

        if end_x is None and end_y is None:
            Fv = load.y_component
            Fh = load.x_component
            B = [Fv,Fh]
            bb = [self.unwrapped_loadpoints[-1]['locals'][0],self.unwrapped_loadpoints[-1]['locals'][0]]
            # T = 1, Fv = 2, Fh = 3, qv = 4, qh = 5
            nn = [2,3]
        else:
            qv = load.y_component
            qh = load.x_component
            B = [qv,-qv,qh,-qh]
            bb = [self.unwrapped_loadpoints[-1]['locals'][0],
                  self.unwrapped_loadpoints[-1]['locals'][1],
                  self.unwrapped_loadpoints[-1]['locals'][0],
                  self.unwrapped_loadpoints[-1]['locals'][1],
                ]
            nn = [4,4,5,5]

        # print(aa,oo,L)
        # qz ############################################################################################
        qz = 0
        x = Symbol('x')

        # bendpoints
        for i in range(len(B)):
            for j in range(len(aa)):
                if bb[i] == aa[-1]:
                    if nn[i] == 1:
                        qz += B[i] * SingularityFunction(x,bb[i],-2)
                        self.beam.apply_load(B[i], bb[i], -2)
                    if nn[i] == 2:
                        qz += B[i] * SingularityFunction(x,bb[i],-1) * cos(oo[-1])
                        self.beam.apply_load(B[i] * cos(oo[-1]), bb[i], -1)
                    if nn[i] == 3:
                        qz += B[i] * SingularityFunction(x,bb[i],-1) * sin(oo[-1])
                        self.beam.apply_load(B[i] * sin(oo[-1]), bb[i], -1)
                    if nn[i] == 4:
                        qz += B[i] * SingularityFunction(x,bb[i],0) * cos(oo[-1])
                        self.beam.apply_load(B[i] * cos(oo[-1]), bb[i], 0)
                    if nn[i] == 5:
                        qz += B[i] * SingularityFunction(x,bb[i],0) * sin(oo[-1])
                        self.beam.apply_load(B[i] * sin(oo[-1]), bb[i], 0)
                    break
                else:
                    if bb[i] < aa[j]:
                        if nn[i] == 1:
                            qz += B[i] * SingularityFunction(x,bb[i],-2)
                            self.beam.apply_load(B[i], bb[i], -2)
                        if nn[i] == 2:
                            qz += B[i] * SingularityFunction(x,bb[i],-1) * cos(oo[j-1])
                            self.beam.apply_load(B[i] * cos(oo[j-1]), bb[i], -1)
                        if nn[i] == 3:
                            qz += B[i] * SingularityFunction(x,bb[i],-1) * sin(oo[j-1])
                            self.beam.apply_load(B[i] * sin(oo[j-1]), bb[i], -1)
                        if nn[i] == 4:
                            qz += B[i] * SingularityFunction(x,bb[i],0) * cos(oo[j-1])
                            self.beam.apply_load(B[i] * cos(oo[j-1]), bb[i], 0)
                        if nn[i] == 5:
                            qz += B[i] * SingularityFunction(x,bb[i],0) * sin(oo[j-1])
                            self.beam.apply_load(B[i] * sin(oo[j-1]), bb[i], 0)
                        break


        for i in range(len(B)):
            for j in range(len(aa)-1):
                if bb[i] < aa[j]:
                    if nn[i] == 2:
                        qz += B[i] * SingularityFunction(x,aa[j],-1) * (cos(oo[j]) - cos(oo[j-1]))
                        self.beam.apply_load(B[i] * (cos(oo[j]) - cos(oo[j-1])), aa[j], -1)
                    if nn[i] == 3:
                        qz += B[i] * SingularityFunction(x,aa[j],-1) * (sin(oo[j]) - sin(oo[j-1]))
                        self.beam.apply_load(B[i] * (sin(oo[j]) - sin(oo[j-1])), aa[j], -1)
                    if nn[i] == 4:
                        # qz += B[i] * ((SingularityFunction(x,aa[j],0) + (SingularityFunction(x,aa[j],-1) * (aa[j] - bb[i]))) * (cos(oo[j]) - cos(oo[j-1])))
                        qz += B[i] * (SingularityFunction(x,aa[j],0)) * (cos(oo[j]) - cos(oo[j-1]))
                        qz += B[i] * ( SingularityFunction(x,aa[j],-1) * (aa[j] - bb[i])) * (cos(oo[j]) - cos(oo[j-1]))

                        self.beam.apply_load(B[i]*(cos(oo[j]) - cos(oo[j-1])), aa[j], 0)
                        self.beam.apply_load(B[i]* (aa[j] - bb[i]) * (cos(oo[j]) - cos(oo[j-1])), aa[j], -1)

                    if nn[i] == 5:
                        # qz += B[i] * ((SingularityFunction(x,aa[j],0) + SingularityFunction(x,aa[j],-1) * (aa[j] - bb[i])) * (sin(oo[j]) - sin(oo[j-1])))
                        qz += B[i] * (SingularityFunction(x,aa[j],0)) * (sin(oo[j]) - sin(oo[j-1]))
                        qz += B[i] * ( SingularityFunction(x,aa[j],-1) * (aa[j] - bb[i])) * (sin(oo[j]) - sin(oo[j-1]))

                        self.beam.apply_load(B[i]*(sin(oo[j]) - sin(oo[j-1])), aa[j], 0)
                        self.beam.apply_load(B[i]* (aa[j] - bb[i]) * (sin(oo[j]) - sin(oo[j-1])), aa[j], -1)


        self.load_qz += qz



    def _unwrap_structure(self):
        unwrapped_len = 0
        unwrapped_bendpoints = []
        unwrapped_loadpoints = []

        for member in self.members:
            unwrapped_bendpoints.append({'bend_point':[unwrapped_len,unwrapped_len+member.length],'angle':member.angle})

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
        # unwrapped_bendpoints.append({'bend_point':[unwrapped_len,unwrapped_len+member.length],'angle':0})

        x = member.x1
        y = member.y1
        for node in self.nodes:
            if node.x == x and node.y == y:
                for load in node.node_loads:
                    unwrapped_loadpoints.append({'l_id':load.load_id+1,'locals':[unwrapped_len]})
                break

        self.unwrapped_bendpoints = unwrapped_bendpoints
        self.unwrapped_loadpoints = sorted(unwrapped_loadpoints,key=lambda x: x['l_id'])
######################################################################
        aa = []
        oo = []

        for i in range(len(unwrapped_bendpoints)):
            aa.append(unwrapped_bendpoints[i]['bend_point'][0])
            oo.append(unwrapped_bendpoints[i]['angle'])
        aa.append(unwrapped_bendpoints[-1]['bend_point'][1])

        L = unwrapped_len
        self.beam.length = L

        return aa,oo,L

    def _find_unwrapped_position(self, x, y):
        aa, oo, L = self._unwrap_structure()

        unwrapped_position_x = 0
        current_length = 0

        for member in self.members:
            if (member.x1 <= x <= member.x2 or member.x2 <= x <= member.x1) and \
            (member.y1 <= y <= member.y2 or member.y2 <= y <= member.y1):

                local_x = sqrt((member.x1 - x)**2 + (member.y1 - y)**2)
                unwrapped_position_x = current_length + local_x
                break

            current_length += member.length

        if unwrapped_position_x == 0:
            unwrapped_position_x = current_length

        return unwrapped_position_x


    def apply_support(self, x, y, type="pin"):
        # support needs to be done manually because the reaction loads also need to be adjusted for by other loads


        support = self.add_or_update_node(x, y, type)
        self.supports.append(support)

        unwarap_x = self._find_unwrapped_position(x, y)

        if type == "pin":
            R1 = self.beam.apply_support(unwarap_x, type)
            return R1
        elif type == "roller":
            R1 = self.beam.apply_support(unwarap_x, type)
            return R1

        elif type == "fixed":
            T1, R1 = self.beam.apply_support(unwarap_x, type)
            return T1, R1

    def solve_for_reaction_loads(self, *args):
        self.beam.solve_for_reaction_loads(*args)

    def shear_force(self):
        return self.beam.shear_force()

    def plot_shear_force(self):
        self.beam.plot_shear_force()

    def bending_moment(self):
        return self.beam.bending_moment()

    def plot_bending_moment(self):
        self.beam.plot_bending_moment()





#################################################################################
    def draw(self):
        fig, ax = plt.subplots()
        # Set aspect ratio to 'equal' so that the members are plotted proportionally
        ax.set_aspect('equal')
        point_load_color = '#38812F'
        distributed_load_color = '#A30000'

        # Draw each member
        for member in self.members:
            x1, y1, x2, y2 = float(member.x1), float(member.y1), float(member.x2), float(member.y2)
            ax.plot([x1, x2], [y1, y2], color='black',linewidth=3)

        # Draw each node
        for node in self.nodes:
            x, y = float(node.x), float(node.y)
            ax.plot(x, y, color='k')

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
                    arrowprops=dict(arrowstyle='->', color=point_load_color,shrinkA=0, shrinkB=0)  # noqa: C408
                ,zorder=200
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
                        arrowprops=dict(arrowstyle='->', color=distributed_load_color,shrinkA=0, shrinkB=0)  # noqa: C408
                    ,zorder=100
                    )
                for i in range(len(arrow_coords)-1):
                    ax.plot([arrow_coords[i][0], arrow_coords[i+1][0]], [arrow_coords[i][1], arrow_coords[i+1][1]], color=distributed_load_color,zorder=100)
                # print(arrow_coords)
                if load.end_x is not None and load.end_y is not None:
                    square_coords = [(x1, y1), arrow_coords[0], arrow_coords[-1], (x2, y2)]
                    square = patches.Polygon(square_coords, closed=True, color=distributed_load_color, alpha=0.15,zorder=10)

                    ax.add_patch(square)

        # Draw each support
        for support in self.supports:
            x, y = float(support.x), float(support.y)
            if support.node_type == "pin":
                ax.plot(x, y, 'o', color='red')
            elif support.node_type == "fixed":
                ax.plot(x, y, 's', color='blue')

        # make y tick the same step as x tick
        x_ticks = plt.xticks()[0]
        x_step = x_ticks[1] - x_ticks[0]
        y_min, y_max = plt.ylim()

        num_ticks_up = int((y_max // x_step) + 1)
        num_ticks_down = int((-y_min // x_step) + 1)
        new_y_ticks = [i * x_step for i in range(-num_ticks_down, num_ticks_up + 1)]

        plt.yticks(new_y_ticks)
        plt.grid()
        plt.show()

    def draw_unwrapped(self):
        fig, ax = plt.subplots()
        # Set aspect ratio to 'equal' for correct proportions
        ax.set_aspect('equal', adjustable='datalim')

        x_start = 0

        # Draw each member in the unwrapped configuration
        for member in self.members:
            plt.plot([x_start, x_start + member.length], [0, 0], 'o-', color='black')
            x_start += member.length

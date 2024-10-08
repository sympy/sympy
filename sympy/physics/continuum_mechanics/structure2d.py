from sympy.physics.continuum_mechanics.beam import Beam
from sympy.core import Symbol, Basic

from sympy.core.numbers import pi
from sympy.functions.elementary.miscellaneous import sqrt
from sympy.functions.elementary.complexes import Abs
from sympy.functions.elementary.trigonometric import sin, cos, atan2, tan
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

np = import_module(
    "numpy",
    import_kwargs={
        "fromlist": [
            "numpy",
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

    def __repr__(self):
        return f"Load(applied_to={self.applied_to},X1={self.start_x}, Y1={self.start_y},X2={self.end_x},Y2={self.end_y} Load={self.value}, Global_Angle={self.global_angle}, Order={self.order}, X-Component={self.x_component}, Y-Component={self.y_component})"

    def _compute_x_component(self):
        return self.value * cos(self.global_angle)

    def _compute_y_component(self):
        return self.value * sin(self.global_angle) *- 1

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
        self.reaction_loads = {}

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
            if simplify(node.x - x) == 0 and simplify(node.y - y) == 0:
                return f'n_{node.node_id}', 0, None

        # Check if point corresponds to a member
        for member in self.members:
            if member.member_eq != 'vertical':
                # Handle non-vertical members
                member_y = member.member_eq.subs({'x': x})
                if (member.x1 <= x <= member.x2 or member.x2 <= x <= member.x1) and simplify(member_y - y) == 0:
                    local_x_start = sqrt((member.x1 - x) ** 2 + (member.y1 - y) ** 2)
                    if x_end is not None and y_end is not None:
                        local_x_start = sqrt((member.x1 - x_start) ** 2 + (member.y1 - y_start) ** 2)
                        local_x_end = member.length - sqrt((member.x2 - x_end) ** 2 + (member.y2 - y_end) ** 2)
                    else:
                        local_x_end = None
                    return f'm_{member.member_id}', local_x_start, local_x_end
            else:
                # Handle vertical members
                if simplify(member.x1 - x) == 0:
                    if (member.y1 <= y <= member.y2 or member.y2 <= y <= member.y1):
                        local_x_start = Abs(member.y1 - y)
                        if x_end is not None and y_end is not None:
                            local_x_start = Abs(member.y1 - y_start)
                            local_x_end = member.length - Abs(member.y2 - y_end)
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
        if load.order != -2:
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
                # if self.unwrapped_loadpoints[-1]['locals'][1] == L:
                #     B[-1] = 0
                    # B[-2] = -qh

                nn = [4,4,5,5]
        elif load.order == -2:
            qz = load.value
            B = [qz]
            bb = [self.unwrapped_loadpoints[-1]['locals'][0]]
            nn = [1]
            # print('T=',B,'@',bb)



        # print(aa,oo,L,bb,value)
        # print(f'Fv={load.y_component},Fh={load.x_component}')
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
            unwrapped_bendpoints.append({'bend_point': [unwrapped_len, unwrapped_len + member.length], 'angle': member.angle})

            # Process loads on the member
            for load in member.member_loads:
                if load.end_x is not None and load.end_y is not None:
                    unwrapped_loadpoints.append({'l_id': load.load_id, 'locals': [unwrapped_len + load.local_start, unwrapped_len + load.local_end]})
                else:
                    unwrapped_loadpoints.append({'l_id': load.load_id, 'locals': [unwrapped_len + load.local_start]})

            # Process loads at the start node
            x_start, y_start = member.x1, member.y1
            for node in self.nodes:
                if simplify(node.x - x_start) == 0 and simplify(node.y - y_start) == 0:
                    for load in node.node_loads:
                        unwrapped_loadpoints.append({'l_id': load.load_id, 'locals': [unwrapped_len]})
                    break

            unwrapped_len += member.length

            # Process loads at the end node
            x_end, y_end = member.x2, member.y2
            for node in self.nodes:
                if simplify(node.x - x_end) == 0 and simplify(node.y - y_end) == 0:
                    for load in node.node_loads:
                        unwrapped_loadpoints.append({'l_id': load.load_id, 'locals': [unwrapped_len]})
                    break

        self.unwrapped_bendpoints = unwrapped_bendpoints
        self.unwrapped_loadpoints = sorted(unwrapped_loadpoints, key=lambda x: x['l_id'])
        ### emulate Alexes input for his algorithm
        aa = []
        oo = []

        for i in range(len(unwrapped_bendpoints)):
            aa.append(unwrapped_bendpoints[i]['bend_point'][0])
            oo.append(unwrapped_bendpoints[i]['angle'])
        aa.append(unwrapped_bendpoints[-1]['bend_point'][1])

        L = unwrapped_len
        self.beam.length = L

        return aa,oo,L
        ### emulate Alexes input for his algorithm

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
        support = self.add_or_update_node(x, y, type)
        self.supports.append(support)

        unwarap_x = self._find_unwrapped_position(x, y)

        Rh = Symbol(f'R_h__{round(x,2)},__{round(y,2)}')
        Rv = Symbol(f'R_v__{round(x,2)},__{round(y,2)}')
        T = Symbol(f'T__{round(x,2)},__{round(y,2)}')

        if type == "pin":
            load_v = -1 * Rv
            load_h = Rh
            self.apply_load(start_x=x, start_y=y, value=load_v, global_angle=90, order=-1)
            self.apply_load(start_x=x, start_y=y, value=load_h, global_angle=0, order=-1)

            # Mark these loads as support reactions
            self.loads[-1].is_support_reaction = True
            self.loads[-2].is_support_reaction = True

            self.beam.bc_deflection.append((unwarap_x, 0))
            return Rv, Rh

        elif type == "roller":
            load_v = -1 * Rv
            self.apply_load(start_x=x, start_y=y, value=load_v, global_angle=90, order=-1)

            # Mark this load as a support reaction
            self.loads[-1].is_support_reaction = True

            self.beam.bc_deflection.append((unwarap_x, 0))
            return Rv

        elif type == "fixed":
            load_t = T
            load_v = -1 * Rv
            load_h = Rh
            self.apply_load(start_x=x, start_y=y, value=load_t, global_angle=0, order=-2)
            self.apply_load(start_x=x, start_y=y, value=load_v, global_angle=90, order=-1)
            self.apply_load(start_x=x, start_y=y, value=load_h, global_angle=0, order=-1)

            # Mark these loads as support reactions
            self.loads[-1].is_support_reaction = True
            self.loads[-2].is_support_reaction = True
            self.loads[-3].is_support_reaction = True

            self.beam.bc_deflection.append((unwarap_x, 0))
            self.beam.bc_slope.append((unwarap_x, 0))
            return T, Rv, Rh

    def solve_for_reaction_loads(self, *args):

        # Split arguments into vertical and horizontal reaction loads
        reaction_loads_vertical = [arg for arg in args if 'R_v' in str(arg)]
        reaction_loads_horizontal = [arg for arg in args if 'R_h' in str(arg)]
        reaction_moments = [arg for arg in args if 'T_' in str(arg)]

        args_for_beam_solver = tuple(reaction_loads_vertical + reaction_loads_horizontal + reaction_moments)

        # print("Debug - reaction load arguments:", test)

        # Solve for moment and vertical reaction loads using the beam solver
        self.beam.solve_for_reaction_loads(*args_for_beam_solver)

        # Compute the horizontal reaction load by summing up all horizontal forces
        sum_horizontal = 0
        for load in self.loads:
            if isinstance(load.value, (int, float)):
                if load.order == -1:
                    sum_horizontal += load.x_component
                elif load.order == 0:
                    length = sqrt((load.end_y - load.start_y)**2 + (load.end_x - load.start_x)**2)
                    sum_horizontal += length * load.x_component

        # Substitute horizontal reactions with their solved values
        horizontal_key = {arg: float(-1 * sum_horizontal) for arg in args if 'R_h' in str(arg)}

        # Update solved reaction loads dictionary with horizontal reaction values
        for key in self.beam._reaction_loads.items():
            key = key[0]
            self.beam._reaction_loads[key] = self.beam._reaction_loads[key].subs(horizontal_key)

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
            vertical_key = {arg: self.beam._reaction_loads[arg] for arg in args if 'R_v' in str(arg)}
            bending_moment_key = {arg: self.beam._reaction_loads[arg] for arg in args if 'T_' in str(arg)}

            for load in self.loads:
                if isinstance(load.value, Basic) and load.order != -2:
                    # Substitute horizontal and vertical reaction loads
                    load.value = load.value.subs(vertical_key).subs(horizontal_key)
                    load.y_component = load.y_component.subs(vertical_key).subs(horizontal_key)
                    load.x_component = load.x_component.subs(vertical_key).subs(horizontal_key)

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


    def shear_force(self,x=None,y=None):
        """
        Compute the shear force at a given point on the beam.

        If no arguments are provided, it returns the shear force equation.

        if x is provided, it returns the shear force at the UNWRAPPED LOC.

        if x AND y are provided, it returns the shear force at the given x and y coordinates.
        """

        if x is None:
            self.beam.shear_force()
        else:
            if y is None:
                x_symbol = Symbol('x')
                return self.beam.shear_force().subs(x_symbol,x)
            else:
                unwarap_x = self._find_unwrapped_position(x, y)
                return self.beam.shear_force().subs(Symbol('x'),unwarap_x)

    def plot_shear_force(self):
        self.beam.plot_shear_force()

    def bending_moment(self,x=None,y=None):
        """
        Compute the bending moment at a given point on the beam.

        If no arguments are provided, it returns the bending moment equation.

        if x is provided, it returns the bending moment at the UNWRAPPED LOC.

        if x AND y are provided, it returns the bending moment at the given x and y coordinates.
        """

        if x is None:
            return self.beam.bending_moment()
        else:
            if y is None:
                x_symbol = Symbol('x')
                return self.beam.bending_moment().subs(x_symbol,x)
            else:
                unwarap_x = self._find_unwrapped_position(x, y)
                return self.beam.bending_moment().subs(Symbol('x'),unwarap_x)

    def plot_bending_moment(self):
        self.beam.plot_bending_moment()

    def summary(self, verbose=True, round_digits=None):
        title = 'Structure Summary'
        line_length = 50

        print(f'{"=" * ((line_length - len(title)) // 2)} {title} {"=" * ((line_length - len(title)) // 2)}')

        if len(self.reaction_loads) == 0:
            print('\nPlease solve for reaction loads first')
        elif verbose and len(self.reaction_loads) != 0:
            dx = 1e-5

            print('\nReaction Loads:')
            for key, value in self.reaction_loads.items():
                support_name = str(key).split('__')[0]
                support_x_loc = float(str(key).split('__')[1].strip(','))
                support_y_loc = float(str(key).split('__')[2].strip(','))
                unwrapped_xl_loc = self._find_unwrapped_position(support_x_loc, support_y_loc)

                # Check if the value contains any symbols. If not, round it.
                if round_digits is not None and not value.has(Symbol):
                    value = round(float(value), round_digits)

                support_str = f'{support_name:<5} [{support_x_loc:.2f},{support_y_loc:.2f}]  ({unwrapped_xl_loc:.2f})'
                print(f'{support_str:<40} = {value}')

            print('\nPoints of Interest - Bending Moment:')
            bend_points = [float(point) for item in self.unwrapped_bendpoints for point in item['bend_point']]
            bend_points = set(bend_points)
            bend_points = sorted(bend_points)
            for point in bend_points:
                bending_moment_value = self.bending_moment(point)
                if round_digits is not None and not bending_moment_value.has(Symbol):
                    bending_moment_value = round(float(bending_moment_value), round_digits)
                string_text = f'bending_moment at [x.xx,y.yy]  ({point:.02f})'
                print(f'{string_text:<40} = {bending_moment_value}')

            print('\nPoints of Interest - Shear Force:')
            load_points = [float(local) for point in self.unwrapped_loadpoints for local in point['locals']]
            load_points = set(load_points + bend_points)
            load_points = sorted(load_points)

            for point in load_points:
                shear_force_value = self.shear_force(point)
                if round_digits is not None and not shear_force_value.has(Symbol):
                    shear_force_value = round(float(shear_force_value), round_digits)

                if point == load_points[0]:
                    string_text = f'shear_force at [x.xx,y.yy]  ({point:.02f}+)'
                    print(f'{string_text:<40} = {shear_force_value}')
                elif point == load_points[-1]:
                    string_text = f'shear_force at [x.xx,y.yy]  ({point:.02f}-)'
                    print(f'{string_text:<40} = {shear_force_value}')
                else:
                    string_text = f'shear_force at [x.xx,y.yy]  ({point:.02f}-)'
                    print(f'{string_text:<40} = {shear_force_value}')
                    string_text = f'shear_force at [x.xx,y.yy]  ({point:.02f}+)'
                    print(f'{string_text:<40} = {shear_force_value}')



#################################################################################
    def draw(self,forced_load_size=None,show_load_values=False):
        fig, ax = plt.subplots()
        ax.set_aspect('equal')

        # Define colors
        point_load_color = 'blue'          # Point load color
        distributed_load_color = '#A30000'    # Distributed load color
        support_reaction_color = 'darkorange'    # Support reaction load color (blue)
        uncomputed_color = '#696969'
        scale = 0.75


        # Draw members
        for member in self.members:
            x1, y1 = float(member.x1), float(member.y1)
            x2, y2 = float(member.x2), float(member.y2)
            ax.plot([x1, x2], [y1, y2], color='black', linewidth=3, zorder=10)

        # Draw loads
        for load in self.loads:
            x, y = float(load.start_x), float(load.start_y)
            is_support_reaction = getattr(load, 'is_support_reaction', False)

            # Set color based on whether it's a support reaction
            if is_support_reaction:
                load_color = support_reaction_color
            else:
                if load.order < 0 or load.order == 1:
                    load_color = point_load_color
                else:
                    load_color = distributed_load_color

            # Set transparency based on whether the load has unsolved symbols
            alpha_value = 1.0
            if isinstance(load.value, Basic) and load.value.has(Symbol):
                # alpha_value = 0.5
                # uncomputed_color = '#696969'
                load_color = uncomputed_color

            if load.order == -2:
                radius = 0.5

                if load.value.is_negative:
                    angle_start = 180+30
                    angle_end = 270-10
                    angle_moment_icon = angle_start
                    angle2moment_icon = angle_start -90

                else:
                    angle_start = 180+10
                    angle_end = 270-30
                    angle_moment_icon = angle_end
                    angle2moment_icon = angle_end +90


                arc = patches.Arc(
                    (x, y),
                    width=2.5*radius,
                    height=2.5*radius,
                    theta1=angle_start,
                    theta2=angle_end,
                    color=load_color,
                    alpha=alpha_value,
                    linewidth=2,
                )
                ax.add_patch(arc)
                # Add arrowhead at staty of arc
                arrowhead = patches.FancyArrow(
                    x + radius*1.25 * np.cos(np.radians(angle_moment_icon)),
                    y + radius*1.25 * np.sin(np.radians(angle_moment_icon)),
                    0.1 * np.cos(np.radians(angle2moment_icon)),
                    0.1 * np.sin(np.radians(angle2moment_icon)),
                    width=0.02,
                    head_width=0.1,
                    head_length=0.2,
                    color=load_color,
                    alpha=alpha_value,
                )
                ax.add_patch(arrowhead)
                if show_load_values:
                    if load_color == uncomputed_color:
                        color_text = 'black'
                    else:
                        color_text = load_color

                    if isinstance(load.value, Basic):
                        if not load.value.has(Symbol):
                            plt.text(x - radius*2, y-radius*2, f'{float(load.value):0.2f}', fontsize=8, color=color_text, zorder=100,
                                    ha="center", va="center", bbox={"facecolor": "lightgrey", "alpha": 0.85, "edgecolor": "none"})
                        else:
                            plt.text(x - radius*2, y-radius*2, f'{str(load.value).split("__")[0]}', fontsize=8, color=color_text, zorder=100,
                                    ha="center", va="center", bbox={"facecolor": "lightgrey", "alpha": 0.85, "edgecolor": "none"})
                    else:
                        plt.text(x - radius*2, y-radius*2, f'{float(load.value):0.2f}', fontsize=8, color=color_text, zorder=100,
                                ha="center", va="center", bbox={"facecolor": "lightgrey", "alpha": 0.85, "edgecolor": "none"})



            elif load.order < 0 or load.order == 1:
                # Point load
                if forced_load_size is not None:
                    load_length = forced_load_size
                else:
                    if isinstance(load.value, (int, float)):
                        load_length = abs(float(load.value)) / 10.0
                    else:
                        if is_support_reaction and isinstance(load.value, Basic) and load.value.has(Symbol):
                            # If the load is an unsolved support reaction, set the length to 0.5
                            load_length = 1.0
                        else:
                            # Default length for other symbolic or unsolved loads
                            load_length = 1.0

                # Only apply the offset if it is a support reaction
                if is_support_reaction:
                    support_offset = scale * 0.65
                    if load.y_component.has(Symbol):
                        sign = load.y_component
                        if sign.is_negative:
                            y = y + support_offset
                        else:
                            y = y - support_offset
                    elif load.x_component.has(Symbol):
                        sign2 = load.x_component
                        if sign2.is_negative:
                            x = x + support_offset
                        else:
                            x = x - support_offset

                    elif load.y_component.has(Symbol) == False and load.x_component.has(Symbol) == False:
                        if is_support_reaction and float(load.y_component) == 0:
                            if float(load.x_component) > 0:
                                x = x -support_offset
                            else:
                                x = x + support_offset
                        if is_support_reaction and float(load.x_component) == 0:
                            if float(load.y_component) > 0:
                                y = y  + support_offset
                            else:
                                y = y - support_offset


                angle = float(load.global_angle)
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
                    alpha=alpha_value,
                    zorder=100,
                )
                ax.add_patch(arrow_patch)
                if load_color == uncomputed_color:
                    color_text = 'black'
                else:
                    color_text = load_color

                if show_load_values:
                    if isinstance(load.value, Basic):
                        if not load.value.has(Symbol):
                            plt.text(x - dx, y - dy, f'{float(load.value):0.2f}', fontsize=8, color=color_text, zorder=100,
                                    ha="center",
                                    va="center",
                                    bbox={"facecolor": "lightgrey", "alpha": 0.850, "edgecolor": "none"},)
                        else:
                            plt.text(x - dx, y - dy, f'{str(load.value).split("__")[0]}', fontsize=8, color=color_text, zorder=100,
                                    ha="center",
                                    va="center",
                                    bbox={"facecolor": "lightgrey", "alpha": 0.850, "edgecolor": "none"},)
                    else:
                        plt.text(x - dx, y - dy, f'{float(load.value):0.2f}', fontsize=8, color=color_text, zorder=100,
                                ha="center",
                                va="center",
                                bbox={"facecolor": "lightgrey", "alpha": 0.850, "edgecolor": "none"},)


            else:
                # Distributed load
                if forced_load_size is not None:
                    load_length = forced_load_size
                else:
                    if isinstance(load.value, (int, float)):
                        load_length = abs(float(load.value)) / 10.0
                    else:
                        # Default length for other symbolic or unsolved loads
                        load_length = 1.0

                x1, y1 = float(load.start_x), float(load.start_y)
                if load.end_x is not None and load.end_y is not None:
                    x2, y2 = float(load.end_x), float(load.end_y)
                else:
                    member = next((m for m in self.members if
                                   (float(m.x1) == x1 and float(m.y1) == y1) or
                                   (float(m.x2) == x1 and float(m.y2) == y1)), None)
                    if member is not None:
                        x2, y2 = float(member.x2), float(member.y2)
                    else:
                        continue  # Skip if member not found

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

                    # load_length = 0.5
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
                        alpha=alpha_value,
                        zorder=90,
                    )
                    ax.add_patch(arrow_patch)
                    arrow_tails_x.append(x_point - dx)
                    arrow_tails_y.append(y_point - dy)

                    arrow_heads_x.append(x_point)
                    arrow_heads_y.append(y_point)

                plt.plot(arrow_tails_x, arrow_tails_y, color=load_color, alpha=alpha_value)

                arrow_coords = [
                    (arrow_tails_x[0], arrow_tails_y[0]),   # First tail point
                    (arrow_tails_x[-1], arrow_tails_y[-1]), # Last tail point
                    (arrow_heads_x[-1], arrow_heads_y[-1]), # Last head point
                    (arrow_heads_x[0], arrow_heads_y[0])    # First head point
]
                # make grey when unsolved
                if isinstance(load.value, Basic) and load.value.has(Symbol):
                    distributed_load_color = '#696969'

                polygon = plt.Polygon(arrow_coords, closed=True, fill=True, color=distributed_load_color, edgecolor=None, alpha=0.25,zorder=89)
                ax.add_patch(polygon)


                if show_load_values:
                    avrage_x = (arrow_tails_x[0] + arrow_tails_x[-1]) / 2
                    avrage_y = (arrow_tails_y[0] + arrow_tails_y[-1]) / 2

                    if load_color == uncomputed_color:
                        color_text = 'black'
                    else:
                        color_text = load_color
                    if isinstance(load.value, Basic):
                        if not load.value.has(Symbol):
                            plt.text(avrage_x, avrage_y, f'{float(load.value):0.2f}', fontsize=8, color=color_text, zorder=100,
                                    ha="center",
                                    va="center",
                                    bbox={"facecolor": "lightgrey", "alpha": 0.850, "edgecolor": "none"},)
                        else:
                            plt.text(avrage_x, avrage_y, f'{str(load.value).split("__")[0]}', fontsize=8, color=color_text, zorder=100,
                                    ha="center",
                                    va="center",
                                    bbox={"facecolor": "lightgrey", "alpha": 0.850, "edgecolor": "none"},)
                    else:
                        plt.text(avrage_x, avrage_y, f'{float(load.value):0.2f}', fontsize=8, color=color_text, zorder=100,
                                ha="center",
                                va="center",
                                bbox={"facecolor": "lightgrey", "alpha": 0.850, "edgecolor": "none"},)




        # Draw supports

        for support in self.supports:
            x, y = float(support.x), float(support.y)
            # x, y = float(support.x), float(support.y)

            if support.node_type == "pin":
                # scale = 0.5  # Adjust this to make the polygon larger or smaller
                triangle_height = 0.5 * scale  # The height of the triangle
                angle = 34 * (pi / 180)  # Convert degrees to radians

                # Calculate the half-width of the triangle using tan(angle)
                half_width = triangle_height * float(tan(angle))

                # Define triangle vertices for the pin support (top point is at x, y)
                triangle_vertices = [
                    [x, y],  # Top point (pivot)
                    [x - half_width, y - triangle_height],  # Bottom-left corner
                    [x + half_width, y - triangle_height]   # Bottom-right corner
                ]


                # Create the triangle with scaled vertices
                triangle = patches.Polygon(
                    triangle_vertices,
                    edgecolor='k',
                    facecolor='none',
                    linewidth=1.5,
                    zorder = 20
                )
                ax.add_patch(triangle)

                # Draw a small circle on top to represent the pivot point
                pivot_circle = patches.Circle(
                    (x, y), 0.1 * scale, edgecolor='k', facecolor='white', linewidth=1.5,zorder=20,
                )
                ax.add_patch(pivot_circle)


            elif support.node_type == "roller":
                # scale = 1  # Adjust this to make the polygon larger or smaller
                triangle_height = 0.5 * scale  # The height of the triangle
                angle = 34 * (pi / 180)  # Convert degrees to radians

                # Calculate the half-width of the triangle using tan(angle)
                half_width = triangle_height * float(tan(angle))

                # Define triangle vertices for the pin support (top point is at x, y)
                triangle_vertices = [
                    [x, y],  # Top point (pivot)
                    [x - half_width, y - triangle_height],  # Bottom-left corner
                    [x + half_width, y - triangle_height]   # Bottom-right corner
                ]

                # Create the triangle with scaled vertices
                triangle = patches.Polygon(
                    triangle_vertices,
                    edgecolor='k',
                    facecolor='none',
                    linewidth=1.5,
                    zorder = 20,
                )
                ax.add_patch(triangle)

                # Draw a small circle on top to represent the pivot point
                pivot_circle = patches.Circle(
                    (x, y), 0.1 * scale, edgecolor='k', facecolor='white', linewidth=1.5,zorder=20,
                )
                ax.add_patch(pivot_circle)

                plt.plot([x - half_width,x + half_width], [y - triangle_height - (scale*0.1),y - triangle_height - (scale*0.1)], color='k', zorder=20)

            elif support.node_type == "fixed":
                box_size =  scale*0.5
                x_box = x - box_size / 2
                y_box = y - box_size / 2


                scaled_vertices = [
                    [x_box, y_box],                      # Bottom-left corner
                    [(x_box + box_size), y_box],          # Bottom-right corner
                    [(x_box + box_size), (y_box + box_size)],  # Top-right corner
                    [x_box, (y_box + box_size)]           # Top-left corner
                ]

                # Create the polygon with scaled vertices
                polygon = patches.Polygon(
                    scaled_vertices,     # Vertices scaled by the scale factor
                    edgecolor='k',
                    facecolor="none",
                    linewidth=1.5       # Optional: set the thickness of the border
                )
                ax.add_patch(polygon)

                # Add the short lines inside the bottom edge of the box
                num_lines = 4  # Number of diagonal lines (can adjust for more or fewer lines)
                for i in range(num_lines):
                    # Start on the bottom edge
                    x_start = x_box + i * (box_size / num_lines)  # Divide the bottom edge into equal segments
                    y_start = y_box  # Start from the bottom edge of the box

                    # End inside the box with a diagonal angle
                    x_end = x_box + (i + 1) * (box_size / num_lines)  # Slight shift in x
                    y_end = y_box + (box_size / num_lines)  # Slight shift in y to create a diagonal

                    ax.plot([x_start, x_end], [y_start, y_end], color='black', linewidth=1)  # Diagonal lines

        x_ticks = plt.xticks()[0]
        y_ticks = plt.yticks()[0]

        if len(x_ticks) > 1 and len(y_ticks) > 1:
            x_step = x_ticks[1] - x_ticks[0]
            y_step = y_ticks[1] - y_ticks[0]

            step = min(x_step, y_step)

            x_min, x_max = plt.xlim()
            y_min, y_max = plt.ylim()

            new_x_ticks = [i * step for i in range(int(x_min // step), int(x_max // step) + 1)]
            new_y_ticks = [i * step for i in range(int(y_min // step), int(y_max // step) + 1)]

            plt.xticks(new_x_ticks)
            plt.yticks(new_y_ticks)

        ax.grid(True,zorder=10)
        plt.show()



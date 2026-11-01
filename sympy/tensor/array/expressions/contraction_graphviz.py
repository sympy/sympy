from dataclasses import dataclass, field
from functools import singledispatchmethod

from sympy import MatrixSymbol, Expr
from sympy.tensor.array.expressions import ArraySymbol, ArrayContraction, ArrayDiagonal, PermuteDims, ArrayTensorProduct


@dataclass
class _GraphUnit:
    name_unique: str
    name: str | None
    indices: list[int | None]
    shape: tuple[Expr | int]

@dataclass
class _GraphUnits:
    units: list[_GraphUnit]
    contractions: list[tuple[str, int]] = field(default_factory=list)
    diagonals: list[tuple[str, int]] = field(default_factory=list)


class ArrayDotGraphPrinter:

    def __init__(self):
        self._counting: int = 0

    @singledispatchmethod
    def _s(self, expr) -> _GraphUnits:
        raise NotImplementedError

    @_s.register
    def _(self, expr: ArrayContraction) -> _GraphUnits:
        sup_graph_units = self._s(expr.expr)
        new_contractions = []
        for contr in expr.contraction_indices:
            new_contraction = []
            counter = 0
            for graph_unit in sup_graph_units.units:
                for i, e in enumerate(graph_unit.indices):
                    if e is None:
                        continue
                    if counter in contr:
                        graph_unit.indices[i] = -1
                        new_contraction.append((graph_unit.name_unique, i))
                    counter += 1
            new_contractions.append(tuple(new_contraction))
        counter = 0
        for graph_unit in sup_graph_units.units:
            for i, e in enumerate(graph_unit.indices):
                if e is None:
                    continue
                elif e == -1:
                    graph_unit.indices[i] = None
                    counter += 1
                else:
                    graph_unit.indices[i] -= counter
        sup_graph_units.contractions.extend(new_contractions)
        return sup_graph_units

    @_s.register
    def _(self, expr: ArrayDiagonal) -> _GraphUnits:
        graph_units = self._s(expr.expr)
        dim_non_diag = len(expr.shape) - len(expr.diagonal_indices)
        for j, diag_ind in enumerate(expr.diagonal_indices):
            for graph_unit in graph_units.units:
                for i, e in enumerate(graph_unit.indices):
                    if e is None:
                        continue
                    if e in diag_ind:
                        graph_unit.indices[i] = dim_non_diag + j
        return graph_units

    @_s.register
    def _(self, expr: PermuteDims) -> _GraphUnits:
        graph_units = self._s(expr.expr)
        inverse_perm = expr.permutation**(-1)
        for graph_unit in graph_units.units:
            for i, e in enumerate(graph_unit.indices):
                if e is None:
                    continue
                graph_unit.indices[i] = inverse_perm(e)  # expr.permutation(e)
        return graph_units

    @_s.register
    def _(self, expr: ArrayTensorProduct) -> _GraphUnits:
        graph_units: list[_GraphUnit] = [i for arg in expr.args for i in self._s(arg).units]
        counter = 0
        for unit in graph_units:
            local_count = 0
            for i, e in enumerate(unit.indices):
                if e is None:
                    continue
                unit.indices[i] = counter + e
                local_count += 1
            counter += local_count
        return _GraphUnits(graph_units)

    @_s.register
    def _(self, expr: ArraySymbol) -> _GraphUnits:
        return _GraphUnits([
            _GraphUnit(
                self._get_name(expr.name),
                expr.name,
                [i for i, e in enumerate(expr.shape)], expr.shape)
        ])

    @_s.register
    def _(self, expr: MatrixSymbol) -> _GraphUnits:
        return _GraphUnits([
            _GraphUnit(
                self._get_name(expr.name),
                expr.name,
                [i for i, e in enumerate(expr.shape)], expr.shape)
        ])

    def _get_name(self, base_name: str | None) -> str:
        if base_name is None:
            base_name = "X"
        self._counting += 1
        return f"{base_name}{self._counting}"

    def _graph_units_to_dot(self, graph_units: _GraphUnits, graph_name="ArrayExpr") -> str:

        lines = [
            f"graph {graph_name} {{",
            "    rankdir = LR;",
            "    splines = true;",
            "",
            "    node [shape=record, style=rounded, fontsize=12, height=0.4, width=1.2];",
            "",
        ]

        for unit in graph_units.units:
            labels = " | ".join([f"<p{i}> {'' if idx is None else str(idx) + ' '}[{shape}]" for i, (shape, idx) in enumerate(zip(unit.shape, unit.indices))])
            line = f"""    {unit.name_unique} [label="{{ {labels} }}", xlabel="{unit.name}"];"""
            lines.append(line)

        lines.append("")

        for contraction in graph_units.contractions:
            line = " "*4 + " -- ".join([f"{name}:p{pindex}" for name, pindex in contraction]) + ";"
            lines.append(line)

        lines.append("}")

        return "\n".join(lines)

    def doprint(self, expr):
        self._counting = 0
        gu = self._s(expr)
        return self._graph_units_to_dot(gu)

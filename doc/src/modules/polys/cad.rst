.. _polys-cad:

=====================================
Cylindrical algebraic decomposition
=====================================

.. currentmodule:: sympy.polys.cad

A cylindrical algebraic decomposition (CAD) of `\mathbb{R}^n` adapted to a
set of polynomials in `x_1, \ldots, x_n` is a partition of `\mathbb{R}^n`
into finitely many connected *cells* on each of which every polynomial has a
constant sign. The decomposition is *cylindrical*: the projections of the
cells onto the first `k` coordinates form a decomposition of `\mathbb{R}^k`
for every `k`. Each cell comes with an exact sample point, so any property
which only depends on the signs of the polynomials can be decided by looking
at finitely many points. This is the basis of the decision procedure and of
quantifier elimination for the first order theory of the real numbers
introduced by Collins.

The implementation follows the classical two phases:

1. **Projection**. Starting from a squarefree basis of the input, a
   projection operator computes polynomials in `x_1, \ldots, x_{n-1}` whose
   sign invariance guarantees that the real roots of the input polynomials
   in `x_n` are delineable, i.e. given by finitely many continuous non
   crossing functions over every cell. This is repeated down to univariate
   polynomials in `x_1`. Two operators are available, McCallum's (the
   default, small but requiring the polynomials to be *well-oriented*) and
   Hong's (always valid).

2. **Lifting**. The real line is decomposed at the real roots of the
   univariate polynomials into points (*sections*) and open intervals
   (*sectors*). Over the sample point of every cell the polynomials of the
   next level become univariate; their real roots split the cylinder over
   the cell into sections and sectors again. Sample points are exact: they
   are rational numbers for sectors and real algebraic numbers, kept in a
   single algebraic number field, for sections.

Examples
========

The circle `x^2 + y^2 = 1` gives 13 cells: the two sides of the vertical
lines `x = \pm 1`, the lines themselves split at the circle, and the strip
in between split by the two arcs of the circle.

>>> from sympy.abc import x, y
>>> from sympy.polys.cad import cylindrical_algebraic_decomposition
>>> cad = cylindrical_algebraic_decomposition([x**2 + y**2 - 1], [x, y])
>>> cad
CAD(13 cells, x, y)
>>> [cell.point for cell in cad if cell.signs == (0,)]
[(-1, 0), (0, -1), (0, 1), (1, 0)]
>>> [cell.point for cell in cad if cell.signs == (-1,)]
[(0, 0)]

Cells are indexed like in Collins' work: the `k`-th entry of the index is
odd for a sector and even for a section of the cylinder over the parent
cell. The dimension of a cell is the number of odd entries.

>>> cell = cad.cells[7]
>>> cell.index, cell.point, cell.dimension
((3, 4), (0, 1), 1)

The signs of the input polynomials on each cell are read off the sample
point, so a system of polynomial equations and inequalities is satisfiable
if and only if it holds at one of the sample points:

>>> from sympy.polys.cad import sample_points
>>> sample_points((x**2 + y**2 < 1) & (x > y), [x, y])
[{x: 0, y: -1/2}, {x: CRootOf(2*x**2 - 1, 1), y: 0}, {x: 3/4, y: 0}]
>>> sample_points((x**2 + y**2 < 1) & (x > 1), [x, y])
[]

Quantifiers are eliminated by propagating the truth of a formula from the
cells of the full decomposition down to the cells of the space of the free
variables. With a single free variable the result is a union of intervals
with exact endpoints; with more free variables it is a formula in the
polynomials computed by the projection.

>>> from sympy import Eq
>>> from sympy.abc import a, b, c
>>> from sympy.polys.cad import quantifier_elimination, decide, solution_set
>>> quantifier_elimination(x**2 + b*x + c > 0, [('forall', x)])
b**2 - 4*c < 0
>>> quantifier_elimination(Eq(x**2 + a*x + b, 0), [('exists', x)])
a**2 - 4*b >= 0
>>> solution_set(Eq(x**2 + y**2, 1) & (y > x), x, [('exists', y)])
Interval.Ropen(-1, CRootOf(2*x**2 - 1, 1))
>>> decide(Eq(y, x**2), [('forall', x), ('exists', y)])
True
>>> decide(Eq(x, y**2), [('forall', x), ('exists', y)])
False

The number of cells grows quickly with the number of variables and the
degrees: this implementation is meant for problems with a few variables
and moderate degrees.

Reference
=========

.. autofunction:: cylindrical_algebraic_decomposition

.. autoclass:: CAD
   :members:

.. autoclass:: CADCell
   :members:

.. autoexception:: NotWellOriented

.. autofunction:: quantifier_elimination

.. autofunction:: decide

.. autofunction:: solution_set

.. autofunction:: sample_points

Projection
----------

.. automodule:: sympy.polys.cad.projection
   :members: projection_sets, mccallum_projection, hong_projection, squarefree_basis

Sample points
-------------

.. automodule:: sympy.polys.cad.samplepoints
   :members: SamplePoint, compare_real, rational_between, rational_below, rational_above, simplest_between

Principal subresultant coefficients
-----------------------------------

.. autofunction:: sympy.polys.euclidtools.dup_psc

.. autofunction:: sympy.polys.euclidtools.dmp_psc

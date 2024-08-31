(active-deprecations)=
# List of active deprecations

This pages lists all active deprecations in the SymPy codebase. See the
{ref}`deprecation-policy` page for a description of SymPy's deprecation
policy, as well as instructions for contributors on how to deprecate things.

In particular, the deprecation policy for SymPy is for deprecations to last at
least **1 year** after the first major release that includes the deprecation.
After that period, the deprecated functionality may be removed from SymPy, and
code will need to be updated to use the replacement feature to continue
working.

During the deprecation period, a `SymPyDeprecationWarning` message will be
printed whenever the deprecated functionality is used. It is recommended for
users to update their code so that it does not use deprecated functionality,
as described below for each given deprecation.

(silencing-sympy-deprecation-warnings)=
## Silencing SymPy Deprecation Warnings

To silence SymPy deprecation warnings, add a filter using the
[`warnings`](https://docs.python.org/3/library/warnings.html) module. For
example:

```py
import warnings
from sympy.utilities.exceptions import SymPyDeprecationWarning

warnings.filterwarnings(
    # replace "ignore" with "error" to make the warning raise an exception.
    # This useful if you want to test you aren't using deprecated code.
    "ignore",

    # message may be omitted to filter all SymPyDeprecationWarnings
    message=r"(?s).*<regex matching the warning message>",

    category=SymPyDeprecationWarning,
    module=r"<regex matching your module>"
)
```

Here `(?s).*<regex matching the warning message>` is a regular expression
matching the warning message. For example, to filter a warning about
`sympy.printing`, you might use `message=r"(?s).*sympy\.printing"`. The
leading `(?s).*` is there because the warnings module matches `message`
against the start of the warning message, and because typical warning messages
span multiple lines.

`<regex matching your module>` should be a regular expression matching your
module that uses the deprecated code. It is recommended to include this so
that you don't also silence the same warning for unrelated modules.

This same pattern may be used to instead turn `SymPyDeprecationWarning` into
an error so that you can test that you aren't using deprecated code. To do
this, replace `"ignore"` with `"error"` in the above example. You may
also omit `message` to make this apply to all `SymPyDeprecationWarning`
warnings.

If you are using pytest, you can use the [pytest warnings filtering
capabilities](https://docs.pytest.org/en/latest/how-to/capture-warnings.html)
to either ignore `SymPyDeprecationWarning` or turn them into errors.

```{note}
The Python [`-W`
flag](https://docs.python.org/3/using/cmdline.html#cmdoption-W) and
[`PYTHONWARNINGS` environment
variable](https://docs.python.org/3/using/cmdline.html#envvar-PYTHONWARNINGS)
will NOT work to filter SymPy deprecation warnings (see [this blog
post](https://nedbatchelder.com/blog/201810/why_warnings_is_mysterious.html)
by Ned Batchelder and [this SymPy
issue](https://github.com/sympy/sympy/issues/15130) for details on why). You
will need to either add a `warnings` filter as above or use pytest to filter
SymPy deprecation warnings.
```

## Version 1.14

There are no deprecations yet for SymPy 1.14.

## Version 1.13

(deprecated-mechanics-body-class)=
### Deprecated mechanics Body class

The ``Body`` class in the ``sympy.physics.mechanics`` module has been
deprecated. It was introduced to support the joints framework. However, it
causes several problems because it represents both rigid bodies and particles.
``Body`` has now been fully replaced by ``RigidBody`` and ``Particle``.
Previously, one could create a simple rigid body or particle using only the
``Body`` class:

```py
>>> from sympy import symbols
>>> from sympy.physics.mechanics import Body
>>> Body("rigid_body")  # doctest: +SKIP
rigid_body
>>> Body("particle", mass=symbols("m"))  # doctest: +SKIP
particle
```

Now they should be created using the ``RigidBody`` and ``Particle`` class:

```py
>>> from sympy.physics.mechanics import RigidBody, Particle
>>> RigidBody("rigid_body")
rigid_body
>>> Particle("particle")
particle
```

(deprecated-mechanics-jointsmethod)=
### Deprecated mechanics JointsMethod

The ``JointsMethod`` class in the ``sympy.physics.mechanics`` module has been
deprecated. It was introduced to support the joints framework, but it has been
fully replaced due to limitations in its design. Previously, one could construct
as system solely consisting out of bodies and joints, which were then parsed by
``JointsMethod`` to a backend, like ``KanesMethod`` to form the equations of
motion.

```py
>>> from sympy import symbols
>>> from sympy.physics.mechanics import (
...   Body, JointsMethod, PinJoint, PrismaticJoint)
>>> g, l = symbols("g l")
>>> wall = Body("wall")
>>> cart = Body("cart")
>>> pendulum = Body("Pendulum")
>>> slider = PrismaticJoint("s", wall, cart, joint_axis=wall.x)
>>> pin = PinJoint("j", cart, pendulum, joint_axis=cart.z,
...                child_point=l * pendulum.y)
>>> pendulum.masscenter.set_vel(pendulum.frame, 0)
>>> cart.apply_force(-g * cart.mass * wall.y)
>>> pendulum.apply_force(-g * pendulum.mass * wall.y)
>>> method = JointsMethod(wall, slider, pin)  # doctest: +SKIP
>>> method.form_eoms()  # doctest: +SKIP
Matrix([
[ Pendulum_mass*l*u_j(t)**2*sin(q_j(t)) - Pendulum_mass*l*cos(q_j(t))*Derivative(u_j(t), t) - (Pendulum_mass + cart_mass)*Derivative(u_s(t), t)],
[-Pendulum_mass*g*l*sin(q_j(t)) - Pendulum_mass*l*cos(q_j(t))*Derivative(u_s(t), t) - (Pendulum_izz + Pendulum_mass*l**2)*Derivative(u_j(t), t)]])
```

The replacement of ``JointsMethod`` is ``System``, which can be used to form the
equations of motion of the same cart pole as follows:

```py
>>> from sympy import symbols
>>> from sympy.physics.mechanics import (
...   Particle, PinJoint, PrismaticJoint, RigidBody, System)
>>> g, l = symbols("g l")
>>> wall = RigidBody("wall")
>>> cart = RigidBody("cart")
>>> pendulum = RigidBody("Pendulum")
>>> slider = PrismaticJoint("s", wall, cart, joint_axis=wall.x)
>>> pin = PinJoint("j", cart, pendulum, joint_axis=cart.z,
...                child_point=l * pendulum.y)
>>> system = System.from_newtonian(wall)
>>> system.add_joints(slider, pin)
>>> system.apply_uniform_gravity(-g * wall.y)
>>> system.form_eoms()
Matrix([
[ Pendulum_mass*l*u_j(t)**2*sin(q_j(t)) - Pendulum_mass*l*cos(q_j(t))*Derivative(u_j(t), t) - (Pendulum_mass + cart_mass)*Derivative(u_s(t), t)],
[-Pendulum_mass*g*l*sin(q_j(t)) - Pendulum_mass*l*cos(q_j(t))*Derivative(u_s(t), t) - (Pendulum_izz + Pendulum_mass*l**2)*Derivative(u_j(t), t)]])
```

(deprecated-matrix-mixins)=
### Deprecated matrix mixin classes

The matrix mixin classes are deprecated. Previously the ``Matrix`` class (aka
``MutableDenseMatrix``) was created through an inheritance hierarchy that
looked like:
```py
class MatrixRequired:
class MatrixShaping(MatrixRequired):
class MatrixSpecial(MatrixRequired):
class MatrixProperties(MatrixRequired):
class MatrixOperations(MatrixRequired):
class MatrixArithmetic(MatrixRequired):
class MatrixCommon(
    MatrixArithmetic,
    MatrixOperations,
    MatrixProperties,
    MatrixSpecial,
    MatrixShaping):
class MatrixDeterminant(MatrixCommon):
class MatrixReductions(MatrixDeterminant):
class MatrixSubspaces(MatrixReductions):
class MatrixEigen(MatrixSubspaces)
class MatrixCalculus(MatrixCommon):
class MatrixDeprecated(MatrixCommon):
class MatrixBase(MatrixDeprecated,
   MatrixCalculus,
   MatrixEigen,
   MatrixCommon,
   Printable):
class RepMatrix(MatrixBase):
class DenseMatrix(RepMatrix):
class MutableRepMatrix(RepMatrix):
class MutableDenseMatrix(DenseMatrix, MutableRepMatrix):
```
As of SymPy 1.13 this has been simplified and all classes above
``MatrixBase``are merged together so the hierarchy looks like:
```py
class MatrixBase(Printable):
class RepMatrix(MatrixBase):
class DenseMatrix(RepMatrix):
class MutableRepMatrix(RepMatrix):
class MutableDenseMatrix(DenseMatrix, MutableRepMatrix):
```
The matrix mixin classes like ``MatrixRequired`` etc are still available
because downstream code might be subclassing these classes but these are all
deprecated and will be removed in a future version of SymPy. Subclassing these
classes is deprecated and anycode that does that should be changed to not
subclass them.

It is also deprecated to use these classes with ``isinstance`` like
``isinstance(M, MatrixCommon)``. Any code doing this should be changed to use
``isinstance(M, Matrixbase)`` instead which will also work with previous SymPy
versions.

More generally importing anything from the ``sympy.matrices.common`` or
``sympy.matrices.matrices`` modules in which these classes are defined is
deprecated. These modules will be removed in a future release of SymPy.

The reason for this change is that the convoluted inheritance hierarchy made it
difficult to improve ``Matrix`` for the majority of users while still providing
all of these classes that could be subclassed. Since these mixin classes are no
longer used as part of ``Matrix`` they no longer serve any function within
SymPy and the removal of this now unused code will simplify the codebase.

(deprecated-sympify-string-fallback)=
### The string fallback in `sympify()`

The `sympify` function would previously convert an unrecognized object to a
string and retry sympification. This was deprecated in SymPy 1.6 and was
removed in SymPy 1.13.

The behavior of {func}`~.sympify` is that `sympify(expr)` tries various methods
to try to convert `expr` into a SymPy objects. Previously if all these methods
would fail, it would take `str(expr)` and try to parse it using
{func}`~.parse_expr`. This string fallback feature was deprecated in SymPy 1.6
and was removed in SymPy 1.13.

This behaviour was problematic for a few reasons:

- It could affect performance in major ways. See for instance issues
  [#18056](https://github.com/sympy/sympy/issues/18056) and
  [#15416](https://github.com/sympy/sympy/issues/15416) where it caused up to
  100x slowdowns. The issue is that SymPy functions automatically call
  `sympify` on their arguments. Whenever a function is passed something that
  `sympify` doesn't know how to convert to a SymPy object, for instance, a
  Python function type, it passes the string to {func}`~.parse_expr`. This is
  significantly slower than the direct conversions that happen by default.
  This occurs specifically whenever `sympify()` is used in library code
  instead of `_sympify()` (or equivalently `sympify(strict=True)`), but
  presently this is done a lot. Using `strict=True` will at some point be the
  default for all library code, but this is a [harder change to
  make](https://github.com/sympy/sympy/issues/11003).

- It can cause security issues, since strings are evaled, and objects can
  return whatever string they want in their `__repr__`. See also
  https://github.com/sympy/sympy/pull/12524.

- It really wasn't very useful to begin with. Just because an object's string
  form can be parsed into a SymPy expression doesn't mean it should be parsed
  that way. This is usually correct for custom numeric types, but an object's
  repr could be anything. For instance, if the string form of an object looks
  like a valid Python identifier, it would parse as a `Symbol`.

There are plenty of ways to make custom objects work inside of
{func}`~.sympify`.

- Firstly, if an object is intended to work alongside other SymPy expressions,
  it should subclass from {class}`~.Basic` (or {class}`~.Expr`). If it does,
  {func}`~.sympify` will just return it unchanged because it will already be a
  valid SymPy object.

- For objects that you control, you can add the `_sympy_` method. The [sympify
  docstring](sympy.core.sympify.sympify) has an example of this.

- For objects that you don't control, you can add a custom converter to the
  `sympy.core.sympify.converter` dictionary. The {func}`~.sympify` docstring
  also has an example of this.

(dmp-rep)=
### Deprecate the DMP.rep attribute.

The internal type of ``Poly`` is the ``DMP`` class which previously could be
used to access the coefficients of a polynomial as a list like:
```pycon
>>> from sympy import symbols, Poly
>>> x = symbols('x')
>>> p = Poly(x**2 + 2*x + 3)
>>> p
Poly(x**2 + 2*x + 3, x, domain='ZZ')
>>> p.rep  # doctest: +SKIP
DMP([1, 2, 3], ZZ)
>>> p.rep.rep  # doctest: +SKIP
[1, 2, 3]
```

As of SymPy 1.13 the ``DMP`` type may be implemented by one of two subclasses:

- ``DMP_Python`` which is like the previous ``DMP`` type and has a list as its
  internal representation.
- ``DUP_Flint`` which wraps a Flint polynomial from python-flint.

The ``DUP_Flint`` type does not have an attribute that is analogous to the list
that ``DMP_Python`` has. Accessing ``.rep`` will still generate a list but now
gives a deprecation warning.

Instead of ``.rep`` use the ``DMP.to_list()`` method which returns an
equivalent list:
```pycon
>>> p.rep.to_list()
[1, 2, 3]
```

The ``.to_list()`` method is also available in previous versions of SymPy and
its behaviour is unchanged.

(pkgdata)=
### Deprecate the pkgdata module

The ``sympy.utilities.pkdata`` module is deprecated and will be removed. It is
no longer used anywhere in SymPy and is unsuitable for use by any downstream
code. Use the stdlib ``importlib.resources`` module instead.

(eq-rewrite-Add)=
### Deprecate Eq.rewrite(Add)
The ability to rewrite ``eq = Eq(x, y)`` like ``eq.rewrite(Add)`` to give ``x - y``
has been deprecated in favor of writing ``eq.lhs - eq.rhs``. A replacement
property/method was not deemed necessary given the clarity of the explicit
use of ``lhs`` and ``rhs``, and the inclusion of this functionality in the
rewrite apparatus leads to failures when a node expecting a Boolean is re-
written as an Expr.


(deprecated-markers-annotations-fill-rectangles)=
### Deprecate markers, annotations, fill, rectangles of the Plot class
The properties ``markers, annotations, fill, rectangles`` (containing
user-provided numerical data to be added on a plot) are deprecated.
The new implementation saves user-provided numerical data into appropriate
data series, which can easily be processed by ``MatplotlibBackend``.
Instead of setting those properties directly, users should pass the homonym
keyword arguments to the plotting functions.

The supported behavior is to pass keyword arguments to the plotting functions,
which works fine for all versions of SymPy (before and after 1.13):

```py
p = plot(x,
  markers=[{"args":[[0, 1], [0, 1]], "marker": "*", "linestyle": "none"}],
  annotations=[{"text": "test", "xy": (0, 0)}],
  fill={"x": [0, 1, 2, 3], "y1": [0, 1, 2, 3]},
  rectangles=[{"xy": (0, 0), "width": 5, "height": 1}])
```

Setting attributes on the plot object is deprecated and will raise warnings:

```py
p = plot(x, show=False)
p.markers = [{"args":[[0, 1], [0, 1]], "marker": "*", "linestyle": "none"}]
p.annotations = [{"text": "test", "xy": (0, 0)}]
p.fill = {"x": [0, 1, 2, 3], "y1": [0, 1, 2, 3]}
p.rectangles = [{"xy": (0, 0), "width": 5, "height": 1}]
p.show()
```

Motivation for this deprecation: the implementation of the ``Plot`` class
suggests that it is ok to add attributes and hard-coded if-statements in the
``MatplotlibBackend`` class to provide more and more functionalities for
user-provided numerical data (e.g. adding horizontal lines, or vertical
lines, or bar plots, etc). However, in doing so one would reinvent the wheel:
plotting libraries already implements the necessary API. There is no need to
hard code these things. The plotting module should facilitate the visualization
of symbolic expressions. The best way to add custom numerical data is to
retrieve the figure created by the plotting module and use the API of a
particular plotting library. For example:

```py
# plot symbolic expression
p = plot(cos(x))
# retrieve Matplotlib's figure and axes object
fig, ax = p._backend.fig, p._backend.ax[0]
# add the desired numerical data using Matplotlib's API
ax.plot([0, 1, 2], [0, 1, -1], "*")
ax.axhline(0.5)
# visualize the figure
fig
```


(moved-mechanics-functions)=
### Moved mechanics functions
With the introduction of some new objects like the ``Inertia`` and load objects
in the ``sympy.physics.mechanics`` module, some functions from
``sympy.physics.mechanics.functions`` have been moved to new modules. This
removes some circular import errors and makes it easier to navigate through the
source code, due to the parity between function names and module names. The
following functions were moved:
- ``inertia`` has been moved to ``sympy.physics.mechanics.inertia``
- ``inertia_of_point_mass`` has been moved to ``sympy.physics.mechanics.inertia``
- ``gravity`` has been moved to ``sympy.physics.mechanics.loads``

Previously you could import the functions from
``sympy.physics.mechanics.functions``:

```py
>>> from sympy.physics.mechanics.functions import inertia, inertia_of_point_mass, gravity
```

Now they should be imported from ``sympy.physics.mechanics``:

```py
>>> from sympy.physics.mechanics import inertia, inertia_of_point_mass
>>> from sympy.physics.mechanics.loads import gravity
```

(modularinteger-compare)=
### Ordered comparisons like ``a < b`` with modular integers

SymPy's ``GF`` domains represent modular integers. Previously it was possible
to compare these with ordered comparisons like ``a < b``:
```py
>>> from sympy import GF
>>> F5 = GF(5)
>>> F5(2) < F5(3) # doctest: +SKIP
True
```
This will now fail with ``TypeError`` when the ground types are set to
``flint``. When the ground types are not ``flint`` these comparisons are now
deprecated: they will still work but will give a deprecation warning when used.

Ordered comparisons of modular integer or finite fields do not make sense
because these are not ordered fields:
```
>>> e = F5(4)
>>> e + 1 > e # doctest: +SKIP
False
```

(modularinteger-to-int)=
### The ``ModularInteger.to_int()`` method

SymPy's ``GF`` domains are for modular integers e.g. ``GF(n)`` is for the
integers modulo ``n`` and can be used like:
```py
>>> from sympy import GF
>>> K = GF(5)
>>> a = K(7)
>>> a
2 mod 5
```

The elements of a modular integer domain have a ``to_int()`` method that is
deprecated since SymPy 1.13:
```py
>>> # this is deprecated:
>>> a.to_int()  # doctest: +SKIP
2
```

Instead the preferred way to achieve equivalent behavior is to use the method
on the domain (added in SymPy 1.13) or alternatively calling ``int`` might be
better:
```py
>>> K.to_int(a)
2
>>> int(a)
2
```

These two ways of converting to an ``int`` are not equivalent. The domain
``GF(p)`` can be defined with ``symmetric=True`` or ``symmetric=False``. This
difference affects the behavior of the ``to_int`` method:
```py
>>> KS = GF(5, symmetric=True)
>>> KU = GF(5, symmetric=False)
>>> [KS.to_int(KS(n)) for n in range(10)]
[0, 1, 2, -2, -1, 0, 1, 2, -2, -1]
>>> [KU.to_int(KU(n)) for n in range(10)]
[0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
>>> [int(KS(n)) for n in range(10)]
[0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
>>> [int(KU(n)) for n in range(10)]
[0, 1, 2, 3, 4, 0, 1, 2, 3, 4]
```

So if ``symmetric=True`` (which is the default) then the ``to_int`` method will
sometimes return negative integers. If ``symmetric=False`` or if the ``int(a)``
method is used the returned result is always a nonnegative integer. Note also
that the behaviour of ``int(a)`` was changed in SymPy 1.13: in previous
versions it was equivalent to ``a.to_int()``. To write code that behaves the
same way in all SymPy versions you can:

1. Use ``symmetric=False`` and use ``int(a)``.
2. Define a function like
    ```py
    def to_int(K, a):
        if hasattr(K, 'to_int'):
            return K.to_int(a)
        else:
            return a.to_int()
    ```

The reason for this change is that it makes it possible to use python-flint's
``nmod`` as an alternative (much faster) implementation for the elements of
``GF(p)``. It is not possible to add a ``to_int`` method to python-flint's
``nmod`` type or to capture the equivalent of ``symmetric=True/False`` by
storing data in the ``nmod`` instance. Deprecating and removing the ``to_int``
method and changing the behavior of the ``int`` method means that the element
instances do not have any behavior that depends on whether the domain is
considered to be "symmetric" or not. Instead the notion of "symmetric" is now
purely a property of the domain object itself rather than of the elements and
so the ``to_int`` method that depends on this must be a domain method rather
than an element method.

(deprecated-ntheory-symbolic-functions)=
### Relocate symbolic functions from ``ntheory`` to ``functions``

The following symbolic functions in ``ntheory`` have been moved to
``functions``:

* ``sympy.ntheory.factor_.divisor_sigma``
* ``sympy.ntheory.factor_.primenu``
* ``sympy.ntheory.factor_.primeomega``
* ``sympy.ntheory.factor_.reduce_totient``
* ``sympy.ntheory.factor_.totient``
* ``sympy.ntheory.generate.primepi``
* ``sympy.partitions_.npartitions``
* ``sympy.ntheory.residue_ntheory.jacobi_symbol``
* ``sympy.ntheory.residue_ntheory.legendre_symbol``
* ``sympy.ntheory.residue_ntheory.mobius``

Code that imports these functions from top-level like ``from sympy import
mobius`` will continue to work fine. However code that imports these from the
fully qualified module like ``from sympy.ntheory import mobius`` or ``from
sympy.ntheory.residue_ntheory import mobius`` will now see a deprecation
warning. The new location for these functions is in ``sympy.functions`` but the
intended way to import them is still from top-level like ``from sympy import
mobius``.

The following symbolic functions in ``ntheory`` have been moved to
``functions``, but cannot be imported at top-level.

* ``sympy.ntheory.factor_.udivisor_sigma``

The following functions have been moved from ``functions`` to ``ntheory``
because they are numeric functions.

* ``sympy.functions.combinatorial.numbers.carmichael.is_carmichael``
* ``sympy.functions.combinatorial.numbers.carmichael.find_carmichael_numbers_in_range``
* ``sympy.functions.combinatorial.numbers.carmichael.find_first_n_carmichaels``

If you are using these functions, change from

```py
>>> from sympy import carmichael
>>> carmichael.is_carmichael(561)
True
```

to

```py
>>> from sympy import is_carmichael
>>> is_carmichael(561)
True
```


## Version 1.12

(managedproperties)=
### The ``ManagedProperties`` metaclass

The ``ManagedProperties`` metaclass was previously the metaclass for ``Basic``.
Now ``Basic`` does not use metaclasses and so its metaclass is just ``type``.
Any code that previously subclassed ``Basic`` and wanted to do anything with
metaclasses would have needed to subclass ``ManagedProperties`` to make the
relevant metaclass. The only relevant method of ``ManagedProperties`` has been
moved to ``Basic.__init_subclass__``. Since ``ManagedProperties`` is not used
as the metaclass for ``Basic`` any more and no longer does anything useful it
should be possible for such code to just subclass ``type`` instead for any
metaclass.


(deprecated-mechanics-joint-coordinate-format)=
### New Joint coordinate format

The format, i.e. type and auto generated name, of the generalized coordinates
and generalized speeds of the joints in the ``sympy.physics.mechanics`` module
has changed. The data type has changed from ``list`` to ``Matrix``, which is the
same as the type for the generalized coordinates within the ``KanesMethod``.
The auto naming of the generalized coordinates and generalized speeds of the
``PinJoint`` and ``PrismaticJoint`` have also changed to ``q_<joint.name>`` and
``u_<joint.name>``. Previously each of those joints had an unique template for
auto generating these names.

(deprecated-mechanics-joint-axis)=
### New Joint intermediate frames

The definition of the joint axis in the ``sympy.physics.mechanics`` module has
changed. Instead of using the arguments ``parent_axis`` and ``child_axis`` to
automatically determine the joint axis and an intermediate reference frame, the
joints now use an intermediate frame argument for both the parent and the child
body, i.e. ``parent_interframe`` and ``child_interframe``. This means that you
can now fully define the joint attachment, consisting of a point and frame, for
both bodies. Furthermore, if a joint like the ``PinJoint`` has a specific joint
axis, e.g. the axis about which the rotation occurs, then this axis can be
specified using the ``joint_axis`` argument. An advantage of this setup is that
one can more accurately define the transformation from the parent body to the
child body.

For example, suppose you want a ``PinJoint`` that rotates the child body about
the ``parent.z`` axis and ``-child.z`` axis. The previous way to specify this
joint was:

```py
>>> from sympy.physics.mechanics import Body, PinJoint
>>> parent, child = Body('parent'), Body('child')
>>> pin = PinJoint('pin', parent, child, parent_axis=parent.z,
...                child_axis=-child.z)   # doctest: +SKIP
>>> parent.dcm(child)   # doctest: +SKIP
Matrix([
[-cos(q_pin(t)), -sin(q_pin(t)),  0],
[-sin(q_pin(t)),  cos(q_pin(t)),  0],
[             0,              0, -1]])
```

When inspecting this matrix you will notice that for ``theta_pin = 0`` the child
body is rotated $\pi$ rad about the ``parent.y`` axis. In the new definition
you can see that we get the same result, but this time we have also specified
this exact rotation:

```py
>>> from sympy import pi
>>> from sympy.physics.mechanics import Body, PinJoint, ReferenceFrame
>>> parent, child, = Body('parent'), Body('child')
>>> int_frame = ReferenceFrame('int_frame')
>>> int_frame.orient_axis(child.frame, child.y, pi)
>>> pin = PinJoint('pin', parent, child, joint_axis=parent.z,
...                child_interframe=int_frame)
>>> parent.dcm(child)
Matrix([
[-cos(q_pin(t)), -sin(q_pin(t)),  0],
[-sin(q_pin(t)),  cos(q_pin(t)),  0],
[             0,              0, -1]])
```

However if you liked the fact that the deprecated arguments aligned the frames
for you, then you can still make use of this feature by providing vectors to
``parent_interframe`` and ``child_interframe``, which are then oriented such
that the joint axis expressed in the intermediate frame is aligned with the
given vector:

```py
>>> from sympy.physics.mechanics import Body, PinJoint
>>> parent, child = Body('parent'), Body('child')
>>> pin = PinJoint('pin', parent, child, parent_interframe=parent.z,
...                child_interframe=-child.z)
>>> parent.dcm(child)
Matrix([
[-cos(q_pin(t)), -sin(q_pin(t)),  0],
[-sin(q_pin(t)),  cos(q_pin(t)),  0],
[             0,              0, -1]])
```

(deprecated-mechanics-joint-pos)=
### Change in joint attachment point argument

The argument names for specifying the attachment points of a joint in
``sympy.physics.mechanics`` , i.e. ``parent_joint_pos`` and ``child_joint_pos``,
have been changed to ``parent_point`` and ``child_point``. This is because these
arguments can now also be ``Point`` objects, so they can be exactly the same as
the ``parent_point`` and ``child_point`` attributes.

For example, suppose you want a ``PinJoint`` in the parent to be positioned at
``parent.frame.x`` with respect to the mass center, and in the child at
``-child.frame.x``. The previous way to specify this was:

```py
>>> from sympy.physics.mechanics import Body, PinJoint
>>> parent, child = Body('parent'), Body('child')
>>> pin = PinJoint('pin', parent, child, parent_joint_pos=parent.frame.x,
...                child_joint_pos=-child.frame.x)   # doctest: +SKIP
>>> pin.parent_point.pos_from(parent.masscenter)   # doctest: +SKIP
parent_frame.x
>>> pin.child_point.pos_from(child.masscenter)   # doctest: +SKIP
- child_frame.x
```

Now you can do the same with either

```py
>>> from sympy.physics.mechanics import Body, PinJoint
>>> parent, child = Body('parent'), Body('child')
>>> pin = PinJoint('pin', parent, child, parent_point=parent.frame.x,
...                child_point=-child.frame.x)
>>> pin.parent_point.pos_from(parent.masscenter)
parent_frame.x
>>> pin.child_point.pos_from(child.masscenter)
- child_frame.x
```

Or

```py
>>> from sympy.physics.mechanics import Body, PinJoint, Point
>>> parent, child = Body('parent'), Body('child')
>>> parent_point = parent.masscenter.locatenew('parent_point', parent.frame.x)
>>> child_point = child.masscenter.locatenew('child_point', -child.frame.x)
>>> pin = PinJoint('pin', parent, child, parent_point=parent_point,
...                child_point=child_point)
>>> pin.parent_point.pos_from(parent.masscenter)
parent_frame.x
>>> pin.child_point.pos_from(child.masscenter)
- child_frame.x
```

## Version 1.11

(deprecated-conv-array-expr-module-names)=
### Modules `sympy.tensor.array.expressions.conv_*` renamed to `sympy.tensor.array.expressions.from_*`

In order to avoid possible naming and tab-completion conflicts with
functions with similar names to the names of the modules, all modules whose
name starts with `conv_*` in `sympy.tensor.array.expressions` have been renamed
to `from_*`.

(mathematica-parser-new)=
### New Mathematica code parser

The old mathematica code parser defined in the module ``sympy.parsing.mathematica``
in the function ``mathematica`` is deprecated. The function ``parse_mathematica``
with a new and more comprehensive parser should be used instead.

The ``additional_translations`` parameter for the Mathematica parser is not available
in ``parse_mathematica``.
Additional translation rules to convert Mathematica expressions into SymPy ones
should be specified after the conversion using SymPy's ``.replace( )`` or ``.subs( )``
methods on the output expression. If the translator fails to recognize the logical
meaning of a Mathematica expression, a form similar to Mathematica's full form
will be returned, using SymPy's ``Function`` object to encode the nodes of the
syntax tree.

For example, suppose you want ``F`` to be a function that returns the maximum
value multiplied by the minimum value, the previous way to
specify this conversion was:

```py
>>> from sympy.parsing.mathematica import mathematica
>>> mathematica('F[7,5,3]', {'F[*x]': 'Max(*x)*Min(*x)'})   # doctest: +SKIP
21
```

Now you can do the same with

```py
>>> from sympy.parsing.mathematica import parse_mathematica
>>> from sympy import Function, Max, Min
>>> parse_mathematica("F[7,5,3]").replace(Function("F"), lambda *x: Max(*x)*Min(*x))
21
```

(deprecated-carmichael-static-methods)=
### Redundant static methods in `carmichael`

A number of static methods in `~.carmichael` are just wrappers around other
functions. Instead of ``carmichael.is_perfect_square`` use
`sympy.ntheory.primetest.is_square` and instead of ``carmichael.is_prime`` use
`~.isprime`. Finally, ``carmichael.divides`` can be replaced by instead checking

```py
n % p == 0
```

(remove-check-argument-from-matrix-operations)=
### The `check` argument to `HadamardProduct`, `MatAdd` and `MatMul`

This argument can be used to pass incorrect values to `~.HadamardProduct`,
`~.MatAdd`, and `~.MatMul` leading to later problems. The `check` argument
will be removed and the arguments will always be checked for correctness, i.e.,
the arguments are matrices or matrix symbols.

## Version 1.10

(deprecated-traversal-functions-moved)=
### Some traversal functions have been moved

Some traversal functions have moved. Specifically, the functions

- `bottom_up`
- `interactive_traversal`
- `postorder_traversal`
- `preorder_traversal`
- `use`

have moved to different SymPy submodules.

These functions should be used from the top-level `sympy` namespace, like

```py
sympy.preorder_traversal
```

or

```py
from sympy import preorder_traversal
```

In general, end-users should use the top-level `sympy` namespace for any
functions present there. If a name is in the top-level namespace, its specific
SymPy submodule should not be relied on, as functions may move around due to
internal refactorings.

(sympy-core-trace-deprecated)=
### `sympy.core.trace`

The trace object `sympy.core.trace.Tr()` was moved to
`sympy.physics.quantum.trace.Tr()`.  This was because it was only used in the
`sympy.physics.quantum` submodule, so it was better to have it there than in
the core.

(deprecated-sympy-core-compatibility)=
### The `sympy.core.compatibility` submodule

The `sympy.core.compatibility` submodule is deprecated.

This submodule was only ever intended for internal use. Now that SymPy no
longer supports Python 2, this module is no longer necessary, and the
remaining helper functions have been moved to more convenient places in the
SymPy codebase.

Some of the functions that were in this module are available from the
top-level SymPy namespace, i.e.,

```
sympy.ordered
sympy.default_sort_key
```

or

```py
from sympy import ordered, default_sort_key
```

In general, end-users should use the top-level `sympy` namespace for any
functions present there. If a name is in the top-level namespace, its specific
SymPy submodule should not be relied on, as functions may move around due to
internal refactorings.

The remaining functions in `sympy.core.compatibility` were only intended for
internal SymPy use and should not be used by user code.

Additionally, these two functions, `ordered` and `default_sort_key`, also used
to be in `sympy.utilities.iterables` but have been moved from there as well.

## Version 1.9

(deprecated-expr-free-symbols)=
### `expr_free_symbols`

The `expr_free_symbols` attribute of various SymPy objects is deprecated.

`expr_free_symbols` was meant to represent indexed objects such as
`MatrixElement` and {class}`~.Indexed` as free symbols. This was
intended to make derivatives of free symbols work. However, this now works
without making use of the method:

```py
>>> from sympy import Indexed, MatrixSymbol, diff
>>> a = Indexed("A", 0)
>>> diff(a**2, a)
2*A[0]
>>> X = MatrixSymbol("X", 3, 3)
>>> diff(X[0, 0]**2, X[0, 0])
2*X[0, 0]
```

This was a general property that was added to solve a very specific problem
but it added a layer of abstraction that is not necessary in general.

1. objects that have structural "non-expression" nodes already allow one to
   focus on the expression node if desired, e.g.

   ```python
   >>> from sympy import Derivative, symbols, Function
   >>> x = symbols('x')
   >>> f = Function('f')
   >>> Derivative(f(x), x).expr
   f(x)
   ```

   introduction of this property encourages imprecise thinking when requesting
   free_symbols since it allows one to get symbols from a specific node of an
   object without specifying the node

2. the property was incorrectly added to `AtomicExpr` so numbers are returned
   as `expr_free_symbols`:

   ```python
   >>> S(2).expr_free_symbols # doctest: +SKIP
   2
   ```

3. the application of the concept was misapplied to define
   `Subs.expr_free_symbols`: it added in `expr_free_symbols` of the point but
   the point is a `Tuple` so nothing was added

4. it was not used anywhere else in the codebase except in the context of
   differentiating a `Subs` object, which suggested that it was not something
   of general use, this is also confirmed by the fact that,

5. it was added without specific tests except for test of the derivatives of
   the Subs object for which it was introduced

See issue [#21494](https://github.com/sympy/sympy/issues/21494) for more
discussion.

(deprecated-sympy-stats-numsamples)=
### `sympy.stats.sample(numsamples=n)`

The `numsamples` parameter to {func}`sympy.stats.sample` is deprecated.

`numsamples` makes `sample()` return a list of size `numsamples`, like

```py
>>> from sympy.stats import Die, sample
>>> X = Die('X', 6)
>>> sample(X, numsamples=3) # doctest: +SKIP
[3, 2, 3]
```

However, this functionality can be easily implemented by the user with a list
comprehension

```py
>>> [sample(X) for i in range(3)] # doctest: +SKIP
[5, 4, 3]
```

Additionally, it is redundant with the `size` parameter, which makes `sample`
return a NumPy array with the given shape.

```py
>>> sample(X, size=(3,)) # doctest: +SKIP
array([6, 6, 1])
```

Historically, `sample` was changed in SymPy 1.7 so it returned an iterator
instead of sample value. Since an iterator was returned, a numsamples
parameter was added to specify the length of the iterator.

However, this new behavior was considered confusing, as discussed in issue
[#21563](https://github.com/sympy/sympy/issues/21563), so it was reverted.
Now, `sample_iter` should be used if a iterator is needed. Consequently, the
`numsamples` parameter is no longer needed for `sample()`.

(deprecated-rawmatrix)=
### `sympy.polys.solvers.RawMatrix`

The `RawMatrix` class is deprecated. The `RawMatrix` class was a subclass
of `Matrix` that used domain elements instead of `Expr` as the elements of
the matrix. This breaks a key internal invariant of `Matrix` and this kind
of subclassing limits improvements to the `Matrix` class.

The only part of SymPy that documented the use of the `RawMatrix` class was
the Smith normal form code, and that has now been changed to use
`DomainMatrix` instead. It is recommended that anyone using `RawMatrix` with
the previous Smith Normal Form code should switch to using `DomainMatrix` as
shown in issue [#21402](https://github.com/sympy/sympy/pull/21402). A better
API for the Smith normal form will be added later.

(deprecated-non-expr-in-matrix)=
### Non-`Expr` objects in a Matrix

In SymPy 1.8 and earlier versions it was possible to put non-{class}`~.Expr`
elements in a [`Matrix`](sympy.matrices.dense.Matrix) and the matrix elements could be any arbitrary
Python object:

```python
>>> M = Matrix([[(1, 2), {}]]) # doctest: +SKIP
```

This is not useful and does not really work, e.g.:

```python
>>> M + M # doctest: +SKIP
Traceback (most recent call last):
...
TypeError: unsupported operand type(s) for +: 'Dict' and 'Dict'
```

The main reason for making this possible was that there were a number of
`Matrix` subclasses in the SymPy codebase that wanted to work with objects
from the polys module, e.g.

1. `RawMatrix` (see [above](deprecated-rawmatrix)) was used in `solve_lin_sys`
   which was part of `heurisch` and was also used by `smith_normal_form`. The
   `NewMatrix` class used domain elements as the elements of the Matrix rather
   than `Expr`.

2. `NewMatrix` was used in the `holonomic` module and also used domain
   elements as matrix elements

3. `PolyMatrix` used a mix of `Poly` and `Expr` as the matrix elements and was
   used by `risch`.

All of these matrix subclasses were broken in different ways and the
introduction of {class}`~.DomainMatrix`
([#20780](https://github.com/sympy/sympy/pull/20780),
[#20759](https://github.com/sympy/sympy/pull/20759),
[#20621](https://github.com/sympy/sympy/pull/20621),
[#19882](https://github.com/sympy/sympy/pull/19882),
[#18844](https://github.com/sympy/sympy/pull/18844)) provides a better
solution for all cases. Previous PRs have removed the dependence of these
other use cases on Matrix
([#21441](https://github.com/sympy/sympy/pull/21441),
[#21427](https://github.com/sympy/sympy/pull/21427),
[#21402](https://github.com/sympy/sympy/pull/21402)) and now
[#21496](https://github.com/sympy/sympy/pull/21496) has deprecated having
non-`Expr` in a `Matrix`.

This change makes it possible to improve the internals of the Matrix class but
it potentially impacts on some downstream use cases that might be similar to
the uses of `Matrix` with non-`Expr` elements that were in the SymPy codebase.
A potential replacement for code that used `Matrix` with non-`Expr` elements
is {class}`~.DomainMatrix` if the elements are something like domain elements
and a domain object can be provided for them. Alternatively if the goal is
just printing support then perhaps `TableForm` can be used.

It isn't clear what to advise as a replacement here without knowing more about
the usecase. If you are unclear how to update your code, please [open an
issue](https://github.com/sympy/sympy/issues/new) or [write to our mailing
list](https://groups.google.com/g/sympy) so we can discuss it.

(deprecated-get-segments)=
### The `get_segments` attribute of plotting objects

The `get_segments` method implemented in {class}`~.Line2DBaseSeries` is used
to convert two list of coordinates, `x` and `y`, into a list of segments used
by Matplotlib's `LineCollection` to plot a line.

Since the list of segments is only required by Matplotlib (for example, Bokeh,
Plotly, Mayavi, K3D only require lists of coordinates), this has been moved
inside the `MatplotlibBackend` class.

Note that previously, the method
{meth}`~sympy.plotting.series.LineOver1DRangeSeries.get_points` always returned
uniformly sampled points, which meant that some functions were not plotted
correctly when using `get_points()` to plot with Matplotlib.

To avoid this problem, the method `get_segments()` could be used, which used
adaptive sampling and which could be used with Matplotlib's `LineCollection`.
However, this has been changed, and now `get_points()` can also use adaptive
sampling. The {meth}`~sympy.plotting.series.Line2DBaseSeries.get_data()` method
can also be used.


(deprecated-physics-mdft)=
### The `mdft` function in `sympy.physics.matrices`

The `sympy.physics.matrices.mdft()` function is deprecated. It can be replaced
with the `DFT` class in `sympy.matrices.expressions.fourier`.

In particular, replace `mdft(n)` with `DFT(n).as_explicit()`. For example:

```py
>>> from sympy.physics.matrices import mdft
>>> mdft(3) # DEPRECATED # doctest: +SKIP
Matrix([
[sqrt(3)/3,                sqrt(3)/3,                sqrt(3)/3],
[sqrt(3)/3, sqrt(3)*exp(-2*I*pi/3)/3,  sqrt(3)*exp(2*I*pi/3)/3],
[sqrt(3)/3,  sqrt(3)*exp(2*I*pi/3)/3, sqrt(3)*exp(-2*I*pi/3)/3]])
```

```py
>>> from sympy.matrices.expressions.fourier import DFT
>>> DFT(3)
DFT(3)
>>> DFT(3).as_explicit()
Matrix([
[sqrt(3)/3,                sqrt(3)/3,                sqrt(3)/3],
[sqrt(3)/3, sqrt(3)*exp(-2*I*pi/3)/3,  sqrt(3)*exp(2*I*pi/3)/3],
[sqrt(3)/3,  sqrt(3)*exp(2*I*pi/3)/3, sqrt(3)*exp(-2*I*pi/3)/3]])
```

This was changed because the `sympy.physics` submodule is supposed to only
contain things that are specific to physics, but the discrete Fourier
transform matrix is a more general mathematical concept, so it is better
located in the `sympy.matrices` module. Furthermore, the `DFT` class is a
[matrix expression](sympy.matrices.expressions), meaning it can be unevaluated
and support symbolic shape.

(deprecated-private-matrix-attributes)=
### The private `SparseMatrix._smat` and `DenseMatrix._mat` attributes

The `._mat` attribute of [`Matrix`](sympy.matrices.dense.Matrix) and the
`._smat` attribute of [`SparseMatrix`](sympy.matrices.sparse.SparseMatrix) are
deprecated.

The internal representation of `Matrix` and `SparseMatrix` was changed to be a
{class}`~.DomainMatrix` in [#21626](https://github.com/sympy/sympy/pull/21626)
so that it is no longer possible to expose a mutable list/dict as a way of
mutating a `Matrix`. Instead of `._mat` the new `.flat()` method can be used,
which returns a new list that cannot be used to mutate the `Matrix` itself.
Instead of `._smat` the `.todok()` method can be used which returns a new
dict.

Note that these attributes are already changed in SymPy 1.9 to return
read-only copies, so that any code that relied on mutating them will be
broken. Also these attributes were technically always private (they started
with an underscore), so user code should not really have been using them in
the first place.

(deprecated-laplace-transform-matrix)=
### laplace_transform of a Matrix with noconds=False

Prior to version 1.9, calling {func}`~.laplace_transform` on a [`Matrix`](sympy.matrices.dense.Matrix) with
`noconds=False` (which is the default), resulted in a Matrix of tuples:

```py
>>> from sympy import laplace_transform, symbols, eye
>>> t, z = symbols('t z')
>>> laplace_transform(eye(2), t, z) # doctest: +SKIP
Matrix([
[(1/z, 0, True),   (0, 0, True)],
[  (0, 0, True), (1/z, 0, True)]])
```

However, `Matrix` is only designed to work with `Expr` objects (see
{ref}`deprecated-non-expr-in-matrix` above).

To avoid this, either use `noconds=True` to remove the convergence conditions

```py
>>> laplace_transform(eye(2), t, z, noconds=True)
Matrix([
[1/z,   0],
[  0, 1/z]])
```

or use `legacy_matrix=False` to return the new behavior, which will be to
return a single tuple with the Matrix in the first argument and the
convergence conditions combined into a single condition for the whole matrix.

```
>>> laplace_transform(eye(2), t, z, legacy_matrix=False)
(Matrix([
[1/z,   0],
[  0, 1/z]]), 0, True)
```

When this deprecation is removed the `legacy_matrix=False` behavior will
become the default, but the flag will be left intact for compatibility.

## Version 1.8

(theanocode-deprecated)=
### `sympy.printing.theanocode`

[Theano](https://github.com/Theano/Theano) has been discontinued, and forked
into a new project called [Aesara](https://github.com/aesara-devs/aesara). The
`sympy.printing.theanocode` module has been renamed to
{mod}`sympy.printing.aesaracode`, and all the corresponding functions have
been renamed (e.g., `theano_code` has been renamed to {func}`~.aesara_code`,
`TheanoPrinter` has been renamed to {class}`~.AesaraPrinter`, and so on).

(deprecated-askhandler)=
### `sympy.assumptions.handlers.AskHandler` and related methods

`Predicate` has experienced a big design change. Previously, its handler was a
list of `AskHandler` classes and registration was done by `add_handler()` and
`remove_handler()` functions. Now, its handler is a multipledispatch instance
and registration is done by `register()` or `register_many()` methods. Users
must define a predicate class to introduce a new one.

Previously, handlers were defined and registered this way:

```python
class AskPrimeHandler(AskHandler):
    @staticmethod
    def Integer(expr, assumptions):
        return expr.is_prime

register_handler('prime', AskPrimeHandler)
```

It should be changed to this:

```python
# Predicate definition.
# Not needed if you are registering the handler to existing predicate.
class PrimePredicate(Predicate):
    name = 'prime'
Q.prime = PrimePredicate()

# Handler registration
@Q.prime.register(Integer)
def _(expr, assumptions):
    return expr.is_prime
```

See GitHub issue [#20209](https://github.com/sympy/sympy/issues/20209).

## Version 1.7.1

(deprecated-distribution-randomindexedsymbol)=
### Calling `sympy.stats.StochasticProcess.distribution` with `RandomIndexedSymbol`

The `distribution` method of `sympy.stats` [stochastic
processes](sympy-stats-stochastic-processes) used to accept a
`RandomIndexedSymbol` (that is, a stochastic process indexed with a
timestamp), but should now only be called with the timestamp.

For example, if you have

```py
>>> from sympy import symbols
>>> from sympy.stats import WienerProcess
>>> W = WienerProcess('W')
>>> t = symbols('t', positive=True)
```

Previously this would work

```py
W.distribution(W(t)) # DEPRECATED
```

It should now be called like

```py
>>> W.distribution(t)
NormalDistribution(0, sqrt(t))
```

This was change was made as part of a change to store only `Basic` objects in
`sympy.stats` `.args`. See issue
[#20078](https://github.com/sympy/sympy/issues/20078) for details.

## Version 1.7

(deprecated-absorbing_probabilites)=
### `sympy.stats.DiscreteMarkovChain.absorbing_probabilites()`

The `absorbing_probabilites` method name was misspelled. The correct spelling
{meth}`~.absorbing_probabilities` ("absorbing probabilit*i*es") should be used
instead.

(deprecated-find-executable)=
### `sympy.utilities.misc.find_executable()`

The function `sympy.utilities.misc.find_executable()` is deprecated. Instead
use the standard library
[`shutil.which()`](https://docs.python.org/3/library/shutil.html#shutil.which)
function, which has been in the standard library since Python 3.3 and is more
powerful.

(deprecated-diffgeom-mutable)=
### Mutable attributes in `sympy.diffgeom`

Several parts of {mod}`sympy.diffgeom` have been updated to no longer be
mutable, which better matches the immutable design used in the rest of SymPy.

- Passing strings for symbol names in {class}`~.CoordSystem` is deprecated.
  Instead you should be explicit and pass symbols with the appropriate
  assumptions, for instance, instead of

  ```py
  CoordSystem(name, patch, ['x', 'y']) # DEPRECATED
  ```

  use

  ```py
  CoordSystem(name, patch, symbols('x y', real=True))
  ```

- Similarly, the `names` keyword argument has been renamed to `symbols`, which
  should be a list of symbols.

- The `Manifold.patches` attribute is deprecated. Patches should be tracked
  separately.

- The `Patch.coord_systems` attribute is deprecated. Coordinate systems should
  be tracked separately.

- The `CoordSystem.transforms` attribute, `CoordSystem.connect_to()` method,
  and `CoordSystem.coord_tuple_transform_to()` method are deprecated. Instead,
  use the `relations` keyword to the `CoordSystem` class constructor and the
  {meth}`.CoordSystem.transformation()` and {meth}`.CoordSystem.transform()`
  methods (see the docstring of {class}`~.CoordSystem` for examples).

(deprecated-pretty-printing-functions)=
### The `unicode` argument and attribute to `sympy.printing.pretty.stringpict.prettyForm` and the `sympy.printing.pretty.pretty_symbology.xstr` function

The `sympy.printing.pretty.pretty_symbology.xstr` function, and the `unicode`
argument and attribute to {class}`sympy.printing.pretty.stringpict.prettyForm`
were both present to support the Unicode behavior of Python 2. Since Unicode
strings are the default in Python 3, these are not needed any more. `xstr()`
should be replaced with just `str()`, the `unicode` argument to `prettyForm`
should be omitted, and the `prettyForm.unicode` attribute should be replaced
with the `prettyForm.s` attribute.

(deprecated-lambdify-arguments-set)=
### Passing the arguments to `lambdify` as a `set`

Passing the function arguments to lambdify as a set is deprecated. Instead
pass them as a list or tuple. For
example, instead of

```py
lambdify({x, y}, x + 2*y) # WRONG
```

use

```py
lambdify((x, y), x + 2*y) # RIGHT
```

This is because sets are unordered. For instance, in the above example it
would be impossible for `lambidfy` to know if it was called with `{x, y}` or
`{y, x}`. Thus, when passed the arguments as a set `lambdify` would have to
guess their order, which would lead to an incorrect function if it guessed
incorrectly.

(non-expr-args-deprecated)=
### Core operators no longer accept non-Expr args

The core operator classes {class}`~.Add`, {class}`~.Mul`, and {class}`~.Pow`
can no longer be constructed directly with objects that are not subclasses of
{class}`~.Expr`.

{class}`~.Expr` is the superclass of all SymPy classes that represent scalar
numeric quantities. For example, {class}`~.sin`, {class}`~.Symbol`, and
{class}`~.Add` are all subclasses of {class}`~.Expr`. However, may objects in
SymPy are not {class}`~.Expr` because they represent some other type of
mathematical object. For example, {class}`~.Set`, {class}`~.Poly`, and
{class}`~.Boolean` are all non-`Expr`. These do not make mathematical sense
inside of `Add`, `Mul`, and `Pow`, which are designed specifically to
represent the addition, multiplication, and exponentiation of scalar complex
numbers.

Manually constructing one of these classes with such an object is possible,
but it will generally create something that will then break. For example

```py
Mul(1, Tuple(2)) # This is deprecated
```

works and creates `Tuple(2)`, but only because `Mul` is "tricked" by always
treating $1 \cdot x = x$. If instead you try

```py
Mul(2, Tuple(2)) # This is deprecated
```

it fails with an exception

```pytb
AttributeError: 'Tuple' object has no attribute 'as_coeff_Mul'
```

because it tries to call a method of `Expr` on the `Tuple` object, which does
not have all the `Expr` methods (because it is not a subclass of `Expr`).

If you want to use the `+`, `*`, or `**` operation on a non-`Expr` object, use
the operator directly rather than using `Mul`, `Add` or `Pow`. If functional
versions of these are desired, you can use a `lambda` or the
[`operator`](https://docs.python.org/3/library/operator.html) module.

## Version 1.6

(deprecated-sympy-utilities-submodules)=
### Various `sympy.utilities` submodules have moved

The following submodules have been renamed.

- `sympy.utilities.benchmarking` → `sympy.testing.benchmarking`
- `sympy.utilities.pytest` → `sympy.testing.pytest`
- `sympy.utilities.randtests` → `sympy.core.random`
- `sympy.utilities.runtests` → `sympy.testing.runtests`
- `sympy.utilities.tmpfiles` → `sympy.testing.tmpfiles`

(deprecated-sympy-testing-randtest)=
### `sympy.testing.randtest`

`sympy.testing.randtest` is deprecated. The functions in it have been moved to
`sympy.core.random`. The following functions have been moved.

- `sympy.testing.randtest.random_complex_number` → `sympy.core.random.random_complex_number`
- `sympy.testing.randtest.verify_numerically` `sympy.core.random.verify_numerically`
- `sympy.testing.randtest.test_derivative_numerically` → `sympy.core.random.test_derivative_numerically`
- `sympy.testing.randtest._randrange` → `sympy.core.random._randrange`
- `sympy.testing.randtest._randint` → `sympy.core.random._randint`

(deprecated-poly-nonpoly-binary-operations)=
### Mixing `Poly` and non-polynomial expressions in binary operations

In previous versions of SymPy, {class}`~.Poly` was a subclass of
{class}`~.Expr`, but it has been changed to only be a subclass of
{class}`~.Basic`. This means that some things that used to work with `Poly`
are now deprecated because they are only designed to work with {class}`~.Expr`
objects.

This includes combining `Poly` with `Expr` objects using binary operations,
for example

```py
Poly(x)*sin(x) # DEPRECATED
```

To do this, either explicitly convert the non-`Poly` operand to a `Poly` using
{meth}`.Expr.as_poly` or convert the `Poly` operand to an {class}`~.Expr`
using {meth}`.Poly.as_expr`, depending on which type you want the result to
be.

(deprecated-permutation-print_cyclic)=
### The `print_cyclic` flag of `sympy.combinatorics.Permutation`

The `print_cyclic` attribute of
[`sympy.combintorics.Permutation`](sympy.combinatorics.permutations.Permutation)
controls whether permutations print as cycles or arrays. This would be done by
setting `Permutation.print_cyclic = True` or `Permutation.print_cyclic =
False`. However, this method of controlling printing is bad because it is a
global flag, but printing should not depend on global behavior.

Instead, users should use the `perm_cyclic` flag of the corresponding printer.
The easiest way to configure this is to set the flag when calling
{func}`~.init_printing`, like

<!-- doctests are skipped here so that it doesn't make the rest of the file -->
<!-- use unicode pretty printing -->

```py
>>> from sympy import init_printing
>>> init_printing(perm_cyclic=False) # Makes Permutation print in array form # doctest: +SKIP
>>> from sympy.combinatorics import Permutation
>>> Permutation(1, 2)(3, 4) # doctest: +SKIP
⎛0 1 2 3 4⎞
⎝0 2 1 4 3⎠
```

The {class}`~.Permutation` docstring contains more details on the
`perm_cyclic` flag.

(deprecated-integrate-poly)=
### Using `integrate` with `Poly`

In previous versions of SymPy, {class}`~.Poly` was a subclass of
{class}`~.Expr`, but it has been changed to only be a subclass of
{class}`~.Basic`. This means that some things that used to work with `Poly`
are now deprecated because they are only designed to work with {class}`~.Expr`
objects.

This includes calling {func}`~.integrate` or {class}`~.Integral` with `Poly`.

To integrate a `Poly`, use the {meth}`.Poly.integrate` method. To compute the
integral as an {class}`~.Expr` object, call the {meth}`.Poly.as_expr` method
first.

See also {ref}`deprecated-poly-nonpoly-binary-operations` above.

(deprecated-indefinite-integral-eq)=
### Creating an indefinite `Integral` with an `Eq` argument

Passing an [`Eq()`](sympy.core.relational.Equality) object to
{func}`~.integrate` is deprecated in the case where the integral is
indefinite. This is because if $f(x) = g(x)$, then $\int f(x)\,dx = \int
g(x)\,dx$ is not true in general, due to the arbitrary constants (which
`integrate` does not include).

If you want to make an equality of indefinite integrals, use
`Eq(integrate(f(x), x), integrate(g(x), x))` explicitly.

If you already have an equality object `eq`, you can use `Eq(integrate(eq.lhs,
x), integrate(eq.rhs, x))`.

## Version 1.5

(deprecated-tensor-fun-eval)=
### `Tensor.fun_eval` and `Tensor.__call__`

`TensExpr.fun_eval` and `Tensor.__call__` (i.e., calling a tensor to evaluate
it) are deprecated. The `Tensor.substitute_indices()` method should be used.
This was changed because `fun_eval` was considered a confusing name and using
function evaluation was considered both confusing and dangerous.

(deprecated-tensortype)=
### `TensorType`

The `TensorType` class is deprecated. Use {func}`~.tensor_heads` instead. The
`TensorType` class had no purpose except shorter creation of
{class}`~.TensorHead` objects.

See also {ref}`deprecated-tensorhead` below.

(deprecated-tensorindextype-dummy-fmt)=
### The `dummy_fmt` argument to `TensorIndexType`

The `dummy_fmt` keyword argument to {class}`~.TensorIndexType` is deprecated.
Setting `dummy_fmt='L'` leads to `_dummy_fmt='L_%d'`, which is confusing and
uses obsolete string formatting. `dummy_name` should be used instead. This
change was made because `dummy_name` is a clearer name.

(deprecated-tensorindextype-metric)=
### The `metric` argument to `TensorIndexType`

The `metric` keyword argument to {class}`~.TensorIndexType` is deprecated.
The name "metric" was ambiguous because it meant "metric symmetry" in some
places and "metric tensor" in others.

Either the `metric_symmetry` keyword or the `TensorIndexType.set_metric()`
method should be used instead.

(deprecated-tensorindextype-methods)=
### The `get_kronecker_delta()` and `get_epsilon()` methods of `TensorIndexType`

The `get_kronecker_delta()` and `get_epsilon()` methods of
{class}`~.TensorIndexType` are deprecated. Use the `TensorIndexType.delta` and
`TensorIndexType.epsilon` properties instead, respectively.

(deprecated-tensorsymmetry)=
### The `tensorsymmetry()` function

The `tensorsymmetry()` function in `sympy.tensor` is deprecated. Use the
{class}`~.TensorSymmetry` class constructor instead.

`TensorSymmetry` is preferred over `tensorsymmetry()` because the latter

1. Does not have any extra functionality
2. Involves obscure Young tableau
3. Is not a member of the `TensorSymmetry` class

(deprecated-tensorhead)=
### The `tensorhead()` function

The `tensorhead()` function is deprecated in favor of {func}`~.tensor_heads`.
`tensor_heads()` is more consistent with other SymPy names (i.e., `Symbol` and
`symbols()` or `TensorIndex` and `tensor_indices()`). It also does not use
Young tableau to denote symmetries.

(deprecated-is-emptyset)=
### The `is_EmptySet` attribute of sets

The `is_EmptySet` attribute of [Set](sets-module) objects is deprecated.
Instead either use

```
from sympy import S
s is S.EmptySet
```

or

```
s.is_empty
```

The difference is that `s.is_empty` may return `None` if it is unknown if the
set is empty.

(deprecated-productset-iterable)=
### `ProductSet(iterable)`

Passing a single iterable as the first argument to {class}`~.ProductSet` is
deprecated. Creating a product set from an iterable should be done using
`ProductSet(*iterable)`, or as each individual argument. For example

```py
>>> from sympy import ProductSet
>>> sets = [{i} for i in range(3)]
>>> ProductSet(*sets)
ProductSet({0}, {1}, {2})
>>> ProductSet({1, 2}, {1})
ProductSet({1, 2}, {1})
```

This is done because sets themselves can be iterables, and sets of sets are
allowed. But the product set of a single set should mathematically be that set
itself (or more exactly, the set of 1-tuples of elements of that set).
Automatically denesting a single iterable makes it impossible to represent
this object and makes `ProductSet` not generalize correctly when passed 1
argument. On the other hand, treating the first argument differently if it is
a set than if it is another type of iterable (which is what is currently done
in the deprecated code path) is confusing behavior.

(deprecated-set-potential-energy)=
### The `set_potential_energy` method in `sympy.physics.mechanics`

The `set_potential_energy()` methods of {class}`sympy.physics.mechanics.particle.Particle`
and {class}`sympy.physics.mechanics.rigidbody.RigidBody` are deprecated.

Instead one should set the {attr}`.Particle.potential_energy` and
{attr}`.RigidBody.potential_energy` attributes to set the potential energy,
like

```py
P.potential_energy = scalar
```

This change was made to be more Pythonic, by using setters and getters of a
`@property` method rather than an explicit `set_` method.

(deprecated-conditionset-set)=
### Using a set for the condition in `ConditionSet`

Using a set for the condition in ConditionSet is deprecated. A boolean should
be used instead. This is because the condition is mathematically a boolean,
and it is ambiguous what a set should mean in this context.

To fix this deprecation, replace

```py
ConditionSet(symbol, set_condition)
```

with

```py
ConditionSet(symbol, And(*[Eq(lhs, 0) for lhs in set_condition]))
```

For example,

```py
ConditionSet((x, y), {x + 1, x + y}, S.Reals) # DEPRECATED
```

would become

```py
ConditionSet((x, y), Eq(x + 1, 0) & Eq(x + y, 0), S.Reals)
```

(deprecated-dixonresultant-properties)=
### The `max_degree` and `get_upper_degree` properties of `sympy.polys.multivariate_resultants.DixonResultant`

The `max_degree` property and `get_upper_degree()` methods of `DixonResultant`
are deprecated. See issue [#17749](https://github.com/sympy/sympy/pull/17749)
for details.

(deprecated-non-tuple-lambda)=
### Non-tuple iterable for the first argument to `Lambda`

Using a non-tuple as the first argument to {class}`~.Lambda` is deprecated. If
you have a non-tuple, convert it to a tuple first, like `Lambda(tuple(args),
expr)`.

This was done so that `Lambda` could support general tuple unpacking, like

```py
>>> from sympy import Lambda, symbols
>>> x, y, z = symbols('x y z')
>>> f = Lambda((x, (y, z)), x + y + z)
>>> f(1, (2, 3))
6
```

(deprecated-differentiate_finite-evaluate)=
### The `evaluate` flag to `differentiate_finite`

The `evaluate` flag to {func}`~.differentiate_finite` is deprecated.

`differentiate_finite(expr, x, evaluate=True)` expands the intermediate
derivatives before computing differences. But this usually not what you want,
as it does not satisfy the product rule.

If you really do want this behavior, you can emulate it with

```py
diff(expr, x).replace(
    lambda arg: arg.is_Derivative,
    lambda arg: arg.as_finite_difference())
```

See the discussion on issue [#17881](https://github.com/sympy/sympy/pull/17881).

## Version 1.4

(deprecated-tensorindextype-attrs)=
### `TensorIndexType.data` and related methods

The `TensorIndexType.data` property is deprecated, as well as several methods
which made use of it including the `get_matrix()`, the `__getitem__()`
(indexing), `__iter__()` (iteration), `_components_data_full_destroy()`, and `__pow__()` (`**`) methods. Storing data on tensor objects was a
design flaw and not consistent with how the rest of SymPy works.

Instead, the {meth}`.TensExpr.replace_with_arrays` method should be
used.

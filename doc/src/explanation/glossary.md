# Glossary

This page is a glossary for various terms used throughout the SymPy
documentation. This glossary is primarily for terms that are specific to
SymPy. For more general Python terms, refer to the [Python
glossary](https://docs.python.org/3/glossary.html). Mathematical terms are
only included here if they have a specific meaning in SymPy. For general
mathematical definitions, refer to other sources such as
[Wikipedia](https://www.wikipedia.org/) or
[MathWorld](https://mathworld.wolfram.com/), as well as the references in the
documentation for the specific SymPy functions.

<!-- Please follow the following style for the glossary:

- Capitalize the names, unless they refer to code.
- Each definition should be short (no more than a paragraph).
- Keep the terms in alphabetical order.
- Use cross-references to refer to full documentation for functions, etc.
- Classes/functions should only be listed here if they are core concepts for
  understanding how to use SymPy. Classes and functions themselves can be
  documented in their respective docstrings in the modules reference.

-->

<!-- To cross reference a glossary item elsewhere in the documentation use
{term}`termname` in Markdown or :term:`termname` in RST. Use
{term}`Text <termname>` or :term:`text <termname>` to make a link with
different text. -->

```{glossary}

Antiderivative

    An *antiderivative* of a function $f(x)$ with respect to $x$ is an
    function $F(x)$ such that $\frac{d}{dx}F(x) = f(x).$ It is also sometimes
    called the "indefinite integral" of $f(x)$, and written as $\int
    f(x)\,dx.$ Antiderivatives in SymPy can be computed with
    {func}`~.integrate`. Any given function has an infinite family of
    antiderivatives, typically obtained by adding an arbitrary constant, but
    {func}`~.integrate` only returns one
    such antiderivative, which is sometimes informally refereed to as "the
    antiderivative" of $f(x)$. Note some sources call this the "primitive" of
    $f(x)$, but that terminology is not used in SymPy because it is not as
    universally used as "antiderivative", and because "primitive" has other
    meanings.

`args`

    The *`args`* of a SymPy {term}`expression` `expr` is a tuple of the
    top-level subexpressions used to create `expr`. They are the arguments to
    the class used to create the expression. The args of any expression can be
    obtained by the `.args` attribute. For example, `(1 + x*y).args` is `(1,
    x*y)`, because it equals `Add(1, x*y)`. The `args` together with
    {term}`func` completely define an expression. It is always possible to
    walk the expression tree and extract any subexpression of a SymPy
    expression by repeated use of `.args`. Every SymPy expression can be
    rebuilt exactly with `func` and `args`, that is, `expr.func(*expr.args) ==
    expr` will always be true of any SymPy expression `expr`. The args of an
    expression may be the empty tuple, for instance when an expression is an
    {term}`atom`.

Assumptions

    *Assumptions* are a set of predicates on a {term}`symbol` or
    {term}`expression` that define the set of possible values it can take.
    Some examples of assumptions are `positive`, `real`, and `integer`.
    Assumptions are related to one another logically, for example, an
    assumption of `integer` will always automatically imply `real`.
    Assumptions use a {term}`three-valued logic` system where predicates are
    either `True`, `False`, or `None`.

    Assumptions are either *assumed* or *queried*. For example, a symbol `x`
    might be assumed to be positive by defining it as `x = symbols('x',
    positive=True)`. Then an assumption might be queried on an expression
    containing this symbol, like `(x + 1).is_real`.

    If no assumptions are assumed on a symbol, then by default symbols are
    assumed to be general complex numbers. Setting assumptions is important
    because certain simplifications are only mathematically true in a
    restricted domain, for example, $\sqrt{x^2} = x$ is not true for general
    complex $x$ but it is true when $x$ is positive. SymPy functions will
    never perform an operation on an expression unless it is true for all
    values allowed by its assumptions.

    SymPy has two separate assumptions systems, which are closely related to
    one another. The first, which is sometimes called the "old assumptions"
    because it is older, defines assumptions on {term}`Symbol` objects, is
    queried by {term}`is_*` attributes. The second, which is sometimes called
    the "new assumptions", defines assumptions using separate predicate
    objects like `Q.positive` and does queries using the {func}`~.ask` function.
    The newer assumptions system is able to support more complex queries, but
    is also not as built out as the older one. Most users of SymPy should
    prefer the older assumptions system at this time.

Atom

    An *atom* is an expression whose {term}`args` is the empty tuple `()`.
    Atoms are the leaves of the expression tree. For example, if a function
    uses recursion to walk an expression tree using `args`, the atomic
    expressions will be the base case of the recursion.

    Note that the class {class}`~.Atom` is sometimes used as the base class of
    atomic expressions, but it is not a requirement for atomic expressions to
    subclass this class. The only requirement for an expression to be atomic
    is for its {term}`args` to be empty.

{class}`~.Basic`

    *`Basic`* is the superclass of all SymPy expressions. It defines the basic
    methods required for an expression, such as {term}`args`, {term}`func`,
    equality, {term}`immutability <immutable>`, and some useful expression
    manipulation functions such as {term}`substitution <substitute>`. Most SymPy classes
    will subclass a more specific `Basic` subclass such as {term}`Boolean`,
    {term}`Expr`, {term}`Function <Function (class)>`, or {term}`Matrix`. An
    object that is not an instance `Basic` typically cannot be used in SymPy
    functions, unless it can be turned into a `Basic` class via
    {term}`sympify()`.

{class}`~.Boolean`

    *`Boolean`* is the base class for the classes in the {mod}`~.logic`
    module. An object that is an instance `Boolean` is a logical predicate
    that is an element of a boolean algebra and can be thought of as having a
    "true" or "false" value (note that `Boolean` objects do not use the
    {term}`three-valued logic` used by the {term}`assumptions`).

Canonicalize

    Often expressions can be written in multiple, mathematically equivalent
    ways. A *canonical form* is a single way of writing an expression, which
    all equivalent expressions can be transformed to. An expression that is
    put into a canonical form is said to be *canonicalized*. Typically
    canonical forms are unique and have properties that make them easier to
    work with. For example, a common canonical form used for rational function
    is $p/q$, where $p$ and $q$ are expanded polynomials with no common
    factors.

Code Generation

    *Code generation* refers to the process of taking a SymPy expression and
    converting it into code for a language or library so that it can be
    evaluated numerically. SymPy supports [code generation](codegen_module)
    for dozens of languages and libraries including C, C++, Fortran, and NumPy.

Core

    The [*core*](core_module) is the submodule that contains the important
    functionality used by all SymPy objects. This includes the {term}`Basic`
    and {term}`Expr` base classes, classes like {class}`~.Add`,
    {class}`~.Mul`, and {class}`~.Pow`, and the {term}`assumptions`.

{class}`~.Expr`

    *`Expr`* is the superclass of all algebraic SymPy expressions. It is
    itself a subclass of {term}`Basic`. SymPy expressions that can be in an
    {class}`~.Add`, {class}`~.Mul`, or {class}`~.Pow` should be `Expr`
    subclasses. Not all SymPy classes are subclasses of `Expr`, for example,
    {term}`Boolean` objects are {term}`Basic` but not `Expr`, because boolean
    expressions do not make mathematical sense in classes like {class}`~.Add`
    or {class}`~.Mul`.

`_eval_*`
    TODO

`evalf`
    TODO

Evaluated
    See {term}`Unevaluated`.

Expression
    TODO

`func`

    *`func` is the function of an {term}`expression`, which can be obtained by
    `expr.func`. This is usually the same as `type(expr)`, but may differ in
    some cases, so it should be preferred to `type(expr)` when rebuilding
    expressions with {term}`args`. Every SymPy expression can be rebuilt
    exactly with `func` and `args`, that is, `expr.func(*expr.args) == expr`
    will always be true of any SymPy expression `expr`.

Function
    TODO

{class}`~.Function` (class)
    TODO

Immutable
    TODO

Interactive
    TODO

`is_*`
    TODO

Kind
    TODO

lamda
    "*Lamda*" is just an alternate spelling of "lambda". It is used sometimes
    in SymPy because `lambda` is a reserved keyword in Python, so a symbol
    named λ must be named something else.

Matrix
    TODO

mpmath
    TODO

Number
    TODO

{class}`oo <sympy.core.numbers.Infinity>`
    *`oo`* is the SymPy object representing positive infinity. It is spelled
    this way, as two lower case letter Os because it resembles the symbol
    $\infty$ and is easy to type.

Polys
    TODO

Printing
    TODO

Relational
    TODO

{class}`S <sympy.core.singleton.Singleton>`
    TODO

Simplification
    TODO

Solve
    TODO

Substitute
    TODO

{class}`~.Symbol`
    TODO

{func}`~.sympify`
    TODO

Three-valued logic

    *Three-valued logic* is a logic with three values, `True`, `False`, and
    `None`. It is also sometimes called *fuzzy logic*, although this term also
    has different meanings in the mathematical literature so "three-valued
    logic" is preferred. `True` and `False` work the same as in the usual
    two-valued predicate logic. `None` is an additional term that represents
    "unknown", "noncomputable", or "could be either True or False"
    (philosophically these are distinct concepts, but logically they all
    function identically). The semantics of `None` are that it absorbs other
    terms in logical operations whenever the result would differ if it were
    replaced with `True` or `False`. For example, `None OR False` is `None`.
    However, `None OR True` would be `True` because the expression is `True`
    whether the `None` "really" represents a value of `True` or `False`. One
    must be careful to not use the usual Python logical operators like `and`,
    `or` and `not` on three-valued logic, since `None` is falsy. See
    [](booleans) for more details on how to code with three-valued logic.

    Three-valued logic is used by the {term}`assumptions` system to represent
    assumptions that are not known. For instance, `x.is_positive` might be
    `None` if `x` could be positive or negative under its given assumptions.
    Note that the predicate logic defined by the {term}`Boolean` class is a
    standard two-valued logic, not three-valued.

Unevaluated
    TODO

{class}`zoo <sympy.core.numbers.ComplexInfinity>`

    *`zoo`* represents complex infinity, i.e., the north pole of the Riemann
    sphere. The reason it is spelled this way is that it is "z-oo", where "z"
    is the symbol commonly used for complex variables, and {term}`oo` is the
    symbol SymPy uses for real positive infinity.

```

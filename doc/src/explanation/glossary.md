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
    such antiderivative, which is sometimes informally refered to as "the
    antiderivative" of $f(x)$. Note some sources call this the "primitive" of
    $f(x)$, but that termonology is not used in SymPy because it is not as
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

Canonical
    TODO

Code Generation
    TODO

Core
    TODO

{class}`~.Expr`
    TODO

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
    TODO

Matrix
    TODO

mpmath
    TODO

Number
    TODO

{class}`oo <sympy.core.numbers.Infinity>`
    TODO

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
    TODO

Unevaluated
    TODO

{class}`zoo <sympy.core.numbers.ComplexInfinity>`
    TODO

```

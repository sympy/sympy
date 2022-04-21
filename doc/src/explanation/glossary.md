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
- Vitalize the first occurrence of the term in the definition.
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

Automatic Simplification

    *Automatic Simplification* refers to any simplification that happens
    automatically inside of a class constructor. For example, `x + x` is
    automatically simplified to `2*x` in the {class}`~.Add` constructor.
    Unlike manual {term}`simplification`, the automatic simplification can
    only be disabled by setting `evaluate=False` (see {term}`Unevaluated`).
    Automatic simplification is often done so that expressions become
    {term}`canonicalized <canonicalize>`. Excessive automatic simplification
    is discouraged as makes it impossible to represent the non-simplified form
    of the expression without using tricks like `evaluate=False`, and it can
    often be an expensive thing to do in a class constructor. Instead, manual
    {term}`simplification` and {term}`canonicalization <canonicalize>` is
    generally preferred.

{class}`~.Basic`

    *`Basic`* is the superclass of all SymPy expressions. It defines the basic
    methods required for an expression, such as {term}`args`, {term}`func`,
    equality, {term}`immutability <immutable>`, and some useful expression
    manipulation functions such as {term}`substitution <substitute>`. Most
    SymPy classes will subclass a more specific `Basic` subclass such as
    {term}`Boolean`, {term}`Expr`, {term}`Function <Function (class)>`, or
    {term}`Matrix`. An object that is not an instance `Basic` typically cannot
    be used in SymPy functions, unless it can be turned into a `Basic` class
    via {term}`sympify()`.

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

Equation

    An *equation* is an {term}`expression` that has an equality. Equations in
    SymPy are represented using the {class}`Eq
    <sympy.core.relational.Equality>` class. Equations are **not** created
    using the `==` operator. The `==` operator does a direct {term}`structural
    equality` check against two expressions, and always returns `True` or
    `False`. To contrast, a symbolic equation may be unevaluated. Equations
    are considered {term}`booleans <boolean>` since they mathematically
    represent a predicate value that is either true or false.

`_eval_*`

    Various methods on {term}`Basic` and {term}`Expr` can be defined on
    subclasses using special *`_eval_*`* methods. For example, an object can
    define how it will be processed by the {func}`~.simplify` function by
    defining a `_eval_simplify` method. `_eval_*` methods are instead of
    overriding the method itself so that the method defined on the base class
    can do pre-processing before dispatching to the `_eval_*` method.

`evalf`

    [*`evalf`*](sympy.core.evalf.EvalfMixin.evalf) is the method present on
    all {term}`Expr` objects that evaluates it to a floating-point numerical
    value, or converts the constant parts of the expression to a numerical
    value if it contains {term}`symbols <symbol>`. "evalf" stands for
    "evaluate floating-point". `evalf` uses {term}`mpmath` under the hood to
    evaluate expressions to arbitrary precision.

Evaluate

    *Evaluate* can refer to:

    - The process of converting an expression into a
      numerical value (see {term}`evalf`).

    - The process of {term}`automatic simplification` that occurs when
      creating an expression (see {term}`Unevaluated`).

{class}`~.Expr`

    *`Expr`* is the superclass of all algebraic SymPy expressions. It is
    itself a subclass of {term}`Basic`. SymPy expressions that can be in an
    {class}`~.Add`, {class}`~.Mul`, or {class}`~.Pow` should be `Expr`
    subclasses. Not all SymPy classes are subclasses of `Expr`, for example,
    {term}`Boolean` objects are {term}`Basic` but not `Expr`, because boolean
    expressions do not make mathematical sense in classes like {class}`~.Add`
    or {class}`~.Mul`.

Expression

    Any SymPy object, that is, any instance of {term}`Basic`, may be called an
    *expression*. Sometimes, the term "expression" is reserved for
    {term}`Expr` objects, which are algebraic expressions. Expressions are not
    to be confused with {term}`equations <equation>`, which are a specific
    type of expression with that represents a mathematical equality.

`func`

    *`func` is the function of an {term}`expression`, which can be obtained by
    `expr.func`. This is usually the same as `type(expr)`, but may differ in
    some cases, so it should be preferred to `type(expr)` when rebuilding
    expressions with {term}`args`. Every SymPy expression can be rebuilt
    exactly with `func` and `args`, that is, `expr.func(*expr.args) == expr`
    will always be true of any SymPy expression `expr`.

Function

    *Function* may refer to:

    - A mathematical function, that is, something which maps values from a
      domain to a range. Sometimes an {term}`expression` containing a
      {term}`symbol` is colloquially called a "function" because the symbol be
      replaced with a value using {term}`substitution <substitute>`,
      evaluating the expression. This usage is colloquial because one must use
      the {meth}`subs <sympy.core.basic.Basic.subs>` method to do this rather
      than the typical Python function calling syntax, and because it is not
      specific about what variable(s) the expression is a function of. An
      expression can be converted into a function object that can be called
      using Python `f(x)` syntax using {class}`~.Lambda`.

    - An instance of the SymPy {term}`Function <Function (class)>` class.

    - A Python function, i.e., a function defined using the `def` keyword.
      Python functions are not {term}`symbolic`, since they must always return
      a value and cannot be {term}`unevaluated`.

{class}`~.Function` (class)

    *`Function`* is the base class of algebraic functions in SymPy. This
    includes common functions like {class}`~.sin` and {class}`~.exp`, special
    functions like {class}`~.zeta` and {class}`~.hyper`, and integral
    functions like {func}`~.primepi` and {class}`~.divisor_sigma`. Function
    classes are always {term}`symbolic`, meaning that they typically remain
    {term}`unevaluated` when passed a {term}`symbol`, like `f(x)`.

    `Function` may also be used to create an {term}`undefined function` by
    passing it a string name for the function, like `Function('f')`.

    Not every function in SymPy is a symbolic `Function` class; some are just
    Python functions which always return a value. For example, most
    simplification functions like {term}`simplify <simplification>` cannot be
    represented symbolically.


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

    *Matrix* refers to the set of classes used by SymPy to represent matrices.
    SymPy actually has several internal classes to represent matrices,
    depending on whether the matrix is symbolic ({class}`~.MatrixExpr`) or
    explicit, mutable or immutable, and what the underlying elements are.

mpmath

    [*mpmath*](https://mpmath.org/) is a pure Python library for arbitrary
    precision numerics. It is a [hard dependency](dependencies-mpmath) of
    SymPy. mpmath is used under the hood whenever SymPy evaluates an
    expression numerically, such as when using {term}`evalf`.

Numeric

    A *numeric* representation or algorithm is one that operates directly on
    numeric inputs. It is in contrast with a {term}`symbolic` representation
    or algorithm, which can work with objects in an unevaluated form. Often a
    numerical algorithm is quite different from a symbolic one. For example,
    numerically solving an ODE typically means evaluating the ODE using an
    algorithm like
    [Runge–Kutta](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods)
    to find a set of numeric points given an initial condition, whereas
    symbolically solving an ODE (such as with SymPy's {func}`~.dsolve`) means
    mathematically manipulating the ODE to produce a {term}`symbolic`
    {term}`equation` that represents the solution. A symbolic ODE solution may
    including symbolic constants which can represent any numerical value.
    Numeric algorithms also need to be concerned with things issues caused by
    floating-point numbers such as loss of precision and numerical stability,
    whereas symbolic algorithms are not concerned with these things because
    they compute things exactly.

    Most scientific libraries other than SymPy, including NumPy or SciPy, are
    strictly numerical, meaning the functions in those libraries can only
    operate on specific numeric inputs. They will not work with SymPy
    expressions, because their algorithms are not designed to work with
    symbolic inputs. SymPy focuses on symbolic functions, leaving purely
    numerical code to other tools like NumPy.

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

Structural Equality

    Two {term}`expressions <expression>` are *structurally equal* if they are
    the same SymPy object. Two structurally equal expressions are considered
    to be identical by SymPy, since all SymPy expressions are
    {term}`immutable`. Structural equality can be checked with the `==`
    operator, which always returns True or False.Symbolic {term}`equality
    <equation>` can be represented using {class}`Eq
    <sympy.core.relational.Equality>`.

    Typically, two expressions are structurally equal if they are the same
    class and (recursively) have the same {term}`args`. Two expressions may be
    mathematically identical but not structurally equal. For example, `(x +
    1)**2` and `x**2 + 2*x + 1` are mathematically equal, but they are not
    structurally equal, because the first is a {class}`~.Pow` of an
    {class}`~.Add` and `Integer(2)`, and the second is an {class}`~.Add` with
    three arguments.

Substitute

    TODO

Symbolic

    A *symbolic* representation of a mathematical object is a representation
    that is partially or completely unevaluated at runtime. It may include
    named {term}`symbols <symbol>` in place of explicit numeric values.
    A symbolic representation is often contrasted with a {term}`numeric` one.
    Symbolic representations are mathematically exact, to contrast with
    numeric representations which are typically rounded so they can fit within
    a floating-point value.
    Symbolic {term}`expression <expression>` representing mathematical objects
    may be aware of mathematical properties of these objects and be able to
    simplify. The goal of SymPy is to represent and manipulate symbolic
    expressions representing mathematical objects.

{class}`~.Symbol`

    A *`Symbol`* is an object that represents a single mathematical variable
    in an expression. {class}`~.Symbol` is a subclass of {term}`Expr`, is
    {term}`atomic <atom>`. A `Symbol` contains a name, which is any string,
    and {term}`assumptions`. Symbols are typically defined with the `Symbol`
    constructor of the {func}`~.symbols` function. Two Symbols with the same
    name and assumptions are considered equal.

{func}`~.sympify`

    *`sympify()`* (not to be confused with *{term}`simplify()
    <simplification>`*) is a function that converts non-SymPy objects into
    SymPy objects. The result of `sympify()` will be an instance of
    {term}`Basic`. Objects that can be *sympified* include native Python
    numeric types such as `int` and `float`, strings that are parsable as
    SymPy {term}`expressions <expression>`, and iterables containing
    *sympifiable* objects (see the documentation for {func}`~.sympify` for
    more information).

    Since all SymPy {term}`expressions <expression>` must be instances of
    {term}`Basic`, all SymPy functions and operations will implicitly call
    `sympify()` on their inputs. For example, `x + 1` implicitly calls
    `sympify(1)` to convert the `1` that is a Python `int` into a SymPy
    {class}`~.Integer`. Functions that accept SymPy expressions should
    typically call `sympify()` on their arguments so that they work even when
    the input is not a SymPy type.

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

Undefined Function

    An *undefined function* is a {term}`Function <Function (class)>` that has
    no mathematical properties defined on it. It always remains
    {term}`unevaluated`. An undefined function can be created by passing a
    string name of the function to `Function`, like `f = Function('f')`.
    Undefined functions are commonly used when working with [ODEs](ode-docs).
    Undefined functions are also the best way to make {term}`symbols <symbol>`
    that depend on other symbols. For example, if `f = Function('f')` and `x =
    Symbol('x')`, then SymPy will know that `f(x)` depends on `x`, meaning,
    for instance, that the derivative `diff(f(x), x)` will not be evaluated to
    0. Otherwise, by default, SymPy assumes that all symbols do not depend on
    other symbols, that is, that they are mathematically "constant" or
    "independent variables" .

Unevaluated

    An expression is *unevaluated* if the {term}`automatic simplification`
    that typically occurs when the expression is created is disabled. This is
    typically done by setting `evaluate=False`, using `with evaluate(False)`,
    or using {class}`~.UnevaluatedExpr`. While unevaluated expressions are
    supported, they can sometimes lead to surprising behavior because the
    expressions are not properly {term}`canonicalized <canonicalize>`.

{class}`zoo <sympy.core.numbers.ComplexInfinity>`

    *`zoo`* represents complex infinity, i.e., the north pole of the Riemann
    sphere. The reason it is spelled this way is that it is "z-oo", where "z"
    is the symbol commonly used for complex variables, and {term}`oo` is the
    symbol SymPy uses for real positive infinity.

```

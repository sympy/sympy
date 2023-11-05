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
- Italicize the first occurrence of the term in the definition.
- Use cross-references to refer to full documentation for functions, etc.
- Cross-reference the first instance of other terms in definitions (see
  below).
- Classes and functions should only be listed here if they are core concepts
  for understanding how to use SymPy. Classes and functions themselves should
  be documented in their respective docstrings in the modules reference.

-->

<!-- To cross reference a glossary item elsewhere in the documentation use
{term}`termname` in Markdown or :term:`termname` in RST. Use
{term}`Text <termname>` or :term:`text <termname>` to make a link with
different text. -->

```{glossary}

Antiderivative

    An *antiderivative* of a function $f(x)$ with respect to $x$ is a function
    $F(x)$ such that $\frac{d}{dx}F(x) = f(x).$ It is also sometimes called an
    "indefinite integral" of $f(x)$, and written as $\int f(x)\,dx.$
    Antiderivatives in SymPy can be computed with {func}`~.integrate`. Note
    some sources call this the "primitive" of $f(x)$, but this terminology is
    not used in SymPy because it is not as universally used as
    "antiderivative", and because "primitive" has other meanings in
    mathematics and in {func}`SymPy <sympy.polys.polytools.primitive>`.

`args`

    The *`args`* property of a SymPy {term}`expression` is a tuple of the
    top-level {term}`subexpressions <subexpression>` used to create it. They
    are the arguments to the class used to create the expression. The args of
    any expression can be obtained by the `.args` attribute. For example,
    `(1 + x*y).args` is `(1, x*y)`, because it equals `Add(1, x*y)`. The `args`
    together with {term}`func` completely define an expression. It is always
    possible to walk the {term}`expression tree` and extract any subexpression
    of a SymPy expression by repeated use of `.args`. Every SymPy expression
    can be rebuilt exactly with `func` and `args`, that is,
    `expr.func(*expr.args) == expr` will always be true of any SymPy
    expression `expr`. The args of an expression may be the empty tuple `()`,
    meaning the expression is an {term}`atom`.

Assumptions

    *Assumptions* are a set of predicates on a {term}`symbol` or
    {term}`expression` that define the set of possible values it can take.
    Some examples of assumptions are `positive`, `real`, and `integer`.
    Assumptions are related to one another logically, for example, an
    assumption of `integer` automatically implies `real`. Assumptions use a
    {term}`three-valued logic` system where predicates are either `True`,
    `False`, or `None`.

    Assumptions are either *assumed* or *queried*. For example, a symbol `x`
    might be *assumed* to be positive by defining it as `x = symbols('x',
    positive=True)`. Then an assumption might be *queried* on the expression
    containing this symbol, like `(x + 1).is_real`, which in this case would
    return `True`.

    If no assumptions are assumed on a symbol, then by default symbols are
    assumed to be general complex numbers. Setting assumptions is important
    because certain simplifications are only mathematically true in a
    restricted domain, for example, $\sqrt{x^2} = x$ is not true for general
    complex $x$ but it is true when $x$ is positive. SymPy functions will
    never perform an operation on an expression unless it is true for all
    values allowed by its assumptions.

    SymPy has two separate assumptions systems, which are closely related to
    one another. In the first, which is sometimes called the "old assumptions"
    because it is older, assumptions are assumed on {term}`Symbol`
    objects and queried with {term}`is_*` attributes. In the second, which is
    sometimes called the "new assumptions",  assumptions are assumed
    using separate predicate objects like `Q.positive` and queried using
    the {func}`~.ask` function. The newer assumptions system is able to
    support more complex queries, but is also not as well developed as the
    older one. Most users of SymPy should prefer the older assumptions system
    at this time.

    See the {ref}`assumptions guide <assumptions-guide>` for more details on
    assumptions.

Atom

    An *atom* is an expression whose {term}`args` is the empty tuple `()`.
    Atoms are the leaves of the {term}`expression tree`. For example, if a
    function uses recursion to walk an expression tree using `args`, the
    atomic expressions will be the base case of the recursion.

    Note that the class {class}`~.Atom` is sometimes used as the base class of
    atomic expressions, but it is not a requirement for atomic expressions to
    subclass this class. The only requirement for an expression to be atomic
    is for its {term}`args` to be empty.

Automatic Simplification

    *Automatic Simplification* refers to any simplification that happens
    automatically inside of a class constructor. For example, `x + x` is
    automatically simplified to `2*x` in the {class}`~.Add` constructor.
    Unlike manual {term}`simplification`, automatic simplification can only be
    disabled by setting `evaluate=False` (see {term}`Unevaluated`). Automatic
    simplification is often done so that expressions become
    {term}`canonicalized <canonicalize>`. Excessive automatic simplification
    is discouraged, as it makes it impossible to represent the non-simplified
    form of the expression without using tricks like `evaluate=False`, and it
    can often be an expensive thing to do in a class constructor. Instead,
    manual {term}`simplification`/{term}`canonicalization <canonicalize>`
    is generally preferred.

{class}`~.Basic`

    *{class}`~.Basic`* is the superclass of all SymPy expressions. It defines
    the basic methods required for a SymPy expression, such as {term}`args`,
    {term}`func`, {term}`equality <structural equality>`, {term}`immutability
    <immutable>`, and some useful expression manipulation functions such as
    {term}`substitution`. Most SymPy classes will subclass a more specific
    `Basic` subclass such as {term}`Boolean`, {term}`Expr`, {term}`Function
    <Function (class)>`, or {term}`Matrix`. An object that is not a `Basic`
    instance typically cannot be used in SymPy functions, unless it can be
    turned into one via {term}`sympify()`.

{class}`~.Boolean`

    *{class}`~.Boolean`* is the base class for the classes in the
    {mod}`~.logic` module. `Boolean` instances represent logical predicates
    that are elements of a [boolean
    algebra](https://en.wikipedia.org/wiki/Boolean_algebra) and can be thought
    of as having a "true" or "false" value (note that `Boolean` objects do not
    use the {term}`three-valued logic` used by the {term}`assumptions`).

Bound symbols

    A {term}`symbol` in an expression is *bound* if it is not {term}`free
    <free symbols>`. A bound symbol can be replaced everywhere with new symbol
    and the resulting expression will still be mathematically equivalent.
    Examples of bound symbols are integration variables in definite integrals
    and substituted variables in a {class}`~.Subs`. Bound symbols are
    sometimes represented by {term}`dummy` symbols, but the are not always
    {class}`~.Dummy` objects, and {class}`~.Dummy` objects are not always
    bound symbols.

Canonical Form
Canonicalize

    Often expressions can be written in multiple, mathematically equivalent
    ways. A *canonical form* is a single way of writing an expression, which
    all equivalent expressions can be transformed to. An expression that is
    put into a canonical form is said to be *canonicalized*. Often canonical
    forms are unique and have properties that make them easier to work with.
    For example, a common canonical form used for rational functions is
    $\frac{p}{q}$, where $p$ and $q$ are expanded polynomials with no common
    factors.

Code Generation

    *Code generation* refers to the process of taking a SymPy expression and
    converting it into code for a language or library so that it can be
    evaluated numerically. SymPy supports code generation for [dozens of
    languages](codegen_module) and libraries including C, C++, Fortran, and
    NumPy.

Core

    The [*core*](core_module) is the submodule that contains the important
    functionality used by all SymPy objects. This includes the {term}`Basic`
    and {term}`Expr` base classes, classes like {class}`~.Add`,
    {class}`~.Mul`, and {class}`~.Pow`, and the {term}`assumptions`.

Dummy

    A *dummy* {term}`symbol` is a symbol that is automatically unequal to any
    other dummy symbol other than itself, even if it has the same name. Dummy
    symbols are used when a function needs to return an expression with a new
    symbol, so that it cannot accidentally clash with a {term}`symbol` of the
    same name. Dummy symbols can be created with {class}`~.Dummy`.

Equation

    An *equation* is an {term}`expression` that has an equals sign $=$.
    Equations in SymPy are represented using the {class}`Eq
    <sympy.core.relational.Equality>` class. Equations are **not** created
    using the `==` operator. The `==` operator does a {term}`structural
    equality` check between two expressions, and always returns `True` or
    `False`. To contrast, a symbolic equation may be {term}`unevaluated`.
    Equations are considered {term}`booleans <boolean>` since they
    mathematically represent a predicate value that is either true or false.

`_eval_*`

    Various methods on {term}`Basic` and {term}`Expr` can be defined on
    subclasses using special *`_eval_*`* methods. For example, an object can
    define how it will be processed by the {func}`~.diff` function by defining
    a `_eval_derivative` method. `_eval_*` methods used are instead of
    overriding the method itself so that the method defined on the base class
    can do pre-processing before dispatching to the `_eval_*` method.

`evalf`

    [*`evalf`*](sympy.core.evalf.EvalfMixin.evalf) is the method present on
    every {term}`Expr` object that evaluates it to a floating-point numerical
    value, or converts the constant parts of the expression to a numerical
    value if it contains {term}`symbols <symbol>`. The {meth}`.n()
    <sympy.core.evalf.EvalfMixin.n>` method and {func}`~.N` function are both
    shorthands for `evalf`. "evalf" stands for "evaluate floating-point".
    `evalf` uses {term}`mpmath` under the hood to evaluate expressions to
    arbitrary precision.

Evaluate

    *Evaluate* can refer to:

    - The process of converting an {term}`expression` into a
      numerical value (see {term}`evalf`).

    - The process of {term}`automatic simplification` that occurs when
      creating an expression (see {term}`Unevaluated`).

    - The process of replacing one or more {term}`symbols <symbol>` in an
      expression with numeric values or with other expressions using
      {term}`substitution`.

{class}`~.Expr`

    *{class}`~.Expr`* is the superclass of all algebraic SymPy expressions. It
    is itself a subclass of {term}`Basic`. SymPy expressions that can be in an
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
    types of expressions that represents mathematical equalities.

Expression Tree

    An *expression tree* is a
    [tree](https://en.wikipedia.org/wiki/Tree_(data_structure)) of
    {term}`expressions <expression>`. Every expression is built up from
    smaller expressions as a tree. The nodes of an expression tree are
    expressions and the children of each node are the direct
    {term}`subexpressions <subexpression>` that constitute that expression.
    Alternatively, one can view an expression tree as a tree where the
    non-leaf nodes are {term}`funcs <func>` and the leaf nodes are
    {term}`atoms <atom>`. An example expression tree is shown in the
    [tutorial](tutorial-expression-trees). The expression tree of any SymPy
    expression can be obtained by recursing through {term}`args`. Note that
    because SymPy expressions are {term}`immutable` and are treated equal
    strictly by {term}`structural equality`, one may also think of an
    expression tree as being a
    [DAG](https://en.wikipedia.org/wiki/Directed_acyclic_graph), where
    identical subexpressions are only represented in the graph once.

Free symbols

    A {term}`symbol` in an expression is *free* if the expression
    mathematically depends on the value of that symbol. That is, if the symbol
    were replaced with a new symbol, the result would be a different
    expression. Symbols that are not free are {term}`bound <bound symbols>`.
    The free symbols of an expression can be accessed with the
    {attr}`free_symbols <sympy.core.basic.Basic.free_symbols>` attribute.

`func`

    The *`func`* property is the function of an {term}`expression`, which can
    be obtained by `expr.func`. This is usually the same as `type(expr)`, but
    may differ in some cases, so it should be preferred to use `expr.func`
    instead of` type(expr)` when rebuilding expressions with {term}`args`.
    Every SymPy expression can be rebuilt exactly with `func` and `args`, that
    is, `expr.func(*expr.args) == expr` will always be true of any SymPy
    expression `expr`.

Function

    *Function* may refer to:

    - A mathematical function, that is, something which maps values from a
      domain to a range. Sometimes an {term}`expression` containing a
      {term}`symbol` is colloquially called a "function" because the symbol
      can be replaced with a value using {term}`substitution`,
      {term}`evaluating <evaluate>` the expression. This usage is colloquial
      because one must use the {meth}`subs <sympy.core.basic.Basic.subs>`
      method to do this rather than the typical Python function calling
      syntax, and because it is not specific about what variable(s) the
      expression is a function of, so generally the term "expression" should
      be preferred unless something is an actual function. An expression can
      be converted into a function object that can be called using the Python
      `f(x)` syntax using {class}`~.Lambda`.

    - An instance of the SymPy {term}`Function <Function (class)>` class.

    - A Python function, i.e., a function defined using the `def` keyword.
      Python functions are not {term}`symbolic`, since they must always return
      a value and thus cannot be {term}`unevaluated`.

{class}`~.Function` (class)

    *{class}`~.Function`* is the base class of symbolic functions in SymPy.
    This includes common functions like {class}`~.sin()` and {class}`~.exp()`,
    special functions like {class}`~.zeta()` and {class}`~.hyper()`, and
    integral functions like {func}`~.primepi` and {class}`~.divisor_sigma()`.
    Function classes are always {term}`symbolic`, meaning that they typically
    remain {term}`unevaluated` when passed a {term}`symbol`, like `f(x)`. Not
    every symbolic {term}`expression` class is a `Function` subclass, for
    example, {term}`core` classes like `Add` and `Mul` are not `Function`
    subclasses.

    `Function` may also be used to create an {term}`undefined function` by
    passing it a string name for the function, like `Function('f')`.

    Not every function in SymPy is a symbolic `Function` class; some are just
    Python functions which always return a value. For example, most
    simplification functions like {term}`simplify() <simplification>` cannot
    be represented symbolically.


Immutable

    In Python, objects are *immutable* if they can not be modified in-place.
    In order to change an immutable object, a new object must be created. In
    SymPy, all {term}`Basic` objects are immutable. This means that all
    functions that operate on {term}`expressions <expression>` will return a
    new expression and leave the original unchanged. Performing an operation
    on an expression will never change other objects or expressions that
    reference that expression. This also means that any two objects that are
    {term}`equal <structural equality>` are completely interchangeable and may
    be thought of as being the same object, even if they happen to be two
    different objects in memory. Immutability makes it easier to maintain a
    mental model of code, because there is no hidden state. SymPy objects
    being immutable also means that they are hashable, which allows them to be
    used as dictionary keys.

Interactive

    *Interactive* usage refers to using SymPy in an interactive REPL
    environment such as the Python prompt, {term}`isympy`,
    [IPython](https://ipython.org/), or the [Jupyter
    notebook](https://jupyter.org/). When using SymPy interactively, all
    commands are typed in real time by the user and all intermediate results
    are shown. *Interactive* use is in contrast with *programmatic* use, which
    is where the code is written in a file which is either executed as a
    script or is part of a larger Python library. Some SymPy idioms are only
    recommended for interactive use and are considered anti-patterns when used
    programmatically. For example, running `from sympy import *` is convenient
    when using SymPy interactively, but is generally frowned upon for
    programmatic usage, where importing names explicitly just using `import
    sympy` is preferred.

`is_*`

    Attributes in SymPy that start with *`is_`* and use a *lowercase* name
    query the given {term}`assumption <assumptions>` on that object (note:
    there are a few properties that are an exception to this because they do
    not use the assumptions system, see {ref}`the assumptions guide
    <assumptions-guide-other-is-properties>`). For example, `x.is_integer`
    will query the `integer` assumption on `x`. `is_*` attributes that use a
    *Capitalized* name test if an object is an instance of the given class.
    Sometimes the same name will exist for both the lowercase and Capitalized
    property, but they are different things. For example, `x.is_Integer` is
    only `True` if `x` is an instance of {class}`~.Integer`, whereas
    `x.is_integer` is `True` if `x` is `integer` in the assumptions system,
    such as `x = symbols('x', integer=True)`. In general, it is recommended to
    not use `is_Capitalized` properties. They exist for historical purposes,
    but they are unneeded because the same thing can be achieved with
    `isinstance()`. See also {term}`Number`.

`isympy`

    *`isympy`* is a command that ships with SymPy that starts an
    {term}`interactive` session on the command line with all SymPy names
    imported and {term}`printing` enabled. It uses
    [IPython](https://ipython.org/) by default when it is installed.

Kind

    The *kind* of a SymPy object represents what sort of mathematical object
    it represents. The kind of an object can be accessed with the `kind`
    attribute. Example kinds are {any}`NumberKind`, which represents complex
    numbers, {any}`MatrixKind`, which represents matrices of some other kind,
    and {any}`BooleanKind`, which represents boolean predicates. The kind of a
    SymPy object is distinct from its Python type, since sometimes a single
    Python type may represent many different kinds of objects. For example,
    `Matrix` could be a matrix of complex numbers or a matrix of objects from
    some other ring of values. See [the classification of SymPy
    objects](kind_classification) page for more details about kinds in SymPy.

lamda

    "*Lamda*" is just an alternate spelling of the Greek letter "lambda". It
    is used sometimes in SymPy because `lambda` is a reserved keyword in
    Python, so a symbol representing λ must be named something else.

{func}`~.lambdify`

    *{func}`~.lambdify`* is a function that converts a SymPy
    expression into a Python function that can be evaluated numerically,
    typically making use of a {term}`numeric` library such as NumPy.

Matrix

    *Matrix* refers to the set of classes used by SymPy to represent matrices.
    SymPy has several internal classes to represent matrices, depending on
    whether the matrix is symbolic ({class}`~.MatrixExpr`) or explicit,
    mutable or immutable, dense or sparse, and what type the underlying
    elements are, but these are often all just called "Matrix".

mpmath

    [*mpmath*](https://mpmath.org/) is a pure Python library for arbitrary
    precision numerics. It is a [hard dependency](dependencies-mpmath) of
    SymPy. mpmath is capable of computing {term}`numerical <numeric>`
    functions to any given number of digits. mpmath is used under the hood
    whenever SymPy evaluates an expression numerically, such as when using
    {term}`evalf`.

Numeric

    A *numeric* representation or algorithm is one that operates directly on
    numeric inputs. It is in contrast with a *{term}`symbolic`* representation
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
    Numeric algorithms are typically designed around issues caused by
    floating-point numbers such as loss of precision and numerical stability,
    whereas symbolic algorithms are not concerned with these things because
    they compute things exactly.

    Most scientific libraries other than SymPy, such as NumPy or SciPy, are
    strictly numerical, meaning the functions in those libraries can only
    operate on specific numeric inputs. They will not work with SymPy
    expressions, because their algorithms are not designed to work with
    symbolic inputs. SymPy focuses on symbolic functions, leaving purely
    numerical code to other tools like NumPy. However, SymPy does interface
    with numerical libraries via tools like {term}`code generation` and
    {term}`lambdify()`.

Number

    *Number* can refer to two things in SymPy:

    - The class {class}`~.Number`, which is the base class for explicit
      numbers ({class}`~.Integer`, {class}`~.Rational`, and {class}`~.Float`).
      Symbolic numeric constants like {class}`pi <sympy.core.numbers.Pi>` are
      not instances of `Number`.

    - Lowercase "*number*", as in the `is_number` property, refers to any
      {term}`expression` that can be {term}`evalfed <evalf>` into an explicit
      `Number`. This includes symbolic constants like {class}`pi
      <sympy.core.numbers.Pi>`. Note that `is_number` is not part of the
      {term}`assumptions` system.

    This distinction is important for the `is_Number` and `is_number`
    properties. `x.is_Number` will check if `x` is an instance of the class
    {class}`~.Number`.

{class}`oo <sympy.core.numbers.Infinity>`

    *{class}`oo <sympy.core.numbers.Infinity>`* is the SymPy object
    representing positive infinity. It is spelled this way, as two lower case
    letter Os, because it resembles the symbol $\infty$ and is easy to type.
    See also {term}`zoo`.

Polys

    The *polys* refers to the {mod}`sympy.polys` submodule, which implements
    the basic data structures and algorithms for polynomial manipulation. The
    polys are a key part of SymPy (though not typically considered part of the
    {term}`core`), because many basic symbolic manipulations can be
    represented as manipulations on polynomials. Many algorithms in SymPy make
    use of the polys under the hood. For example, {func}`~.factor` is a
    wrapper around the polynomial factorization algorithms that are
    implemented in the polys. The classes in the polys are implemented using
    efficient data structures, and are not subclasses of {term}`Basic` like
    the other classes in SymPy.

Printing

    *Printing* refers to the act of taking an {term}`expression` and
    converting it into a form that can be viewed on screen. Printing is also
    often used to refer to {term}`code generation`. SymPy has several printers
    which represent expressions using different formats. Some of the more
    common printers are the string printer (`str()`), the pretty printer
    ({func}`pprint() <sympy.printing.pretty.pretty.pretty_print>`) the LaTeX
    printer ({func}`~.latex`), and code printers.

Relational

    A *relational* is an {term}`expression` that is a {term}`symbolic`
    {term}`equality <equation>` (like $a=b$), or a symbolic inequality like
    "less than" ($a<b$). Equality ($=$) and non-equality ($\neq$) relationals
    are created with {class}`Eq <sympy.core.relational.Equality>` and
    {class}`Ne <sympy.core.relational.Unequality>`, respectively. For example,
    `Eq(x, 0)` represents $x=0$. These should be used instead of `==` or `!=`,
    as these are used for {term}`structural <structural equality>` rather than
    symbolic equality. Inequality relationals can be created directly using
    `<`, `<=`, `>`, and `>=`, like `x < 0`.

{class}`S <sympy.core.singleton.Singleton>`

    The *{class}`S <sympy.core.singleton.Singleton>`* object in SymPy has two
    purposes:

    - It holds all singleton classes as attributes. Some special classes in
      SymPy are singletonized, meaning that there is always exactly one
      instance of them. This is an optimization that allows saving memory. For
      instance, there is only ever one instance of `Integer(0)`, which is
      available as `S.Zero`.

    - It serves as a shorthand for {term}`sympify()`, that is `S(a)` is the
      same as `sympify(a)`. This is useful for converting integers to SymPy
      Integers in expressions to avoid dividing Python ints (see [the gotchas
      section of the tutorial](tutorial-gotchas-final-notes)).

Simplification

    *Simplification* (not to be confused with {term}`sympify <sympify()>`)
    refers to the process of taking an {term}`expression` and transforming it
    into another expression that is mathematically equivalent but which is
    somehow "simpler". The adjective "simple" is actually not very
    well-defined. What counts as simpler depends on the specific use-case and
    personal aesthetics.

    The SymPy function {func}`~.simplify` heuristically tries various
    simplification algorithms to try to find a "simpler" form of an
    expression. If you aren't particular about what you want from "simplify",
    it may be a good fit. But if you have an idea about what sort of
    simplification you want to apply, it is generally better to use one or
    more of targeted [simplification functions](simplify-docs) which apply
    very specific mathematical manipulations to an expression.

Solve
Solvers

    To *solve* an {term}`equation` or system of equations means to find a set
    of {term}`expressions <expression>` that make the equation(s) true when
    the given {term}`symbol(s) <symbol>` are {term}`substituted
    <substitution>` with them. For example, the solution to the equation $x^2
    = 1$ with respect to $x$ would be the set $\{-1, 1\}$. Different types of
    equations can be solved by SymPy using different
    [*solvers*](solving-guide) functions. For instance, algebraic equations
    can be solved with {func}`~.solve`, differential equations can be solved
    with {func}`~.dsolve`, and so on.

    SymPy generally uses the word "solve" and "solvers" to mean equation
    solving in this sense. It is not used in the sense of "solving a problem".
    For instance, one would generally prefer to say "compute an integral" or
    "evaluate an integral" rather than "solve an integral" to refer to
    symbolic integration using the function {func}`~.integrate`.

Structural Equality

    Two SymPy objects are *structurally equal* if they are equal as
    {term}`expressions <expression>`, that is, they have the same
    {term}`expression trees <expression tree>`. Two structurally equal
    expressions are considered to be identical by SymPy, since all SymPy
    expressions are {term}`immutable`. Structural equality can be checked with
    the `==` operator, which always returns `True` or `False`. Symbolic
    {term}`equality <equation>` can be represented using {class}`Eq
    <sympy.core.relational.Equality>`.

    Typically, two expressions are structurally equal if they are the same
    class and (recursively) have the same {term}`args`. Two expressions may be
    mathematically identical but not structurally equal. For example, `(x +
    1)**2` and `x**2 + 2*x + 1` are mathematically equal, but they are not
    structurally equal, because the first is a {class}`~.Pow` whose
    {term}`args` consist of an {class}`~.Add` and an {class}`~.Integer`, and
    the second is an {class}`~.Add` whose {term}`args` consist of a
    {class}`~.Pow`, a {class}`~.Mul`, and an {class}`~.Integer`.

    Two apparently different expressions may be structurally equal if they are
    {term}`canonicalized <canonicalize>` to the same thing by {term}`automatic
    simplification`. For example, `x + y` and `y + x` are structurally equal
    because the {class}`~.Add` constructor automatically sorts its arguments,
    making them both the same.

Subexpression

    A *subexpression* is an {term}`expression` that is contained within a
    larger expression. A subexpression appears somewhere in the
    {term}`expression tree`. For `Add` and `Mul` terms, commutative and
    associative laws may be taken into account when determining what is a
    subexpression. For instance, `x + y` may sometimes be considered a
    subexpression of `x + y + z`, even though the expression tree for `Add(x,
    y)` is not a direct child of the expression tree for `Add(x, y, z)`.

Substitution

    *Substitution* refers to the act of replacing a {term}`symbol` or
    {term}`subexpression` inside of an {term}`expression` with another
    expression. There are different methods in SymPy for performing
    substitution, including {meth}`subs <sympy.core.basic.Basic.subs>`,
    {meth}`replace <sympy.core.basic.Basic.replace>`, and {meth}`xreplace
    <sympy.core.basic.Basic.xreplace>`. The methods may differ depending on
    whether they perform substitution using only strict {term}`structural
    equality` or by making use of mathematical knowledge when determining
    where a subexpression appears in an expression. Substitution is the
    standard way to treat an expression as a mathematical {term}`function` and
    evaluate it at a point.

Symbolic

    A *symbolic* representation of a mathematical object is a representation
    that is partially or completely unevaluated at runtime. It may include
    named {term}`symbolic constants <symbol>` in place of explicit numeric
    values. A symbolic representation is often contrasted with a
    {term}`numeric` one. Symbolic representations are mathematically exact, to
    contrast with numeric representations which are typically rounded so they
    can fit within a floating-point value. Symbolic {term}`expressions
    <expression>` representing mathematical objects may be aware of
    mathematical properties of these objects and be able to {term}`simplify
    <simplification>` to equivalent symbolic expressions using those
    properties. The goal of SymPy is to represent and manipulate symbolic
    expressions representing various mathematical objects.

    Some sources use the phrases "analytic solution" or "closed-form" to refer
    to the concept of "symbolic", but this terminology is not used in SymPy.
    If used in SymPy, "analytic" would refer to the property of being [an
    analytic function](https://en.wikipedia.org/wiki/Analytic_function), and
    in SymPy {term}`solve` refers only to a certain type of symbolic
    operation. "Closed-form" in SymPy would typically refer to the
    mathematical sense of the term, whereas "symbolic" would generally refer
    to the implementation detail of how a mathematical concept is implemented,
    and be in contrast with a {term}`numeric` implementation of the same
    mathematical concept.

{class}`~.Symbol`

    *{class}`~.Symbol`* is the class for symbol objects. A symbol represents a
    single mathematical variable in an expression. The {class}`~.Symbol` class
    is a subclass of {term}`Expr` and is {term}`atomic <atom>`. A `Symbol`
    contains a name, which is any string, and {term}`assumptions`. Symbols are
    typically defined with the `Symbol` constructor or the {func}`~.symbols`
    function. Two Symbols with the same name and assumptions are considered
    {term}`equal <structural equality>`. Symbols are implicitly assumed to be
    independent or constant with respect to one another. Constants, variables,
    and parameters are all represented by Symbols. The distinction is
    generally made in the way the Symbols are used in a given SymPy function.

{func}`~.sympify`

    *{func}`~.sympify()`* (not to be confused with *{term}`simplify()
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
    has different meanings in the mathematical literature, so "three-valued
    logic" is preferred. `True` and `False` work the same as in the usual
    two-valued predicate logic. `None` is an additional term that represents
    "unknown", "noncomputable", or "could be either True or False"
    (philosophically these are distinct concepts, but logically they all
    function identically). The semantics of `None` are that it absorbs other
    terms in logical operations whenever the result would differ if it were
    replaced with `True` or `False`. For example, `None OR False` is `None`,
    but `None OR True` is `True` because the predicate is `True` whether the
    `None` "really" represents a value of `True` or `False`. One must be
    careful when using the usual Python logical operators like `and`, `or` and
    `not` on three-valued logic, since `None` is false. See [the guide for
    symbolic and fuzzy booleans](booleans-guide) for more details on how to
    code with three-valued logic.

    Three-valued logic is used by the {term}`assumptions` system to represent
    assumptions that are not known. For instance, `x.is_positive` might be
    `None` if `x` could be positive or negative under its given assumptions.
    Note that the predicate logic defined by {term}`Boolean` subclasses
    represents a standard two-valued logic, not three-valued logic.

Undefined Function

    An *undefined function* is a {term}`Function <Function (class)>` that has
    no mathematical properties defined on it. It always remains
    {term}`unevaluated`, like `f(x)`. An undefined function can be created by
    passing a string name of the function to `Function`, like `f =
    Function('f')`. Undefined functions are commonly used when working with
    [ODEs](ode-docs). Undefined functions are also the easiest way to make
    {term}`symbols <symbol>` that mathematically depend on other symbols. For
    example, if `f = Function('f')` and `x = Symbol('x')`, then SymPy will
    know that `f(x)` depends on `x`, meaning, for instance, that the
    derivative `diff(f(x), x)` will not be evaluated to `0`.

Unevaluated

    An expression is *unevaluated* if the {term}`automatic simplification`
    that typically occurs when the expression is created is disabled. This is
    typically done by setting `evaluate=False`, using `with evaluate(False)`,
    or using {class}`~.UnevaluatedExpr`. While unevaluated expressions are
    supported, they can sometimes lead to surprising behavior because the
    expressions are not properly {term}`canonicalized <canonicalize>`.

    The term *unevaluated* is also sometimes used to denote the fact that an
    expression does not {term}`evaluate` to a specific value when its
    arguments are {term}`symbolic`.

{class}`zoo <sympy.core.numbers.ComplexInfinity>`

    *{class}`zoo <sympy.core.numbers.ComplexInfinity>`* represents [complex
    infinity](https://mathworld.wolfram.com/ComplexInfinity.html), i.e., the
    north pole of the [Riemann
    sphere](https://en.wikipedia.org/wiki/Riemann_sphere). The reason it is
    spelled this way is that it is "z-oo", where "z" is the symbol commonly
    used for complex variables, and {term}`oo` is the symbol SymPy uses for
    real positive infinity.

```

# CLAUDE.md - SymPy Project Guidelines

This file provides guidelines for AI assistants (like Claude) contributing to SymPy.
For complete documentation, see the [Contributing to SymPy](https://docs.sympy.org/dev/contributing/index.html) guide.

## Project Overview

SymPy is a Python library for symbolic mathematics. It requires only **mpmath** as a
hard dependency and supports Python 3.10+.

See [Dependencies](https://docs.sympy.org/dev/contributing/dependencies.html) for more information.

## Code Style

See [Development Workflow Process](https://docs.sympy.org/dev/contributing/new-contributors-guide/workflow-process.html) for complete guidelines.

### General Rules
- Use 4 spaces for indentation (no tabs)
- Maximum line length: 88 characters (configured in pyproject.toml)
- Use spaces around `+` and `-`, no spaces around `*`, `**`, and `/`
- Use exact symbolic values (`S(1)/2`) instead of floats (`0.5`) unless testing float behavior
- Run quality checks before committing: `python bin/test quality`, `flake8 sympy/`, `ruff check sympy`

### Naming Conventions
- **Classes**: CamelCase (e.g., `Basic`, `Function`, `SympifyError`)
- **Functions and methods**: snake_case (e.g., `as_real_imag`, `_eval_derivative`)
- **Private methods**: prefix with underscore (e.g., `_sympify`, `_eval_simplify`)

### Imports
- Use `from sympy import ...` rather than `import sympy` in examples and docstrings
- Use `from sympy.abc import x, y, z` for common symbols in examples
- Within SymPy source code:
  - Use relative imports within the same subpackage: `from .basic import Basic`
  - Use absolute imports for other SymPy subpackages: `from sympy.utilities.decorator import deprecated`
- For optional dependencies, use conditional imports:
  ```python
  from sympy.external import import_module
  numpy = import_module('numpy')
  ```

### Symbols and Assumptions
- Define symbols with assumptions inside test functions to prevent unintended reuse
- Use `is True`, `is False`, `is None` for assumption checks—avoid relying on truthiness
- Build expressions directly; avoid string representations of expressions

## SymPy Architecture

See [Custom Functions Guide](https://docs.sympy.org/dev/guides/custom-functions.html) for complete details.

### Core Classes
SymPy expressions inherit from these base classes:
- **Basic**: The root class for all SymPy objects
- **Expr**: For mathematical expressions (numbers, symbols, functions, operations)
- **Function**: For mathematical functions (sin, cos, custom functions)

### Immutability
All SymPy objects are immutable. Operations return new objects rather than modifying existing ones:
```python
expr = x + 1
expr.subs(x, 2)  # Returns Integer(3), does not modify expr
print(expr)       # Still x + 1
```

### The `args` Property
Every SymPy object stores its arguments in the `args` tuple. This is fundamental for traversal and manipulation:
```python
expr = sin(x) + cos(y)
expr.args        # (sin(x), cos(y))
expr.args[0].args  # (x,)
```

### `sympify()` and `S()`
- `sympify(expr)` converts Python objects to SymPy types
- `S(expr)` is shorthand for `sympify(expr)`
- Use `S(1)/2` instead of `1/2` to get a SymPy Rational, not a Python float

See [Basic Operations Tutorial](https://docs.sympy.org/dev/tutorials/intro-tutorial/basic_operations.html) for more details.

### The `_eval_*` Method Pattern
SymPy uses `_eval_*` methods to define behavior for operations on custom classes:
- `_eval_simplify()` - custom simplification logic
- `_eval_derivative(x)` - differentiation
- `_eval_evalf(prec)` - numerical evaluation
- `_eval_rewrite(rule, args, **hints)` - rewriting in terms of other functions
- `_eval_is_<assumption>()` - assumption handlers (e.g., `_eval_is_positive`)

### The `eval()` Classmethod
For `Function` subclasses, implement `eval()` to define automatic evaluation:
```python
class MyFunc(Function):
    @classmethod
    def eval(cls, x):
        if x.is_Number:
            return ...  # Return simplified result
        # Return None to leave unevaluated
```

### The `doit()` Method
Use `doit()` to explicitly evaluate unevaluated expressions:
```python
d = Derivative(sin(x), x)  # Unevaluated
d.doit()                    # Returns cos(x)
```

## Creating Custom SymPy Classes

See [Custom Functions Guide](https://docs.sympy.org/dev/guides/custom-functions.html) for complete details.

### Use `__new__`, Not `__init__`
SymPy objects are immutable, so use `__new__` instead of `__init__`:
```python
class MyExpr(Expr):
    def __new__(cls, arg):
        arg = sympify(arg)
        obj = Expr.__new__(cls, arg)
        return obj
```

**Why**: `__init__` is called after object creation and implies mutability. `__new__` creates and returns the immutable object in one step.

### Always Define `__slots__`
Every SymPy class must define `__slots__` for memory efficiency. The CI runs `slotscheck` and will fail without it:
```python
class MyExpr(Expr):
    __slots__ = ()  # Empty tuple if no new attributes

class MyExprWithAttr(Expr):
    __slots__ = ('name',)  # Tuple of new attribute names

    def __new__(cls, name, arg):
        obj = Expr.__new__(cls, arg)
        obj.name = name
        return obj
```

### Function Subclasses: Don't Override `__new__`
For `Function` subclasses, **do not** override `__new__` or `__init__`. Use `eval()` instead:
```python
class MyFunc(Function):
    __slots__ = ()

    @classmethod
    def eval(cls, x):
        if x.is_zero:
            return S.One
        # Return None to leave unevaluated
```

### Implement `_hashable_content` for Custom Attributes
If your class has attributes beyond `args`, implement `_hashable_content`:
```python
class MyExpr(Expr):
    __slots__ = ('name',)

    def __new__(cls, name, arg):
        obj = Expr.__new__(cls, arg)
        obj.name = name
        return obj

    def _hashable_content(self):
        return (self.name,) + self.args
```

### Implement `__getnewargs__` for Pickling
For classes with custom attributes, implement `__getnewargs__`:
```python
def __getnewargs__(self):
    return (self.name,) + self.args
```

### The `func` Property
Use `self.func` (not `self.__class__`) to reconstruct expressions:
```python
def some_method(self):
    # Correct: use func to create new instance
    return self.func(*new_args)

    # Wrong: don't use __class__ directly
    return self.__class__(*new_args)
```

## Assumptions System

See [Assumptions Guide](https://docs.sympy.org/dev/guides/assumptions.html) for complete details.

### Creating Symbols with Assumptions
```python
x = Symbol('x', positive=True)
x, y = symbols('x y', real=True, finite=True)
n = Symbol('n', integer=True)
```

### Common Assumptions
- **Number types**: `integer`, `rational`, `irrational`, `real`, `complex`
- **Sign predicates**: `positive`, `negative`, `zero`, `nonzero`, `nonnegative`, `nonpositive`
- **Other**: `finite`, `infinite`, `even`, `odd`, `prime`, `commutative`

### Querying Assumptions
Assumptions use three-valued logic (`True`, `False`, or `None` for unknown):
```python
x = Symbol('x', positive=True)
x.is_positive   # True
x.is_negative   # False
x.is_integer    # None (unknown)
```

### Important: Symbols with Different Assumptions are Different Objects
```python
x1 = Symbol('x', positive=True)
x2 = Symbol('x')
x1 == x2  # False - these are different symbols
```

## Printing System

See [Printing Module](https://docs.sympy.org/dev/modules/printing.html) for complete details.

### Custom Printing Methods
Implement these methods for custom classes to control their output:
- `_sympystr(self, printer)` - string representation
- `_latex(self, printer)` - LaTeX output
- `_pretty(self, printer)` - ASCII pretty printing
- `_sympyrepr(self, printer)` - repr output

### Printing Method Pattern
Always use `printer._print()` for nested expressions:
```python
def _latex(self, printer):
    # Correct: use printer._print for subexpressions
    return r'\operatorname{MyFunc}{%s}' % printer._print(self.args[0])

    # Wrong: don't use str() or create new printers
```

## Testing

See [Writing Tests](https://docs.sympy.org/dev/contributing/new-contributors-guide/writing-tests.html) for complete guidelines.

### Requirements
- All new functionality must have tests
- Bug fixes must include regression tests (tests that would fail before the fix)
- All tests must pass before merging (GitHub Actions CI runs automatically)

### Test Organization
- Tests go in `sympy/<submodule>/tests/test_<file>.py`
- Test functions must start with `test_` prefix
- Helper functions should start with underscore

### Running Tests
```bash
python bin/test              # Full test suite
python bin/test solvers      # Specific module
pytest -m 'not slow' sympy/  # Using pytest
python bin/doctest           # Doctests only
```

### Testing Mathematical Equivalence

See [Gotchas - Equals](https://docs.sympy.org/dev/tutorials/intro-tutorial/gotchas.html#equals-signs) for why `==` doesn't test mathematical equality.

The `==` operator tests **structural equality**, not mathematical equality:
```python
(x + 1)**2 == x**2 + 2*x + 1  # False (different structure)
```

Use these approaches to test mathematical equivalence:
```python
# Method 1: simplify the difference
assert simplify(expr1 - expr2) == 0

# Method 2: use equals() for numerical verification
assert expr1.equals(expr2)

# Method 3: expand both sides
assert expand(expr1) == expand(expr2)

# Method 4: trigsimp for trigonometric expressions
assert trigsimp(expr1 - expr2) == 0
```

### Writing Assertions
```python
# Basic assertion
assert function(arguments) == result

# For unchanged expressions
from sympy.core.expr import unchanged
assert unchanged(sin, 1)

# For expressions with Dummy symbols
assert factorial(k).rewrite(Product).dummy_eq(Product(_i, (_i, 1, k)))

# Testing exceptions
from sympy.testing.pytest import raises
raises(TypeError, lambda: cos(x, y))

# Testing warnings
from sympy.testing.pytest import warns
with warns(UserWarning, match=r'warning'):
    function_that_emits_warning()

# For deprecation warnings
from sympy.testing.pytest import warns_deprecated_sympy
```

### Test Markers
```python
from sympy.testing.pytest import slow, XFAIL, skip

@slow  # For tests taking >1 minute
def test_expensive_computation():
    pass

@XFAIL  # For expected failures
def test_known_failing():
    pass

# Conditional skip
if not numpy:
    skip('numpy is not installed')
```

## Docstrings

SymPy uses reStructuredText (reST) format for docstrings as specified in
[PEP 287](https://peps.python.org/pep-0287/), processed by Sphinx with the
[numpydoc](https://numpydoc.readthedocs.io/) extension.

See [SymPy Docstring Style Guide](https://docs.sympy.org/dev/contributing/docstring.html) for complete guidelines.

### Required Elements
Every public function/class docstring must include:
1. Single-sentence summary (required)
2. Examples section with doctests (required)

### Format
- Use triple double quotes: `"""docstring"""`
- Maximum 80 characters per line
- Blank line before closing quotes
- reStructuredText (reST) format with Sphinx extensions

### Section Order
1. Single-sentence summary
2. Explanation (if needed)
3. Examples
4. Parameters
5. See Also
6. References

### Example Docstring
```python
def function_name(x, y):
    """
    Short one-line summary ending with period.

    Longer explanation if needed, providing context
    and clarifying behavior.

    Examples
    ========

    >>> from sympy import symbols
    >>> from sympy.module import function_name
    >>> x, y = symbols('x y')
    >>> function_name(x, y)
    expected_output

    Parameters
    ==========

    x : type
        Description of x.
    y : type
        Description of y.

    See Also
    ========

    related_function : Brief description.

    """
```

### Doctest Guidelines
- Doctests are user-facing examples, not comprehensive tests
- Include all imports and symbol definitions (self-contained)
- Use `...` for continuation lines
- Use `<BLANKLINE>` for blank lines in output
- Use `# doctest: +SKIP` for non-executable examples

## Documentation

See [Documentation Style Guide](https://docs.sympy.org/dev/contributing/documentation-style-guide.html) for complete guidelines.

### Formatting
- Use dollar signs for LaTeX math: `$x^2$`
- Use double backticks for inline code: ``` ``code`` ```
- Use asterisks for parameter references: `*n*`
- American English spelling and punctuation
- Present tense; first-person plural ("we")

### Cross-References
```rst
:obj:`~.function_name()`           # For exported objects
:obj:`sympy.submodule.object()`    # For non-exported objects
```

### Headings (RST)
```rst
Level 1 Heading
===============

Level 2 Heading
---------------

Level 3 Heading
^^^^^^^^^^^^^^^
```

## Git Workflow

See [Development Workflow Process](https://docs.sympy.org/dev/contributing/new-contributors-guide/workflow-process.html) for complete guidelines.

### Branch Naming
- Never commit directly to `master`
- Use short, descriptive branch names: `fix-solve-bug`, `add-trig-identities`
- Avoid issue numbers in branch names

### Commit Messages
- Title: maximum 71 characters, no period at end
- Blank line between title and body
- Body: maximum 78 characters per line
- Provide context: `integrals: Improved speed of heurisch()`
- Reference issues by number in body

### Pull Request Checklist
- [ ] Code quality checks pass
- [ ] Tests added for new functionality
- [ ] Regression tests for bug fixes
- [ ] Public APIs have docstrings with doctests
- [ ] All tests pass locally
- [ ] First-time contributors: add name to `.mailmap`
- [ ] Reference related issues with "fixes #123"

## Deprecation Policy

See [Deprecation Policy](https://docs.sympy.org/dev/contributing/deprecation.html) for complete guidelines.

When deprecating functionality:
- Use `sympy_deprecation_warning()` from `sympy.utilities.exceptions`
- Maintain backwards compatibility for at least one minor version
- Document deprecation in docstring and release notes
- Test that deprecation warning is emitted using `warns_deprecated_sympy()`

## Running Quality Checks

See [Development Workflow Process](https://docs.sympy.org/dev/contributing/new-contributors-guide/workflow-process.html#code-quality-checks) for more information.

```bash
python bin/test quality       # Basic quality tests
flake8 sympy                  # Flake8 linting
ruff check sympy              # Ruff linting
mypy sympy                    # Type checking
```

The project uses flake8 and ruff. Common issues:
- Trailing whitespace (configure editor to auto-strip)
- Unused imports
- Line too long (max 88 chars)

Add `# noqa: <CODE>` for legitimate false positives.

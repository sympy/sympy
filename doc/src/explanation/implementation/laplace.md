# Laplace Transform: Design Decisions and Implementation Details

## Personal Motivation

SymPy used to do all integral transforms by calculating integrals, which may not sound surprising. And yet, when we teach students to use the Laplace transform, we teach them to never integrate, but always use the known rules and tables of known Laplace transform pairs. I found it sad that SymPy was much weaker and slower for most problems than humans using tables. 

This is why, starting in 2021, I re-wrote the Laplace transform almost completely to use rules whenever possible, and only fall back to integration when it cannot solve a problem with rules. Now, in the beginning of 2026, the state of the code is such that integration is almost never used in practice, and when it is used, it also does not provide a solution.

The advantage of such an implementation is speed (it is now massively faster) and power. The disadvantage is that every rule can have its own bugs, and that the implementation is tightly aligned with one possible strategy (here: my strategy) of applying the rules.

This strategy comes from years of using and teaching the Laplace transform in electrical engineering. It would be nice to have a similar implementation for the Fourier transform, but  I do not have enough experience with the Fourier transform to confidently re-write it. So one of the reasons to write this document is to show my design decisions in the hope they can be used for other parts of SymPy where using rules and tables might improve the algoreithms.

## Mathematical Background

The Laplace and inverse Laplace transforms are integral transforms that transform functions from one domain to another domain. A function $f(t)$ can be transformed to a function $F(s)$ by the integral

$$F(s) = \int_{0^{-}}^\infty f(t) e^{-st} \, dt\;,$$

where $s$ is a complex variable and the lower bound $0^{-}$ is to be understood as infinitely little below zero, such that if $f(t)$ contains the Dirac delta function $\delta(t)$, the whole $\delta(t)$ is to be used to calculate the result.

The inverse Laplace transform can be calculated with the Fourier-Mellin integral

$$f(t) = \frac{1}{2 \pi i} \lim_{\omega\to\infty} \int_{\sigma - i \omega}^{\sigma + i \omega} e^{st} F(s)\, ds\;,$$

where both $\sigma$ and $\omega$ are real variables. There is always a condition on $\sigma$ such that this integral converges, and it is mostly very hard to calculate.

This is why practitioners have started to collect known transform pairs and rules since the Laplace transform was introduced. One of the most complete tables was collected and published by Bateman in 1954, made available in SymPy's [Development Docs Repository](https://github.com/sympy/sympy-development-docs/tree/main/integrals/laplace). Bateman did not include anything about the Dirac delta function $\delta(t)$, these rules can be found in other textbooks and in the [Wikipedia Article on the Laplace Transform](https://en.wikipedia.org/wiki/Laplace_transform).

Please note that a function $f(t)$ calculated with the Fourier-Mellin integral must be zero for $t<0$. SymPy covers this by returning a factor `Heaviside(t)` in the result. Other computer algebra tools do not do this, but that is mathematically incorrect.

The most important property of these transforms is that, being integrals, they are linear, so if $f(t)$ and $F(s)$, and $g(t)$ and $G(s)$, are transform pairs, then so are $af(t)+bg(t)$ and $aF(s)+bG(s)$ for any $a$ and $b$. So if any function can be rewritten as a sum, the sum terms can be transformed individually. This is the linearity principle.

## Objects, Functions, and Structure of the Code

The code in `laplace.py` introduces two objects, `LaplaceTransform` and `InverseLaplaceTransform`. They can be used to write unevaluated transforms, which can be evaluated by calling the object's `doit()` method.

```py
>>> from sympy import LaplaceTransform, sin, symbols
>>> t = symbols('t', real=True)
>>> s = symbols('s')
>>> LaplaceTransform(sin(t), t, s)
LaplaceTransform(sin(t), t, s)
>>> LaplaceTransform(sin(t), t, s).doit()
1/(s**2 + 1)
```

If the intention is to calculate a transform, it is better to use the funtions `laplace_transform` and `inverse_laplace_transform`:

```py
>>> from sympy import laplace_transform, sin, symbols
>>> t = symbols('t', real=True)
>>> s = symbols('s')
>>> laplace_transform(sin(t), t, s)
(1/(s**2 + 1), 0, True)
```

This also returns the convergence plane and possible additional conditions. Calling `laplace_transform` with the additional argument `noconds=True` will return `1/(s**2 + 1)`, while using the method `.doit(noconds=False)` further above will return `(1/(s**2 + 1), 0, True)`.  (The inverse Laplace transform does not output conditions.)

The actual work is done for both object and function by the functions `_laplace_transform` and `_inverse_laplace_transform`.

All relevant functions have a debug wrapper, which makes it possible to see what the algorithm is doing, using indents to show recursion levels. Debugging is enabled with

```py
import sympy
sympy.SYMPY_DEBUG=True
```

If the example above is executed with debugging enabled, the following output is generated:

```
[LT doit] (5*sin(t), t, s)

------------------------------------------------------------------------------
-LT- _laplace_transform(5*sin(t), t, s)
-LT-   _laplace_apply_simple_rules(sin(t), t, s)
-LT-     _laplace_build_rules is building rules
-LT-     _laplace_deep_collect(sin(_t), _t)
-LT-       _laplace_deep_collect(_t, _t)
-LT-       ---> _t
-LT-     ---> sin(_t)
-LT-   ---> (1/(s**2 + 1), 0, True)
-LT- ---> (5/(s**2 + 1), 0, True)
------------------------------------------------------------------------------

(5/(s**2 + 1), 0, True)
```

This first shows that the algorithm has decided to transform `5*sin(t)` from the `t` to the `s` domain. 

Reading from outer levels to inner levels: `laplace_transform` passes `(5*sin(t), t, s)` to `_laplace_transform`, which returns `(5/(s**2 + 1), 0, True)`.

It does the calculation by passing `(sin(t), t, s)` to `_laplace_apply_simple_rules`, which returns `(1/(s**2 + 1), 0, True)`.

`_laplace_apply_simple_rules` first has to build the rules with `_laplace_build_rules`. It then uses `_laplace_deep_collect` to recursively simplify `sin(_t)`, which, in this case, is of course pointless, and then applies one of its many rules.

It is very advisable, when reading the texts below, to look at several such examples by yourself.

## Implementation of the Laplace Transform

The front-end function `laplace_transform` is able to do a Laplace transform of a matrix by doing element-wise Laplace transforms and collecting all conditions. For a non-matrix element, it simply hands the arguments to the `LaplaceTransform` object and executes `.doit(noconds=False, simplify=_simplify)`, which we will describe below. Then it returns either only the transformed function, or also the conditions, depending on its own setting of the parameter `noconds`. This was a conscious design decision: it is much easier (and not computationally costly) to write all functions such that they always return conditions than to propagate `noconds` throught the recursion.

The object `LaplaceTransform` uses the method `doit()` to run `_laplace_transform`, giving it the function, both variables, and only one key, `simplify`. It decides by itself, basing on its own setting of `noconds`, whether it should return conditions or not.

### Core of the algorithm: `_laplace_transform`

The core of the algorithm is `_laplace_transform`. This function first splits an expression like `a*f + b*g` into sum terms `a*f` and `b*g` to transform them inependently. It also attempts to rewrite piecewise functions using `Heaviside` with the funstion `_piecewise_to_heaviside`. It will return the sum of all the transforms, the maximum of all convergence planes, and the union of all conditions.

The main part of the algorithm first removes any `Heaviside(t)` factor and then does the following step by step:
- try `_laplace_apply_simple_rules`,
- try `_laplace_apply_prog_rules`,
- try `_laplace_expand`,
- check whether undefined functions are present in the term to transform; if yes, return a `LaplaceTransform` object,
- try `_laplace_transform_integration`.

All these functions (and also the functions described in other sections) either return a result, if they have one, or they return `None`; the walrus operator `:=` is used in many places to do this, e.g.,

```py
if (
        (r := _laplace_apply_simple_rules(ft, t_, s_))
        is not None or
        (r := _laplace_apply_prog_rules(ft, t_, s_))
        is not None or
        (r := _laplace_expand(ft, t_, s_)) is not None):
    pass
elif ...
``` 

If any of them succeeds, the variable `r` will contain the result, if not, the code after `elif` will be executed. If any  result is found, it is returned; if not, a `LaplaceTransform` object is returned. The algorithm is recursive by nature, several of the functions that are executed can again execute `laplace_transform`. Endless loops are prevented only by design, there is no check for the recursion level depth.

There is one major disadvantage in this algorithm: there is a very small number of cases in Bateman's tables where $f(t)+g(t)$ can be transformed but $f(t)$ and $g(t)$ cannot, because only the integral of the sum converges. Splitting sums first makes it hard to deal with these cases. To our knowledge, these cases are not practically relevant, but if they ever need to be covered, it needs to be done by an additional function executed before sums are split into their terms.

### Simple rules: `_laplace_apply_simple_rules`

This function checks whether any rules from tables apply, by using the `.match()` method. For this, it builds a list of rules with `_laplace_build_rules`. This function is cached, so in one session the rules are only built once.

Rules are written as follows:

```py
    (Heaviside(a*t-b), exp(-s*b/a)/s,
     And(a > 0, b > 0), S.Zero, dco),  # 4.4.1
```

This rule will try to match the input with `Heaviside(a*t-b)` and return `exp(-s*b/a)/s` if it finds a match. It will only apply if `And(a > 0, b > 0)` is true, and will return `S.Zero` as the convergence plane. Before the match is attempted, the function `dco` is applied, which is short for `_laplace_deep_collect`. That function traverses through the input expression and collect expressions for `t` to bring them into the form that can be matched. Example:

```py
>>> from sympy import Heaviside, laplace_transform, symbols
>>> s = symbols('s')
>>> t = symbols('t', real=True)
>>> a, b = symbols('a, b', positive=True)
>>> laplace_transform(Heaviside(a*(t - b)), t, s)
(exp(-b*s)/s, 0, True)
```

The debugging output is:

```
[LT doit] (Heaviside(a*(-b + t)), t, s)

------------------------------------------------------------------------------
-LT- _laplace_transform(Heaviside(a*(-b + t)), t, s)
-LT-   _laplace_apply_simple_rules(Heaviside(a*(-b + t)), t, s)
-LT-     _laplace_build_rules is building rules
-LT-     _laplace_deep_collect(Heaviside(a*(_t - b)), _t)
-LT-       _laplace_deep_collect(a*(_t - b), _t)
-LT-       ---> _t*a - a*b
-LT-       _laplace_deep_collect(1/2, _t)
-LT-       ---> 1/2
-LT-     ---> Heaviside(_t*a - a*b)
-LT-   ---> (exp(-b*s)/s, 0, True)
-LT- ---> (exp(-b*s)/s, 0, True)
------------------------------------------------------------------------------
```

Here the value of `dco` becomes apparent; the rule would not match `Heaviside(a*(-b + t))`, but it matches `Heaviside(_t*a - a*b)`.

The [Development Docs Repository](https://github.com/sympy/sympy-development-docs/tree/main/integrals/laplace) contains information on which table rules are already implemented. Several rules are written but commemnted out, because they are not yet tested well enough or cause problems. Extending the rules table would be a good point for fresh contributors to help.


### Algorithmic rules: `_laplace_apply_prog_rules`

Next, algorithmic rules that cannot be resolved by simple pattern matching are applied in the following order:

- `_laplace_rule_heaviside` can also deal with products of time-shifted `Heaviside()` functions and unknown functions `f(t)`, and with rectangular windows on functions.
- `_laplace_rule_delta` deals with `DiracDelta` factors in the time domain by using the Delta distribution's masking property.
- `_laplace_rule_timescale` deals with functions where the time variable has a factor, like `f(a*t)`.
- `_laplace_rule_exp` removes exponential factors `exp(a*t)`, expressing them as a frequency shift `s + a` in the result.
- `_laplace_rule_trig` is a relatively complicated algroithm that converts trigonometric functions to exponential functions, transforms everything, and then makes nice expressions in `s` by collecting complex conjugate poles cominng from `sin` and `cos` and symmetric real poles coming from `sinh` and `cosh`.
- `_laplace_rule_diff` looks for derivatives in the `t` domain to replace them with factors of `s` in the result.
- `_laplace_rule_sdiff` looks for factors `t` in the `t` domain to replace them with derivatives in the `s` domain.

Many of the above functions call `_laplace_transform` recursively, but with `simplify=False`, knowing that simplification is done at the top of the recursion.

### Further expansion: `_laplace_expand`

If all of the above fail, `_laplace_expand` by using, in this order,
- `expand(f, deep=False)`,
- `expand_mul(f)`,
- `expand(f)`,
- `expand(expand_trig(f))`.
If any of these return an `Add` object (i.e., a sum), it calls `_laplace_transform` recursively. If `expand(f)` canges the input (this is one of the design decisions that prevents infinite recursion), it also calls `_laplace_transform` recursively.

### Calculate an integral: `_laplace_transform_integration`

Finally, the algorithm uses `integrate(f*exp(-s*t), (t, S.Zero, S.Infinity))` to attempt a solution. The main part of this code is to re-write the conditions returned by `integrate`; that code pre-dates the re-write of the Laplace transform and depends closely on the implementation of `integrate`.

### Long debugging example

The following example shows the recursive nature of the algorithm nicely:

```py
>>> from sympy import cos, laplace_transform, sinh, symbols
>>> import sympy
>>> sympy.SYMPY_DEBUG = True
>>> s = symbols('s')
>>> t = symbols('t', real=True)
>>> a, b, c = symbols('a, b, c', real=True)
>>> laplace_transform(cos(a*t)*sinh(b*t)*sinh(c*t), t, s)
((-s**3/2 + s*(-a**2/2 + (-b + c)**2/2))/(a**4 + 2*a**2*(-b + c)**2 + s**4 + s**2*(2*a**2 - 2*(-b + c)**2) + (-b + c)**4) + (s**3/2 + s*(a**2/2 - (-b - c)**2/2))/(a**4 + 2*a**2*(-b - c)**2 + s**4 + s**2*(2*a**2 - 2*(-b - c)**2) + (-b - c)**4), Max(Abs(b - c), Abs(b + c)), True)
```

The debugging output shows how many functions attempt to solve it, and `_laplace_rule_trig` finally succeeds.

```
[LT doit] (cos(a*t)*sinh(b*t)*sinh(c*t), t, s)

------------------------------------------------------------------------------
-LT- _laplace_transform(cos(a*t)*sinh(b*t)*sinh(c*t), t, s)
-LT-   _laplace_apply_simple_rules(cos(a*t)*sinh(b*t)*sinh(c*t), t, s)
-LT-     _laplace_deep_collect(cos(_t*a)*sinh(_t*b)*sinh(_t*c), _t)
-LT-       _laplace_deep_collect(cos(_t*a), _t)
-LT-         _laplace_deep_collect(_t*a, _t)
-LT-         ---> _t*a
-LT-       ---> cos(_t*a)
-LT-       _laplace_deep_collect(sinh(_t*b), _t)
-LT-         _laplace_deep_collect(_t*b, _t)
-LT-         ---> _t*b
-LT-       ---> sinh(_t*b)
-LT-       _laplace_deep_collect(sinh(_t*c), _t)
-LT-         _laplace_deep_collect(_t*c, _t)
-LT-         ---> _t*c
-LT-       ---> sinh(_t*c)
-LT-     ---> cos(_t*a)*sinh(_t*b)*sinh(_t*c)
-LT-   ---> None
-LT-   _laplace_apply_prog_rules(cos(a*t)*sinh(b*t)*sinh(c*t), t, s)
-LT-     _laplace_rule_heaviside(cos(a*t)*sinh(b*t)*sinh(c*t), t, s)
-LT-     ---> None
-LT-     _laplace_rule_delta(cos(a*t)*sinh(b*t)*sinh(c*t), t, s)
-LT-     ---> None
-LT-     _laplace_rule_timescale(cos(a*t)*sinh(b*t)*sinh(c*t), t, s)
-LT-     ---> None
-LT-     _laplace_rule_exp(cos(a*t)*sinh(b*t)*sinh(c*t), t, s)
-LT-     ---> None
-LT-     _laplace_rule_trig(cos(a*t)*sinh(b*t)*sinh(c*t), t, s)
-LT-       _laplace_trig_split(cos(_t*a)*sinh(_t*b)*sinh(_t*c),)
-LT-       ---> (cos(_t*a)*sinh(_t*b)*sinh(_t*c), 1)
-LT-       _laplace_trig_expsum(cos(_t*a)*sinh(_t*b)*sinh(_t*c), _t)
-LT-         _laplace_deep_collect(-exp(_t*I*a + _t*b - _t*c)/8, _t)
-LT-           _laplace_deep_collect(-1/8, _t)
-LT-           ---> -1/8
-LT-           _laplace_deep_collect(exp(_t*I*a + _t*b - _t*c), _t)
-LT-             _laplace_deep_collect(_t*I*a + _t*b - _t*c, _t)
-LT-             ---> _t*(I*a + b - c)
-LT-           ---> exp(_t*(I*a + b - c))
-LT-         ---> -exp(_t*(I*a + b - c))/8
-LT-         _laplace_deep_collect(-exp(-_t*I*a + _t*b - _t*c)/8, _t)
-LT-           _laplace_deep_collect(-1/8, _t)
-LT-           ---> -1/8
-LT-           _laplace_deep_collect(exp(-_t*I*a + _t*b - _t*c), _t)
-LT-             _laplace_deep_collect(-_t*I*a + _t*b - _t*c, _t)
-LT-             ---> _t*(-I*a + b - c)
-LT-           ---> exp(_t*(-I*a + b - c))
-LT-         ---> -exp(_t*(-I*a + b - c))/8
-LT-         _laplace_deep_collect(-exp(_t*I*a - _t*b + _t*c)/8, _t)
-LT-           _laplace_deep_collect(-1/8, _t)
-LT-           ---> -1/8
-LT-           _laplace_deep_collect(exp(_t*I*a - _t*b + _t*c), _t)
-LT-             _laplace_deep_collect(_t*I*a - _t*b + _t*c, _t)
-LT-             ---> _t*(I*a - b + c)
-LT-           ---> exp(_t*(I*a - b + c))
-LT-         ---> -exp(_t*(I*a - b + c))/8
-LT-         _laplace_deep_collect(-exp(-_t*I*a - _t*b + _t*c)/8, _t)
-LT-           _laplace_deep_collect(-1/8, _t)
-LT-           ---> -1/8
-LT-           _laplace_deep_collect(exp(-_t*I*a - _t*b + _t*c), _t)
-LT-             _laplace_deep_collect(-_t*I*a - _t*b + _t*c, _t)
-LT-             ---> _t*(-I*a - b + c)
-LT-           ---> exp(_t*(-I*a - b + c))
-LT-         ---> -exp(_t*(-I*a - b + c))/8
-LT-         _laplace_deep_collect(exp(_t*I*a + _t*b + _t*c)/8, _t)
-LT-           _laplace_deep_collect(1/8, _t)
-LT-           ---> 1/8
-LT-           _laplace_deep_collect(exp(_t*I*a + _t*b + _t*c), _t)
-LT-             _laplace_deep_collect(_t*I*a + _t*b + _t*c, _t)
-LT-             ---> _t*(I*a + b + c)
-LT-           ---> exp(_t*(I*a + b + c))
-LT-         ---> exp(_t*(I*a + b + c))/8
-LT-         _laplace_deep_collect(exp(-_t*I*a + _t*b + _t*c)/8, _t)
-LT-           _laplace_deep_collect(1/8, _t)
-LT-           ---> 1/8
-LT-           _laplace_deep_collect(exp(-_t*I*a + _t*b + _t*c), _t)
-LT-             _laplace_deep_collect(-_t*I*a + _t*b + _t*c, _t)
-LT-             ---> _t*(-I*a + b + c)
-LT-           ---> exp(_t*(-I*a + b + c))
-LT-         ---> exp(_t*(-I*a + b + c))/8
-LT-         _laplace_deep_collect(exp(_t*I*a - _t*b - _t*c)/8, _t)
-LT-           _laplace_deep_collect(1/8, _t)
-LT-           ---> 1/8
-LT-           _laplace_deep_collect(exp(_t*I*a - _t*b - _t*c), _t)
-LT-             _laplace_deep_collect(_t*I*a - _t*b - _t*c, _t)
-LT-             ---> _t*(I*a - b - c)
-LT-           ---> exp(_t*(I*a - b - c))
-LT-         ---> exp(_t*(I*a - b - c))/8
-LT-         _laplace_deep_collect(exp(-_t*I*a - _t*b - _t*c)/8, _t)
-LT-           _laplace_deep_collect(1/8, _t)
-LT-           ---> 1/8
-LT-           _laplace_deep_collect(exp(-_t*I*a - _t*b - _t*c), _t)
-LT-             _laplace_deep_collect(-_t*I*a - _t*b - _t*c, _t)
-LT-             ---> _t*(-I*a - b - c)
-LT-           ---> exp(_t*(-I*a - b - c))
-LT-         ---> exp(_t*(-I*a - b - c))/8
-LT-       ---> ([{'k': -1/8, 'a': I*a + b - c, re: b - c, im: a}, {'k': -1/8, 'a': -I*a + b - c, re: b - c, im: -a}, {'k': -1/8, 'a': I*a - b + c, re: -b + c, im: a}, {'k': -1/8, 'a': -I*a - b + c, re: -b + c, im: -a}, {'k': 1/8, 'a': I*a + b + c, re: b + c, im: a}, {'k': 1/8, 'a': -I*a + b + c, re: b + c, im: -a}, {'k': 1/8, 'a': I*a - b - c, re: -b - c, im: a}, {'k': 1/8, 'a': -I*a - b - c, re: -b - c, im: -a}], [])
-LT-       _laplace_trig_ltex([{'k': -1/8, 'a': I*a + b - c, re: b - c, im: a}, {'k': -1/8, 'a': -I*a + b - c, re: b - c, im: -a}, {'k': -1/8, 'a': I*a - b + c, re: -b + c, im: a}, {'k': -1/8, 'a': -I*a - b + c, re: -b + c, im: -a}, {'k': 1/8, 'a': I*a + b + c, re: b + c, im: a}, {'k': 1/8, 'a': -I*a + b + c, re: b + c, im: -a}, {'k': 1/8, 'a': I*a - b - c, re: -b - c, im: a}, {'k': 1/8, 'a': -I*a - b - c, re: -b - c, im: -a}], _t, s)
-LT-       ---> ((-s**3/2 + s*(-a**2/2 + (-b + c)**2/2))/(a**4 + 2*a**2*(-b + c)**2 + s**4 + s**2*(2*a**2 - 2*(-b + c)**2) + (-b + c)**4) + (s**3/2 + s*(a**2/2 - (-b - c)**2/2))/(a**4 + 2*a**2*(-b - c)**2 + s**4 + s**2*(2*a**2 - 2*(-b - c)**2) + (-b - c)**4), Max(Abs(b - c), Abs(b + c)))
-LT-     ---> ((-s**3/2 + s*(-a**2/2 + (-b + c)**2/2))/(a**4 + 2*a**2*(-b + c)**2 + s**4 + s**2*(2*a**2 - 2*(-b + c)**2) + (-b + c)**4) + (s**3/2 + s*(a**2/2 - (-b - c)**2/2))/(a**4 + 2*a**2*(-b - c)**2 + s**4 + s**2*(2*a**2 - 2*(-b - c)**2) + (-b - c)**4), Max(Abs(b - c), Abs(b + c)), True)
-LT-   ---> ((-s**3/2 + s*(-a**2/2 + (-b + c)**2/2))/(a**4 + 2*a**2*(-b + c)**2 + s**4 + s**2*(2*a**2 - 2*(-b + c)**2) + (-b + c)**4) + (s**3/2 + s*(a**2/2 - (-b - c)**2/2))/(a**4 + 2*a**2*(-b - c)**2 + s**4 + s**2*(2*a**2 - 2*(-b - c)**2) + (-b - c)**4), Max(Abs(b - c), Abs(b + c)), True)
-LT- ---> ((-s**3/2 + s*(-a**2/2 + (-b + c)**2/2))/(a**4 + 2*a**2*(-b + c)**2 + s**4 + s**2*(2*a**2 - 2*(-b + c)**2) + (-b + c)**4) + (s**3/2 + s*(a**2/2 - (-b - c)**2/2))/(a**4 + 2*a**2*(-b - c)**2 + s**4 + s**2*(2*a**2 - 2*(-b - c)**2) + (-b - c)**4), Max(Abs(b - c), Abs(b + c)), True)
------------------------------------------------------------------------------
```


## Implementation of the Inverse Transform

to be written ...

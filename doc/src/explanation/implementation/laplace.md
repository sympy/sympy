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

The inverse Laplace transform can be calculated with

$$f(t) = \frac{1}{2 \pi i} \lim_{\omega\to\infty} \int_{\sigma - i \omega}^{\sigma + i \omega} e^{st} F(s)\, ds\;,$$

where both $\sigma$ and $\omega$ are real variables. There is always a condition on $\sigma$ such that this integral converges, and it is mostly very hard to calculate.

This is why practitioners have started to collect known transform pairs and rules since the Laplace transform was introduced. One of the most complete tables was collected and published by Bateman in 1954, made available in SymPy's [Development Docs Repository](https://github.com/sympy/sympy-development-docs/tree/main/integrals/laplace). Bateman did not include anything about the Dirac delta function $\delta(t)$, these rules can be found in other textbooks and in the [Wikipedia Article on the Laplace Transform](https://en.wikipedia.org/wiki/Laplace_transform).

The most important property of these transforms is that, being integrals, they are linear, so if $f(t)$ and $F(s)$, and $g(t)$ and $G(s)$, are transform pairs, then so are $af(t)+bg(t)$ and $aF(s)+bG(s)$ for any $a$ and $b$. So if any function can be rewritten as a sum, the sum terms can be transformed individually.

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

```py>>> from sympy import laplace_transform, sin, symbols
>>> t = symbols('t', real=True)
>>> s = symbols('s')
>>> laplace_transform(sin(t), t, s)
(1/(s**2 + 1), 0, True)
```

This also returns the convergence plane and possible additional conditions. Calling `laplace_transform` with the additional argument `noconds=True` will return `1/(s**2 + 1)`, while using the method `.doit(noconds=False)` further above will return `(1/(s**2 + 1), 0, True)`.  (The inverse Laplace transform does not output conditions.)

The actual work is done for both object and function by the functions `_laplace_transform` and `_inverse_laplace_transform`.

All relevant functions have a debug wrapper, which makes it possible to see what the algorithm is doing, using indents to show recursion levels. For example,

```py
>>> from sympy import laplace_transform, sin, symbols
>>> import sympy
>>> sympy.SYMPY_DEBUG=True
>>> t = symbols('t', real=True)
>>> s = symbols('s')
>>> laplace_transform(5*sin(t), t, s)
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
>>>```

shows that the algorithm has decided to transform `5*sin(t)` from `t` to `s`. 

Reading from outer levels to inner levels: `laplace_transform` passes `(5*sin(t), t, s)` to `_laplace_transform`, which returns `(5/(s**2 + 1), 0, True)`.

It does the calculation by passing `(sin(t), t, s)` to `_laplace_apply_simple_rules`, which returns `(1/(s**2 + 1), 0, True)`.

`_laplace_apply_simple_rules` first has to build the rules with `_laplace_build_rules`. It then uses `_laplace_deep_collect` to recursively simplify `sin(_t)`, which, in this case, is of course pointless, and then applies one of its many rules.

It is very advisable, when reading the texts below, to look at several such examples by yourself.

## Implementation of the Laplace Transform

to be written ...

## Implementation of the Inverse Transform

to be written ...

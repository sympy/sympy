"""Reproducible timing benchmark for ``sympy.cse`` on large inputs.

This is NOT a pytest test (it is slow); run it directly::

    python sympy/simplify/tests/bench_cse_scaling.py

It builds a meta-GGA / PBE-correlation-like scalar expression in 8
variables, differentiates it to a requested order to obtain expression
sets with up to millions of nodes, and times ``cse`` for both
``optimizations=[]`` (the default, fast path) and ``optimizations='basic'``
(the factor_terms preprocessing path).

The benchmark was used to evaluate two performance changes:

* memoizing ``factor_terms`` over shared (DAG) subexpressions
  (``sympy/core/exprtools.py``), and
* backing ``OrderedSet`` with a plain ``dict`` instead of ``OrderedDict``
  (``sympy/core/containers.py``).

Neither changes the output of ``cse``.
"""

import time

from sympy import symbols, sqrt, log, exp, cbrt, Rational, cse, count_ops


def build_functional():
    """A deeply nested rational+power expression in 8 variables."""
    ra, rb, saa, sab, sbb, ta, tb, _lap = symbols(
        'rho_a rho_b sigma_aa sigma_ab sigma_bb tau_a tau_b lapl',
        positive=True)
    rho = ra + rb
    zeta = (ra - rb)/rho
    sig = saa + 2*sab + sbb
    rs = cbrt(Rational(3, 4)/Rational(314159, 100000)/rho)
    t2 = sig/rho**Rational(7, 3)
    phi = (cbrt(1 + zeta) + cbrt(1 - zeta))/2
    A = 1/(exp(-1/phi**3) - 1 + Rational(1, 100))
    H = phi**3*log(1 + t2/(1 + A*t2)*(1 + A*t2)/(1 + A*t2 + A**2*t2**2))
    tau = ta + tb
    tw = sig/(8*rho)
    alpha = (tau - tw)/(rho**Rational(5, 3) + Rational(1, 1000))
    fa = 1/(1 + alpha**2)
    eps = -rs/(1 + rs + rs**2) + H*fa + log(1 + rho)*sqrt(1 + t2)
    return rho*eps, (ra, rb, saa, sab, sbb, ta, tb)


def make_exprs(order, nvars):
    f, dvars = build_functional()
    dvars = dvars[:nvars]
    exprs = [f]
    cur = [f]
    for _ in range(order):
        nxt = [e.diff(v) for e in cur for v in dvars]
        exprs.extend(nxt)
        cur = nxt
    return exprs


def bench(order, nvars, do_basic=True):
    exprs = make_exprs(order, nvars)
    nops = sum(count_ops(e) for e in exprs)
    print(f"order={order} nvars={nvars} nexprs={len(exprs)} ops={nops}")
    t0 = time.perf_counter()
    cse(exprs, optimizations=[])
    print(f"    optimizations=[]      {time.perf_counter()-t0:8.3f} s")
    if do_basic:
        t0 = time.perf_counter()
        cse(exprs, optimizations='basic')
        print(f"    optimizations='basic' {time.perf_counter()-t0:8.3f} s")


if __name__ == '__main__':
    # 'basic' becomes very slow on the larger sizes on stock sympy, so only
    # run it on the smaller inputs.
    bench(1, 4, do_basic=True)
    bench(2, 3, do_basic=True)
    bench(2, 4, do_basic=False)
    bench(3, 2, do_basic=False)

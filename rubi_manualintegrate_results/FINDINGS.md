# manualintegrate vs. the Rubi test suite — key findings

`sympy.integrals.manualintegrate` was run over all **67,883** integrals of
the Rubi (MathematicaSyntaxTestSuite) test suite with a 30 s per-case
timeout; every returned antiderivative was verified by differentiation
(`simplify`, with numeric spot checks as fallback).  Full breakdown in
[REPORT.md](REPORT.md); per-case data in the `*.jsonl` files.

## Headline numbers

| outcome | cases | share |
|---|---|---|
| solved, verified correct | 14,471 | 21.3% |
| solved but wrong | 3 | ~0% |
| not solved (DontKnow / unevaluated Integral) | 48,436 | 71.4% |
| timeout (30 s) | 4,826 | 7.1% |
| exceptions | 147 | 0.2% |

## Findings

* **Reliability is exceptional**: of 14,474 antiderivatives produced,
  14,471 verified correct — 99.98%.  The 3 wrong ones are a
  `TrigSubstitutionRule` branch-cut bug on `sqrt(1/(x**2 - 1))` (the
  returned Piecewise has no branch for `|x| > 1` and the wrong branch of
  the square root inside `(-1, 1)`) plus two parametric-exponent
  Piecewise cases:
  * `x**p*(a*x**n + b*x**(13*n + p + 1))**12` (RewriteRule)
  * `x**(-15*n - 1)*(a + b*x**n)**8` (AlternativeRule)
* **Coverage is the weakness, and it's concentrated**: algebraic
  functions solve at 48.6%, but trig (3.2%), hyperbolic (3.7%) and
  inverse-hyperbolic (1.0%) families are almost entirely `DontKnowRule` —
  Rubi's parametric `(a + b*sin(x))**m`-style families have no
  manualintegrate counterpart.
* **Difficulty scaling is sane**: the solved rate falls monotonically
  with Rubi's own step count after 2 steps — 52% at 2 steps down to 8%
  at 13+.  (Oddly, 1-step cases solve *less* often than 2-step ones —
  many are special-function forms Rubi handles with a single table
  rule.)
* **Timeouts cluster in logarithms** (818 of 3,036 cases, 27%) — likely
  a good hunting ground for performance fixes.
* The report also breaks down results per file, per top-level rule
  (`RewriteRule`, `AlternativeRule` and `URule` dominate the solved
  set), verification method, exception classes and the slowest solved
  cases.

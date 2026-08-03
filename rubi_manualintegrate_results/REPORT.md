# manualintegrate on the Rubi integration test suite

Results of running `sympy.integrals.manualintegrate` over the 67883 integrals of the Rubi (MathematicaSyntaxTestSuite) test suite, with a 30 s per-case timeout and correctness verified by differentiating the returned antiderivative (`simplify`, then numeric spot checks).

## Overall outcome

| outcome | cases | share |
|---|---|---|
| solved, verified correct | 14471 | 21.3% |
| solved, WRONG antiderivative | 3 | 0.0% |
| not solved (DontKnowRule / unevaluated Integral) | 48436 | 71.4% |
| timeout | 4826 | 7.1% |
| exception from integral_steps/manualintegrate | 147 | 0.2% |
| **any antiderivative found** | **14474** | **21.3%** |

Of the 14474 solved integrals, 14471 (100.0%) were verified correct, 0 could not be verified either way, and 3 were wrong.

## By test-suite section

| section | cases | solved | verified | wrong | dont_know | timeout | error | solved % |
|---|---|---|---|---|---|---|---|---|
| t_1_algebraic_functions | 25327 | 12299 | 12296 | 3 | 10046 | 2926 | 56 | 48.6% |
| t_2_exponentials | 963 | 226 | 226 | 0 | 719 | 16 | 2 | 23.5% |
| t_3_logarithms | 3036 | 453 | 453 | 0 | 1741 | 818 | 24 | 14.9% |
| t_4_trig_functions | 22221 | 703 | 703 | 0 | 21142 | 376 | 0 | 3.2% |
| t_5_inverse_trig_functions | 3943 | 411 | 411 | 0 | 2992 | 505 | 35 | 10.4% |
| t_6_hyperbolic_functions | 5053 | 187 | 187 | 0 | 4817 | 49 | 0 | 3.7% |
| t_7_inverse_hyperbolic_functions | 5488 | 55 | 55 | 0 | 5331 | 74 | 28 | 1.0% |
| t_8_special_functions | 1852 | 140 | 140 | 0 | 1648 | 62 | 2 | 7.6% |
| **TOTAL** | **67883** | **14474** | **14471** | **3** | **48436** | **4826** | **147** | **21.3%** |

## Solved rate by Rubi step count

Number of integration steps Rubi itself needed, as a difficulty proxy.

| Rubi steps | cases | solved | verified | wrong | dont_know | timeout | error | solved % |
|---|---|---|---|---|---|---|---|---|
| <=0 | 2549 | 9 | 9 | 0 | 2432 | 104 | 4 | 0.4% |
| 1 | 3008 | 869 | 869 | 0 | 2073 | 46 | 20 | 28.9% |
| 2 | 8841 | 4601 | 4600 | 1 | 4002 | 211 | 27 | 52.0% |
| 3-4 | 18678 | 4902 | 4900 | 2 | 12873 | 886 | 17 | 26.2% |
| 5-7 | 18633 | 2583 | 2583 | 0 | 14385 | 1627 | 38 | 13.9% |
| 8-12 | 10657 | 1068 | 1068 | 0 | 8383 | 1183 | 23 | 10.0% |
| 13+ | 5517 | 442 | 442 | 0 | 4288 | 769 | 18 | 8.0% |
| **TOTAL** | **67883** | **14474** | **14471** | **3** | **48436** | **4826** | **147** | **21.3%** |

## Top-level rule of solved integrals

| rule | solved | verified | wrong |
|---|---|---|---|
| RewriteRule | 5274 | 5273 | 1 |
| AlternativeRule | 3495 | 3494 | 1 |
| URule | 2510 | 2510 | 0 |
| PiecewiseRule | 2231 | 2231 | 0 |
| PartsRule | 233 | 233 | 0 |
| AddRule | 188 | 188 | 0 |
| NestedPowRule | 161 | 161 | 0 |
| RatintRule | 136 | 136 | 0 |
| ConstantTimesRule | 64 | 64 | 0 |
| SqrtQuadraticDenomRule | 49 | 49 | 0 |
| SqrtQuadraticRule | 28 | 28 | 0 |
| PowerRule | 24 | 24 | 0 |
| TrigSubstitutionRule | 19 | 18 | 1 |
| ConstantRule | 11 | 11 | 0 |
| CompleteSquareRule | 11 | 11 | 0 |
| CyclicPartsRule | 10 | 10 | 0 |
| ReciprocalSqrtQuadraticRule | 9 | 9 | 0 |
| UpperGammaRule | 4 | 4 | 0 |
| FresnelCRule | 4 | 4 | 0 |
| FresnelSRule | 4 | 4 | 0 |
| PolylogRule | 3 | 3 | 0 |
| ErfRule | 1 | 1 | 0 |
| ReciprocalRule | 1 | 1 | 0 |
| ArcsinRule | 1 | 1 | 0 |
| ArcsinhRule | 1 | 1 | 0 |
| LiRule | 1 | 1 | 0 |
| EiRule | 1 | 1 | 0 |

## Verification method of correct results

| method | cases |
|---|---|
| simplify | 7567 |
| numeric | 6357 |
| syntactic | 558 |

## Wrong antiderivatives (3)

* `x**p*(a*x**n + b*x**(13*n + p + 1))**12`  
  source: 1 Algebraic functions/1.1 Binomial products/1.1.4 Improper/1.1.4.2 (c x)^m (a x^j+b x^n)^p.m #346, rule `RewriteRule`
* `x**(-15*n - 1)*(a + b*x**n)**8`  
  source: 1 Algebraic functions/1.1 Binomial products/1.1.3 General/1.1.3.2 (c x)^m (a+b x^n)^p.m #2520, rule `AlternativeRule`
* `sqrt(1/(x**2 - 1))`  
  source: 1 Algebraic functions/1.3 Miscellaneous/1.3.2 Algebraic functions.m #680, rule `TrigSubstitutionRule`

## Exceptions raised

| exception | cases |
|---|---|
| TypeError | 80 |
| PolynomialError | 36 |
| PolificationFailed | 21 |
| CoercionFailed | 6 |
| PolynomialDivisionFailed | 4 |

## Timing

Total CPU-case time 89.4 h; per case: mean 4.74 s, median 1.37 s, p95 30.07 s, max 33.1 s.

Slowest solved cases:

* 30.1 s — `(A + B*x)*(d + e*x)**(S(5)/2)/(a**S(2) + S(2)*a*b*x + b**S(2)*x**S(2))**S(2)` (1 Algebraic functions/1.2 Trinomial products/1.2.1 Quadratic/1.2.1.3 (d+e x)^m (f+g x) (a+b x+c x^2)^p.m #1816)
* 30.1 s — `(b + 2*c*x)*(d + e*x)**3/(a + b*x + c*x**2)**2` (1 Algebraic functions/1.2 Trinomial products/1.2.1 Quadratic/1.2.1.3 (d+e x)^m (f+g x) (a+b x+c x^2)^p.m #1531)
* 29.9 s — `(c + d*x)/(a + b*(c + d*x)**3)**3` (1 Algebraic functions/1.1 Binomial products/1.1.3 General/1.1.3.2 (c x)^m (a+b x^n)^p.m #2779)
* 29.9 s — `x**7*(d + e*x**2 + f*x**4)/(a + b*x**2 + c*x**4)**2` (1 Algebraic functions/1.2 Trinomial products/1.2.2 Quartic/1.2.2.6 P(x) (d x)^m (a+b x^2+c x^4)^p.m #60)
* 29.8 s — `(a + b/x)**(S(3)/2)/sqrt(x)` (1 Algebraic functions/1.1 Binomial products/1.1.3 General/1.1.3.2 (c x)^m (a+b x^n)^p.m #1729)

## Appendix: per-file breakdown

| file | cases | solved | verified | wrong | dont_know | timeout | error | solved % |
|---|---|---|---|---|---|---|---|---|
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_1_linear/t_1_1_1_2.jsonl | 1914 | 985 | 985 | 0 | 910 | 19 | 0 | 51.5% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_1_linear/t_1_1_1_3.jsonl | 3173 | 1606 | 1606 | 0 | 1530 | 37 | 0 | 50.6% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_1_linear/t_1_1_1_4.jsonl | 148 | 20 | 20 | 0 | 127 | 1 | 0 | 13.5% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_1_linear/t_1_1_1_5.jsonl | 34 | 20 | 20 | 0 | 10 | 4 | 0 | 58.8% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_1_linear/t_1_1_1_6.jsonl | 71 | 0 | 0 | 0 | 71 | 0 | 0 | 0.0% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_1_linear/t_1_1_1_7.jsonl | 35 | 0 | 0 | 0 | 35 | 0 | 0 | 0.0% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_2_quadratic/t_1_1_2_2.jsonl | 1036 | 559 | 559 | 0 | 443 | 34 | 0 | 54.0% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_2_quadratic/t_1_1_2_3.jsonl | 344 | 98 | 98 | 0 | 229 | 16 | 1 | 28.5% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_2_quadratic/t_1_1_2_4.jsonl | 1156 | 636 | 636 | 0 | 373 | 147 | 0 | 55.0% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_2_quadratic/t_1_1_2_5.jsonl | 112 | 23 | 23 | 0 | 88 | 1 | 0 | 20.5% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_2_quadratic/t_1_1_2_6.jsonl | 51 | 8 | 8 | 0 | 37 | 6 | 0 | 15.7% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_2_quadratic/t_1_1_2_8.jsonl | 172 | 141 | 141 | 0 | 4 | 27 | 0 | 82.0% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_3_general/t_1_1_3_2.jsonl | 2970 | 1737 | 1736 | 1 | 1186 | 43 | 4 | 58.5% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_3_general/t_1_1_3_3.jsonl | 278 | 104 | 104 | 0 | 169 | 5 | 0 | 37.4% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_3_general/t_1_1_3_4.jsonl | 899 | 416 | 416 | 0 | 435 | 45 | 3 | 46.3% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_3_general/t_1_1_3_6.jsonl | 45 | 6 | 6 | 0 | 33 | 6 | 0 | 13.3% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_3_general/t_1_1_3_8.jsonl | 590 | 318 | 318 | 0 | 175 | 97 | 0 | 53.9% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_4_improper/t_1_1_4_2.jsonl | 451 | 116 | 115 | 1 | 321 | 14 | 0 | 25.7% |
| t_1_algebraic_functions/t_1_1_binomial_products/t_1_1_4_improper/t_1_1_4_3.jsonl | 298 | 168 | 168 | 0 | 106 | 24 | 0 | 56.4% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_1_quadratic/t_1_2_1_1.jsonl | 140 | 109 | 109 | 0 | 30 | 1 | 0 | 77.9% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_1_quadratic/t_1_2_1_2.jsonl | 2545 | 1259 | 1259 | 0 | 734 | 549 | 3 | 49.5% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_1_quadratic/t_1_2_1_3.jsonl | 2643 | 1429 | 1429 | 0 | 396 | 816 | 2 | 54.1% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_1_quadratic/t_1_2_1_4.jsonl | 937 | 116 | 116 | 0 | 542 | 270 | 9 | 12.4% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_1_quadratic/t_1_2_1_5.jsonl | 123 | 77 | 77 | 0 | 6 | 31 | 9 | 62.6% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_1_quadratic/t_1_2_1_6.jsonl | 142 | 30 | 30 | 0 | 17 | 95 | 0 | 21.1% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_1_quadratic/t_1_2_1_9.jsonl | 397 | 256 | 256 | 0 | 22 | 112 | 7 | 64.5% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_2_quartic/t_1_2_2_2.jsonl | 1120 | 494 | 494 | 0 | 552 | 74 | 0 | 44.1% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_2_quartic/t_1_2_2_3.jsonl | 404 | 156 | 156 | 0 | 220 | 22 | 6 | 38.6% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_2_quartic/t_1_2_2_4.jsonl | 398 | 168 | 168 | 0 | 125 | 105 | 0 | 42.2% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_2_quartic/t_1_2_2_5.jsonl | 111 | 74 | 74 | 0 | 4 | 33 | 0 | 66.7% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_2_quartic/t_1_2_2_6.jsonl | 135 | 71 | 71 | 0 | 7 | 57 | 0 | 52.6% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_2_quartic/t_1_2_2_7.jsonl | 37 | 0 | 0 | 0 | 36 | 1 | 0 | 0.0% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_2_quartic/t_1_2_2_8.jsonl | 4 | 0 | 0 | 0 | 4 | 0 | 0 | 0.0% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_3_general/t_1_2_3_2.jsonl | 657 | 241 | 241 | 0 | 327 | 89 | 0 | 36.7% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_3_general/t_1_2_3_3.jsonl | 96 | 25 | 25 | 0 | 52 | 16 | 3 | 26.0% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_3_general/t_1_2_3_4.jsonl | 154 | 86 | 86 | 0 | 40 | 28 | 0 | 55.8% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_3_general/t_1_2_3_5.jsonl | 17 | 0 | 0 | 0 | 14 | 3 | 0 | 0.0% |
| t_1_algebraic_functions/t_1_2_trinomial_products/t_1_2_4_improper/t_1_2_4_2.jsonl | 140 | 58 | 58 | 0 | 73 | 9 | 0 | 41.4% |
| t_1_algebraic_functions/t_1_3_miscellaneous/t_1_3_1.jsonl | 484 | 398 | 398 | 0 | 18 | 68 | 0 | 82.2% |
| t_1_algebraic_functions/t_1_3_miscellaneous/t_1_3_2.jsonl | 866 | 291 | 290 | 1 | 545 | 21 | 9 | 33.6% |
| t_2_exponentials/t_2_1.jsonl | 98 | 12 | 12 | 0 | 86 | 0 | 0 | 12.2% |
| t_2_exponentials/t_2_2.jsonl | 93 | 21 | 21 | 0 | 72 | 0 | 0 | 22.6% |
| t_2_exponentials/t_2_3.jsonl | 772 | 193 | 193 | 0 | 561 | 16 | 2 | 25.0% |
| t_3_logarithms/t_3_1_2.jsonl | 193 | 33 | 33 | 0 | 154 | 0 | 6 | 17.1% |
| t_3_logarithms/t_3_1_4.jsonl | 433 | 69 | 69 | 0 | 185 | 179 | 0 | 15.9% |
| t_3_logarithms/t_3_1_5.jsonl | 248 | 19 | 19 | 0 | 72 | 154 | 3 | 7.7% |
| t_3_logarithms/t_3_2_1.jsonl | 309 | 25 | 25 | 0 | 200 | 84 | 0 | 8.1% |
| t_3_logarithms/t_3_2_2.jsonl | 259 | 13 | 13 | 0 | 52 | 194 | 0 | 5.0% |
| t_3_logarithms/t_3_2_3.jsonl | 108 | 9 | 9 | 0 | 45 | 54 | 0 | 8.3% |
| t_3_logarithms/t_3_3.jsonl | 546 | 42 | 42 | 0 | 427 | 76 | 1 | 7.7% |
| t_3_logarithms/t_3_4.jsonl | 630 | 142 | 142 | 0 | 424 | 50 | 14 | 22.5% |
| t_3_logarithms/t_3_5.jsonl | 310 | 101 | 101 | 0 | 182 | 27 | 0 | 32.6% |
| t_4_trig_functions/t_4_1_sine/t_4_1_0.jsonl | 535 | 24 | 24 | 0 | 511 | 0 | 0 | 4.5% |
| t_4_trig_functions/t_4_1_sine/t_4_1_10.jsonl | 348 | 3 | 3 | 0 | 341 | 4 | 0 | 0.9% |
| t_4_trig_functions/t_4_1_sine/t_4_1_11.jsonl | 113 | 1 | 1 | 0 | 112 | 0 | 0 | 0.9% |
| t_4_trig_functions/t_4_1_sine/t_4_1_12.jsonl | 355 | 11 | 11 | 0 | 335 | 9 | 0 | 3.1% |
| t_4_trig_functions/t_4_1_sine/t_4_1_13.jsonl | 36 | 4 | 4 | 0 | 32 | 0 | 0 | 11.1% |
| t_4_trig_functions/t_4_1_sine/t_4_1_1_1.jsonl | 72 | 0 | 0 | 0 | 72 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_1_sine/t_4_1_1_2.jsonl | 650 | 35 | 35 | 0 | 611 | 4 | 0 | 5.4% |
| t_4_trig_functions/t_4_1_sine/t_4_1_1_3.jsonl | 197 | 0 | 0 | 0 | 197 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_1_sine/t_4_1_2_1.jsonl | 814 | 25 | 25 | 0 | 785 | 4 | 0 | 3.1% |
| t_4_trig_functions/t_4_1_sine/t_4_1_2_2.jsonl | 1534 | 76 | 76 | 0 | 1410 | 48 | 0 | 5.0% |
| t_4_trig_functions/t_4_1_sine/t_4_1_2_3.jsonl | 50 | 4 | 4 | 0 | 46 | 0 | 0 | 8.0% |
| t_4_trig_functions/t_4_1_sine/t_4_1_3_1.jsonl | 354 | 13 | 13 | 0 | 286 | 55 | 0 | 3.7% |
| t_4_trig_functions/t_4_1_sine/t_4_1_4_1.jsonl | 19 | 1 | 1 | 0 | 18 | 0 | 0 | 5.3% |
| t_4_trig_functions/t_4_1_sine/t_4_1_4_2.jsonl | 33 | 0 | 0 | 0 | 33 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_1_sine/t_4_1_7.jsonl | 578 | 9 | 9 | 0 | 569 | 0 | 0 | 1.6% |
| t_4_trig_functions/t_4_1_sine/t_4_1_8.jsonl | 9 | 0 | 0 | 0 | 9 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_1_sine/t_4_1_9.jsonl | 19 | 0 | 0 | 0 | 19 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_0.jsonl | 292 | 1 | 1 | 0 | 291 | 0 | 0 | 0.3% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_10.jsonl | 189 | 3 | 3 | 0 | 184 | 2 | 0 | 1.6% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_12.jsonl | 99 | 11 | 11 | 0 | 88 | 0 | 0 | 11.1% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_13.jsonl | 34 | 3 | 3 | 0 | 31 | 0 | 0 | 8.8% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_1_1.jsonl | 62 | 0 | 0 | 0 | 62 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_1_2.jsonl | 86 | 6 | 6 | 0 | 80 | 0 | 0 | 7.0% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_1_3.jsonl | 22 | 0 | 0 | 0 | 22 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_2_1.jsonl | 923 | 18 | 18 | 0 | 904 | 1 | 0 | 2.0% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_2_2.jsonl | 4 | 0 | 0 | 0 | 4 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_2_3.jsonl | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_3_1.jsonl | 643 | 27 | 27 | 0 | 610 | 6 | 0 | 4.2% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_4_1.jsonl | 393 | 1 | 1 | 0 | 392 | 0 | 0 | 0.3% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_4_2.jsonl | 1538 | 46 | 46 | 0 | 1482 | 10 | 0 | 3.0% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_7.jsonl | 97 | 0 | 0 | 0 | 97 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_8.jsonl | 21 | 0 | 0 | 0 | 21 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_2_cosine/t_4_2_9.jsonl | 20 | 0 | 0 | 0 | 20 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_3_tangent/t_4_3_0.jsonl | 386 | 3 | 3 | 0 | 383 | 0 | 0 | 0.8% |
| t_4_trig_functions/t_4_3_tangent/t_4_3_10.jsonl | 63 | 0 | 0 | 0 | 63 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_3_tangent/t_4_3_11.jsonl | 66 | 1 | 1 | 0 | 65 | 0 | 0 | 1.5% |
| t_4_trig_functions/t_4_3_tangent/t_4_3_1_2.jsonl | 682 | 6 | 6 | 0 | 672 | 4 | 0 | 0.9% |
| t_4_trig_functions/t_4_3_tangent/t_4_3_1_3.jsonl | 91 | 0 | 0 | 0 | 91 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_3_tangent/t_4_3_2_1.jsonl | 1316 | 34 | 34 | 0 | 1277 | 5 | 0 | 2.6% |
| t_4_trig_functions/t_4_3_tangent/t_4_3_3_1.jsonl | 855 | 17 | 17 | 0 | 798 | 40 | 0 | 2.0% |
| t_4_trig_functions/t_4_3_tangent/t_4_3_4_2.jsonl | 171 | 2 | 2 | 0 | 168 | 1 | 0 | 1.2% |
| t_4_trig_functions/t_4_3_tangent/t_4_3_7.jsonl | 499 | 18 | 18 | 0 | 469 | 12 | 0 | 3.6% |
| t_4_trig_functions/t_4_3_tangent/t_4_3_9.jsonl | 51 | 0 | 0 | 0 | 47 | 4 | 0 | 0.0% |
| t_4_trig_functions/t_4_4_cotangent/t_4_4_0.jsonl | 52 | 1 | 1 | 0 | 51 | 0 | 0 | 1.9% |
| t_4_trig_functions/t_4_4_cotangent/t_4_4_10.jsonl | 61 | 0 | 0 | 0 | 61 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_4_cotangent/t_4_4_1_2.jsonl | 23 | 1 | 1 | 0 | 22 | 0 | 0 | 4.3% |
| t_4_trig_functions/t_4_4_cotangent/t_4_4_1_3.jsonl | 19 | 0 | 0 | 0 | 19 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_4_cotangent/t_4_4_2_1.jsonl | 106 | 0 | 0 | 0 | 106 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_4_cotangent/t_4_4_7.jsonl | 64 | 5 | 5 | 0 | 59 | 0 | 0 | 7.8% |
| t_4_trig_functions/t_4_4_cotangent/t_4_4_9.jsonl | 32 | 0 | 0 | 0 | 28 | 4 | 0 | 0.0% |
| t_4_trig_functions/t_4_5_secant/t_4_5_0.jsonl | 299 | 1 | 1 | 0 | 298 | 0 | 0 | 0.3% |
| t_4_trig_functions/t_4_5_secant/t_4_5_10.jsonl | 46 | 0 | 0 | 0 | 46 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_5_secant/t_4_5_11.jsonl | 83 | 2 | 2 | 0 | 72 | 9 | 0 | 2.4% |
| t_4_trig_functions/t_4_5_secant/t_4_5_1_2.jsonl | 878 | 4 | 4 | 0 | 874 | 0 | 0 | 0.5% |
| t_4_trig_functions/t_4_5_secant/t_4_5_1_3.jsonl | 306 | 0 | 0 | 0 | 306 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_5_secant/t_4_5_1_4.jsonl | 351 | 12 | 12 | 0 | 339 | 0 | 0 | 3.4% |
| t_4_trig_functions/t_4_5_secant/t_4_5_2_1.jsonl | 235 | 0 | 0 | 0 | 234 | 1 | 0 | 0.0% |
| t_4_trig_functions/t_4_5_secant/t_4_5_2_3.jsonl | 243 | 1 | 1 | 0 | 241 | 1 | 0 | 0.4% |
| t_4_trig_functions/t_4_5_secant/t_4_5_3_1.jsonl | 630 | 1 | 1 | 0 | 627 | 2 | 0 | 0.2% |
| t_4_trig_functions/t_4_5_secant/t_4_5_4_1.jsonl | 70 | 1 | 1 | 0 | 69 | 0 | 0 | 1.4% |
| t_4_trig_functions/t_4_5_secant/t_4_5_4_2.jsonl | 1369 | 2 | 2 | 0 | 1367 | 0 | 0 | 0.1% |
| t_4_trig_functions/t_4_5_secant/t_4_5_7.jsonl | 471 | 8 | 8 | 0 | 462 | 1 | 0 | 1.7% |
| t_4_trig_functions/t_4_6_cosecant/t_4_6_0.jsonl | 70 | 1 | 1 | 0 | 69 | 0 | 0 | 1.4% |
| t_4_trig_functions/t_4_6_cosecant/t_4_6_11.jsonl | 84 | 3 | 3 | 0 | 76 | 5 | 0 | 3.6% |
| t_4_trig_functions/t_4_6_cosecant/t_4_6_1_2.jsonl | 59 | 0 | 0 | 0 | 59 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_6_cosecant/t_4_6_1_3.jsonl | 16 | 0 | 0 | 0 | 16 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_6_cosecant/t_4_6_1_4.jsonl | 21 | 0 | 0 | 0 | 21 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_6_cosecant/t_4_6_3_1.jsonl | 24 | 2 | 2 | 0 | 22 | 0 | 0 | 8.3% |
| t_4_trig_functions/t_4_6_cosecant/t_4_6_4_2.jsonl | 1 | 0 | 0 | 0 | 1 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_6_cosecant/t_4_6_7.jsonl | 26 | 0 | 0 | 0 | 26 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_7_miscellaneous/t_4_7_1.jsonl | 254 | 6 | 6 | 0 | 193 | 55 | 0 | 2.4% |
| t_4_trig_functions/t_4_7_miscellaneous/t_4_7_2.jsonl | 282 | 15 | 15 | 0 | 262 | 5 | 0 | 5.3% |
| t_4_trig_functions/t_4_7_miscellaneous/t_4_7_3.jsonl | 393 | 5 | 5 | 0 | 346 | 42 | 0 | 1.3% |
| t_4_trig_functions/t_4_7_miscellaneous/t_4_7_4.jsonl | 9 | 0 | 0 | 0 | 9 | 0 | 0 | 0.0% |
| t_4_trig_functions/t_4_7_miscellaneous/t_4_7_5.jsonl | 248 | 10 | 10 | 0 | 238 | 0 | 0 | 4.0% |
| t_4_trig_functions/t_4_7_miscellaneous/t_4_7_6.jsonl | 141 | 21 | 21 | 0 | 119 | 1 | 0 | 14.9% |
| t_4_trig_functions/t_4_7_miscellaneous/t_4_7_7.jsonl | 935 | 198 | 198 | 0 | 696 | 41 | 0 | 21.2% |
| t_5_inverse_trig_functions/t_5_1_inverse_sine/t_5_1_2.jsonl | 224 | 17 | 17 | 0 | 195 | 3 | 9 | 7.6% |
| t_5_inverse_trig_functions/t_5_1_inverse_sine/t_5_1_4.jsonl | 108 | 22 | 22 | 0 | 73 | 13 | 0 | 20.4% |
| t_5_inverse_trig_functions/t_5_1_inverse_sine/t_5_1_5.jsonl | 466 | 41 | 41 | 0 | 348 | 73 | 4 | 8.8% |
| t_5_inverse_trig_functions/t_5_2_inverse_cosine/t_5_2_2.jsonl | 227 | 17 | 17 | 0 | 196 | 5 | 9 | 7.5% |
| t_5_inverse_trig_functions/t_5_2_inverse_cosine/t_5_2_5.jsonl | 151 | 26 | 26 | 0 | 91 | 34 | 0 | 17.2% |
| t_5_inverse_trig_functions/t_5_3_inverse_tangent/t_5_3_2.jsonl | 166 | 58 | 58 | 0 | 104 | 4 | 0 | 34.9% |
| t_5_inverse_trig_functions/t_5_3_inverse_tangent/t_5_3_3.jsonl | 31 | 8 | 8 | 0 | 23 | 0 | 0 | 25.8% |
| t_5_inverse_trig_functions/t_5_3_inverse_tangent/t_5_3_4.jsonl | 1286 | 81 | 81 | 0 | 1017 | 187 | 1 | 6.3% |
| t_5_inverse_trig_functions/t_5_3_inverse_tangent/t_5_3_5.jsonl | 70 | 14 | 14 | 0 | 49 | 7 | 0 | 20.0% |
| t_5_inverse_trig_functions/t_5_3_inverse_tangent/t_5_3_6.jsonl | 381 | 5 | 5 | 0 | 376 | 0 | 0 | 1.3% |
| t_5_inverse_trig_functions/t_5_3_inverse_tangent/t_5_3_7.jsonl | 151 | 38 | 38 | 0 | 65 | 42 | 6 | 25.2% |
| t_5_inverse_trig_functions/t_5_4_inverse_cotangent/t_5_4_1.jsonl | 232 | 63 | 63 | 0 | 126 | 37 | 6 | 27.2% |
| t_5_inverse_trig_functions/t_5_4_inverse_cotangent/t_5_4_2.jsonl | 12 | 1 | 1 | 0 | 11 | 0 | 0 | 8.3% |
| t_5_inverse_trig_functions/t_5_5_inverse_secant/t_5_5_1.jsonl | 168 | 4 | 4 | 0 | 129 | 35 | 0 | 2.4% |
| t_5_inverse_trig_functions/t_5_5_inverse_secant/t_5_5_2.jsonl | 50 | 6 | 6 | 0 | 37 | 7 | 0 | 12.0% |
| t_5_inverse_trig_functions/t_5_6_inverse_cosecant/t_5_6_1.jsonl | 171 | 4 | 4 | 0 | 114 | 53 | 0 | 2.3% |
| t_5_inverse_trig_functions/t_5_6_inverse_cosecant/t_5_6_2.jsonl | 49 | 6 | 6 | 0 | 38 | 5 | 0 | 12.2% |
| t_6_hyperbolic_functions/t_6_1_hyperbolic_sine/t_6_1_1.jsonl | 502 | 5 | 5 | 0 | 494 | 3 | 0 | 1.0% |
| t_6_hyperbolic_functions/t_6_1_hyperbolic_sine/t_6_1_3.jsonl | 101 | 6 | 6 | 0 | 95 | 0 | 0 | 5.9% |
| t_6_hyperbolic_functions/t_6_1_hyperbolic_sine/t_6_1_4.jsonl | 33 | 0 | 0 | 0 | 33 | 0 | 0 | 0.0% |
| t_6_hyperbolic_functions/t_6_1_hyperbolic_sine/t_6_1_5.jsonl | 368 | 16 | 16 | 0 | 352 | 0 | 0 | 4.3% |
| t_6_hyperbolic_functions/t_6_1_hyperbolic_sine/t_6_1_7.jsonl | 523 | 5 | 5 | 0 | 518 | 0 | 0 | 1.0% |
| t_6_hyperbolic_functions/t_6_2_hyperbolic_cosine/t_6_2_1.jsonl | 183 | 0 | 0 | 0 | 183 | 0 | 0 | 0.0% |
| t_6_hyperbolic_functions/t_6_2_hyperbolic_cosine/t_6_2_2.jsonl | 111 | 1 | 1 | 0 | 110 | 0 | 0 | 0.9% |
| t_6_hyperbolic_functions/t_6_2_hyperbolic_cosine/t_6_2_3.jsonl | 67 | 6 | 6 | 0 | 61 | 0 | 0 | 9.0% |
| t_6_hyperbolic_functions/t_6_2_hyperbolic_cosine/t_6_2_4.jsonl | 33 | 0 | 0 | 0 | 33 | 0 | 0 | 0.0% |
| t_6_hyperbolic_functions/t_6_2_hyperbolic_cosine/t_6_2_5.jsonl | 336 | 17 | 17 | 0 | 319 | 0 | 0 | 5.1% |
| t_6_hyperbolic_functions/t_6_2_hyperbolic_cosine/t_6_2_7.jsonl | 84 | 0 | 0 | 0 | 84 | 0 | 0 | 0.0% |
| t_6_hyperbolic_functions/t_6_3_hyperbolic_tangent/t_6_3_1.jsonl | 77 | 0 | 0 | 0 | 66 | 11 | 0 | 0.0% |
| t_6_hyperbolic_functions/t_6_3_hyperbolic_tangent/t_6_3_2.jsonl | 202 | 45 | 45 | 0 | 155 | 2 | 0 | 22.3% |
| t_6_hyperbolic_functions/t_6_3_hyperbolic_tangent/t_6_3_7.jsonl | 257 | 23 | 23 | 0 | 203 | 31 | 0 | 8.9% |
| t_6_hyperbolic_functions/t_6_4_hyperbolic_cotangent/t_6_4_1.jsonl | 61 | 0 | 0 | 0 | 61 | 0 | 0 | 0.0% |
| t_6_hyperbolic_functions/t_6_4_hyperbolic_cotangent/t_6_4_2.jsonl | 181 | 1 | 1 | 0 | 180 | 0 | 0 | 0.6% |
| t_6_hyperbolic_functions/t_6_4_hyperbolic_cotangent/t_6_4_7.jsonl | 53 | 0 | 0 | 0 | 53 | 0 | 0 | 0.0% |
| t_6_hyperbolic_functions/t_6_5_hyperbolic_secant/t_6_5_1.jsonl | 16 | 0 | 0 | 0 | 16 | 0 | 0 | 0.0% |
| t_6_hyperbolic_functions/t_6_5_hyperbolic_secant/t_6_5_2.jsonl | 84 | 2 | 2 | 0 | 81 | 1 | 0 | 2.4% |
| t_6_hyperbolic_functions/t_6_5_hyperbolic_secant/t_6_5_3.jsonl | 201 | 3 | 3 | 0 | 198 | 0 | 0 | 1.5% |
| t_6_hyperbolic_functions/t_6_5_hyperbolic_secant/t_6_5_7.jsonl | 217 | 8 | 8 | 0 | 209 | 0 | 0 | 3.7% |
| t_6_hyperbolic_functions/t_6_6_hyperbolic_cosecant/t_6_6_1.jsonl | 29 | 0 | 0 | 0 | 29 | 0 | 0 | 0.0% |
| t_6_hyperbolic_functions/t_6_6_hyperbolic_cosecant/t_6_6_2.jsonl | 83 | 2 | 2 | 0 | 81 | 0 | 0 | 2.4% |
| t_6_hyperbolic_functions/t_6_6_hyperbolic_cosecant/t_6_6_3.jsonl | 174 | 3 | 3 | 0 | 171 | 0 | 0 | 1.7% |
| t_6_hyperbolic_functions/t_6_6_hyperbolic_cosecant/t_6_6_7.jsonl | 26 | 0 | 0 | 0 | 26 | 0 | 0 | 0.0% |
| t_6_hyperbolic_functions/t_6_7_miscellaneous/t_6_7_1.jsonl | 1051 | 44 | 44 | 0 | 1006 | 1 | 0 | 4.2% |
| t_7_inverse_hyperbolic_functions/t_7_1_inverse_hyperbolic_sine/t_7_1_2.jsonl | 156 | 0 | 0 | 0 | 156 | 0 | 0 | 0.0% |
| t_7_inverse_hyperbolic_functions/t_7_1_inverse_hyperbolic_sine/t_7_1_4.jsonl | 58 | 0 | 0 | 0 | 58 | 0 | 0 | 0.0% |
| t_7_inverse_hyperbolic_functions/t_7_1_inverse_hyperbolic_sine/t_7_1_5.jsonl | 371 | 6 | 6 | 0 | 333 | 4 | 28 | 1.6% |
| t_7_inverse_hyperbolic_functions/t_7_2_inverse_hyperbolic_cosine/t_7_2_2.jsonl | 166 | 0 | 0 | 0 | 166 | 0 | 0 | 0.0% |
| t_7_inverse_hyperbolic_functions/t_7_2_inverse_hyperbolic_cosine/t_7_2_4.jsonl | 105 | 0 | 0 | 0 | 105 | 0 | 0 | 0.0% |
| t_7_inverse_hyperbolic_functions/t_7_2_inverse_hyperbolic_cosine/t_7_2_5.jsonl | 291 | 0 | 0 | 0 | 286 | 5 | 0 | 0.0% |
| t_7_inverse_hyperbolic_functions/t_7_3_inverse_hyperbolic_tangent/t_7_3_2.jsonl | 243 | 0 | 0 | 0 | 243 | 0 | 0 | 0.0% |
| t_7_inverse_hyperbolic_functions/t_7_3_inverse_hyperbolic_tangent/t_7_3_3.jsonl | 47 | 0 | 0 | 0 | 47 | 0 | 0 | 0.0% |
| t_7_inverse_hyperbolic_functions/t_7_3_inverse_hyperbolic_tangent/t_7_3_4.jsonl | 526 | 9 | 9 | 0 | 517 | 0 | 0 | 1.7% |
| t_7_inverse_hyperbolic_functions/t_7_3_inverse_hyperbolic_tangent/t_7_3_5.jsonl | 62 | 0 | 0 | 0 | 62 | 0 | 0 | 0.0% |
| t_7_inverse_hyperbolic_functions/t_7_3_inverse_hyperbolic_tangent/t_7_3_6.jsonl | 1370 | 9 | 9 | 0 | 1324 | 37 | 0 | 0.7% |
| t_7_inverse_hyperbolic_functions/t_7_3_inverse_hyperbolic_tangent/t_7_3_7.jsonl | 359 | 13 | 13 | 0 | 346 | 0 | 0 | 3.6% |
| t_7_inverse_hyperbolic_functions/t_7_4_inverse_hyperbolic_cotangent/t_7_4_1.jsonl | 296 | 10 | 10 | 0 | 286 | 0 | 0 | 3.4% |
| t_7_inverse_hyperbolic_functions/t_7_4_inverse_hyperbolic_cotangent/t_7_4_2.jsonl | 921 | 8 | 8 | 0 | 888 | 25 | 0 | 0.9% |
| t_7_inverse_hyperbolic_functions/t_7_5_inverse_hyperbolic_secant/t_7_5_1.jsonl | 187 | 0 | 0 | 0 | 185 | 2 | 0 | 0.0% |
| t_7_inverse_hyperbolic_functions/t_7_5_inverse_hyperbolic_secant/t_7_5_2.jsonl | 90 | 0 | 0 | 0 | 90 | 0 | 0 | 0.0% |
| t_7_inverse_hyperbolic_functions/t_7_6_inverse_hyperbolic_cosecant/t_7_6_1.jsonl | 171 | 0 | 0 | 0 | 170 | 1 | 0 | 0.0% |
| t_7_inverse_hyperbolic_functions/t_7_6_inverse_hyperbolic_cosecant/t_7_6_2.jsonl | 69 | 0 | 0 | 0 | 69 | 0 | 0 | 0.0% |
| t_8_special_functions/t_8_1.jsonl | 311 | 27 | 27 | 0 | 243 | 41 | 0 | 8.7% |
| t_8_special_functions/t_8_2.jsonl | 218 | 16 | 16 | 0 | 198 | 4 | 0 | 7.3% |
| t_8_special_functions/t_8_3.jsonl | 208 | 31 | 31 | 0 | 168 | 8 | 1 | 14.9% |
| t_8_special_functions/t_8_4.jsonl | 136 | 4 | 4 | 0 | 126 | 6 | 0 | 2.9% |
| t_8_special_functions/t_8_5.jsonl | 136 | 4 | 4 | 0 | 132 | 0 | 0 | 2.9% |
| t_8_special_functions/t_8_6.jsonl | 233 | 39 | 39 | 0 | 192 | 2 | 0 | 16.7% |
| t_8_special_functions/t_8_7.jsonl | 14 | 0 | 0 | 0 | 14 | 0 | 0 | 0.0% |
| t_8_special_functions/t_8_8.jsonl | 198 | 4 | 4 | 0 | 192 | 1 | 1 | 2.0% |
| t_8_special_functions/t_8_9.jsonl | 398 | 15 | 15 | 0 | 383 | 0 | 0 | 3.8% |
| **TOTAL** | **67883** | **14474** | **14471** | **3** | **48436** | **4826** | **147** | **21.3%** |

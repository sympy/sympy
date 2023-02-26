from sympy import gcd
def chicken_mc(m, v, check=True):
    r"""Chicken McNugget Theorm (Frobenius Coin Problem).

  for any two relatively prime integers m and n,
  the largest integer that cannot be expressed in the form am + bn
  for non-negative integers a and b is mn − m − n.

  If the a, b are co-prime the correct result will be returned
  else the return value will be None if there is no solution.

  Examples
  ========

  As an example consider two numbers ``M = 5``
  ``V = 9``. Then we have::

     >>> from sympy.ntheory.modular import chicken_mc

     >>> chicken_mc(5, 9)
     31


  If the moduli are not co-prime, you may receive an incorrect result
  if you use ``check=False``:

     >>> chicken_mc(4, 2,check=False)
     2
     >>> chicken_mc(4, 2) is None
     True

  """

    result = m * v - (m + v)
    if gcd(m, v) != 1 and check:
        return None
    return result

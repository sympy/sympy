from __future__ import division, print_function

from sympy.core.function import expand_mul


def _iszero(x):
    """Returns True if x is zero."""
    return getattr(x, 'is_zero', None)


def _is_zero_after_expand_mul(x):
    """Tests by expand_mul only, suitable for polynomials and rational
    functions."""
    return expand_mul(x) == 0



from sympy.simplify.simplify import simplify as _simplify, dotprodsimp as _dotprodsimp

_simpopts_default = {}
_simpopts_none    = {'simplify': False, 'simpfunc': _noplamb, 'dotprodsimp': False}
_simpopts_full    = {'simplify': True, 'simpfunc': _simplify, 'dotprodsimp': True}
_simpopts_dotprod = {'dotprodsimp': True}

_simpopts_dicts   = {
  None:      _simpopts_default,
  False:     _simpopts_none,
  True:      _simpopts_full,
  'default': _simpopts_default,
  'none':    _simpopts_none,
  'full':    _simpopts_full,
  'dotprod': _simpopts_dotprod,
}

def _get_simpopts(simpopts, *args, simplify_as_func=False, dotprodsimp_as_func=False, **kwargs):
  simpopts = _simpopts_dicts.get(simpopts, simpopts) # if 'simpopts' is dict it is used
  ret      = []

  for opt, default in zip(args[::2], args[1::2]):
    val = kwargs.get(opt, ret) # ret used as sentinel

    if val is ret:
      val = simpopts.get(opt, default)

    if opt == 'simplify':
      if simplify_as_func and not callable(_simplify):
        val = (lambda x: x, _simplify)[bool(val)]
    elif opt == 'dotprodsimp':
      if dotprodsimp_as_func and not callable(_dotprodsimp):
        val = (lambda x: x, _dotprodsimp)[bool(val)]

    ret.append(val)

  return ret

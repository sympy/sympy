from core import Symbol

_latin = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
# COSINEQ should not be imported as they clash; gamma, pi and zeta clash, too
_greek = 'alpha beta gamma delta epsilon zeta eta theta iota kappa lamda '\
  'mu nu xi omicron pi rho sigma tau upsilon phi chi psi omega'.split(' ')
# Note: We import lamda since lambda is a reserved keyword in Python

for _s in _latin + _greek:
    exec "%s = Symbol('%s')" % (_s, _s)

del _latin, _greek, _s

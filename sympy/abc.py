from core import Symbol

_latin = list('abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ')
_greek = 'alpha beta gamma delta epsilon zeta eta theta iota kappa '\
  'mu nu xi omicron pi rho sigma tau upsilon phi chi psi omega'.split(' ')

for s in _latin + _greek:
    exec "%s = Symbol('%s')" % (s, s)

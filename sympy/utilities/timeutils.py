"""Simple tools for timing functions execution, when IPython is not available. """

import timeit, math

_scales = [1e0, 1e3, 1e6, 1e9]
_units  = [u's', u'ms', u'\u03bcs', u'ns']

def timed(func):
    """Adaptively measure execution time of a function. """
    timer = timeit.Timer(func)

    repeat, number = 3, 1

    for i in range(1, 10):
        if timer.timeit(number) >= 0.2:
            break
        else:
            number *= 10

    time = min(timer.repeat(repeat, number)) / number

    if time > 0.0:
        order = min(-int(math.floor(math.log10(time)) // 3), 3)
    else:
        order = 3

    return (number, time, time*_scales[order], _units[order])

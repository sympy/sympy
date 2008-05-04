from sympy import *

def textplot(expr, a, b, W=55, H=18):
    """
    Print a crude ASCII art plot of the SymPy expression 'expr' (which
    should contain a single symbol, e.g. x or something else) over the
    interval [a, b].

    Example: textplot(sin(t)*t, 0, 15)
    """

    f = None
    for x in expr.atoms():
        if isinstance(x, Symbol):
            f = lambdify([x], expr)
            break
    assert f is not None
    a = float(a)
    b = float(b)

    # Calculate function values
    y = [0] * W
    for x in range(W):
        try:
            y[x] = f(a + (b-a)/float(W)*x)
        except:
            y[x] = 0

    # Normalize height to screen space
    ma = max(y)
    mi = min(y)
    for x in range(W):
        y[x] = int(float(H) * (y[x] - mi) / (ma - mi))
    margin = 7
    print

    for h in range(H-1, -1, -1):
        s = [' '] * W
        for x in range(W):
            if y[x] == h:
                s[x] = '.'

        # Print y values
        if h == H-1:    prefix = ("%g" % ma).rjust(margin)[:margin]
        elif h == H//2: prefix = ("%g" % ((mi+ma)/2)).rjust(margin)[:margin]
        elif h == 0:    prefix = ("%g" % mi).rjust(margin)[:margin]
        else:           prefix = " "*margin
        s = "".join(s)
        if h == H//2: s = s.replace(" ", "-")
        print prefix + " | " + s

    # Print x values
    bottom = " " * (margin + 3)
    bottom += ("%g" % a).ljust(W//2-4)
    bottom += ("%g" % ((a+b)/2)).ljust(W//2)
    bottom += "%g" % b
    print bottom

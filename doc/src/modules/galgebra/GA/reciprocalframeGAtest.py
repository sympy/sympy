    >>> MV.setup('e1 e2 e3',metric)

    >>> print 'Example: Reciprocal Frames e1, e2, and e3 unit vectors.\n\n'

    Example: Reciprocal Frames e1, e2, and e3 unit vectors.



    >>> E = e1^e2^e3
    >>> Esq = (E*E)()

    >>> print 'E =',E
    E = e1^e2^e3

    >>> print 'E^2 =',Esq
    E^2 = -1 - 2*(e1.e2)*(e1.e3)*(e2.e3) + (e1.e2)**2 + (e1.e3)**2 + (e2.e3)**2

    >>> Esq_inv = 1/Esq

    >>> E1 = (e2^e3)*E
    >>> E2 = (-1)*(e1^e3)*E
    >>> E3 = (e1^e2)*E

    >>> print 'E1 = (e2^e3)*E =',E1
    E1 = (e2^e3)*E = {-1 + (e2.e3)**2}e1+{(e1.e2) - (e1.e3)*(e2.e3)}e2+{(e1.e3) - (e1.e2)*(e2.e3)}e3

    >>> print 'E2 =-(e1^e3)*E =',E2
    E2 =-(e1^e3)*E = {(e1.e2) - (e1.e3)*(e2.e3)}e1+{-1 + (e1.e3)**2}e2+{(e2.e3) - (e1.e2)*(e1.e3)}e3

    >>> print 'E3 = (e1^e2)*E =',E3
    E3 = (e1^e2)*E = {(e1.e3) - (e1.e2)*(e2.e3)}e1+{(e2.e3) - (e1.e2)*(e1.e3)}e2+{-1 + (e1.e2)**2}e3

    >>> w = (E1|e2)
    >>> w.collect(MV.g)
    >>> w = w().expand()

    >>> print 'E1|e2 =',w
    E1|e2 = 0

    >>> w = (E1|e3)
    >>> w.collect(MV.g)
    >>> w = w().expand()

    >>> print 'E1|e3 =',w
    E1|e3 = 0

    >>> w = (E2|e1)
    >>> w.collect(MV.g)
    >>> w = w().expand()

    >>> print 'E2|e1 =',w
    E2|e1 = 0

    >>> w = (E2|e3)
    >>> w.collect(MV.g)
    >>> w = w().expand()

    >>> print 'E2|e3 =',w
    E2|e3 = 0

    >>> w = (E3|e1)
    >>> w.collect(MV.g)
    >>> w = w().expand()

    >>> print 'E3|e1 =',w
    E3|e1 = 0

    >>> w = (E3|e2)
    >>> w.collect(MV.g)
    >>> w = w().expand()

    >>> print 'E3|e2 =',w
    E3|e2 = 0

    >>> w = (E1|e1)
    >>> w = w().expand()
    >>> Esq = Esq.expand()

    >>> print '(E1|e1)/E^2 =',w/Esq
    (E1|e1)/E^2 = 1

    >>> w = (E2|e2)
    >>> w = w().expand()

    >>> print '(E2|e2)/E^2 =',w/Esq
    (E2|e2)/E^2 = 1

    >>> w = (E3|e3)
    >>> w = w().expand()

    >>> print '(E3|e3)/E^2 =',w/Esq
    (E3|e3)/E^2 = 1


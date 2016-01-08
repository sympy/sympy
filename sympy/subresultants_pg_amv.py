# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 13:25:02 2015

@author: alkis
"""

# from __future__ import print_function, division

from sympy import *

def sylvester(f, g, x, method = 1):
    '''
      f, g polys in x, m = deg(f) >= deg(g) = n.

      a. If method = 1 (default), computes sylvester1, Sylvester's matrix of 1840 
          of dimension (m + n) x (m + n). The determinants of properly chosen 
          submatrices of this matrix (a.k.a. subresultants) can be 
          used to compute the coefficients of the Euclidean PRS of f, g.     
 
      b. If method = 2, computes sylvester2, Sylvester's matrix of 1853
          of dimension (2*m) x (2*m). The determinants of properly chosen 
          submatrices of this matrix (a.k.a. ``modified'' subresultants) can be 
          used to compute the coefficients of the Sturmian PRS of f, g.

      Applications of these matrices can be found in the references below.
      Especially, for applications of sylvester2, see the first reference!!

      References: 
      ===========
      1. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``On a Theorem 
      by Van Vleck Regarding Sturm Sequences. Serdica Journal of Computing, 
      Vol. 7, No 4, 101–134, 2013. 

      2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``Sturm Sequences 
      and Modified Subresultant Polynomial Remainder Sequences.'' 
      Serdica Journal of Computing, Vol. 8, No 1, 29–46, 2014. 

    '''
    # obtain degrees of polys
    m, n = degree( Poly(f, x), x), degree( Poly(g, x), x)

    # Special cases:
    # A:: case m = n = -oo (i.e. both polys are 0)
    if m == n and n == -oo:
        return Matrix([])

    # B:: case m = n = 0  (i.e. both polys are constants)
    if m == n and n == 0:
        return Matrix([])

    # C:: m == 0 and n == -oo or m == -oo and n == 0
    # (i.e. one polys is constant and the other is 0)
    if m == 0 and n == -oo:
        return Matrix([])
    elif m == -oo and n == 0:
        return Matrix([])

    # D:: m >= 1 and n == -oo or m == -oo and n >=1  
    # (i.e. one poly is of degree >=1 and the other is 0) 
    if m >= 1 and n == -oo:
        return Matrix([0])
    elif m == -oo and n >= 1:
        return Matrix([0])

    fp = Poly(f, x).all_coeffs()  
    gp = Poly(g, x).all_coeffs()

    # order lists of poly coeffs according to degree
    if m < n:
        fp, gp = gp, fp
        m, n = n, m

    # Sylvester's matrix of 1840 (default; a.k.a. sylvester1)
    if method <= 1:
        M = zeros(m + n)
        k = 0
        for i in range(n):
            j = k
            for coeff in fp:
                M[i, j] = coeff
                j = j + 1
            k = k + 1
        k = 0
        for i in range(n, m + n):
            j = k
            for coeff in gp:
                M[i, j] = coeff
                j = j + 1
            k = k + 1
        return M

    # Sylvester's matrix of 1853 (a.k.a sylvester2)
    if method >= 2:
        dim = 2*m
        M = zeros( dim )
        k = 0
        for i in range( m ):
            j = k
            for coeff in fp:
                if i == 0:
                    M[i, j] = coeff
                else:
                    M[2*i, j] = coeff
                j = j + 1
            j = m - n + k
            for coeff in gp:
                if i == 0:
                    M[i+1, j] = coeff
                else:   
                    M[2*i + 1, j] = coeff
                j = j + 1
            k = k + 1
        return M

def sturm_PG(p,q,x, method=0):
    """
    p, q are polynomials in Z[x] or Q[x].

    Computes the (generalized) Sturm sequence of p and q in Z[x] or Q[x]. 
    If q = diff(p, x, 1) it is the usual Sturm sequence.

    A. If method == 0, default, the remainder coefficients of the sequence 
       are (in absolute value) ``modified'' subresultants, which for non-monic 
       polynomials are greater than the coefficients of the corresponding 
       subresultants by the factor Abs(LC(p)**( deg(p)- deg(q))). 
       

    B. If method == 1, the remainder coefficients of the sequence are (in 
       absolute value) subresultants, which for non-monic polynomials are 
       smaller than the coefficients of the corresponding ``modified'' 
       subresultants by the factor Abs(LC(p)**( deg(p)- deg(q))).

    If the Sturm sequence is complete and method=0 then the coefficients 
    of the polynomials in the sequence are ``modified'' subresultants. 
    That is, they are  determinants of appropriately selected submatrices of
    sylvester2, Sylvester's matrix of 1853. In this case the Sturm sequence 
    coincides with the ``modified'' subresultant prs, of the polynomials 
    p, q.

    If the Sturm sequence is incomplete and method=0 then the signs of the 
    coefficients of the polynomials in the sequence may differ from the signs 
    of the coefficients of the corresponding polynomials in the ``modified'' 
    subresultant prs; however, the absolute values are the same. 

    To compute the coefficients, no determinant evaluation takes place. Instead, 
    polynomial divisions in Q[x] are performed, using the function rem(p, q, x);  
    the coefficients of the remainders computed this way become (``modified'') 
    subresultants with the help of the Pell-Gordon Theorem of 1917.
    See also the function euclid_PG(p, q, x).
    
    References:
    =========== 
    1. Pell A. J., R. L. Gordon. The Modified Remainders Obtained in Finding
    the Highest Common Factor of Two Polynomials. Annals of Mathematics,
    Second Series, 18 (1917), No. 4, 188–193.

    2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``Sturm Sequences 
    and Modified Subresultant Polynomial Remainder Sequences.'' 
    Serdica Journal of Computing, Vol. 8, No 1, 29–46, 2014. 
     
    """
    # make sure neither p nor q is 0
    if p == 0 or q == 0:
        return [p, q]

    FLAG = 0
    d0 =  degree(p, x)
    d1 =  degree(q, x)
    # make sure proper degrees  
    if d0 == 0 and d1 == 0:
        return []
    if d1 > d0:
        d0, d1 = d1, d0
        p, q = q, p
    if d0 > 0 and d1 == 0:
        return [p,q]
        
    # make sure LC(p) > 0
    if  LC(p,x) < 0:
        FLAG = 1
        p = -p
        q = -q
        
    # initialize
    lcf = LC(p, x)**(d0 - d1)               # lcf * subr = modified subr 
    a0, a1 = p, q                           # the input polys
    sturmSeq = [a0, a1]                     # the output list 
    del0 = d0 - d1                          # degree difference
    rho1 =  LC(a1, x)                       # leading coeff of a1
    expDeg = d1 - 1                         # expected degree of a2
    a2 = - rem(a0, a1, domain=QQ)           # first remainder
    rho2 =  LC(a2,x)                        # leading coeff of a2
    d2 =  degree(a2, x)                     # actual degree of a2
    degDiff_NEW = expDeg - d2               # expected - actual degree
    del1 = d1 - d2                          # degree difference

    # mulFac is the factor by which a2 is multiplied to 
    # get integer coefficients
    mulFac_OLD = rho1**(del0 + del1 - degDiff_NEW) 

    # update variable and append
    degDiff_OLD = degDiff_NEW 
    if method == 0:
        sturmSeq.append( simplify(lcf * a2 *  Abs(mulFac_OLD)))
    else:
        sturmSeq.append( simplify( a2 *  Abs(mulFac_OLD)))
    # main loop
    while d2 > 0:     
        a0, a1, d0, d1 = a1, a2, d1, d2      # update polys and degrees
        del0 = del1                          # update degree difference
        expDeg = d1 - 1                      # new expected degree 
        a2 = - rem(a0, a1, domain=QQ)        # new remainder
        rho3 =  LC(a2, x)                    # leading coeff of a2
        d2 =  degree(a2, x)                  # actual degree of a2
        degDiff_NEW = expDeg - d2            # expected - actual degree
        del1 = d1 - d2                       # degree difference

        # take into consideration the power 
        # rho1**degDiff_OLD that was "left out"
        expo_OLD = degDiff_OLD               # rho1 raised to this power
        expo_NEW = del0 + del1 - degDiff_NEW # rho2 raised to this power

        mulFac_NEW = rho2**(expo_NEW) * rho1**(expo_OLD) * mulFac_OLD 

        # update variables and append
        degDiff_OLD, mulFac_OLD = degDiff_NEW, mulFac_NEW
        rho1, rho2 = rho2, rho3
        if method == 0:
            sturmSeq.append( simplify(lcf * a2 *  Abs(mulFac_OLD)))
        else:
            sturmSeq.append( simplify( a2 *  Abs(mulFac_OLD)))

    if FLAG:          # change the sign of the sequence
        sturmSeq = [-i for i in sturmSeq]
    
    # gcd is of degree > 0 ?   
    m = len(sturmSeq)
    if sturmSeq[m - 1] == nan or sturmSeq[m - 1] == 0:
        sturmSeq.pop(m - 1)

    return sturmSeq

def sturm_Q(p, q, x):
    """
    p, q are polynomials in Z[x] or Q[x];

    Computes the (generalized) Sturm sequence of p and q in Q[x]. 
    Polynomial divisions in Q[x] are performed, using the function rem(p, q, x).  

    The coefficients of the polynomials in the Sturm sequence can be uniquely 
    determined from the corresponding coefficients of the polynomials found 
    either in: 

        (a) the ``modified'' subresultant prs, (references 1, 2)

    or in

        (b) the subresultant prs (reference 3).

    References:
    =========== 
    1. Pell A. J., R. L. Gordon. The Modified Remainders Obtained in Finding
    the Highest Common Factor of Two Polynomials. Annals of Mathematics,
    Second Series, 18 (1917), No. 4, 188–193.

    2 Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``Sturm Sequences 
    and Modified Subresultant Polynomial Remainder Sequences.'' 
    Serdica Journal of Computing, Vol. 8, No 1, 29–46, 2014. 

    3. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``A Basic Result 
    on the Theory of Subresultants.'' Submitted for publication.

    """
    # make sure neither p nor q is 0
    if p == 0 or q == 0:
        return [p, q]

    FLAG = 0
    d0 =  degree(p, x)
    d1 =  degree(q, x)

    # make sure proper degrees  
    if d0 == 0 and d1 == 0:
        return []
    if d1 > d0:
        d0, d1 = d1, d0
        p, q = q, p
    if d0 > 0 and d1 == 0:
        return [p,q]
        
    # make sure LC(p) > 0
    if  LC(p,x) < 0:
        FLAG = 1
        p = -p
        q = -q

    # initialize
    a0, a1 = p, q                           # the input polys
    sturmSeq = [a0, a1]                     # the output list 
    a2 = -rem(a0, a1, domain=QQ)            # first remainder
    d2 =  degree(a2, x)                     # degree of a2
    sturmSeq.append( a2 ) 


    # main loop
    while d2 > 0:     
        a0, a1, d0, d1 = a1, a2, d1, d2      # update polys and degrees
        a2 = -rem(a0, a1, domain=QQ)         # new remainder
        d2 =  degree(a2, x)                  # actual degree of a2
        sturmSeq.append( a2 ) 
    
    if FLAG:          # change the sign of the sequence
        sturmSeq = [-i for i in sturmSeq]

    # gcd is of degree > 0 ?   
    m = len(sturmSeq)
    if sturmSeq[m - 1] == nan or sturmSeq[m - 1] == 0:
        sturmSeq.pop(m - 1)

    return sturmSeq

def sturm_AMV(p, q, x, method=0):
    """
    p, q are polynomials in Z[x] or Q[x].

    Computes the (generalized) Sturm sequence of p and q in Z[x] or Q[x]. 
    If q = diff(p, x, 1) it is the usual Sturm sequence.

    A. If method == 0, default, the remainder coefficients of the 
       sequence are (in absolute value) ``modified'' subresultants, which 
       for non-monic polynomials are greater than the coefficients of the 
       corresponding subresultants by the factor Abs(LC(p)**( deg(p)- deg(q))). 

    B. If method == 1, the remainder coefficients of the sequence are (in 
       absolute value) subresultants, which for non-monic polynomials are 
       smaller than the coefficients of the corresponding ``modified''
       subresultants by the factor Abs(LC(p)**( deg(p)- deg(q))).

    If the Sturm sequence is complete and method=0 then the coefficients 
    of the polynomials in the sequence are ``modified'' subresultants. 
    That is, they are  determinants of appropriately selected submatrices of
    sylvester2, Sylvester's matrix of 1853. In this case the Sturm sequence 
    coincides with the ``modified'' subresultant prs, of the polynomials 
    p, q.

    If the Sturm sequence is incomplete and method=0 then the signs of the 
    coefficients of the polynomials in the sequence may differ from the signs 
    of the coefficients of the corresponding polynomials in the ``modified'' 
    subresultant prs; however, the absolute values are the same. 

    To compute the coefficients, no determinant evaluation takes place.    
    Instead, we first compute the euclidean sequence  of p and q using 
    euclid_AMV(p, q, x) and then: (a) change the signs of the remainders in the 
    Euclidean sequence according to the pattern "-, -, +, +, -, -, +, +,..." 
    (see Lemma 1 in the 1st reference or Theorem 3 in the 2nd reference)
    and (b) if method=0, assuming deg(p) > deg(q), we multiply the remainder 
    coefficients of the Euclidean sequence times the factor 
    Abs(LC(p)**( deg(p)- deg(q))) to make them modified subresultants.
    See also the function sturm_PG(p, q, x).

    References:
    ===========
    1. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``A Basic Result 
    on the Theory of Subresultants.'' Submitted for publication.

    2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``On the Remainders 
    Obtained in Finding the Greatest Common Divisor of Two Polynomials.'' Serdica 
    Journal of Computing, to appear.

    3. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``Subresultant Polynomial   
    Remainder Sequences Obtained by Polynomial Divisions in Q[x] or in Z[x].''  
    Submitted for publication.
   
    """
    # compute the euclidean sequence 
    prs = euclid_AMV(p, q, x)

    # defensive
    if prs == [] or len(prs) == 2:
        return prs

    # the coefficients in prs are subresultants and hence are smaller
    # than the corresponding subresultants by the factor 
    # |LC(prs[0])**( deg(prs[0]) - deg(prs[1]))|; Theorem 2, 2nd reference.
    lcf = Abs( LC(prs[0])**( degree(prs[0]) - degree(prs[1]) ) )
    
    # the signs of the first two polys in the sequence stay the same
    sturmSeq = [prs[0], prs[1]]
    flag = 0
    m = len(prs)
    i = 2
    # change the signs according to "-, -, +, +, -, -, +, +,..."
    # and multiply times lcf if needed
    while i <= m-1:
        if  flag == 0:
            sturmSeq.append( - prs[i] )
            i = i + 1
            if i == m: 
                break
            sturmSeq.append( - prs[i] )
            i = i + 1
            flag = 1
        elif flag == 1:
            sturmSeq.append( prs[i] )
            i = i + 1
            if i == m: 
                break
            sturmSeq.append( prs[i] )
            i = i + 1
            flag = 0

    # subresultants or modified subresultants?  
    if method == 0 and lcf > 1:
        auxSeq = [sturmSeq[0], sturmSeq[1]]
        for i in range(2, m):
            auxSeq.append(simplify(sturmSeq[i] * lcf ))
        sturmSeq = auxSeq

    return sturmSeq

def euclid_PG(p, q, x):
    """
    p, q are polynomials in Z[x] or Q[x];

    Computes the Euclidean sequence of p and q in Z[x] or Q[x]. 
 
    If the Euclidean sequence is complete the coefficients of the polynomials 
    in the sequence are subresultants. That is, they are  determinants of 
    appropriately selected submatrices of sylvester1, Sylvester's matrix of 1840. 
    In this case the Euclidean sequence coincides with the subresultant prs 
    of the polynomials p, q.

    If the Euclidean sequence is incomplete the signs of the coefficients of the 
    polynomials in the sequence may differ from the signs of the coefficients of
    the corresponding polynomials in the subresultant prs; however, the absolute 
    values are the same. 

    To compute the Euclidean sequence, no determinant evaluation takes place. 
    We first compute the (generalized) Sturm sequence  of p and q using 
    sturm_PG(p, q, x, 1), in which case the coefficients are (in absolute value) 
    equal to subresultants. Then we change the signs of the remainders in the 
    Sturm sequence according to the pattern "-, -, +, +, -, -, +, +,..." ;
    see Lemma 1 in the 1st reference or Theorem 3 in the 2nd reference as well as
    the function sturm_PG(p, q, x).
    
    References:
    ===========
    1. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``A Basic Result 
    on the Theory of Subresultants.'' Submitted for publication.

    2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``On the Remainders 
    Obtained in Finding the Greatest Common Divisor of Two Polynomials.'' Serdica 
    Journal of Computing, to appear.

    3. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``Subresultant Polynomial   
    Remainder Sequences Obtained by Polynomial Divisions in Q[x] or in Z[x].''  
    Submitted for publication.
   
    """
    
    # compute the sturmian sequence using the Pell-Gordon theorem
    # with the coefficients in the prs being (in absolute value) subresultants
    prs = sturm_PG(p, q, x, 1)  ## or prs = sturm_AMV(p, q, x, 1)

    # defensive
    if prs == [] or len(prs) == 2:
        return prs

    # the signs of the first two polys in the sequence stay the same
    euclidSeq = [prs[0], prs[1]]
    flag = 0
    m = len(prs)
    i = 2

    # change the signs according to "-, -, +, +, -, -, +, +,..."
    while i <= m-1:
        if  flag == 0:
            euclidSeq.append(- prs[i] )
            i = i + 1
            if i == m: 
                break
            euclidSeq.append(- prs[i] )
            i = i + 1
            flag = 1
        elif flag == 1:
            euclidSeq.append(prs[i] )
            i = i + 1
            if i == m: 
                break
            euclidSeq.append(prs[i] )
            i = i + 1
            flag = 0

    return euclidSeq

def euclid_Q(p, q, x):
    """
    p, q are polynomials in Z[x] or Q[x];

    Computes the Euclidean sequence of p and q in Q[x].
    Polynomial divisions in Q[x] are performed, using the function rem(p, q, x).  

    The coefficients of the polynomials in the Euclidean sequence can be uniquely 
    determined from the corresponding coefficients of the polynomials found 
    either in: 

        (a) the ``modified'' subresultant polynomial remainder sequence,
    (references 1, 2) 

    or in

        (b) the subresultant polynomial remainder sequence (references 3).

    References:
    =========== 
    1. Pell A. J., R. L. Gordon. The Modified Remainders Obtained in Finding
    the Highest Common Factor of Two Polynomials. Annals of Mathematics,
    Second Series, 18 (1917), No. 4, 188–193.

    2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``Sturm Sequences 
    and Modified Subresultant Polynomial Remainder Sequences.'' 
    Serdica Journal of Computing, Vol. 8, No 1, 29–46, 2014. 

    3. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``A Basic Result 
    on the Theory of Subresultants.'' Submitted for publication.

    """
    # make sure neither p nor q is 0
    if p == 0 or q == 0:
        return [p, q]

    FLAG = 0
    d0 =  degree(p, x)
    d1 =  degree(q, x)
    # make sure proper degrees  
    if d0 == 0 and d1 == 0:
        return []
    if d1 > d0:
        d0, d1 = d1, d0
        p, q = q, p
    if d0 > 0 and d1 == 0:
        return [p,q]
        
    # make sure LC(p) > 0
    if  LC(p,x) < 0:
        FLAG = 1
        p = -p
        q = -q

    # initialize
    a0, a1 = p, q                           # the input polys
    euclidSeq = [a0, a1]                    # the output list 
    a2 = rem(a0, a1, domain=QQ)             # first remainder
    d2 =  degree(a2, x)                     # degree of a2
    euclidSeq.append( a2 ) 

    # main loop
    while d2 > 0:     
        a0, a1, d0, d1 = a1, a2, d1, d2      # update polys and degrees
        a2 = rem(a0, a1, domain=QQ)          # new remainder
        d2 =  degree(a2, x)                  # actual degree of a2
        euclidSeq.append( a2 ) 

    
    if FLAG:          # change the sign of the sequence
        euclidSeq = [-i for i in euclidSeq]

    # gcd is of degree > 0 ?   
    m = len(euclidSeq)
    if euclidSeq[m - 1] == nan or euclidSeq[m - 1] == 0:
        euclidSeq.pop(m - 1)

    return euclidSeq

def euclid_AMV(f, g, x):
    """
    f, g are polynomials in Z[x] or Q[x].

    Computes the Euclidean sequence of p and q in Z[x] or Q[x]. 
 
    If the Euclidean sequence is complete the coefficients of the polynomials 
    in the sequence are subresultants. That is, they are  determinants of 
    appropriately selected submatrices of sylvester1, Sylvester's matrix of 1840. 
    In this case the Euclidean sequence coincides with the subresultant prs, 
    of the polynomials p, q.

    If the Euclidean sequence is incomplete the signs of the coefficients of the 
    polynomials in the sequence may differ from the signs of the coefficients of
    the corresponding polynomials in the subresultant prs; however, the absolute 
    values are the same. 

    To compute the coefficients, no determinant evaluation takes place. 
    Instead, polynomial divisions in Z[x] or Q[x] are performed, using 
    the function remZ(f, g, x);  the coefficients of the remainders 
    computed this way become subresultants with the help of the 
    Collins-Brown-Traub formula for coefficient reduction.
   
    References:
    ===========  
    1. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``A Basic Result 
    on the Theory of Subresultants.'' Submitted for publication.

    2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``Subresultant Polynomial   
    Remainder Sequences Obtained by Polynomial Divisions in Q[x] or in Z[x].''  
    Submitted for publication.
    
    """
    # make sure neither f nor g is 0
    if f == 0 or g == 0:
        return [f, g]

    d0 =  degree(f, x)
    d1 =  degree(g, x)

        # make sure proper degrees  
    if d0 == 0 and d1 == 0:
        return []
    if d0 > 0 and d1 == 0:
        return [f, g]
    if d1 > d0:
        d0, d1 = d1, d0
        f, g = g, f
    if d0 == 1:
        return [f, g]

    a0 = f
    a1 = g

    euclidSeq = [a0, a1]
    degdifP1, c = degree(a0) - degree(a1) + 1, -1
    
    i = 0                                          # counter for remainders 

        # compute the first polynomial of the prs
    i += 1
    a2 = remZ(a0, a1, x) / Abs( (-1)**degdifP1 )     # first remainder
    euclidSeq.append( a2 )
    d2 =  degree(a2, x)                              # actual degree of a2

    while d2 >= 1:
        a0, a1, d0, d1 = a1, a2, d1, d2       # update polys and degrees
        i += 1
        sigma0 = -LC(a0)
        c = (sigma0**(degdifP1 - 1)) / (c**(degdifP1 - 2))
        degdifP1 = degree(a0, x) - d2 + 1
        a2 = remZ(a0, a1, x) / Abs( ((c**(degdifP1 - 1)) * sigma0) )
        euclidSeq.append( a2 )
        d2 =  degree(a2, x)                   # actual degree of a2
        
    # gcd is of degree > 0 ? 
    m = len(euclidSeq)  
    if euclidSeq[m - 1] == nan or euclidSeq[m - 1] == 0:
        euclidSeq.pop(m - 1)
        
    return euclidSeq


def modified_subresultants_PG(p,q,x):
    """
    p, q are polynomials in Z[x] or Q[x].

    Computes the ``modified'' subresultant prs of p and q in Z[x] or Q[x]; 
    the coefficients of the polynomials in the sequence are 
    ``modified'' subresultants. That is, they are  determinants of appropriately 
    selected submatrices of sylvester2, Sylvester's matrix of 1853. 

    To compute the coefficients, no determinant evaluation takes place. Instead,   
    polynomial divisions in Q[x] are performed, using the function rem(p, q, x);  
    the coefficients of the remainders computed this way become ``modified'' 
    subresultants with the help of the Pell-Gordon Theorem of 1917.
   
    If the ``modified'' subresultant prs is complete, then it coincides with the 
    (generalized) Sturm sequence of the polynomials p, q.

    References:
    ===========  
    1. Pell A. J., R. L. Gordon. The Modified Remainders Obtained in Finding
    the Highest Common Factor of Two Polynomials. Annals of Mathematics,
    Second Series, 18 (1917), No. 4, 188–193.

    2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``Sturm Sequences 
    and Modified Subresultant Polynomial Remainder Sequences.'' 
    Serdica Journal of Computing, Vol. 8, No 1, 29–46, 2014. 
    
    """
    # make sure neither p nor q is 0
    if p == 0 or q == 0:
        return [p, q]

    k =  var('k')                              # index in summation formula
    u_List = []                                # of elements (-1)**u_i 
    subresL = [p, q]                           # mod. subr. prs output list    
    d0 =  degree(p,x)
    d1 =  degree(q,x)

    # make sure proper degrees  
    if d0 == 0 and d1 == 0:
        return []
    if d1 > d0:
        d0, d1 = d1, d0
        p, q = q, p
        subresL = [p, q]                       # mod. subr. prs output list    
    if d0 > 0 and d1 == 0:
        return [p,q]
                
    # initialize
    a0, a1 = p, q                               # the input polys
    del0 = d0 - d1                              # degree difference
    degdif = del0                               # save it
    rho_1 = LC(a0)                              # lead. coeff (a0)

    # Initialize Pell-Gordon variables
    rho_List_Minus_1 =  sign( LC(a0, x))      # sign of LC(a0)
    rho1 =  LC(a1, x)                         # leading coeff of a1
    rho_List = [ sign(rho1)]                  # of signs
    p_List = [del0]                           # of degree differences
    u =  summation(k, (k, 1, p_List[0]))      # value of u
    u_List.append(u)                          # of u values
    v = sum(p_List)                           # v value

    expDeg = d1 - 1                           # expected degree of a2
    a2 = - rem(a0, a1, domain=QQ)             # first remainder
    rho2 =  LC(a2, x)                         # leading coeff of a2
    d2 =  degree(a2, x)                       # actual degree of a2
    degDiff_NEW = expDeg - d2                 # expected - actual degree
    del1 = d1 - d2                            # degree difference

    # mulFac is the factor by which a2 is multiplied to 
    # get integer coefficients
    mulFac_OLD = rho1**(del0 + del1 - degDiff_NEW) 
        
    # update Pell-Gordon variables
    p_List.append(1 + degDiff_NEW)              # degDiff_NEW is 0 for complete seq

    # apply Pell-Gordon formula (7) in second reference
    num = 1                                     # numerator of fraction
    for k in range(len(u_List)):
        num *= (-1)**u_List[k]
    num = num * (-1)**v

    # denominator depends on complete / incomplete seq
    if degDiff_NEW == 0:                        # complete seq
        den = 1
        for k in range(len(rho_List)):
            den *= rho_List[k]**(p_List[k] + p_List[k + 1])
        den = den * rho_List_Minus_1
    else:                                       # incomplete seq
        den = 1
        for k in range(len(rho_List)-1):
            den *= rho_List[k]**(p_List[k] + p_List[k + 1])
        den = den * rho_List_Minus_1 
        expo = (p_List[len(rho_List) - 1] + p_List[len(rho_List)] - degDiff_NEW)
        den = den * rho_List[len(rho_List) - 1]**expo 
    
    # the sign of the determinant depends on sg(num / den)
    if  sign(num / den) > 0:
        subresL.append( simplify(rho_1**degdif*a2* Abs(mulFac_OLD) ) ) 
    else:
        subresL.append(- simplify(rho_1**degdif*a2* Abs(mulFac_OLD) ) )      
  
    # update Pell-Gordon variables
    k =  var('k')
    rho_List.append( sign(rho2))
    u =  summation(k, (k, 1, p_List[len(p_List) - 1]))   
    u_List.append(u)   
    v = sum(p_List)
    
    degDiff_OLD=degDiff_NEW
    
    # main loop
    while d2 > 0:     
        a0, a1, d0, d1 = a1, a2, d1, d2      # update polys and degrees
        del0 = del1                          # update degree difference
        expDeg = d1 - 1                      # new expected degree 
        a2 = - rem(a0, a1, domain=QQ)        # new remainder
        rho3 =  LC(a2, x)                    # leading coeff of a2
        d2 =  degree(a2, x)                  # actual degree of a2
        degDiff_NEW = expDeg - d2            # expected - actual degree
        del1 = d1 - d2                       # degree difference

        # take into consideration the power 
        # rho1**degDiff_OLD that was "left out"
        expo_OLD = degDiff_OLD               # rho1 raised to this power
        expo_NEW = del0 + del1 - degDiff_NEW # rho2 raised to this power

        mulFac_NEW = rho2**(expo_NEW) * rho1**(expo_OLD) * mulFac_OLD 

        # update variables and append
        degDiff_OLD, mulFac_OLD = degDiff_NEW, mulFac_NEW
        rho1, rho2 = rho2, rho3
        
        # update Pell-Gordon variables
        p_List.append(1 + degDiff_NEW)       # degDiff_NEW is 0 for complete seq
            
        # apply Pell-Gordon formula (7) in second reference
        num = 1                              # numerator 
        for k in range(len(u_List)): 
            num *= (-1)**u_List[k]
        num = num * (-1)**v

        # denominator depends on complete / incomplete seq
        if degDiff_NEW == 0:                 # complete seq
            den = 1
            for k in range(len(rho_List)):
                den *= rho_List[k]**(p_List[k] + p_List[k + 1])
            den = den * rho_List_Minus_1
        else:                                # incomplete seq
            den = 1
            for k in range(len(rho_List)-1):
                den *= rho_List[k]**(p_List[k] + p_List[k + 1])
            den = den * rho_List_Minus_1 
            expo = (p_List[len(rho_List) - 1] + p_List[len(rho_List)] - degDiff_NEW)
            den = den * rho_List[len(rho_List) - 1]**expo 

                    
        # the sign of the determinant depends on sg(num / den)          
        if  sign(num / den) > 0:
            subresL.append( simplify(rho_1**degdif*a2* Abs(mulFac_OLD) ) ) 
        else:
            subresL.append(- simplify(rho_1**degdif*a2* Abs(mulFac_OLD) ) )
            
        # update Pell-Gordon variables
        k =  var('k')
        rho_List.append( sign(rho2))
        u =  summation(k, (k, 1, p_List[len(p_List) - 1]))   
        u_List.append(u)   
        v = sum(p_List)

    # gcd is of degree > 0 ?   
    m = len(subresL)
    if subresL[m - 1] == nan or subresL[m - 1] == 0:
        subresL.pop(m - 1)

    return  subresL 

def subresultants_PG(p,q,x):
    """
    p, q are polynomials in Z[x] or Q[x].

    Computes the subresultant prs of p and q in Z[x] or Q[x], from 
    the modified subresultant prs of p and q.
    The coefficients of the polynomials in the two sequences differ only 
    in sign and a power of the absolute value of the leading coefficient 
    of p, as stated in Theorem 2 of the reference.

    The coefficients of the polynomials in the output sequence are 
    subresultants. That is, they are  determinants of appropriately 
    selected submatrices of sylvester1, Sylvester's matrix of 1840. 
  
    If the subresultant prs is complete, then it coincides with the 
    Euclidean sequence of the polynomials p, q.

    References:
    ===========
    1. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ‘‘On the Remainders 
    Obtained in Finding the Greatest Common Divisor of Two Polynomials.'' 
    Serdica Journal of Computing, to appear.
    
    """
    # compute the modified subresultant prs 
    lst = modified_subresultants_PG(p,q,x)  ## or lst = modified_subresultants_AMV(p,q,x)

    # defensive
    if lst == [] or len(lst) == 2:
        return lst

    # the coefficients in lst are modified subresultants and, hence, are 
    # greater than those of the corresponding subresultants by the factor 
    # |LC(lst[0])**( deg(lst[0]) - deg(lst[1]))|; see Theorem 2 in reference.
    lcf = Abs( LC(lst[0])**( degree(lst[0]) - degree(lst[1]) ) )

    # Initialize the subresultant prs list 
    subrSeq = [lst[0], lst[1]]

    # compute the degree sequences m_i and j_i of Theorem 2 in reference.
    degSeq = [degree(Poly(poly, x)) for poly in lst]
    deg = degSeq[0]
    degSeqS = degSeq[1:-1]
    mSeq = [m-1 for m in degSeqS]
    jSeq = [deg - m for m in mSeq]
    
    # compute the AMV factors of Theorem 2 in reference.
    fact = [(-1)**( j*(j-1)/S(2) ) for j in jSeq]

    # shortened list without the first two polys
    lstS = lst[2:]

    # poly lstS[k] is multiplied times fact[k], divided by lcf
    # and appended to the subresultant prs list
    m = len(fact)
    for k in range(m):
        if sign(fact[k]) == -1:
            subrSeq.append(-lstS[k] / lcf)
        else:
            subrSeq.append(lstS[k] / lcf)
    return subrSeq

def subresultants_AMV_Q(p,q,x):
    """
    p, q are polynomials in Z[x] or Q[x].

    Computes the subresultant prs of p and q in Q[x]; 
    the coefficients of the polynomials in the sequence are 
    subresultants. That is, they are  determinants of appropriately 
    selected submatrices of sylvester1, Sylvester's matrix of 1840. 

    To compute the coefficients, no determinant evaluation takes place. 
    Instead, polynomial divisions in Q[x] are performed, using the
    function rem(p, q, x);  the coefficients of the remainders 
    computed this way become subresultants with the help of the 
    Akritas-Malaschonok-Vigklas Theorem of 2015.
   
    If the subresultant prs is complete, then it coincides with the 
    Euclidean sequence of the polynomials p, q.

    References:
    ===========  
    1. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``A Basic Result 
    on the Theory of Subresultants.'' Submitted for publication.

    2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``Subresultant Polynomial   
    Remainder Sequences Obtained by Polynomial Divisions in Q[x] or in Z[x].''  
    Submitted for publication.
    
    """
    # make sure neither p nor q is 0
    if p == 0 or q == 0:
        return [p, q]

    d0 =  degree(p, x)
    d1 =  degree(q, x)

    # make sure proper degrees  
    if d0 == 0 and d1 == 0:
        return []
    if d0 > 0 and d1 == 0:
        return [p, q]
    if d1 > d0:
        d0, d1 = d1, d0
        p, q = q, p
        
        
    # initialize
    i, s = 0, 0                              # counters for remainders & odd elements
    p_odd_index_sum = 0                      # contains the sum of p_1, p_3, etc
    subresL = [p, q]                         # subresultant prs output list
    a0, a1 = p, q                            # the input polys
    sigma1 =  LC(a1, x)                      # leading coeff of a1
    p0 = d0 - d1                             # degree difference
    if p0 % 2 == 1:
        s += 1       
    phi = floor( (s + 1) / 2 )     
    mulFac = 1
    d2 = d1

    # main loop
    while d2 > 0:   
        i += 1
        a2 = rem(a0, a1, domain= QQ)          # new remainder
        if i == 1:
            sigma2 =  LC(a2, x)                       
        else:
            sigma3 =  LC(a2, x)                       
            sigma1, sigma2 = sigma2, sigma3

            
        d2 =  degree(a2, x)                         
        p1 = d1 - d2                                       
        psi = i + phi + p_odd_index_sum

        # new mulFac
        mulFac = sigma1**(p0 + 1) * mulFac
            
        ## compute the sign of the first fraction in formula (9) of the paper
        # numerator 
        num = (-1)**psi     
        # denominator 
        den = sign(mulFac)
         
        # the sign of the determinant depends on sign( num / den ) != 0          
        if  sign(num / den) > 0:
            subresL.append( simplify(expand(a2* Abs(mulFac))))
        else:
            subresL.append(- simplify(expand(a2* Abs(mulFac))))
        

        ## bring into mulFac the missing power of sigma if there was a degree gap
        if p1 - 1 > 0:       
            mulFac = mulFac * sigma1**(p1 - 1)

       # update AMV variables
        a0, a1, d0, d1 = a1, a2, d1, d2       
        p0 = p1                                           
        if p0 % 2 ==1:
            s += 1        
        phi = floor( (s + 1) / 2 )
        if i%2 == 1: 
            p_odd_index_sum += p0             # p_i has odd index
       
    # gcd is of degree > 0 ?   
    m = len(subresL)
    if subresL[m - 1] == nan or subresL[m - 1] == 0:
        subresL.pop(m - 1)

    return  subresL 


def computeSign(base, expo):
    '''
    base != 0 and expo >= 0 are integers;
    
    returns the sign of base**expo without 
    evaluating the power itself!
    '''
    sb = sign(base)
    if sb == 1:
        return 1
    pe = expo % 2
    if pe == 0:
        return -sb
    else:
        return sb

def remZ(p, q, x):
    '''
    Intended mainly for p, q polynomials in Z[x] so that,
    on dividing p by q, the remainder will also be in Z[x]. (However,
    it also works fine for polynomials in Q[x].) 

    It premultiplies p by the _absolute_ value of the leading coefficient 
    of q, raised to the power deg(p) - deg(q) + 1 and then performs 
    polynomial division in Q[x], using the function rem(p, q, x). 

    By contrast the function prem(p, q, x) does _not_ use the absolute 
    value of the leading coefficient of q. 
    This results not only in ``messing up the signs'' of the Euclidean and 
    Sturmian prs's as mentioned in the second reference, 
    but also in violation of the main results of the first and third 
    references --- Theorem 4 and Theorem 1 respectively. Theorems 4 and 1 
    establish a one-to-one correspondence between the Euclidean and the 
    Sturmian prs of p, q, on one hand, and the subresultant prs of p, q, 
    on the other.

    References:
    ===========
    1. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``On the Remainders 
    Obtained in Finding the Greatest Common Divisor of Two Polynomials.'' 
    Serdica Journal of Computing, to appear.

    2. http://planetmath.org/sturmstheorem

    3. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``A Basic Result on 
    the Theory of Subresultants.'' Submitted for publication.

    '''
    delta = (degree(p, x) - degree(q, x) + 1)
    return rem(Abs(LC(q, x))**delta  *  p, q, x) 

 
def subresultants_AMV(f, g, x):
    """
    p, q are polynomials in Z[x] or Q[x].

    Computes the subresultant prs of p and q in Z[x] or Q[x]; 
    the coefficients of the polynomials in the sequence are 
    subresultants. That is, they are  determinants of appropriately 
    selected submatrices of sylvester1, Sylvester's matrix of 1840. 

    To compute the coefficients, no determinant evaluation takes place. 
    Instead, polynomial divisions in Z[x] or Q[x] are performed, using 
    the function remZ(p, q, x);  the coefficients of the remainders 
    computed this way become subresultants with the help of the 
    Akritas-Malaschonok-Vigklas Theorem of 2015 and the Collins-Brown-
    Traub formula for coefficient reduction.
   
    If the subresultant prs is complete, then it coincides with the 
    Euclidean sequence of the polynomials p, q.

    References:
    ===========  
    1. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``A Basic Result 
    on the Theory of Subresultants.'' Submitted for publication.

    2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``Subresultant Polynomial   
    Remainder Sequences Obtained by Polynomial Divisions in Q[x] or in Z[x].''  
    Submitted for publication.
    
    """
    # make sure neither f nor g is 0
    if f == 0 or g == 0:
        return [f, g]

    d0 =  degree(f, x)
    d1 =  degree(g, x)

        # make sure proper degrees  
    if d0 == 0 and d1 == 0:
        return []
    if d0 > 0 and d1 == 0:
        return [f, g]
    if d1 > d0:
        d0, d1 = d1, d0
        f, g = g, f
    if d0 == 1:
        return [f, g]

    a0 = f
    a1 = g

    subresL = [a0, a1]
    degdifP1, c = degree(a0) - degree(a1) + 1, -1
    
        # initialize AMV variables
    sigma1 =  LC(a1, x)                      # leading coeff of a1
    i, s = 0, 0                              # counters for remainders & odd elements
    p_odd_index_sum = 0                      # contains the sum of p_1, p_3, etc
    p0 = degdifP1 - 1
    if p0 % 2 == 1:
        s += 1       
    phi = floor( (s + 1) / 2 )
        # compute the first polynomial of the prs
    i += 1
    a2 = remZ(a0, a1, x) / Abs( (-1)**degdifP1 )     # first remainder
    sigma2 =  LC(a2, x)                       # leading coeff of a2
    d2 =  degree(a2, x)                       # actual degree of a2
    p1 = d1 - d2                              # degree difference
        # sgnDen is the factor, the denominator 1st fraction of (9),  
        # by which a2 is multiplied to get integer coefficients
    sgnDen = computeSign( sigma1, p0 + 1 )
        ## compute sign of the 1st fraction in formula (9) of the paper 
    # numerator
    psi = i + phi + p_odd_index_sum 
    num = (-1)**psi  
    # denominator
    den = sgnDen
        # the sign of the determinant depends on sign(num / den) != 0
    if  sign(num / den) > 0:
        subresL.append( a2 )
    else:
        subresL.append( -a2 )
        # update AMV variables
    if p1 % 2 == 1:
        s += 1       
    ## bring in the missing power of sigma if there was gap
    if p1 - 1 > 0:
        sgnDen = sgnDen * computeSign( sigma1, p1 - 1 )

    while d2 >= 1:
        phi = floor( (s + 1) / 2 )
        if i%2 == 1: 
            p_odd_index_sum += p1             # p_i has odd index
        a0, a1, d0, d1 = a1, a2, d1, d2       # update polys and degrees
        p0 = p1                               # update degree difference
        i += 1

        sigma0 = -LC(a0)
        c = (sigma0**(degdifP1 - 1)) / (c**(degdifP1 - 2))
        degdifP1 = degree(a0, x) - d2 + 1
        a2 = remZ(a0, a1, x) / Abs( ((c**(degdifP1 - 1)) * sigma0) )
        sigma3 =  LC(a2, x)                   # leading coeff of a2
        d2 =  degree(a2, x)                   # actual degree of a2
        p1 = d1 - d2                          # degree difference
        psi = i + phi + p_odd_index_sum
        # update variables 
        sigma1, sigma2 = sigma2, sigma3       
        # new sgnDen
        sgnDen = computeSign( sigma1, p0 + 1 ) * sgnDen            
        ## compute the sign of the first fraction in formula (9) of the paper
        # numerator 
        num = (-1)**psi     
        # denominator 
        den = sgnDen
         
        # the sign of the determinant depends on sign( num / den ) != 0          
        if  sign(num / den) > 0:
            subresL.append( a2 )
        else:
            subresL.append( -a2 )
        
        # update AMV variables
        if p1 % 2 ==1:
            s += 1        
        ## bring in the missing power of sigma if there was gap
        if p1 - 1 > 0:       
            sgnDen = sgnDen * computeSign( sigma1, p1 - 1 )
        
    # gcd is of degree > 0 ?   
    m = len(subresL)
    if subresL[m - 1] == nan or subresL[m - 1] == 0:
        subresL.pop(m - 1)

    return subresL

def modified_subresultants_AMV(p,q,x):
    """
    p, q are polynomials in Z[x] or Q[x].

    Computes the modified subresultant prs of p and q in Z[x] or Q[x], 
    from the subresultant prs of p and q.
    The coefficients of the polynomials in the two sequences differ only 
    in sign and a power of the absolute value of the leading coefficient 
    of p, as stated in Theorem 2 of the reference.

    The coefficients of the polynomials in the output sequence are 
    modified subresultants. That is, they are  determinants of appropriately 
    selected submatrices of sylvester2, Sylvester's matrix of 1853. 
  
    If the modified subresultant prs is complete, then it coincides with the 
    (generalized) Sturm's sequence of the polynomials p, q.

    References:
    ===========
    1. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ‘‘On the Remainders 
    Obtained in Finding the Greatest Common Divisor of Two Polynomials.'' 
    Serdica Journal of Computing, to appear.
    
    """
    # compute the subresultant prs 
    lst = subresultants_AMV(p,q,x)     ## or lst = subresultants_AMV_Q(p, q, x)

    # defensive
    if lst == [] or len(lst) == 2:
        return lst

    # the coefficients in lst are subresultants and, hence, smaller than those 
    # of the corresponding modified subresultants by the factor 
    # Abs(LC(lst[0])**( deg(lst[0]) - deg(lst[1]))); see Theorem 2.
    lcf = Abs( LC(lst[0])**( degree(lst[0]) - degree(lst[1]) ) )

    # Initialize the modified subresultant prs list 
    subrSeq = [lst[0], lst[1]]

    # compute the degree sequences m_i and j_i of Theorem 2
    degSeq = [degree(Poly(poly, x)) for poly in lst]
    deg = degSeq[0]
    degSeqS = degSeq[1:-1]
    mSeq = [m-1 for m in degSeqS]
    jSeq = [deg - m for m in mSeq]
    
    # compute the AMV factors of Theorem 2
    fact = [(-1)**( j*(j-1)/S(2) ) for j in jSeq]

    # shortened list without the first two polys
    lstS = lst[2:]

    # poly lstS[k] is multiplied times fact[k] and times lcf
    # and appended to the subresultant prs list
    m = len(fact)
    for k in range(m):
        if sign(fact[k]) == -1:
            subrSeq.append( simplify(-lstS[k] * lcf) )
        else:
            subrSeq.append( simplify(lstS[k] * lcf) )
    return subrSeq

def correctSign(degF, degG, S1, rdel, cdel):
    """
    Used in various subresultant prs algorithms.

    Evaluates the determinant, (a.k.a. subresultant) of a properly selected
    submatrix of S1, Sylvester's matrix of 1840, to get the correct sign 
    and value of the leading coefficient of a given polynomial remainder.

    degF, degG are the degrees of the original polynomials p, q for which the
    matrix S1 = sylvester(p, q, x, 1) was constructed.
    
    rdel denotes the expected degree of the remainder; it is the number of 
    rows to be deleted from each group of rows in S1 as described in the 
    reference below.

    cdel denotes the expected degree minus the actual degree of the remainder; 
    it is the number of columns to be deleted --- starting with the last column 
    forming the square matrix --- from the matrix resulting after the row deletions.

    References:
    ===========
    Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``Sturm Sequences 
    and Modified Subresultant Polynomial Remainder Sequences.'' 
    Serdica Journal of Computing, Vol. 8, No 1, 29–46, 2014. 

    """
    M = S1[:, :]   # copy of matrix S1

    # eliminate rdel rows from the first degG rows
    for i in range(M.rows - degF - 1, M.rows - degF - rdel - 1, -1):
        M.row_del(i)
    
    # eliminate rdel rows from the last degF rows
    for i in range(M.rows - 1, M.rows - rdel - 1, -1):
        M.row_del(i)

    # eliminate cdel columns
    for i in range(cdel):
        M.col_del(M.rows - 1)
    
    # define submatrix        
    Md = M[:, 0: M.rows]

    return Md.det()

def subresultants_rem(p, q, x):
    """
    p, q are polynomials in Z[x] or Q[x].

    Computes the subresultant prs of p and q in Z[x] or Q[x]; 
    the coefficients of the polynomials in the sequence are 
    subresultants. That is, they are  determinants of appropriately 
    selected submatrices of sylvester1, Sylvester's matrix of 1840. 

    To compute the coefficients polynomial divisions in Q[x] are 
    performed, using the function rem(p, q, x). The coefficients 
    of the remainders computed this way become subresultants by evaluating
    one subresultant per remainder --- that of the leading coefficient. 
    This way we obtain the correct sign and value of the leading coefficient 
    of the remainder and we easily ``force'' the rest of the coefficients 
    to become subresultants.
   
    If the subresultant prs is complete, then it coincides with the 
    Euclidean sequence of the polynomials p, q.

    References:
    ===========  
    1. Akritas, A. G.:``Three New Methods for Computing Subresultant 
    Polynomial Remainder Sequences (PRS’s).'' Serdica Journal of Computing, 
    to appear.

    """
    # make sure neither p nor q is 0
    if p == 0 or q == 0:
        return [p, q]

    f, g = p, q

    n = degF = degree(f, x)
    m = degG = degree(g, x)

    # make sure proper degrees
    if n == 0 and m == 0:
        return []
    if n > 0 and m == 0:
        return [f, g]
    if n < m:
        n, m, degF, degG, f, g = m, n, degG, degF, g, f

    S1 = sylvester(f, g, x, 1)
    SR_list = [f, g]      # subresultant list

    
    while degG > 0:
        # compute the remainder in Z[x] by doing division in Z[x]
        r = rem(p, q, x)        
        # actual degree
        d = degree(r, x)
        # deg(0) == -oo
        if d == -oo:
            return SR_list
        # expected degree
        expDeg = degG - 1
        # make coefficients subresultants
        # by evaluating ONE determinant
        signValue = correctSign(n, m, S1, expDeg, expDeg - d)
        r = simplify((r / LC(r, x)) * signValue)

        # append poly with subresultant coeffs
        SR_list.append(r)
                    
        # update degrees and polys
        degF, degG = degG, d
        p, q = q, r
        
    # gcd is of degree > 0 ? 
    m = len(SR_list)
    if SR_list[m - 1] == nan or SR_list[m - 1] == 0:
        SR_list.pop(m - 1)

    return SR_list

def pivot(M, i, j):
    ''' 
    M is a matrix, and M[i, j] specifies the pivot element.
    
    All elements below M[i, j], in the j-th column, will       
    be zeroed, if they are not already 0, according to       
    Dodgson-Bareiss' integer preserving transormation. 

    References:
    ===========
    1. Akritas, A. G.: ``A new method for computing polynomial greatest 
    common divisors and polynomial remainder sequences.''  
    Numerische Mathematik 52, 119-127, 1988.

    2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``On a Theorem 
    by Van Vleck Regarding Sturm Sequences.'' 
    Serdica Journal of Computing, 7, No 4, 101–134, 2013. 

    '''
    Ma = M[:, :] # copy of matrix M
    rs = Ma.rows # No. of rows
    cs = Ma.cols # No. of cols
    for r in range(i+1, rs):
        if Ma[r, j] != 0:
            for c in range(j + 1, cs):
                Ma[r, c] = Ma[i, j] * Ma[r, c] - Ma[i, c] * Ma[r, j]
            Ma[r, j] = 0
    return Ma

def rotateR(L, k):
    ''' 
    Rotates right by k. L is a row of a matrix or a list.
 
    '''
    LL = list(L)
    if LL == []:
        return []
    for i in range(k):
        el = LL.pop(len(LL) - 1)
        LL.insert(0, el)
    return LL if type(L) is list else Matrix([LL])

def rotateL(L, k):
    ''' 
    Rotates left by k. L is a row of a matrix or a list.

    '''
    LL = list(L)
    if LL == []:
        return []
    for i in range(k):
        el = LL.pop(0)
        LL.insert(len(LL) - 1, el)
    return LL if type(L) is list else Matrix([LL])


def row2poly(row, deg, x):
    '''
    Converts the row of a matrix to a poly of degree deg and variable x.    
    Some entries at the beginning and/or at the end of the row may be zero.

    '''
    k = 0
    poly = []
    leng = len(row)
    
    # find the beginning of the poly ; i.e. the first non-    
    # zero element of the row
    while row[k] == 0:
        k = k + 1
    
    # append the next deg + 1 elements to poly
    for j in range( deg + 1):
        if k + j <= leng:
            poly.append(row[k + j])

    return Poly(poly, x)        
    
def createM(degF, degG, row1, row2, colNum):
    '''
    Creates a ``small'' matrix M to be triangularized.
    
    degF, degG are the degrees of the divident and of the 
    divisor polynomials respectively, degG > degF.
    
    The coefficients of the divident poly are the elements  
    in row2 and those of the divisor poly are the elements 
    in row1.
    
    colNum indicates the number of columns of the matrix M.

    '''
    if degG - degF >= 1:
        return "ERROR! Reverse input degrees"
        # print "ERROR! Reverse input degrees"
        # return

    M = zeros(degF - degG + 2, colNum)
    
    for i in range(degF - degG + 1):
        M[i, :] = rotateR(row1, i)
    M[degF - degG + 1, :] = row2

    return M


def findDegree(M, degF):
    '''
    Finds the degree of the poly corresponding 
    to the _last_ row of the ``small'' matrix M.

    degF is the degree of the divident poly.

    If row is all 0's returns None.

    '''
    j = degF
    for i in range(0, M.cols):
        if M[M.rows - 1, i] == 0:
            j = j - 1
        else:
            return j if j >= 0 else 0

def subresultants_VanVleck(p, q, x, method = 0):
    """
    p, q are polynomials in Z[x] (intended) or Q[x].

    Computes the subresultant prs of p, q by triangularizing,
    in Z[x] or in Q[x], sylvester2, Sylvester's matrix of 1853.
    See references 1 and 2 for Van Vleck's method. 

    Use this version if sylvester2 has small dimensions; otherwise,
    use the second version, subresultants_VanVleck2(p, q, x),
    where sylvester2 is used implicitly.

    Sylvester's matrix sylvester1  is also used to compute one 
    subresultant per remainder; namely, that of the leading 
    coefficient, in order to obtain the correct sign and to  
    force the remainder coefficients to become subresultants.
    
    If the dimensions of sylvester2 are big use method=0 to not 
    print the final, triangularized matrix sylvester2. 
    Use method=1 to print the final matrix.

    If the subresultant prs is complete, then it coincides with the 
    Euclidean sequence of the polynomials p, q.

    References:
    ===========
    1. Akritas, A. G.: ``A new method for computing polynomial greatest 
    common divisors and polynomial remainder sequences.''  
    Numerische Mathematik 52, 119-127, 1988.

    2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``On a Theorem 
    by Van Vleck Regarding Sturm Sequences.'' 
    Serdica Journal of Computing, 7, No 4, 101–134, 2013. 

    3. Akritas, A. G.:``Three New Methods for Computing Subresultant 
    Polynomial Remainder Sequences (PRS’s).'' Serdica Journal of Computing, 
    to appear.


    """
    # make sure neither p nor q is 0
    if p == 0 or q == 0:
        return [p, q]

    f, g = p, q

    n = degF = degree(f, x)
    m = degG = degree(g, x)

    # make sure proper degrees
    if n == 0 and m == 0:
        return []
    if n > 0 and m == 0:
        return [f, g]
    if n < m:
        n, m, degF, degG, f, g = m, n, degG, degF, g, f

    S1 = sylvester(f, g, x, 1)
    S2 = sylvester(f, g, x, 2)

    SR_list = [f, g]      # subresultant list
    colNum = 2 * n
                
    row0 = Poly(f, x, domain = QQ).all_coeffs()
    leng0 = len(row0)
    for i in range(colNum - leng0):
        row0.append(0)
    row0 = Matrix([row0])

    row1 = Poly(g,x, domain = QQ).all_coeffs()
    leng1 = len(row1)
    for i in range(colNum - leng1):
        row1.append(0)
    row1 = Matrix([row1]) 
    
    r = 2    
    if degF - degG > 1:
        # Several last rows in S2 will remain unprocessed
        # Same thing happens when gcd is of deg > 0        
        r = 1
        # insert first poly shifted degF - degG - 1 times
        for i in range(degF - degG - 1):
            S2 = S2.row_insert(i + 1, rotateR(row0, i + 1) )
            r = r + 1
        # insert second poly shifted degF - degG times
        for i in range(degF - degG):
            S2 = S2.row_insert(r + i, rotateR(row1, r + i) )
        r = r + degF - degG

    if degF - degG == 0:
        # the first poly will disappear from S2
        r = 0
   
    while degG > 0:
        
        M = createM(degF, degG, row1, row0, colNum)
        
        # pivot as many times as needed
        for i in range(degF - degG + 1):
            M1 = pivot(M, i, i)
            M = M1[:, :]
            
        d = findDegree(M, degF)
        if d == None:
            if method != 0:
                pprint(S2)
            return SR_list
        expDeg = degG - 1

        # compute subresultant & make poly coeffs subresultants
        poly = row2poly(M[M.rows - 1, :], d, x)
        signValue = correctSign(n, m, S1, expDeg, expDeg - d)
        temp2 = LC(poly, x)
        poly = simplify((poly / LC(poly, x)) * signValue)
        
        # update S2
        for i in range(degG - d):
            # insert first row of M
            S2[r + i, :] = rotateR(M[0, :], r + i)
        r = r + degG - d
        row0 = M[0, :]   

        row1 = rotateL(M[M.rows - 1, :], degF - d)
        row1 = (row1 / temp2) * signValue
        for i in range(degG - d):
            # insert last row of M
            S2[r + i, :] = rotateR(row1, r + i)
        r = r + degG - d

        # update degrees 
        degF, degG = degG, d
        # append poly with subresultant coeffs
        SR_list.append(poly)
                    
    if method != 0:
        pprint(S2)
    return SR_list

def subresultants_VanVleck2(p, q, x):
    """
    p, q are polynomials in Z[x] (intended) or Q[x].

    Computes the subresultant prs of p, q by triangularizing, 
    in Z[x] or in Q[x], all the smaller matrices encountered in the 
    process of triangularizing sylvester2, Sylvester's matrix of 1853. 
    See references 1 and 2 for Van Vleck's method. 

    If the sylvester2 matrix has big dimensions use this version,
    where sylvester2 is used implicitly. If you want to see the final,
    triangularized matrix sylvester2, then use the first version,
    subresultants_VanVleck(p, q, x).

    sylvester1, Sylvester's matrix of 1840, is also used to compute 
    one subresultant per remainder; namely, that of the leading 
    coefficient, in order to obtain the correct sign and to  
    ``force'' the remainder coefficients to become subresultants.
    
    If the subresultant prs is complete, then it coincides with the 
    Euclidean sequence of the polynomials p, q.

    References:
    ===========
    1. Akritas, A. G.: ``A new method for computing polynomial greatest 
    common divisors and polynomial remainder sequences.''  
    Numerische Mathematik 52, 119-127, 1988.

    2. Akritas, A. G., G.I. Malaschonok and P.S. Vigklas: ``On a Theorem 
    by Van Vleck Regarding Sturm Sequences.'' 
    Serdica Journal of Computing, 7, No 4, 101–134, 2013. 

    3. Akritas, A. G.:``Three New Methods for Computing Subresultant 
    Polynomial Remainder Sequences (PRS’s).'' Serdica Journal of Computing, 
    to appear.

    """
    # make sure neither p nor q is 0
    if p == 0 or q == 0:
        return [p, q]

    f, g = p, q

    n = degF = degree(f, x)
    m = degG = degree(g, x)

    # make sure proper degrees
    if n == 0 and m == 0:
        return []
    if n > 0 and m == 0:
        return [f, g]
    if n < m:
        n, m, degF, degG, f, g = m, n, degG, degF, g, f

    S1 = sylvester(f, g, x, 1)

    SR_list = [f, g]      # subresultant list
    # cols of sylvester2
    colNum = 2 * n

    row0 = Poly(f, x, domain = QQ).all_coeffs()
    leng0 = len(row0)
    for i in range(colNum - leng0):
        row0.append(0)
    row0 = Matrix([row0])

    row1 = Poly(g,x, domain = QQ).all_coeffs()
    leng1 = len(row1)
    for i in range(colNum - leng1):
        row1.append(0)
    row1 = Matrix([row1])
    
    while degG > 0:
        
        M = createM(degF, degG, row1, row0, colNum)
            
        for i in range(degF - degG + 1):
            M1 = pivot(M, i, i)
            M = M1[:, :]
            
        d = findDegree(M, degF)
        if d == None:
            return SR_list
        expDeg = degG - 1
        # make entries of S2[r + 1, :] subresultants
        poly = row2poly(M[M.rows - 1, :], d, x)
        signValue = correctSign(n, m, S1, expDeg, expDeg - d)
        poly = simplify((poly / LC(poly, x)) * signValue)

        # append poly with subresultant coeffs
        SR_list.append(poly)
                    
        # update degrees and rows
        degF, degG = degG, d
        row0 = row1
        row1 = Poly(poly, x, domain = QQ).all_coeffs()
        leng1 = len(row1)
        for i in range(colNum - leng1):
            row1.append(0)
        row1 = Matrix([row1])

    return SR_list
 


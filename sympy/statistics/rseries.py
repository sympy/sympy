from sympy import Poly, Integer, factorial, Mul, Subs


def rseries(fct, x, pt, n):
    """
    Develop inverse series of a function.
    fct: function
    x: function variable
    pt: point of the approximation.
    n: Degree of resulting polynomial.
    fct must be differentiable n times.
    The first derivative at point pt may not be zero.
    """

    def split_fct(fct, var):
        """Split Mul by variable, if possible"""
        constant = Integer(1)
        factor = Integer(1)
        if isinstance(fct,Mul):
            for f in fct.args:
                if var in f:
                    factor*= f
                else:
                    constant*= f
        else:
            factor = fct
        return factor, constant

    def allsolutions(n):
        """
        Return all possible solutions for positive integers s1,s2,s3,...,sn
        that satisfy the equation:
        s1 + 2*s2 + 3*s3 + ... + n*sn = n
        for a given positive integer n.
        """
        if n == 0:
            return [[0],]
        comb = []
        n = int(n)
        assert(n>=0)
        factors = [i+1 for i in range(n)]

        def combinations(n, fpos, path=[]): #recursively append combinations to comb
            if fpos>n-1:
                return
            #print "Called with: ", n, fpos
            fctr = factors[fpos]
            #print "NEW FCT: ", "Fpos: ", fpos
            #print "n: ", n
            current = n
            while(current>=0):
                if fctr*current == n:
                    #print "Append: ", "fctr:",fctr, "current: ",current
                    comb.append(path+[current])
                if fctr*current < n:
                    combinations(n-fctr*current,fpos+1,path+[current])
                current-=1

        combinations(n,0)
        return comb

    def make_coefficient( a, k):
        """Get coefficient degree k"""
        #print "_"*50
        sol = allsolutions(k-1)

        c0 = 1/(k*a[0]**k)
        ret = Integer(0)
        for s in sol:
            c1 = Integer(0)
            c2 = Integer(1)
            c3 = Integer(1)
            c4 = Integer(1)
            for _s in s:
                c1 +=_s
            for _s in s:
                c2*=factorial(_s)
            for i in range(k,k+c1):
                c3*=i
            for i in range(len(s)):
                if (s[i] != 0):
                    c4 *= (a[i+1]/a[0])**Integer(s[i])
            
            #print "c0: ", c0, "c1: ", c1 , "c2: ", c2 , "c3: ", c3 , "c4:  ", c4
            ret += (-1)**(c1)*(c3)/(c2)*c4

        ret = c0*ret
        print "Coefficient: %s: %s" % (k, ret)
        return ret

    factor, constant = split_fct(fct,x)
    a=[]
    for i in range(1,n+1):
        a.append( Integer(1)/factorial(i)*factor.diff(x,i).subs(x,pt) )

    koeffs = []
    for i in range(1,n+1):
        koeffs.append( make_coefficient( a, i ) )
    koeffs.reverse()
    koeffs+=[0]
    # Poly translated to be correct inverse
    #return Poly(koeffs,x).as_expr().subs({x,x/constant}).subs(x,x-fct.subs(x,pt))+pt
    return Subs(Poly(koeffs,x).as_expr(),x,x/constant-fct.subs(x,pt)) + pt



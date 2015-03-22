###Phase II

Implement multivariate Series class in Sympy.

Series expansion in Sympy suffers from the following issues:   
1. All types of series are grouped under the same name  
2. Series expansion is slow.  
```
%time s = (1/cos(x/log(x))).series(x, 0, 10)
CPU times: user 2.75 s, sys: 50.4 ms, total: 2.8 s
Wall time: 4.1 s
```
It takes 4 seconds!
```
In [89]: t1 = clock()
    ...: h = exp(x).series(x, 0, 10)
    ...: t2 =clock()
    ...: t2-t1
    ...: 
Out:[89]: 0.5608868598937988
```
Now, 0.56s might not seem much, but when used in an algorithm that needs repeated calls for series expansion,
for example PyDy, the total time taken is exorbitant.
3. Series expansions are stored as sum of core objects
```
In [48]: from sympy.polys.ring_series import rs_exp, rs_mul
In [49]: q = rs_exp(x**2, x, 101)  #Series of exp(x**2) stored in a sparse dictionary
In [50]: t1 = clock()
    ...: b = rs_mul(q,q,y,101)     #Sparse multiplication
    ...: t2 =clock()
    ...: t2-t1
    ...: 
Out[50]: 0.0010459423065185547
In [56]: p = exp(y**2).series(y, 0, 100)  #Series of exp(y**2) stored as a core object
In [57]: t1 = clock()
    ...: a = expand(p*p) 
    ...: t2 =clock()
    ...: t2-t1
    ...: 
Out[53]: 0.0587308406829834    
```
Manipulation over core objects is slow by a factor of 50. This, of course is a slightly sparse series, 
however sparse multiplication will still be faster with dense series.
4. Core has to deal with order

I propose the following solutions to these issues:

1. Implement a class based representation:
     FormalSeries
     FormalPowerSeries
     FormalLaurentSeries
     FormalPuiseuxSeries
A semi-working example involving univariate series on QQ:  
Note: This example covers only the univariate case. However, the final
    implementation will support multivariate series as well.
```
    class FormalSeries(object):
    def __init__(cls, data):
        pass

class FormalPowerSeries(FormalSeries):
    def __init__(self, data):
        if isinstance(data, PolyElement):
            self.series = data
            self.ring = data.ring
        else:
             self.ring, self.series = sring(data, domain=QQ)

    def __repr__(self):
        from sympy.printing import sstr
        return sstr(self.series, order=None)

    def __str__(self):
        from sympy.printing import sstr
        return sstr(self.series, order=None)

    def __add__(self, other):
        return FormalPowerSeries(self.series + other.series)

    def __mul__(self, other):
        if(self.ring.ngens == 1 and other.ring.ngens == 1):
            x = ((self.ring).gens)[0]
            prec = min((self.series).degree(), (other.series).degree()) + 1
            return FormalPowerSeries(rs_mul(self.series, other.series, x, prec))

    def __pow__(self, n):
        if(self.ring.ngens == 1):
            x = ((self.ring).gens)[0]
            prec = (self.series).degree() + 1
            return FormalPowerSeries(rs_pow(self.series, n, x, prec))

    def __div__(self, other):
        if(self.ring.ngens == 1 and other.ring.ngens == 1):
            x = ((self.ring).gens)[0]
            return self * other**(-1)

class FormalLaurentSeries(FormalSeries):
    def __init__(self, data):          #Converts a laurent series into a power series by multiplying the 
        x = (data.atoms(Symbol)).pop() #  whole series by the lowest negative exponent.
        lead = data.extract_leading_order(x)
        if lead < 0:
            self.ring, self.series = sring(data * (x**(-lead)), domain=QQ)
        else FormalPowerSeries.__init__(data)
```
Ultimately we should be able to initialise the series as
```
x = Symbol('x')
expr = x**2 + 2*x + 1
f = FormalPowerSeries(expr, gens=x, prec=3, domain=QQ) 
f.ring
Out[78]: Polynomial ring in x over QQ with lex order
or

R, x = Ring("x", QQ)
f = FormalPowerSeries(x**2 + 2*x + 1, prec=3)
```
2. Implement sparse representation for series, on top of polynomial rings in
   ring_series using a dictionary. For univariate monomials the key is just the
   exponent of the monomial.(Do we need Kronecker's trick for use as dict key in
   Sympy?) Kronecker's substitution trick will further speed it up.  
3. A class based approach is superior as it provides more flexibility and type
   safety. The user knows beforehand what type of series to expect. This also
   solves [4] as the core need not bother about printing and handling the order
   of the series.
4. Core need not bother about the order of the series. The series instance
   stores its order. When required to be printed, the series can be printed
   along with its order term in the proper format.

####About ring_series

Sparse power series based polynomial rings is partially implemented in 
polys/ring_series.py based on the work of Mario Pernici. However,
 a greater part of his work is unmerged (https://github.com/sympy/sympy/pull/609/files).
* The polynomial uses distributed format, i.e , it is stored as sum of
monomials.
* It supports multivariate series as it uses tuples to store the exponenets. So
`x**2*y**3` is stored as (2,3)

 I will build on his representation and complete it. Such an approach will give
 the best speed[2]. Thus, I will

1. Add the series of the following functions: sin, cos, tan
sinh, cosh, tanh lambert
2. Add more user functions.

####Return Type

There are two ways to return any expansion. We either return a specified number
of terms or else we return a generator that can return the next term as many
times as the user wants. The former is implemented in Sympy as nseries while the
latter, called lazy expansion is implemented as lseries. Sympy's lseries is
rather slow as by default, it internally uses nseries repeatedly.

One issue with nseries is that given a nested function, you never know how many
terms you need to expand in the innermost function to finally get the desired
order. There are two ways to tackle this:

Start with 'n' terms in the innermost part and if the order falls short, keep
increasing 'n' till you get the desired order.  Start with 'n' terms in the
innermost part and output whatever ordered expansion you get. This method would
be much faster than the 1st one but there is no guarantee you'd get the
requested order.

####References
1. http://www.cs.berkeley.edu/~fateman/papers/polysbyGMP.pdf
2. https://groups.google.com/forum/#!searchin/sympy/series/sympy/hiRuIHa8ImA/JLOBsMr9yUcJ  
3. http://link.springer.com/chapter/10.1007%2F978-3-319-02297-0_8#page-1: Contains definitions and algorithm for multiplication of DMP (Distributed Multivariate Polynomials)
4. http://www.cybertester.com/data/gruntz.pdf
5. https://mattpap.github.io/masters-thesis/html/src/conclusions.html#polynomial-arithmetics
6. https://groups.google.com/forum/#!topic/sympy/TVEp72mZ3Uo
7. https://mattpap.github.io/masters-thesis/html/src/internals.html
8. https://groups.google.com/forum/#!topic/sage-gsoc/WbmAJAaGlhs
9. http://arxiv.org/pdf/0712.4046.pdf?origin=publication_detail

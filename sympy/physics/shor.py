from sympy import Expr, sympify, Add, Mul, Pow, I, Function, Integer, S, sympify, Matrix, elementary
from sympy.core.numbers import *
from sympy.core.basic import S, sympify
from sympy.core.function import Function
from sympy.functions.elementary.exponential import *
from sympy.functions.elementary.miscellaneous import *
from sympy.matrices.matrices import *
from sympy.simplify import *
from sympy.core.symbol import *
from sympy.physics.qbit import *
import math

#This uses Euclid's Algorithm to find the gcd of a and b  
def gcd(a,b):
    remainder = a
    while remainder >= b:
        remainder = remainder - b
    if remainder == 0:
        return b
    else:
        return gcd(b, remainder)

def shor(N):
    a = random.randrange(N-2)+2
    if gcd(N,a) != 1:
        print "got lucky with rand"
        return gcd(N,a)
    print "a= ",a
    print "N= ",N
    r = periodfind(a,N)
    print "r= ",r
    if r%2 == 1:
        print "r not even, begin again"
        shor(N)
    if a**(r/2)%N == -1:
        print "r not meet prereques, begin again"
        shor(N)
    answer = (gcd(a**(r/2)-1, N), gcd(a**(r/2)+1, N))
    return answer

def arr(num, t):
    car = []
    for i in reversed(range(t)):
        car.append((num>>i)&1)
    return car

def getr(x, y, N):
    fraction = continuedFraction(x,y)
    #now convert into r
    total = ratioize(fraction, N)
    return total

def ratioize(list, N):
    if list[0] > N:
        return 0
    if len(list) == 1:
        return list[0]
    return list[0] + ratioize(list[1:], N)
     
def continuedFraction(x, y):
    x = int(x)
    y = int(y)
    temp = x/y
    if temp*y == x:
        return [temp,]

    list = continuedFraction(y, x-temp*y)
    list.insert(0, temp)
    return list

    
# This is quantum part of Shor's algorithm
# Takes two registers, puts first in superposition of states with Hadamards so: |k>|0> with k being all possible choices 
def periodfind(a, N):
    epsilon = .5
    #picks out t's such that maintains accuracy within epsilon
    t = int(2*math.ceil(log(N,2)))
    # make the first half of register be 0's |000...000>
    start = [0 for x in range(t)]
    #Put second half into superposition of states so we have |1>x|0> + |2>x|0> + ... |k>x>|0> + ... + |2**n-1>x|0>
    factor = 1/sqrt(2**t)
    qbits = 0
    for i in range(2**t):
        qbitArray = arr(i, t) + start
        qbits = qbits + Qbit(*qbitArray)
    circuit = (factor*qbits).expand()
    #Controlled second half of register so that we have:
    # |1>x|a**1 %N> + |2>x|a**2 %N> + ... + |k>x|a**k %N >+ ... + |2**n-1=k>x|a**k % n>
    circuit = controlledMod(t,a,N)*circuit
    #will measure first half of register giving one of the a**k%N's 
    circuit = apply_gates(circuit)
    print circuit
    print "controlled Mod'd"
    for i in range(t):
        circuit = measure(i)*circuit
    circuit = apply_gates(circuit)
    print circuit
    print "measured 1"
    #Now apply Inverse Quantum Fourier Transform on the second half of the register
    circuit = apply_gates(QFT(t, t*2).decompose()*circuit, floatingPoint = True)
    print circuit
    print "QFT'd"
    for i in range(t):
        circuit = measure(i+t)*circuit
    circuit = apply_gates(circuit)
    print circuit
    if isinstance(circuit, Qbit):
        register = circuit
    elif isinstance(circuit, Mul):
        register = circuit.args[-1]
    else:
        register = circuit.args[-1].args[-1]
    
    print register
    n = 1
    answer = 0
    for i in range(len(register)/2):
        answer += n*register[i+t]
        n = n<<1
    if answer == 0:
        raise OrderFindingException("Order finder returned 0. Happens with chance %f" % epsilon)
    #turn answer into r using continued fractions
    g = getr(answer, 2**t, N)
    print g
    return g
        
class OrderFindingException(Exception):
    pass


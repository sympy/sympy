from sympy import mpmath as mp
from sympy.functions.elementary.trigonometric import asin,atan
from sympy.functions.elementary.complexes import re,im
from sympy import N,I

'''
Inverse Jacobi functions
sn-1
ns-1
dn-1
nd-1
cs-1
sc-1
sd-1
ds-1
cn-1
nc-1
cd-1
dc-1
 are defined here only for the real arguements.Though the Code wirtten for sn-1 works for
 complex/real numbers the other functions which are based on sn-1 doesnot yeild result.
 
 sn-1 = F(psi,m)
 and entire Documentation is available at
 
REFERENCE : https://researchspace.auckland.ac.nz/bitstream/handle/2292/5042/390.pdf?sequence=1

USAGE
    >>>from sympy import inversej
    >>>inversej('sn',0.3,0.4)
    0.306574557462381
    >>>inversej('ds','4','0.5')
    3.45812493060054
To call Inverse jacobian 'sn' for complex numbers
    >>>from sympy import sni
    >>>sni(4+3*I,0.4)
    0.250044686201085 + 1.75297022418506*I

TODO :
Analyse the function behaviour for Complex inputs

'''

def sni(z=None,l= None):
    
    if 0<abs(z)<1:
        
        a = N(asin(z))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(mp.ellipf(r+i*mp.j,l))
    
    elif 1<abs(z)<(1.0/(l**0.5)):
        
        a = N(asin(((1-z**-2)/(1-l))**0.5))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(mp.ellipk(l) + mp.j*mp.ellipf(r+i*mp.j,1-l))
    
    elif (1.0/(l**0.5))<abs(z):
        
        a = N(asin(1.0/(z*(l)**0.5)))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(mp.ellipf(r+i*mp.j,l) + mp.j*mp.ellipk(1-l))
        
def nsi(z=None,l=None):

    if z>0:
        return sni(1.0/z,l)
    elif z==0:
        return N(I*(mp.ellipk(1-l)))
    else:
        return N(sni(1.0/z,l) + 2*I*(mp.ellipk(1-l)))
        
def cni(z=None,l=None):
    
    if 0<=z<1:
        a = N(asin((1-z**2)**0.5))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(mp.ellipf(r+i*mp.j,l))
        
    elif -1<=z<0:
        a = N(asin((1-z**2)**0.5))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(2*mp.ellipk(l)-mp.ellipf(a,l))
        
    elif z>=1:
        a = N(asin((1-z**-2)**0.5))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(I*mp.ellipf(a,1-l))
        
    else:
        a = N(asin((1-z**-2)**0.5))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(2*mp.ellipk(l)-I*mp.ellipf(a,1-l))
        
def nci(z=None,l=None):

    if z>0:
        return cni(1.0/z,l)
    if z == 0:
        return N(I*mp.ellipk(1-l))
    if z<0:
        return N(cni(1.0/z,l)-2*mp.ellipk(l) + I*2*mp.ellipk(1-l))
        
def dni(z=None,l=None):
    
    if (1-l)**0.5 <= z < 1:
        a = ((1-z**2)/l)**0.5
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(mp.ellipf(r+i*mp.j,l))
        
    elif z >=1:
        a = N(asin(((1-z**-2)/(1-(1-l)*(z**-2)))**0.5))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(I*mp.ellipf(r+i*mp.j,1-l))
        
    elif 0<=z< (1-l)**0.5 :
        a = N(asin(((1-(z**2)/(1-l))/(1-z**2))**0.5))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(mp.ellipk(l)+I*mp.ellipf(r+i*mp.j,1-l))
        
    elif -(1-l)**0.5<=z<0:
        a = N(asin(((1-(z**2)/(1-l))/(1-z**2))**0.5))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(mp.ellipk(l)+I*(2*mp.ellipk(1-l)-mp.ellipf(r+i*mp.j,1-l)))
        
    elif -1<=z<-(1-l)**0.5:
        a = N(asin(((1-z**2)/l)**0.5))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(2*mp.ellipk(l)-mp.ellipf(r+i*mp.j,l)+I*mp.ellipk(1-l))
        
    else:
        a = N(asin(((1-z**-2)/(1-(1-l)*(z**-2)))**0.5))
        i = mp.mpf(im(a))
        r = mp.mpf(re(a))
        return N(2*mp.ellipk(l)+I*(2*mp.ellipk(1-l)-mp.ellipf(r+i*mp.j,1-l)))
        
def ndi(z=None,l=None):
    
    if z>0:
        return dni(1.0/z,l)
    if z==0:
        return N(I*mp.ellipk(1-l))
    if z<0:
        return N(dni(1.0/z,l)-2*mp.ellipk(l))
        
def sci(z=None,l=None):
    
    a = N(atan(z))
    i = mp.mpf(im(a))
    r = mp.mpf(re(a))
    return N(mp.ellipf(r+i*mp.j,l))
    
def csi(z=None,l=None):
    
    if z>0:
        return N(mp.ellipf(atan(1.0/z),l))
    if z==0:
        return N(mp.ellipk(l))
    if z<0:
        return N(mp.ellipf(atan(1.0/z),l)+2*mp.ellipk(l))
        
def cdi(z=None,l=None):
    return N(sni(z,l)-mp.ellipk(l))
    
def dci(z=None,l=None):
    return N(nsi(z,l)-mp.ellipk(l))
    
def sdi(z=None,l=None):
    return N(cni(z*(1-l)**0.5,l)+mp.ellipk(l))
    
def dsi(z=None,l=None):
    return N(nci(z/(1-l)**0.5,l)+mp.ellipk(l))
    
inverse_jacobi = {
    'sn': sni,
    'ns': nsi,
    'dn': dni,
    'nd': ndi,
    'cn': cni,
    'nc': nci,
    'sc': sci,
    'cs': csi,
    'cd': cdi,
    'dc': dci,
    'sd': sdi,
    'ds': dsi
}
def inversej(key,z=None,l=None):
    if im(z)==0 and im(l)==0:
        value = inverse_jacobi[key]
        if (value == dni or value == ndi) and l<1:
            return value(z,l)
        else:
            return value(z,l)
    else:
        raise ValueError("Enter only real values")
# File: SineIntegrationRules.py
#
# This file is a Python code generated from the publicly available rule-based
# integrator (Rubi) system created by Albert Rich available at
# http://www.apmaths.uwo.ca/~arich/. This file is under the modified BSD
# license (see the LICENSE file for more details).

from sympy import Integer, sin, cos
from .rubi_definitions import eq, gt, lt, ge, le, integer, subst, integrate

def intsin5(a,b,c,n,x):
    if eq(n,0) or eq(b,0) or eq(c,0):
        return x*(c*sin(a+b*x))**n
    else:
        if ge(n,1):
            if eq(n,1):
                return -c*cos(a+b*x)/b
            else:
                if integer(Integer(1)/Integer(2)*(Integer(-1)+n)):
                    return -c**n*subst(integrate((Integer(1)-x**2)**(Integer(1)/Integer(2)*(Integer(-1)+n)),x),x,cos(a+b*x))/b
                else:
                    return c**2*(Integer(-1)+n)*intsin5(a,b,c,Integer(-2)+n,x)/n-c*cos(a+b*x)*(c*sin(a+b*x))**(Integer(-1)+n)/(b*n)
        else:
            if le(n,-1):
                if eq(n,-1):
                    return intcsc5(a,b,1/c,Integer(1),x)
                else:
                    return (Integer(2)+n)*intsin5(a,b,c,Integer(2)+n,x)/(c**2*(Integer(1)+n))+cos(a+b*x)*(c*sin(a+b*x))**(Integer(1)+n)/(b*c*(Integer(1)+n))
            else:
                if eq(n**2,Integer(1)/Integer(4)):
                    if eq(c,1):
                        if eq(n,Integer(1)/Integer(2)):
                            return -Integer(2)*ellipe(Integer(1)/Integer(4)*Pi+Integer(1)/Integer(2)*(-a-b*x),Integer(2))/b
                        else:
                            return -Integer(2)*ellipf(Integer(1)/Integer(4)*Pi+Integer(1)/Integer(2)*(-a-b*x),Integer(2))/b
                    else:
                        return intsin5(a,b,Integer(1),n,x)*(c*sin(a+b*x))**n/sin(a+b*x)**n
                else:
                    if gt(n,0):
                        return -c*cos(a+b*x)*hyp2f1(Integer(1)/Integer(2),Integer(1)/Integer(2)*(Integer(1)-n),Integer(3)/Integer(2),cos(a+b*x)**2)*(c*sin(a+b*x))**(Integer(-1)+n)*(sin(a+b*x)**2)**(Integer(1)/Integer(2)*(Integer(1)-n))/b
                    else:
                        return -cos(a+b*x)*hyp2f1(Integer(1)/Integer(2),Integer(1)/Integer(2)*(Integer(1)-n),Integer(3)/Integer(2),cos(a+b*x)**2)*(c*sin(a+b*x))**(Integer(1)+n)*(sin(a+b*x)**2)**(Integer(1)/Integer(2)*(Integer(-1)-n))/(b*c)


def intsin6(a,b,c,d,n,x):
    if eq(n,0) or eq(b,0) or eq(d,0):
        return x*(a+b*sin(c+d*x))**n
    else:
        if eq(a,0):
            return intsin5(c,d,b,n,x)
        else:
            if eq(n,1):
                return a*x+intsin5(c,d,b,n,x)
            else:
                if eq(n,2):
                    return intsin6(a**2+Integer(1)/Integer(2)*b**2,Integer(2)*a*b,c,d,Integer(1),x)-Integer(1)/Integer(2)*b**2*cos(c+d*x)*sin(c+d*x)/d
                else:
                    if eq(a**2-b**2,0):
                        if integer(Integer(2)*n):
                            if gt(n,0):
                                if eq(n,Integer(1)/Integer(2)):
                                    return -Integer(2)*b*cos(c+d*x)/(d*sqrt(a+b*sin(c+d*x)))
                                else:
                                    return a*(Integer(-1)+Integer(2)*n)*intsin6(a,b,c,d,Integer(-1)+n,x)/n-b*cos(c+d*x)*(a+b*sin(c+d*x))**(Integer(-1)+n)/(d*n)
                            else:
                                if eq(n,-1):
                                    return -cos(c+d*x)/(d*(b+a*sin(c+d*x)))
                                else:
                                    if eq(n,Integer(-1)/Integer(2)):
                                        return -Integer(2)*subst(integrate(1/(Integer(2)*a-x**2),x),x,b*cos(c+d*x)/sqrt(a+b*sin(c+d*x)))/d
                                    else:
                                        return (Integer(1)+n)*intsin6(a,b,c,d,Integer(1)+n,x)/(a*(Integer(1)+Integer(2)*n))+b*cos(c+d*x)*(a+b*sin(c+d*x))**n/(a*d*(Integer(1)+Integer(2)*n))
                        else:
                            return a*cos(c+d*x)*hyp2f1(Integer(1)/Integer(2),Integer(1)/Integer(2)+n,Integer(3)/Integer(2)+n,Integer(1)/Integer(2)*(a+b*sin(c+d*x))/a)*(a+b*sin(c+d*x))**n*sqrt(Integer(2))/(b*d*(Integer(1)+Integer(2)*n)*sqrt((a-b*sin(c+d*x))/a))
                    else:
                        if integer(Integer(2)*n):
                            if gt(n,0):
                                if eq(n,Integer(1)/Integer(2)):
                                    if gt(a+b,0):
                                        return -Integer(2)*ellipe(Integer(1)/Integer(4)*Pi+Integer(1)/Integer(2)*(-c-d*x),Integer(2)*b/(a+b))*sqrt(a+b)/d
                                    else:
                                        if gt(a-b,0):
                                            return Integer(2)*ellipe(Integer(1)/Integer(4)*Pi+Integer(1)/Integer(2)*(c+d*x),-Integer(2)*b/(a-b))*sqrt(a-b)/d
                                        else:
                                            return intsin6(a/(a+b),b/(a+b),c,d,Integer(1)/Integer(2),x)*sqrt(a+b*sin(c+d*x))/sqrt((a+b*sin(c+d*x))/(a+b))
                                else:
                                    return intsin9(a,b,b**2*(Integer(-1)+n)+a**2*n,a*b*(Integer(-1)+Integer(2)*n),c,d,Integer(-2)+n,Integer(1),x)/n-b*cos(c+d*x)*(a+b*sin(c+d*x))**(Integer(-1)+n)/(d*n)
                            else:
                                if eq(n,-1):
                                    if False:
                                        return Integer(2)*subst(integrate(1/(a+b+(a-b)*x**2),x),x,tan(Integer(1)/Integer(2)*(c-Integer(1)/Integer(2)*Pi+d*x)))/d
                                    else:
                                        return Integer(2)*subst(integrate(1/(a+Integer(2)*b*x+a*x**2),x),x,tan(Integer(1)/Integer(2)*(c+d*x)))/d
                                else:
                                    if eq(n,Integer(-1)/Integer(2)):
                                        if gt(a+b,0):
                                            return -Integer(2)*ellipf(Integer(1)/Integer(4)*Pi+Integer(1)/Integer(2)*(-c-d*x),Integer(2)*b/(a+b))/(d*sqrt(a+b))
                                        else:
                                            if gt(a-b,0):
                                                return Integer(2)*ellipf(Integer(1)/Integer(4)*Pi+Integer(1)/Integer(2)*(c+d*x),-Integer(2)*b/(a-b))/(d*sqrt(a-b))
                                            else:
                                                return intsin6(a/(a+b),b/(a+b),c,d,Integer(-1)/Integer(2),x)*sqrt((a+b*sin(c+d*x))/(a+b))/sqrt(a+b*sin(c+d*x))
                                    else:
                                        return intsin9(a,b,a*(Integer(1)+n),-b*(Integer(2)+n),c,d,Integer(1)+n,Integer(1),x)/((a**2-b**2)*(Integer(1)+n))-b*cos(c+d*x)*(a+b*sin(c+d*x))**(Integer(1)+n)/((a**2-b**2)*d*(Integer(1)+n))
                        else:
                            return appellf1(Integer(1)+n,Integer(1)/Integer(2),Integer(1)/Integer(2),Integer(2)+n,(a+b*sin(c+d*x))/(a-b),(a+b*sin(c+d*x))/(a+b))*(a+b*sin(c+d*x))**(Integer(1)+n)*sqrt(b*(Integer(1)-sin(c+d*x))/(a+b))*sqrt(-b*(Integer(1)+sin(c+d*x))/(a-b))/(b*d*(Integer(1)+n)*cos(c+d*x))


def intsin9(a,b,c,d,e,f,m,n,x):
    if eq(f,0):
        return x*(a+b*sin(e))**m*(c+d*sin(e))**n
    else:
        if eq(m,0) or eq(b,0):
            return intsin6(c,d,e,f,n,x)*(a+b*sin(e+f*x))**m
        else:
            if eq(n,0) or eq(d,0):
                return intsin6(a,b,e,f,m,x)*(c+d*sin(e+f*x))**n
            else:
                if eq(b*c-a*d,0):
                    if integer(m) or gt(b/d,0):
                        return (b/d)**m*intsin6(c,d,e,f,m+n,x)
                    else:
                        return intsin6(c,d,e,f,m+n,x)*(a+b*sin(e+f*x))**m/(c+d*sin(e+f*x))**m
                else:
                    if eq(m,1) and eq(n,1):
                        return Integer(1)/Integer(2)*(Integer(2)*a*c+b*d)*x-(b*c+a*d)*cos(e+f*x)/f-Integer(1)/Integer(2)*b*d*cos(e+f*x)*sin(e+f*x)/f
                    else:
                        if eq(m,1) and eq(n,-1):
                            return b*x/d-(b*c-a*d)*intsin6(c,d,e,f,Integer(-1),x)/d
                        else:
                            if eq(m,-1) and eq(n,1):
                                return intsin9(c,d,a,b,e,f,n,m,x)
                            else:
                                if eq(a**2-b**2,0) and eq(c**2-d**2,0):
                                    if integer(m):
                                        if integer(Integer(1)/Integer(2)+n):
                                            if gt(n,0):
                                                if eq(n,Integer(1)/Integer(2)):
                                                    return -Integer(2)*a**m*c**m*d*cos(e+f*x)**(Integer(1)+Integer(2)*m)*(c+d*sin(e+f*x))**(Integer(-1)/Integer(2)-m)/(f*(Integer(1)+Integer(2)*m))
                                                else:
                                                    return c*(Integer(-1)+Integer(2)*n)*intsin9(a,b,c,d,e,f,m,Integer(-1)+n,x)/(m+n)-a**m*c**m*d*cos(e+f*x)**(Integer(1)+Integer(2)*m)*(c+d*sin(e+f*x))**(Integer(-1)-m+n)/(f*(m+n))
                                            else:
                                                if gt(m,0):
                                                    if eq(n,Integer(-1)/Integer(2)):
                                                        return Integer(2)*a*intsin9(a,b,c,d,e,f,Integer(-1)+m,Integer(-1)/Integer(2),x)-Integer(2)*a**(Integer(-1)+m)*b*c**(Integer(-1)+m)*cos(e+f*x)**(Integer(-1)+Integer(2)*m)*(c+d*sin(e+f*x))**(Integer(1)/Integer(2)-m)/(f*(Integer(-1)+Integer(2)*m))
                                                    else:
                                                        return -b*(Integer(-1)+Integer(2)*m)*intsin9(a,b,c,d,e,f,Integer(-1)+m,Integer(1)+n,x)/(d*(Integer(1)+Integer(2)*n))-Integer(2)*a**(Integer(-1)+m)*b*c**(Integer(-1)+m)*cos(e+f*x)**(Integer(-1)+Integer(2)*m)*(c+d*sin(e+f*x))**(Integer(1)-m+n)/(f*(Integer(1)+Integer(2)*n))
                                                else:
                                                    return (Integer(1)+m+n)*intsin9(a,b,c,d,e,f,Integer(1)+m,n,x)/(a*(Integer(1)+Integer(2)*m))+a**(Integer(-1)+m)*b*c**m*cos(e+f*x)**(Integer(1)+Integer(2)*m)*(c+d*sin(e+f*x))**(-m+n)/(f*(Integer(1)+Integer(2)*m))
                                        else:
                                            return a**m*c**m*AppellF1(Integer(1)-m+n,Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m),Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m),Integer(2)-m+n,(c+d*Sin(e+f*x))/(c+d),(c+d*Sin(e+f*x))/(c-d))*Cos(e+f*x)**(Integer(-1)+Integer(2)*m)*(c+d*Sin(e+f*x))**(Integer(1)-m+n)*(Integer(1)+(-c-d*Sin(e+f*x))/(c-d))**(Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m))*(Integer(1)+(-c-d*Sin(e+f*x))/(c+d))**(Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m))/(d*f*(Integer(1)-m+n))
                                    else:
                                        if integer(n):
                                            return intsin9(c,d,a,b,e,f,n,m,x)
                                        else:
                                            if eq(Integer(1)+m+n,0) and not eq(m,Integer(-1)/Integer(2)):
                                                return b*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n/(a*f*(Integer(1)+Integer(2)*m))
                                            else:
                                                if eq(Integer(1)+m+n,0) and not eq(n,Integer(-1)/Integer(2)):
                                                    return intsin9(c,d,a,b,e,f,n,m,x)
                                                else:
                                                    if integer(Integer(2)*m):
                                                        if integer(Integer(2)*n):
                                                            if eq(m,Integer(1)/Integer(2)) and not eq(n,Integer(-1)/Integer(2)):
                                                                if eq(n,Integer(1)/Integer(2)) and False:
                                                                    return intsin9(c,d,a,b,e,f,n,m,x)
                                                                else:
                                                                    return -Integer(2)*b*cos(e+f*x)*(c+d*sin(e+f*x))**n/(f*(Integer(1)+Integer(2)*n)*sqrt(a+b*sin(e+f*x)))
                                                            else:
                                                                if eq(n,Integer(1)/Integer(2)) and not eq(m,Integer(-1)/Integer(2)):
                                                                    return intsin9(c,d,a,b,e,f,n,m,x)
                                                                else:
                                                                    if gt(m,1) and lt(n,-1) and not integer(m+n) and gt(Integer(-1)-m-n,0) and le(Integer(-1)-m-n,Integer(-1)/Integer(2)+m):
                                                                        return -b*(Integer(-1)+Integer(2)*m)*intsin9(a,b,c,d,e,f,Integer(-1)+m,Integer(1)+n,x)/(d*(Integer(1)+Integer(2)*n))-Integer(2)*b*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(-1)+m)*(c+d*sin(e+f*x))**n/(f*(Integer(1)+Integer(2)*n))
                                                                    else:
                                                                        if gt(n,1) and lt(m,-1) and not integer(m+n) and gt(Integer(-1)-m-n,0) and le(Integer(-1)-m-n,Integer(-1)/Integer(2)+n):
                                                                            return intsin9(c,d,a,b,e,f,n,m,x)
                                                                        else:
                                                                            if integer(Integer(-1)/Integer(2)+m) and gt(m,1) and not eq(m+n,0) and not integer(Integer(-1)/Integer(2)+n) and gt(n,1) and lt(-m+n,0) and not integer(m+n) and gt(Integer(-1)-m-n,0) and le(Integer(-1)-m-n,Integer(-1)/Integer(2)+m):
                                                                                return a*(Integer(-1)+Integer(2)*m)*intsin9(a,b,c,d,e,f,Integer(-1)+m,n,x)/(m+n)-b*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(-1)+m)*(c+d*sin(e+f*x))**n/(f*(m+n))
                                                                            else:
                                                                                if integer(Integer(-1)/Integer(2)+n) and gt(n,1) and not eq(m+n,0) and not integer(Integer(-1)/Integer(2)+m) and gt(m,1) and lt(m-n,0) and not integer(m+n) and gt(Integer(-1)-m-n,0) and le(Integer(-1)-m-n,Integer(-1)/Integer(2)+n):
                                                                                    return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                else:
                                                                                    if lt(m,-1) and not gt(-m+n,0) and lt(n,-1):
                                                                                        return (Integer(1)+m+n)*intsin9(a,b,c,d,e,f,Integer(1)+m,n,x)/(a*(Integer(1)+Integer(2)*m))+b*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n/(a*f*(Integer(1)+Integer(2)*m))
                                                                                    else:
                                                                                        if lt(n,-1) and not gt(m-n,0) and lt(m,-1):
                                                                                            return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                        else:
                                                                                            if integer(normal(Integer(1)+m+n)) and lt(normal(Integer(1)+m+n),0) and not SumSimplerQ(n,Integer(1)) and not eq(m,Integer(-1)/Integer(2)):
                                                                                                return (Integer(1)+m+n)*intsin9(a,b,c,d,e,f,Integer(1)+m,n,x)/(a*(Integer(1)+Integer(2)*m))+b*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n/(a*f*(Integer(1)+Integer(2)*m))
                                                                                            else:
                                                                                                if integer(normal(Integer(1)+m+n)) and lt(normal(Integer(1)+m+n),0) and not SumSimplerQ(m,Integer(1)) and not eq(n,Integer(-1)/Integer(2)):
                                                                                                    return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                                else:
                                                                                                    return a**(Integer(1)/Integer(2)+m)*c**(Integer(1)/Integer(2)+m)*AppellF1(Integer(1)-m+n,Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m),Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m),Integer(2)-m+n,(c+d*Sin(e+f*x))/(c+d),(c+d*Sin(e+f*x))/(c-d))*cos(e+f*x)*Cos(e+f*x)**(Integer(-1)+Integer(2)*m)*(c+d*Sin(e+f*x))**(Integer(1)-m+n)*(Integer(1)+(-c-d*Sin(e+f*x))/(c-d))**(Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m))*(Integer(1)+(-c-d*Sin(e+f*x))/(c+d))**(Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m))/(d*f*(Integer(1)-m+n)*sqrt(a+b*sin(e+f*x))*sqrt(c+d*sin(e+f*x)))
                                                        else:
                                                            if eq(m,Integer(1)/Integer(2)):
                                                                return -Integer(2)*b*cos(e+f*x)*(c+d*sin(e+f*x))**n/(f*(Integer(1)+Integer(2)*n)*sqrt(a+b*sin(e+f*x)))
                                                            else:
                                                                return AppellF1(Integer(1)+m-n,Integer(1)/Integer(2)*(Integer(1)-Integer(2)*n),Integer(1)/Integer(2)*(Integer(1)-Integer(2)*n),Integer(2)+m-n,(a+b*Sin(e+f*x))/(a+b),(a+b*Sin(e+f*x))/(a-b))*Cos(e+f*x)**(Integer(-1)+Integer(2)*n)*(a+b*sin(e+f*x))**n*(c+d*sin(e+f*x))**n*(a+b*Sin(e+f*x))**(Integer(1)+m-n)*(Integer(1)+(-a-b*Sin(e+f*x))/(a-b))**(Integer(1)/Integer(2)*(Integer(1)-Integer(2)*n))*(Integer(1)+(-a-b*Sin(e+f*x))/(a+b))**(Integer(1)/Integer(2)*(Integer(1)-Integer(2)*n))/(b*f*(Integer(1)+m-n)*cos(e+f*x)**(Integer(2)*n))
                                                    else:
                                                        if integer(Integer(2)*n):
                                                            return intsin9(c,d,a,b,e,f,n,m,x)
                                                        else:
                                                            if True:
                                                                return AppellF1(Integer(1)-m+n,Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m),Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m),Integer(2)-m+n,(c+d*Sin(e+f*x))/(c+d),(c+d*Sin(e+f*x))/(c-d))*Cos(e+f*x)**(Integer(-1)+Integer(2)*m)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**m*(c+d*Sin(e+f*x))**(Integer(1)-m+n)*(Integer(1)+(-c-d*Sin(e+f*x))/(c-d))**(Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m))*(Integer(1)+(-c-d*Sin(e+f*x))/(c+d))**(Integer(1)/Integer(2)*(Integer(1)-Integer(2)*m))/(d*f*(Integer(1)-m+n)*cos(e+f*x)**(Integer(2)*m))
                                                            else:
                                                                return intsin9(c,d,a,b,e,f,n,m,x)
                                else:
                                    if eq(m,2) and eq(n,-1):
                                        return -b**2*cos(e+f*x)/(d*f)+intsin9(a**2*d,-b*(b*c-Integer(2)*a*d),c,d,e,f,Integer(1),Integer(-1),x)/d
                                    else:
                                        if eq(m,-1) and eq(n,2):
                                            return intsin9(c,d,a,b,e,f,n,m,x)
                                        else:
                                            if eq(m,-1) and eq(n,-1):
                                                return b*intsin6(a,b,e,f,Integer(-1),x)/(b*c-a*d)-d*intsin6(c,d,e,f,Integer(-1),x)/(b*c-a*d)
                                            else:
                                                if eq(n,1):
                                                    if eq(a,0):
                                                        return c*intsin5(e,f,b,m,x)+d*intsin5(e,f,b,Integer(1)+m,x)/b
                                                    else:
                                                        if eq(a**2-b**2,0):
                                                            if eq(a*d*m+b*c*(Integer(1)+m),0):
                                                                return -d*cos(e+f*x)*(a+b*sin(e+f*x))**m/(f*(Integer(1)+m))
                                                            else:
                                                                if lt(m,Integer(-1)/Integer(2)):
                                                                    return (a*d*m+b*c*(Integer(1)+m))*intsin6(a,b,e,f,Integer(1)+m,x)/(a*b*(Integer(1)+Integer(2)*m))+(b*c-a*d)*cos(e+f*x)*(a+b*sin(e+f*x))**m/(a*f*(Integer(1)+Integer(2)*m))
                                                                else:
                                                                    return (a*d*m+b*c*(Integer(1)+m))*intsin6(a,b,e,f,m,x)/(b*(Integer(1)+m))-d*cos(e+f*x)*(a+b*sin(e+f*x))**m/(f*(Integer(1)+m))
                                                        else:
                                                            if gt(m,0):
                                                                return intsin9(a,b,b*d*m+a*c*(Integer(1)+m),a*d*m+b*c*(Integer(1)+m),e,f,Integer(-1)+m,Integer(1),x)/(Integer(1)+m)-d*cos(e+f*x)*(a+b*sin(e+f*x))**m/(f*(Integer(1)+m))
                                                            else:
                                                                if lt(m,-1):
                                                                    return intsin9(a,b,(a*c-b*d)*(Integer(1)+m),-(b*c-a*d)*(Integer(2)+m),e,f,Integer(1)+m,Integer(1),x)/((a**2-b**2)*(Integer(1)+m))-(b*c-a*d)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/((a**2-b**2)*f*(Integer(1)+m))
                                                                else:
                                                                    if eq(c**2-d**2,0) and not integer(Integer(2)*m):
                                                                        return -Integer(2)*c*appellf1(Integer(1)/Integer(2),Integer(-1)/Integer(2),-m,Integer(3)/Integer(2),Integer(1)/Integer(2)*(c-d*sin(e+f*x))/c,b*(c-d*sin(e+f*x))/(b*c+a*d))*(a+b*sin(e+f*x))**m*(c-d*sin(e+f*x))*sqrt(Integer(2))*sqrt((c+d*sin(e+f*x))/c)/(d*f*cos(e+f*x)*(c*(a+b*sin(e+f*x))/(a*c+b*d))**m)
                                                                    else:
                                                                        return (b*c-a*d)*intsin6(a,b,e,f,m,x)/b+d*intsin6(a,b,e,f,Integer(1)+m,x)/b
                                                else:
                                                    if eq(m,1):
                                                        return intsin9(c,d,a,b,e,f,n,m,x)
                                                    else:
                                                        if eq(n,2):
                                                            if eq(a,0):
                                                                return Integer(2)*c*d*intsin6(Integer(0),b,e,f,Integer(1)+m,x)/b+intsin12(Integer(0),b,Integer(1),Integer(0),e,f,c**2,Integer(0),d**2,m,Integer(0),x)
                                                            else:
                                                                if eq(a**2-b**2,0):
                                                                    if eq(c,0):
                                                                        if lt(m,Integer(-1)/Integer(2)):
                                                                            return -d**2*intsin9(a,b,a*m,-b*(Integer(1)+Integer(2)*m),e,f,Integer(1)+m,Integer(1),x)/(a**2*(Integer(1)+Integer(2)*m))+b*d**2*cos(e+f*x)*(a+b*sin(e+f*x))**m/(a*f*(Integer(1)+Integer(2)*m))
                                                                        else:
                                                                            return d**2*intsin9(a,b,b*(Integer(1)+m),-a,e,f,m,Integer(1),x)/(b*(Integer(2)+m))-d**2*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b*f*(Integer(2)+m))
                                                                    else:
                                                                        if lt(m,-1):
                                                                            return intsin9(a,b,a*c*d*(Integer(-1)+m)+b*(d**2+c**2*(Integer(1)+m)),d*(a*d*(Integer(-1)+m)+b*c*(Integer(2)+m)),e,f,Integer(1)+m,Integer(1),x)/(a*b*(Integer(1)+Integer(2)*m))+(b*c-a*d)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))/(a*f*(Integer(1)+Integer(2)*m))
                                                                        else:
                                                                            return intsin9(a,b,b*(d**2*(Integer(1)+m)+c**2*(Integer(2)+m)),-d*(a*d-Integer(2)*b*c*(Integer(2)+m)),e,f,m,Integer(1),x)/(b*(Integer(2)+m))-d**2*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b*f*(Integer(2)+m))
                                                                else:
                                                                    if integer(Integer(2)*m) or not eq(c**2-d**2,0):
                                                                        if lt(m,-1):
                                                                            return -intsin9(a,b,b*(Integer(2)*b*c*d-a*(c**2+d**2))*(Integer(1)+m),a**2*d**2-Integer(2)*a*b*c*d*(Integer(2)+m)+b**2*(d**2*(Integer(1)+m)+c**2*(Integer(2)+m)),e,f,Integer(1)+m,Integer(1),x)/(b*(a**2-b**2)*(Integer(1)+m))-(b**2*c**2-Integer(2)*a*b*c*d+a**2*d**2)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b*(a**2-b**2)*f*(Integer(1)+m))
                                                                        else:
                                                                            return intsin9(a,b,b*(d**2*(Integer(1)+m)+c**2*(Integer(2)+m)),-d*(a*d-Integer(2)*b*c*(Integer(2)+m)),e,f,m,Integer(1),x)/(b*(Integer(2)+m))-d**2*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b*f*(Integer(2)+m))
                                                                    else:
                                                                        if eq(c**2-d**2,0):
                                                                            return -2**(Integer(1)/Integer(2)+n)*appellf1(Integer(1)/Integer(2),Integer(1)/Integer(2)-n,-m,Integer(3)/Integer(2),Integer(1)/Integer(2)*(c-d*sin(e+f*x))/c,b*(c-d*sin(e+f*x))/(b*c+a*d))*sec(e+f*x)*(a+b*sin(e+f*x))**m*(c-d*sin(e+f*x))*(c+d*sin(e+f*x))**n*((c+d*sin(e+f*x))/c)**(Integer(1)/Integer(2)-n)/(d*f*(c*(a+b*sin(e+f*x))/(a*c+b*d))**m)
                                                                        else:
                                                                            return integrate((a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n,x)
                                                        else:
                                                            if eq(m,2):
                                                                return intsin9(c,d,a,b,e,f,n,m,x)
                                                            else:
                                                                if eq(a**2-b**2,0):
                                                                    if gt(m,1) and integer(Integer(2)*m) and integer(Integer(2)*n) or integer(Integer(1)/Integer(2)+m) or integer(m) and eq(c,0):
                                                                        if lt(n,-1):
                                                                            return b**2*intsin(a,b,c,d,e,f,a*c*(Integer(-2)+m)-b*d*(Integer(-4)+m-Integer(2)*n),-b*c*(Integer(-1)+m)+a*d*(Integer(1)+m+Integer(2)*n),Integer(-2)+m,Integer(1)+n,x)/(d*(b*c+a*d)*(Integer(1)+n))-b**2*(b*c-a*d)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(-2)+m)*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*(b*c+a*d)*f*(Integer(1)+n))
                                                                        else:
                                                                            return intsin(a,b,c,d,e,f,a*b*c*(Integer(-2)+m)+b**2*d*(Integer(1)+n)+a**2*d*(m+n),-b*(b*c*(Integer(-1)+m)-a*d*(Integer(-2)+Integer(3)*m+Integer(2)*n)),Integer(-2)+m,n,x)/(d*(m+n))-b**2*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(-2)+m)*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*f*(m+n))
                                                                    else:
                                                                        if lt(m,-1) and integer(Integer(2)*m) and integer(Integer(2)*n) or integer(m) and eq(c,0):
                                                                            if gt(n,0):
                                                                                if lt(n,1):
                                                                                    return -intsin(a,b,c,d,e,f,-b*c*(Integer(1)+m)+a*d*n,-b*d*(Integer(1)+m+n),Integer(1)+m,Integer(-1)+n,x)/(a*b*(Integer(1)+Integer(2)*m))+b*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n/(a*f*(Integer(1)+Integer(2)*m))
                                                                                else:
                                                                                    return intsin(a,b,c,d,e,f,b*(c**2*(Integer(1)+m)+d**2*(Integer(-1)+n))+a*c*d*(Integer(1)+m-n),d*(a*d*(Integer(1)+m-n)+b*c*(m+n)),Integer(1)+m,Integer(-2)+n,x)/(a*b*(Integer(1)+Integer(2)*m))+(b*c-a*d)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(-1)+n)/(a*f*(Integer(1)+Integer(2)*m))
                                                                            else:
                                                                                return intsin(a,b,c,d,e,f,b*c*(Integer(1)+m)-a*d*(Integer(2)+Integer(2)*m+n),b*d*(Integer(2)+m+n),Integer(1)+m,n,x)/(a*(b*c-a*d)*(Integer(1)+Integer(2)*m))+b**2*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(1)+n)/(a*(b*c-a*d)*f*(Integer(1)+Integer(2)*m))
                                                                        else:
                                                                            if eq(m,-1) and integer(Integer(2)*n) or eq(c,0):
                                                                                if gt(n,1):
                                                                                    return -d*intsin9(b*d*(Integer(-1)+n)-a*c*n,b*c*(Integer(-1)+n)-a*d*n,c,d,e,f,Integer(1),Integer(-2)+n,x)/(a*b)-(b*c-a*d)*cos(e+f*x)*(c+d*sin(e+f*x))**(Integer(-1)+n)/(a*f*(a+b*sin(e+f*x)))
                                                                                else:
                                                                                    if lt(n,0):
                                                                                        return d*intsin9(a*n,-b*(Integer(1)+n),c,d,e,f,Integer(1),n,x)/(a*(b*c-a*d))-b**2*cos(e+f*x)*(c+d*sin(e+f*x))**(Integer(1)+n)/(a*(b*c-a*d)*f*(a+b*sin(e+f*x)))
                                                                                    else:
                                                                                        return d*n*intsin9(a,-b,c,d,e,f,Integer(1),Integer(-1)+n,x)/(a*b)-b*cos(e+f*x)*(c+d*sin(e+f*x))**n/(a*f*(a+b*sin(e+f*x)))
                                                                            else:
                                                                                if eq(m,Integer(1)/Integer(2)):
                                                                                    if integer(Integer(2)*n):
                                                                                        if gt(n,0):
                                                                                            return Integer(2)*(b*c+a*d)*n*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(-1)+n,x)/(b*(Integer(1)+Integer(2)*n))-Integer(2)*b*cos(e+f*x)*(c+d*sin(e+f*x))**n/(f*(Integer(1)+Integer(2)*n)*sqrt(a+b*sin(e+f*x)))
                                                                                        else:
                                                                                            if eq(n,-1):
                                                                                                return -Integer(2)*b*subst(integrate(1/(b*c+a*d-d*x**2),x),x,b*cos(e+f*x)/sqrt(a+b*sin(e+f*x)))/f
                                                                                            else:
                                                                                                if eq(n,Integer(-1)/Integer(2)):
                                                                                                    if eq(c,0) and eq(-a/b+d,0):
                                                                                                        return -Integer(2)*subst(integrate(1/sqrt(Integer(1)-x**2/a),x),x,b*cos(e+f*x)/sqrt(a+b*sin(e+f*x)))/f
                                                                                                    else:
                                                                                                        return -Integer(2)*b*subst(integrate(1/(b+d*x**2),x),x,b*cos(e+f*x)/(sqrt(a+b*sin(e+f*x))*sqrt(c+d*sin(e+f*x))))/f
                                                                                                else:
                                                                                                    if eq(n,Integer(-3)/Integer(2)):
                                                                                                        return -Integer(2)*b**2*cos(e+f*x)/((b*c+a*d)*f*sqrt(a+b*sin(e+f*x))*sqrt(c+d*sin(e+f*x)))
                                                                                                    else:
                                                                                                        return Integer(1)/Integer(2)*(b*c-a*d)*(Integer(3)+Integer(2)*n)*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(1)+n,x)/(b*(c**2-d**2)*(Integer(1)+n))+(b*c-a*d)*cos(e+f*x)*(c+d*sin(e+f*x))**(Integer(1)+n)/((c**2-d**2)*f*(Integer(1)+n)*sqrt(a+b*sin(e+f*x)))
                                                                                    else:
                                                                                        return -Integer(2)*b*cos(e+f*x)*hyp2f1(Integer(1)/Integer(2),-n,Integer(3)/Integer(2),d*(a-b*sin(e+f*x))/(b*c+a*d))*(c+d*sin(e+f*x))**n/(f*(a*(c+d*sin(e+f*x))/(a*c+b*d))**n*sqrt(a+b*sin(e+f*x)))
                                                                                else:
                                                                                    if eq(m,Integer(-1)/Integer(2)) and integer(Integer(2)*n):
                                                                                        if gt(n,0):
                                                                                            if eq(n,Integer(1)/Integer(2)):
                                                                                                return (b*c-a*d)*intsin9(a,b,c,d,e,f,Integer(-1)/Integer(2),Integer(-1)/Integer(2),x)/b+d*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(-1)/Integer(2),x)/b
                                                                                            else:
                                                                                                return -intsin(a,b,c,d,e,f,a*c*d-b*(Integer(2)*d**2*(Integer(-1)+n)+c**2*(Integer(-1)+Integer(2)*n)),d*(a*d-b*c*(Integer(-3)+Integer(4)*n)),Integer(-1)/Integer(2),Integer(-2)+n,x)/(b*(Integer(-1)+Integer(2)*n))-Integer(2)*d*cos(e+f*x)*(c+d*sin(e+f*x))**(Integer(-1)+n)/(f*(Integer(-1)+Integer(2)*n)*sqrt(a+b*sin(e+f*x)))
                                                                                        else:
                                                                                            if eq(n,-1):
                                                                                                return b*intsin6(a,b,e,f,Integer(-1)/Integer(2),x)/(b*c-a*d)-d*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(-1),x)/(b*c-a*d)
                                                                                            else:
                                                                                                if eq(n,Integer(-1)/Integer(2)):
                                                                                                    if eq(c,0) and eq(-a/b+d,0) and gt(a,0):
                                                                                                        return -sqrt(Integer(2))*subst(integrate(1/sqrt(Integer(1)-x**2),x),x,b*cos(e+f*x)/(a+b*sin(e+f*x)))/(f*sqrt(a))
                                                                                                    else:
                                                                                                        return -Integer(2)*a*subst(integrate(1/(Integer(2)*b**2-(a*c-b*d)*x**2),x),x,b*cos(e+f*x)/(sqrt(a+b*sin(e+f*x))*sqrt(c+d*sin(e+f*x))))/f
                                                                                                else:
                                                                                                    return -Integer(1)/Integer(2)*intsin(a,b,c,d,e,f,a*d-Integer(2)*b*c*(Integer(1)+n),b*d*(Integer(3)+Integer(2)*n),Integer(-1)/Integer(2),Integer(1)+n,x)/(b*(c**2-d**2)*(Integer(1)+n))-d*cos(e+f*x)*(c+d*sin(e+f*x))**(Integer(1)+n)/((c**2-d**2)*f*(Integer(1)+n)*sqrt(a+b*sin(e+f*x)))
                                                                                    else:
                                                                                        if gt(n,1) and integer(n):
                                                                                            return intsin(a,b,c,d,e,f,d*(a*c*m+b*d*(Integer(-1)+n))+b*c**2*(m+n),d*(a*d*m+b*c*(Integer(-1)+m+Integer(2)*n)),m,Integer(-2)+n,x)/(b*(m+n))-d*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(-1)+n)/(f*(m+n))
                                                                                        else:
                                                                                            if eq(Integer(1)+m+n,0):
                                                                                                if gt(c**2-d**2,0) and gt(b*c/a,0) or lt(c**2-d**2,0) and gt(b*d/a,0):
                                                                                                    return -2**(Integer(1)+m)*b*((a*c+b*d)/a)**(Integer(-1)-m)*cos(e+f*x)*hyp2f1(Integer(1)/Integer(2),Integer(1)+m,Integer(3)/Integer(2),-(a*c-b*d)*(a-b*sin(e+f*x))/((a*c+b*d)*(a+b*sin(e+f*x))))*(a+b*sin(e+f*x))**m*((a+b*sin(e+f*x))/a)**(Integer(-1)-m)/(a*f)
                                                                                                else:
                                                                                                    if gt(c**2-d**2,0) and lt(b*c/a,0) or lt(c**2-d**2,0) and lt(b*d/a,0):
                                                                                                        return -intsin9(a,b,-c,-d,e,f,m,n,x)*(-c-d*sin(e+f*x))**(Integer(-1)-n)*(c+d*sin(e+f*x))**(Integer(1)+n)
                                                                                                    else:
                                                                                                        return -2**m*a**2*cos(e+f*x)*hyp2f1(Integer(1)/Integer(2),Integer(1)/Integer(2)-m,Integer(3)/Integer(2),Integer(1)/Integer(2)*(b*c-a*d)*(a-b*sin(e+f*x))/(a*b*(c+d*sin(e+f*x))))*(a+b*sin(e+f*x))**(Integer(-1)+m)*((b*c+a*d)*(a+b*sin(e+f*x))/(a*b*(c+d*sin(e+f*x))))**(Integer(1)/Integer(2)-m)*sqrt(Integer(2))/((b*c+a*d)*f*(c+d*sin(e+f*x))**m)
                                                                                            else:
                                                                                                if integer(m):
                                                                                                    return -2**m*a**m*appellf1(Integer(1)/Integer(2),Integer(1)/Integer(2)-m,-n,Integer(3)/Integer(2),Integer(1)/Integer(2)*(a-b*sin(e+f*x))/a,d*(a-b*sin(e+f*x))/(b*c+a*d))*(a-b*sin(e+f*x))*(c+d*sin(e+f*x))**n*sqrt(Integer(2))*sqrt((a+b*sin(e+f*x))/a)/(b*f*cos(e+f*x)*(a*(c+d*sin(e+f*x))/(a*c+b*d))**n)
                                                                                                else:
                                                                                                    return -2**(Integer(1)/Integer(2)+m)*appellf1(Integer(1)/Integer(2),Integer(1)/Integer(2)-m,-n,Integer(3)/Integer(2),Integer(1)/Integer(2)*(a-b*sin(e+f*x))/a,d*(a-b*sin(e+f*x))/(b*c+a*d))*(a-b*sin(e+f*x))*(a+b*sin(e+f*x))**m*((a+b*sin(e+f*x))/a)**(Integer(1)/Integer(2)-m)*(c+d*sin(e+f*x))**n/(b*f*cos(e+f*x)*(a*(c+d*sin(e+f*x))/(a*c+b*d))**n)
                                                                else:
                                                                    if eq(c**2-d**2,0):
                                                                        return intsin9(c,d,a,b,e,f,n,m,x)
                                                                    else:
                                                                        if integer(m) or integer(Integer(2)*m) and integer(Integer(2)*n):
                                                                            if gt(m,2):
                                                                                if lt(n,-1):
                                                                                    return intsin12(a,b,c,d,e,f,b*(b*c-a*d)**2*(Integer(-2)+m)+a*d*((a**2+b**2)*c-Integer(2)*a*b*d)*(Integer(1)+n),b*(a*b*c**2+(a**2+b**2)*c*d-Integer(3)*a*b*d**2)*(Integer(1)+n)-a*(b*c-a*d)**2*(Integer(2)+n),b*(b**2*(c**2-d**2)-(b*c-a*d)**2*m+d*(Integer(2)*a*b*c-(a**2+b**2)*d)*n),Integer(-3)+m,Integer(1)+n,x)/(d*(c**2-d**2)*(Integer(1)+n))-(b**2*c**2-Integer(2)*a*b*c*d+a**2*d**2)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(-2)+m)*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*(c**2-d**2)*f*(Integer(1)+n))
                                                                                else:
                                                                                    return intsin12(a,b,c,d,e,f,a**3*d*(m+n)+b**2*(b*c*(Integer(-2)+m)+a*d*(Integer(1)+n)),-b*(a*b*c-b**2*d*(Integer(-1)+m+n)-Integer(3)*a**2*d*(m+n)),-b**2*(b*c*(Integer(-1)+m)-a*d*(Integer(-2)+Integer(3)*m+Integer(2)*n)),Integer(-3)+m,n,x)/(d*(m+n))-b**2*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(-2)+m)*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*f*(m+n))
                                                                            else:
                                                                                if gt(n,2):
                                                                                    return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                else:
                                                                                    if lt(m,-1):
                                                                                        if gt(n,0):
                                                                                            if lt(n,1):
                                                                                                if eq(m,Integer(-3)/Integer(2)) and eq(n,Integer(1)/Integer(2)):
                                                                                                    if eq(c,0):
                                                                                                        return -d**2*intsin9(a,b,Integer(0),d,e,f,Integer(1)/Integer(2),Integer(-3)/Integer(2),x)/(a**2-b**2)-Integer(2)*a*d*cos(e+f*x)/((a**2-b**2)*f*sqrt(d*sin(e+f*x))*sqrt(a+b*sin(e+f*x)))
                                                                                                    else:
                                                                                                        return (c-d)*intsin9(a,b,c,d,e,f,Integer(-1)/Integer(2),Integer(-1)/Integer(2),x)/(a-b)-(b*c-a*d)*intsin(a,b,c,d,e,f,Integer(1),Integer(1),Integer(-3)/Integer(2),Integer(-1)/Integer(2),x)/(a-b)
                                                                                                else:
                                                                                                    return intsin12(a,b,c,d,e,f,a*c*(Integer(1)+m)+b*d*n,a*d*(Integer(1)+m)-b*c*(Integer(2)+m),-b*d*(Integer(2)+m+n),Integer(1)+m,Integer(-1)+n,x)/((a**2-b**2)*(Integer(1)+m))-b*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)*(c+d*sin(e+f*x))**n/((a**2-b**2)*f*(Integer(1)+m))
                                                                                            else:
                                                                                                if eq(m,Integer(-3)/Integer(2)) and eq(n,Integer(3)/Integer(2)):
                                                                                                    if eq(c,0):
                                                                                                        return -a*d*intsin9(a,b,Integer(0),d,e,f,Integer(-3)/Integer(2),Integer(1)/Integer(2),x)/b+d*intsin9(a,b,Integer(0),d,e,f,Integer(-1)/Integer(2),Integer(1)/Integer(2),x)/b
                                                                                                    else:
                                                                                                        return d**2*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(-1)/Integer(2),x)/b**2+(b*c-a*d)*intsin(a,b,c,d,e,f,b*c+a*d,Integer(2)*b*d,Integer(-3)/Integer(2),Integer(-1)/Integer(2),x)/b**2
                                                                                                else:
                                                                                                    return intsin12(a,b,c,d,e,f,c*(a*c-b*d)*(Integer(1)+m)+d*(b*c-a*d)*(Integer(-1)+n),d*(a*c-b*d)*(Integer(1)+m)-c*(b*c-a*d)*(Integer(2)+m),-d*(b*c-a*d)*(Integer(1)+m+n),Integer(1)+m,Integer(-2)+n,x)/((a**2-b**2)*(Integer(1)+m))-(b*c-a*d)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)*(c+d*sin(e+f*x))**(Integer(-1)+n)/((a**2-b**2)*f*(Integer(1)+m))
                                                                                        else:
                                                                                            if eq(m,Integer(-3)/Integer(2)) and eq(n,Integer(-1)/Integer(2)):
                                                                                                if eq(c,0):
                                                                                                    return d*intsin(a,b,Integer(0),d,e,f,b,a,Integer(-1)/Integer(2),Integer(-3)/Integer(2),x)/(a**2-b**2)+Integer(2)*b*cos(e+f*x)/((a**2-b**2)*f*sqrt(d*sin(e+f*x))*sqrt(a+b*sin(e+f*x)))
                                                                                                else:
                                                                                                    return intsin9(a,b,c,d,e,f,Integer(-1)/Integer(2),Integer(-1)/Integer(2),x)/(a-b)-b*intsin(a,b,c,d,e,f,Integer(1),Integer(1),Integer(-3)/Integer(2),Integer(-1)/Integer(2),x)/(a-b)
                                                                                            else:
                                                                                                if eq(a,0) and integer(m) and not integer(n) or not integer(Integer(2)*n) and lt(n,-1) and integer(n) and not integer(m) or eq(a,0):
                                                                                                    return intsin12(a,b,c,d,e,f,a*(b*c-a*d)*(Integer(1)+m)+b**2*d*(Integer(2)+m+n),-b**2*c-b*(b*c-a*d)*(Integer(1)+m),-b**2*d*(Integer(3)+m+n),Integer(1)+m,n,x)/((a**2-b**2)*(b*c-a*d)*(Integer(1)+m))-b**2*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)*(c+d*sin(e+f*x))**(Integer(1)+n)/((a**2-b**2)*(b*c-a*d)*f*(Integer(1)+m))
                                                                                                else:
                                                                                                    return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                    else:
                                                                                        if lt(n,-1):
                                                                                            return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                        else:
                                                                                            if eq(m,-1) and eq(n,Integer(1)/Integer(2)):
                                                                                                return d*intsin6(c,d,e,f,Integer(-1)/Integer(2),x)/b+(b*c-a*d)*intsin9(a,b,c,d,e,f,Integer(-1),Integer(-1)/Integer(2),x)/b
                                                                                            else:
                                                                                                if eq(m,Integer(1)/Integer(2)) and eq(n,-1):
                                                                                                    return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                                else:
                                                                                                    if eq(m,Integer(3)/Integer(2)) and eq(n,-1):
                                                                                                        return b*intsin6(a,b,e,f,Integer(1)/Integer(2),x)/d-(b*c-a*d)*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(-1),x)/d
                                                                                                    else:
                                                                                                        if eq(m,-1) and eq(n,Integer(3)/Integer(2)):
                                                                                                            return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                                        else:
                                                                                                            if eq(m,-1) and eq(n,Integer(-1)/Integer(2)):
                                                                                                                if gt(c+d,0):
                                                                                                                    return Integer(2)*ellippi(Integer(2)*b/(a+b),-Integer(1)/Integer(4)*Pi+Integer(1)/Integer(2)*(e+f*x),Integer(2)*d/(c+d))/((a+b)*f*sqrt(c+d))
                                                                                                                else:
                                                                                                                    if gt(c-d,0):
                                                                                                                        return Integer(2)*ellippi(-Integer(2)*b/(a-b),Integer(1)/Integer(4)*Pi+Integer(1)/Integer(2)*(e+f*x),-Integer(2)*d/(c-d))/((a-b)*f*sqrt(c-d))
                                                                                                                    else:
                                                                                                                        return intsin9(a,b,c/(c+d),d/(c+d),e,f,Integer(-1),Integer(-1)/Integer(2),x)*sqrt((c+d*sin(e+f*x))/(c+d))/sqrt(c+d*sin(e+f*x))
                                                                                                            else:
                                                                                                                if eq(m,Integer(-1)/Integer(2)) and eq(n,-1):
                                                                                                                    return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                                                else:
                                                                                                                    if eq(m,Integer(1)/Integer(2)) and eq(n,Integer(-1)/Integer(2)):
                                                                                                                        if eq(a,0):
                                                                                                                            if pos((c+d)/b):
                                                                                                                                if gt(c**2-d**2,0) and gt(c**2,0):
                                                                                                                                    return Integer(2)*b**(Integer(1)/Integer(2))*c*(c+d)**(Integer(1)/Integer(2))*ellippi((c+d)/d,asin(b**(Integer(1)/Integer(2))*sqrt(c+d*sin(e+f*x))/((c+d)**(Integer(1)/Integer(2))*sqrt(b*sin(e+f*x)))),(-c-d)/(c-d))*sqrt(Integer(1)-csc(e+f*x))*sqrt(Integer(1)+csc(e+f*x))*tan(e+f*x)/(d*f*sqrt(c**2-d**2))
                                                                                                                                else:
                                                                                                                                    return Integer(2)*b**(Integer(1)/Integer(2))*(c+d)**(Integer(1)/Integer(2))*ellippi((c+d)/d,asin(b**(Integer(1)/Integer(2))*sqrt(c+d*sin(e+f*x))/((c+d)**(Integer(1)/Integer(2))*sqrt(b*sin(e+f*x)))),(-c-d)/(c-d))*sqrt(c*(Integer(1)-csc(e+f*x))/(c+d))*sqrt(c*(Integer(1)+csc(e+f*x))/(c-d))*tan(e+f*x)/(d*f)
                                                                                                                            else:
                                                                                                                                return intsin9(Integer(0),-b,c,d,e,f,Integer(1)/Integer(2),Integer(-1)/Integer(2),x)*sqrt(b*sin(e+f*x))/sqrt(-b*sin(e+f*x))
                                                                                                                        else:
                                                                                                                            if eq(c,0):
                                                                                                                                return a*intsin9(a,b,Integer(0),d,e,f,Integer(-1)/Integer(2),Integer(-1)/Integer(2),x)+b*intsin9(a,b,Integer(0),d,e,f,Integer(-1)/Integer(2),Integer(1)/Integer(2),x)/d
                                                                                                                            else:
                                                                                                                                if pos((a+b)/(c+d)):
                                                                                                                                    return Integer(2)*(c+d)**(Integer(1)/Integer(2))*ellippi(b*(c+d)/((a+b)*d),asin((a+b)**(Integer(1)/Integer(2))*sqrt(c+d*sin(e+f*x))/((c+d)**(Integer(1)/Integer(2))*sqrt(a+b*sin(e+f*x)))),(a-b)*(c+d)/((a+b)*(c-d)))*(a+b*sin(e+f*x))*sqrt(-(b*c-a*d)*(Integer(1)-sin(e+f*x))/((c+d)*(a+b*sin(e+f*x))))*sqrt((b*c-a*d)*(Integer(1)+sin(e+f*x))/((c-d)*(a+b*sin(e+f*x))))/((a+b)**(Integer(1)/Integer(2))*d*f*cos(e+f*x))
                                                                                                                                else:
                                                                                                                                    return intsin9(a,b,-c,-d,e,f,Integer(1)/Integer(2),Integer(-1)/Integer(2),x)*sqrt(-c-d*sin(e+f*x))/sqrt(c+d*sin(e+f*x))
                                                                                                                    else:
                                                                                                                        if eq(m,Integer(-1)/Integer(2)) and eq(n,Integer(1)/Integer(2)):
                                                                                                                            return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                                                        else:
                                                                                                                            if eq(m,Integer(-1)/Integer(2)) and eq(n,Integer(-1)/Integer(2)):
                                                                                                                                if eq(c,0):
                                                                                                                                    if lt(a**2-b**2,0) and gt(b**2,0):
                                                                                                                                        if eq(Integer(-1)+d**2,0) and gt(b*d,0):
                                                                                                                                            return -Integer(2)*d*ellipf(asin(cos(e+f*x)/(Integer(1)+d*sin(e+f*x))),(-a+b*d)/(a+b*d))/(f*sqrt(a+b*d))
                                                                                                                                        else:
                                                                                                                                            return intsin9(a,b,Integer(0),Sign(b),e,f,Integer(-1)/Integer(2),Integer(-1)/Integer(2),x)*sqrt(Sign(b)*sin(e+f*x))/sqrt(d*sin(e+f*x))
                                                                                                                                    else:
                                                                                                                                        if pos((a+b)/d):
                                                                                                                                            if gt(a**2-b**2,0) and gt(a**2,0):
                                                                                                                                                return -Integer(2)*(a+b)**(Integer(1)/Integer(2))*ellipf(asin(d**(Integer(1)/Integer(2))*sqrt(a+b*sin(e+f*x))/((a+b)**(Integer(1)/Integer(2))*sqrt(d*sin(e+f*x)))),(-a-b)/(a-b))*sqrt(a**2)*sqrt(-cot(e+f*x)**2)/(a*d**(Integer(1)/Integer(2))*f*cot(e+f*x)*sqrt(a**2-b**2))
                                                                                                                                            else:
                                                                                                                                                return -Integer(2)*(a+b)**(Integer(1)/Integer(2))*ellipf(asin(d**(Integer(1)/Integer(2))*sqrt(a+b*sin(e+f*x))/((a+b)**(Integer(1)/Integer(2))*sqrt(d*sin(e+f*x)))),(-a-b)/(a-b))*sqrt(a*(Integer(1)-csc(e+f*x))/(a+b))*sqrt(a*(Integer(1)+csc(e+f*x))/(a-b))*tan(e+f*x)/(a*d**(Integer(1)/Integer(2))*f)
                                                                                                                                        else:
                                                                                                                                            return intsin9(a,b,Integer(0),-d,e,f,Integer(-1)/Integer(2),Integer(-1)/Integer(2),x)*sqrt(-d*sin(e+f*x))/sqrt(d*sin(e+f*x))
                                                                                                                                else:
                                                                                                                                    if eq(a,0):
                                                                                                                                        return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                                                                    else:
                                                                                                                                        if pos((c+d)/(a+b)):
                                                                                                                                            if True:
                                                                                                                                                return Integer(2)*(a+b)**(Integer(1)/Integer(2))*ellipf(asin((c+d)**(Integer(1)/Integer(2))*sqrt(a+b*sin(e+f*x))/((a+b)**(Integer(1)/Integer(2))*sqrt(c+d*sin(e+f*x)))),(a+b)*(c-d)/((a-b)*(c+d)))*(c+d*sin(e+f*x))*sqrt((b*c-a*d)*(Integer(1)-sin(e+f*x))/((a+b)*(c+d*sin(e+f*x))))*sqrt(-(b*c-a*d)*(Integer(1)+sin(e+f*x))/((a-b)*(c+d*sin(e+f*x))))/((c+d)**(Integer(1)/Integer(2))*(b*c-a*d)*f*cos(e+f*x))
                                                                                                                                            else:
                                                                                                                                                return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                                                                        else:
                                                                                                                                            return intsin9(-a,-b,c,d,e,f,Integer(-1)/Integer(2),Integer(-1)/Integer(2),x)*sqrt(-a-b*sin(e+f*x))/sqrt(a+b*sin(e+f*x))
                                                                                                                            else:
                                                                                                                                if eq(m,Integer(-1)/Integer(2)) and eq(n,Integer(3)/Integer(2)) and eq(c,0):
                                                                                                                                    return -Integer(1)/Integer(2)*a*d*intsin9(a,b,Integer(0),d,e,f,Integer(-1)/Integer(2),Integer(1)/Integer(2),x)/b+Integer(1)/Integer(2)*d*intsin(a,b,Integer(0),d,e,f,a,Integer(2)*b,Integer(-1)/Integer(2),Integer(1)/Integer(2),x)/b
                                                                                                                                else:
                                                                                                                                    if eq(n,Integer(-1)/Integer(2)) and eq(m,Integer(3)/Integer(2)) and eq(a,0):
                                                                                                                                        return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                                                                    else:
                                                                                                                                        if gt(m,0) and lt(m,2) and gt(n,-1) and lt(n,2):
                                                                                                                                            if eq(c,0) and eq(m,Integer(1)/Integer(2)) and eq(n,Integer(1)/Integer(2)):
                                                                                                                                                return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                                                                            else:
                                                                                                                                                return intsin12(a,b,c,d,e,f,a**2*c*d*(m+n)+b*d*(b*c*(Integer(-1)+m)+a*d*n),a*d*(Integer(2)*b*c+a*d)*(m+n)-b*d*(a*c-b*d*(Integer(-1)+m+n)),b*d*(b*c*n+a*d*(Integer(-1)+Integer(2)*m+n)),Integer(-2)+m,Integer(-1)+n,x)/(d*(m+n))-b*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(-1)+m)*(c+d*sin(e+f*x))**n/(f*(m+n))
                                                                                                                                        else:
                                                                                                                                            if gt(n,0) and lt(n,2) and gt(m,-1) and lt(m,2):
                                                                                                                                                return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                                                                            else:
                                                                                                                                                if integer(m) and gt(m,0):
                                                                                                                                                    return -(b*c-a*d)*intsin9(a,b,c,d,e,f,Integer(-1)+m,n,x)/d+b*intsin9(a,b,c,d,e,f,Integer(-1)+m,Integer(1)+n,x)/d
                                                                                                                                                else:
                                                                                                                                                    if integer(n) and gt(n,0):
                                                                                                                                                        return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                                                                                    else:
                                                                                                                                                        return integrate((a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n,x)
                                                                        else:
                                                                            if integer(m) and gt(m,0):
                                                                                return -(b*c-a*d)*intsin9(a,b,c,d,e,f,Integer(-1)+m,n,x)/d+b*intsin9(a,b,c,d,e,f,Integer(-1)+m,Integer(1)+n,x)/d
                                                                            else:
                                                                                if integer(n) and gt(n,0):
                                                                                    return intsin9(c,d,a,b,e,f,n,m,x)
                                                                                else:
                                                                                    return integrate((a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n,x)


def intsin(a,b,c,d,e,f,A,B,m,n,x):
    if eq(f,0):
        return x*(a+b*sin(e))**m*(A+B*sin(e))*(c+d*sin(e))**n
    else:
        if eq(m,0) or eq(b,0):
            return intsin9(c,d,A,B,e,f,n,Integer(1),x)*(a+b*sin(e+f*x))**m
        else:
            if eq(n,0) or eq(d,0):
                return intsin9(a,b,A,B,e,f,m,Integer(1),x)*(c+d*sin(e+f*x))**n
            else:
                if eq(B,0):
                    if eq(A,0):
                        return Integer(0)
                    else:
                        return A*intsin9(a,b,c,d,e,f,m,n,x)
                else:
                    if eq(A*b-a*B,0):
                        return B*intsin9(a,b,c,d,e,f,Integer(1)+m,n,x)/b
                    else:
                        if eq(B*c-A*d,0):
                            return B*intsin9(a,b,c,d,e,f,m,Integer(1)+n,x)/d
                        else:
                            if eq(b*c-a*d,0):
                                if integer(m) or gt(b/d,0):
                                    return (b/d)**m*intsin9(c,d,A,B,e,f,m+n,Integer(1),x)
                                else:
                                    if integer(n) or gt(d/b,0):
                                        return (d/b)**n*intsin9(a,b,A,B,e,f,m+n,Integer(1),x)
                                    else:
                                        return intsin9(c,d,A,B,e,f,m+n,Integer(1),x)*(a+b*sin(e+f*x))**m/(c+d*sin(e+f*x))**m
                            else:
                                if eq(n,1):
                                    if eq(A,0):
                                        if eq(c,0):
                                            return intsin9(a,b,Integer(0),B*d,e,f,m,Integer(2),x)
                                        else:
                                            return intsin(a,b,Integer(0),B,e,f,c,d,m,Integer(1),x)
                                    else:
                                        if eq(a**2-b**2,0):
                                            if lt(m,Integer(-1)/Integer(2)):
                                                return intsin9(a,b,B*(b*c-a*d)*m+A*(b*d*m+a*c*(Integer(1)+m)),b*B*d*(Integer(1)+Integer(2)*m),e,f,Integer(1)+m,Integer(1),x)/(b**2*(Integer(1)+Integer(2)*m))+(A*b-a*B)*(b*c-a*d)*cos(e+f*x)*(a+b*sin(e+f*x))**m/(a*b*f*(Integer(1)+Integer(2)*m))
                                            else:
                                                return intsin9(a,b,b*(B*d*(Integer(1)+m)+A*c*(Integer(2)+m)),A*b*d*(Integer(2)+m)-B*(a*d-b*c*(Integer(2)+m)),e,f,m,Integer(1),x)/(b*(Integer(2)+m))-B*d*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b*f*(Integer(2)+m))
                                        else:
                                            if lt(m,-1):
                                                return intsin9(a,b,b*(-B*(b*c-a*d)+A*(a*c-b*d))*(Integer(1)+m),-A*b*(b*c-a*d)*(Integer(2)+m)+B*(a*b*c*(Integer(2)+m)-d*(a**2+b**2*(Integer(1)+m))),e,f,Integer(1)+m,Integer(1),x)/(b*(a**2-b**2)*(Integer(1)+m))-(A*b-a*B)*(b*c-a*d)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b*(a**2-b**2)*f*(Integer(1)+m))
                                            else:
                                                return intsin9(a,b,b*(B*d*(Integer(1)+m)+A*c*(Integer(2)+m)),A*b*d*(Integer(2)+m)-B*(a*d-b*c*(Integer(2)+m)),e,f,m,Integer(1),x)/(b*(Integer(2)+m))-B*d*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b*f*(Integer(2)+m))
                                else:
                                    if eq(m,1):
                                        return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                    else:
                                        if eq(a**2-b**2,0):
                                            if eq(c**2-d**2,0):
                                                if integer(m):
                                                    if eq(a*B*(m-n)+A*b*(Integer(1)+m+n),0):
                                                        return -a**m*B*c**m*cos(e+f*x)**(Integer(1)+Integer(2)*m)*(c+d*sin(e+f*x))**(-m+n)/(f*(Integer(1)+m+n))
                                                    else:
                                                        if integer(Integer(1)/Integer(2)+n):
                                                            if gt(m,0):
                                                                if lt(Integer(1)/Integer(2)+n,0):
                                                                    return -(B*c*(m-n)-A*d*(Integer(1)+m+n))*intsin9(a,b,c,d,e,f,m,Integer(1)+n,x)/(c*d*(Integer(1)+Integer(2)*n))-a**m*c**(Integer(-1)+m)*(B*c-A*d)*cos(e+f*x)**(Integer(1)+Integer(2)*m)*(c+d*sin(e+f*x))**(-m+n)/(f*(Integer(1)+Integer(2)*n))
                                                                else:
                                                                    return -(B*c*(m-n)-A*d*(Integer(1)+m+n))*intsin9(a,b,c,d,e,f,m,n,x)/(d*(Integer(1)+m+n))-a**m*B*c**m*cos(e+f*x)**(Integer(1)+Integer(2)*m)*(c+d*sin(e+f*x))**(-m+n)/(f*(Integer(1)+m+n))
                                                            else:
                                                                return (a*B*(m-n)+A*b*(Integer(1)+m+n))*intsin9(a,b,c,d,e,f,Integer(1)+m,n,x)/(a*b*(Integer(1)+Integer(2)*m))+a**(Integer(-1)+m)*(A*b-a*B)*c**m*cos(e+f*x)**(Integer(1)+Integer(2)*m)*(c+d*sin(e+f*x))**(-m+n)/(f*(Integer(1)+Integer(2)*m))
                                                        else:
                                                            return a**m*c**m*Integrate(Cos(e+f*x)**(Integer(2)*m)*(A+B*Sin(e+f*x))*(c+d*Sin(e+f*x))**(-m+n),x)
                                                else:
                                                    if integer(n):
                                                        return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                    else:
                                                        if eq(a*B*(m-n)+A*b*(Integer(1)+m+n),0):
                                                            if eq(m,Integer(-1)/Integer(2)) and eq(n,Integer(-1)/Integer(2)):
                                                                return Integer(1)/Integer(2)*(B*c+A*d)*intsin9(a,b,c,d,e,f,Integer(-1)/Integer(2),Integer(1)/Integer(2),x)/(c*d)+Integer(1)/Integer(2)*(A*b+a*B)*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(-1)/Integer(2),x)/(a*b)
                                                            else:
                                                                return -B*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n/(f*(Integer(1)+m+n))
                                                        else:
                                                            if eq(a*B*(-m+n)+A*b*(Integer(1)+m+n),0):
                                                                return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                            else:
                                                                if eq(m,Integer(1)/Integer(2)):
                                                                    return -(B*c-A*d)*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),n,x)/d+B*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(1)+n,x)/d
                                                                else:
                                                                    if eq(n,Integer(1)/Integer(2)):
                                                                        return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                                    else:
                                                                        if lt(m,Integer(-1)/Integer(2)):
                                                                            return (a*B*(m-n)+A*b*(Integer(1)+m+n))*intsin9(a,b,c,d,e,f,Integer(1)+m,n,x)/(a*b*(Integer(1)+Integer(2)*m))+(A*b-a*B)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n/(a*f*(Integer(1)+Integer(2)*m))
                                                                        else:
                                                                            if lt(n,Integer(-1)/Integer(2)):
                                                                                return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                                            else:
                                                                                if integer(m+n) and lt(m+n,0) and not SumSimplerQ(n,Integer(1)) and not eq(m,Integer(-1)/Integer(2)):
                                                                                    return (a*B*(m-n)+A*b*(Integer(1)+m+n))*intsin9(a,b,c,d,e,f,Integer(1)+m,n,x)/(a*b*(Integer(1)+Integer(2)*m))+(A*b-a*B)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n/(a*f*(Integer(1)+Integer(2)*m))
                                                                                else:
                                                                                    if integer(m+n) and lt(m+n,0) and not SumSimplerQ(m,Integer(1)) and not eq(n,Integer(-1)/Integer(2)):
                                                                                        return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                                                    else:
                                                                                        return -(B*c*(m-n)-A*d*(Integer(1)+m+n))*intsin9(a,b,c,d,e,f,m,n,x)/(d*(Integer(1)+m+n))-B*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n/(f*(Integer(1)+m+n))
                                            else:
                                                if eq(Integer(2)+m+n,0) and eq(A*(a*d*m+b*c*(Integer(1)+n))-B*(a*c*m+b*d*(Integer(1)+n)),0):
                                                    return (B*c-A*d)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(1)+n)/((c**2-d**2)*f*(Integer(1)+n))
                                                else:
                                                    if gt(m,0) and integer(Integer(2)*m):
                                                        if eq(m,Integer(1)/Integer(2)):
                                                            if eq(A*b*d*(Integer(3)+Integer(2)*n)-B*(b*c-Integer(2)*a*d*(Integer(1)+n)),0):
                                                                return -Integer(2)*b*B*cos(e+f*x)*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*f*(Integer(3)+Integer(2)*n)*sqrt(a+b*sin(e+f*x)))
                                                            else:
                                                                if lt(n,-1):
                                                                    return Integer(1)/Integer(2)*(A*b*d*(Integer(3)+Integer(2)*n)-B*(b*c-Integer(2)*a*d*(Integer(1)+n)))*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(1)+n,x)/(d*(b*c+a*d)*(Integer(1)+n))-b**2*(B*c-A*d)*cos(e+f*x)*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*(b*c+a*d)*f*(Integer(1)+n)*sqrt(a+b*sin(e+f*x)))
                                                                else:
                                                                    return (A*b*d*(Integer(3)+Integer(2)*n)-B*(b*c-Integer(2)*a*d*(Integer(1)+n)))*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),n,x)/(b*d*(Integer(3)+Integer(2)*n))-Integer(2)*b*B*cos(e+f*x)*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*f*(Integer(3)+Integer(2)*n)*sqrt(a+b*sin(e+f*x)))
                                                        else:
                                                            if lt(n,-1):
                                                                return -b*intsin(a,b,c,d,e,f,a*A*d*(Integer(-2)+m-n)-B*(a*c*(Integer(-1)+m)+b*d*(Integer(1)+n)),-A*b*d*(Integer(1)+m+n)+B*(b*c*m-a*d*(Integer(1)+n)),Integer(-1)+m,Integer(1)+n,x)/(d*(b*c+a*d)*(Integer(1)+n))-b**2*(B*c-A*d)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(-1)+m)*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*(b*c+a*d)*f*(Integer(1)+n))
                                                            else:
                                                                return intsin(a,b,c,d,e,f,a*A*d*(Integer(1)+m+n)+B*(a*c*(Integer(-1)+m)+b*d*(Integer(1)+n)),A*b*d*(Integer(1)+m+n)-B*(b*c*m-a*d*(Integer(2)*m+n)),Integer(-1)+m,n,x)/(d*(Integer(1)+m+n))-b*B*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(-1)+m)*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*f*(Integer(1)+m+n))
                                                    else:
                                                        if lt(m,Integer(-1)/Integer(2)) and integer(Integer(2)*m):
                                                            if gt(n,0):
                                                                return -intsin(a,b,c,d,e,f,A*(-b*c*(Integer(1)+m)+a*d*n)-B*(a*c*m+b*d*n),-d*(a*B*(m-n)+A*b*(Integer(1)+m+n)),Integer(1)+m,Integer(-1)+n,x)/(a*b*(Integer(1)+Integer(2)*m))+(A*b-a*B)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n/(a*f*(Integer(1)+Integer(2)*m))
                                                            else:
                                                                return intsin(a,b,c,d,e,f,B*(a*c*m+b*d*(Integer(1)+n))+A*(b*c*(Integer(1)+m)-a*d*(Integer(2)+Integer(2)*m+n)),(A*b-a*B)*d*(Integer(2)+m+n),Integer(1)+m,n,x)/(a*(b*c-a*d)*(Integer(1)+Integer(2)*m))+b*(A*b-a*B)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(1)+n)/(a*(b*c-a*d)*f*(Integer(1)+Integer(2)*m))
                                                        else:
                                                            if gt(n,0) and integer(n) or eq(Integer(1)/Integer(2)+m,0):
                                                                return intsin(a,b,c,d,e,f,A*b*c*(Integer(1)+m+n)+B*(a*c*m+b*d*n),A*b*d*(Integer(1)+m+n)+B*(a*d*m+b*c*n),m,Integer(-1)+n,x)/(b*(Integer(1)+m+n))-B*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n/(f*(Integer(1)+m+n))
                                                            else:
                                                                if le(n,-1) and integer(n) or eq(Integer(1)/Integer(2)+m,0):
                                                                    if eq(n,-1):
                                                                        if eq(m,Integer(-1)/Integer(2)):
                                                                            return (A*b-a*B)*intsin6(a,b,e,f,Integer(-1)/Integer(2),x)/(b*c-a*d)+(B*c-A*d)*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(-1),x)/(b*c-a*d)
                                                                        else:
                                                                            return B*intsin6(a,b,e,f,m,x)/d-(B*c-A*d)*intsin9(a,b,c,d,e,f,m,Integer(-1),x)/d
                                                                    else:
                                                                        return intsin(a,b,c,d,e,f,A*(a*d*m+b*c*(Integer(1)+n))-B*(a*c*m+b*d*(Integer(1)+n)),b*(B*c-A*d)*(Integer(2)+m+n),m,Integer(1)+n,x)/(b*(c**2-d**2)*(Integer(1)+n))+(B*c-A*d)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(1)+n)/((c**2-d**2)*f*(Integer(1)+n))
                                                                else:
                                                                    if eq(m,Integer(-1)/Integer(2)) and eq(n,Integer(-1)/Integer(2)):
                                                                        return (A*b-a*B)*intsin9(a,b,c,d,e,f,Integer(-1)/Integer(2),Integer(-1)/Integer(2),x)/b+B*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(-1)/Integer(2),x)/b
                                                                    else:
                                                                        return (A*b-a*B)*intsin9(a,b,c,d,e,f,m,n,x)/b+B*intsin9(a,b,c,d,e,f,Integer(1)+m,n,x)/b
                                        else:
                                            if eq(c**2-d**2,0):
                                                return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                            else:
                                                if gt(m,1):
                                                    if lt(n,-1):
                                                        return intsin12(a,b,c,d,e,f,b*(b*c-a*d)*(B*c-A*d)*(Integer(-1)+m)+a*d*(a*A*c+b*B*c-(A*b+a*B)*d)*(Integer(1)+n),b*(b*d*(B*c-A*d)+a*(A*c*d+B*(c**2-Integer(2)*d**2)))*(Integer(1)+n)-a*(b*c-a*d)*(B*c-A*d)*(Integer(2)+n),b*(d*(A*b*c+a*B*c-a*A*d)*(Integer(1)+m+n)-b*B*(c**2*m+d**2*(Integer(1)+n))),Integer(-2)+m,Integer(1)+n,x)/(d*(c**2-d**2)*(Integer(1)+n))-(b*c-a*d)*(B*c-A*d)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(-1)+m)*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*(c**2-d**2)*f*(Integer(1)+n))
                                                    else:
                                                        if integer(n) and gt(n,1) and not integer(m) or eq(a,0) and not eq(c,0):
                                                            return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                        else:
                                                            return intsin12(a,b,c,d,e,f,a**2*A*d*(Integer(1)+m+n)+b*B*(b*c*(Integer(-1)+m)+a*d*(Integer(1)+n)),a*(Integer(2)*A*b+a*B)*d*(Integer(1)+m+n)-b*B*(a*c-b*d*(m+n)),b*(A*b*d*(Integer(1)+m+n)-B*(b*c*m-a*d*(Integer(2)*m+n))),Integer(-2)+m,n,x)/(d*(Integer(1)+m+n))-b*B*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(-1)+m)*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*f*(Integer(1)+m+n))
                                                else:
                                                    if gt(n,1):
                                                        return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                    else:
                                                        if lt(m,-1):
                                                            if eq(m,Integer(-3)/Integer(2)) and eq(n,Integer(-1)/Integer(2)):
                                                                if eq(c,0):
                                                                    return d*intsin(a,b,Integer(0),d,e,f,A*b-a*B,a*A-b*B,Integer(-1)/Integer(2),Integer(-3)/Integer(2),x)/(a**2-b**2)+Integer(2)*(A*b-a*B)*cos(e+f*x)/((a**2-b**2)*f*sqrt(d*sin(e+f*x))*sqrt(a+b*sin(e+f*x)))
                                                                else:
                                                                    if eq(A-B,0):
                                                                        if eq(a,0):
                                                                            if pos((c+d)/b):
                                                                                return -Integer(2)*A*(c-d)*(c+d)**(Integer(1)/Integer(2))*ellipe(asin(b**(Integer(1)/Integer(2))*sqrt(c+d*sin(e+f*x))/((c+d)**(Integer(1)/Integer(2))*sqrt(b*sin(e+f*x)))),(-c-d)/(c-d))*sqrt(c*(Integer(1)-csc(e+f*x))/(c+d))*sqrt(c*(Integer(1)+csc(e+f*x))/(c-d))*tan(e+f*x)/(b**(Integer(3)/Integer(2))*c**2*f)
                                                                            else:
                                                                                return -intsin(Integer(0),-b,c,d,e,f,A,B,Integer(-3)/Integer(2),Integer(-1)/Integer(2),x)*sqrt(-b*sin(e+f*x))/sqrt(b*sin(e+f*x))
                                                                        else:
                                                                            if pos((a+b)/(c+d)):
                                                                                return -Integer(2)*A*(c-d)*(c+d)**(Integer(1)/Integer(2))*ellipe(asin((a+b)**(Integer(1)/Integer(2))*sqrt(c+d*sin(e+f*x))/((c+d)**(Integer(1)/Integer(2))*sqrt(a+b*sin(e+f*x)))),(a-b)*(c+d)/((a+b)*(c-d)))*(a+b*sin(e+f*x))*sqrt(-(b*c-a*d)*(Integer(1)-sin(e+f*x))/((c+d)*(a+b*sin(e+f*x))))*sqrt((b*c-a*d)*(Integer(1)+sin(e+f*x))/((c-d)*(a+b*sin(e+f*x))))/((a+b)**(Integer(1)/Integer(2))*(b*c-a*d)**2*f*cos(e+f*x))
                                                                            else:
                                                                                return intsin(a,b,-c,-d,e,f,A,B,Integer(-3)/Integer(2),Integer(-1)/Integer(2),x)*sqrt(-c-d*sin(e+f*x))/sqrt(c+d*sin(e+f*x))
                                                                    else:
                                                                        return (A-B)*intsin9(a,b,c,d,e,f,Integer(-1)/Integer(2),Integer(-1)/Integer(2),x)/(a-b)-(A*b-a*B)*intsin(a,b,c,d,e,f,Integer(1),Integer(1),Integer(-3)/Integer(2),Integer(-1)/Integer(2),x)/(a-b)
                                                            else:
                                                                if gt(n,0):
                                                                    return intsin12(a,b,c,d,e,f,(a*A-b*B)*c*(Integer(1)+m)+(A*b-a*B)*d*n,(a*A-b*B)*d*(Integer(1)+m)-(A*b-a*B)*c*(Integer(2)+m),-(A*b-a*B)*d*(Integer(2)+m+n),Integer(1)+m,Integer(-1)+n,x)/((a**2-b**2)*(Integer(1)+m))+(-A*b+a*B)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)*(c+d*sin(e+f*x))**n/((a**2-b**2)*f*(Integer(1)+m))
                                                                else:
                                                                    if eq(a,0) and integer(m) and not integer(n) or not integer(Integer(2)*n) and lt(n,-1) and integer(n) and not integer(m) or eq(a,0):
                                                                        return intsin12(a,b,c,d,e,f,(a*A-b*B)*(b*c-a*d)*(Integer(1)+m)+b*(A*b-a*B)*d*(Integer(2)+m+n),(A*b-a*B)*(a*d*(Integer(1)+m)-b*c*(Integer(2)+m)),-b*(A*b-a*B)*d*(Integer(3)+m+n),Integer(1)+m,n,x)/((a**2-b**2)*(b*c-a*d)*(Integer(1)+m))-b*(A*b-a*B)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)*(c+d*sin(e+f*x))**(Integer(1)+n)/((a**2-b**2)*(b*c-a*d)*f*(Integer(1)+m))
                                                                    else:
                                                                        return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                        else:
                                                            if lt(n,-1):
                                                                return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                            else:
                                                                if gt(m,0) and lt(m,1) and gt(n,-1):
                                                                    return intsin12(a,b,c,d,e,f,a*A*c*(Integer(1)+m+n)+B*(b*c*m+a*d*n),B*(a*c+b*d)*(m+n)+A*(b*c+a*d)*(Integer(1)+m+n),A*b*d*(Integer(1)+m+n)+B*(a*d*m+b*c*n),Integer(-1)+m,Integer(-1)+n,x)/(Integer(1)+m+n)-B*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n/(f*(Integer(1)+m+n))
                                                                else:
                                                                    if gt(n,0) and lt(n,1) and gt(m,-1):
                                                                        return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                                    else:
                                                                        if eq(m,-1) and eq(n,-1):
                                                                            return (A*b-a*B)*intsin6(a,b,e,f,Integer(-1),x)/(b*c-a*d)+(B*c-A*d)*intsin6(c,d,e,f,Integer(-1),x)/(b*c-a*d)
                                                                        else:
                                                                            if eq(m,Integer(-1)/Integer(2)) and eq(n,Integer(-1)/Integer(2)):
                                                                                if eq(c,0) and gt(b,0) and gt(-a**2+b**2,0) and eq(A-B,0):
                                                                                    if eq(d,1):
                                                                                        return Integer(4)*A*ellippi(Integer(-1),-asin(cos(e+f*x)/(Integer(1)+sin(e+f*x))),(-a+b)/(a+b))/(f*sqrt(a+b))
                                                                                    else:
                                                                                        return intsin(a,b,Integer(0),Integer(1),e,f,A,B,Integer(-1)/Integer(2),Integer(-1)/Integer(2),x)*sqrt(sin(e+f*x))/sqrt(d*sin(e+f*x))
                                                                                else:
                                                                                    if eq(c,0) or True:
                                                                                        return -(B*c-A*d)*intsin9(a,b,c,d,e,f,Integer(-1)/Integer(2),Integer(-1)/Integer(2),x)/d+B*intsin9(a,b,c,d,e,f,Integer(-1)/Integer(2),Integer(1)/Integer(2),x)/d
                                                                                    else:
                                                                                        return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                                            else:
                                                                                if eq(m,Integer(1)/Integer(2)) and eq(n,Integer(-3)/Integer(2)):
                                                                                    return -(B*c-A*d)*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(-3)/Integer(2),x)/d+B*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(-1)/Integer(2),x)/d
                                                                                else:
                                                                                    if eq(n,Integer(1)/Integer(2)) and eq(m,Integer(-3)/Integer(2)):
                                                                                        return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                                                    else:
                                                                                        if eq(n,-1):
                                                                                            return B*intsin6(a,b,e,f,m,x)/d-(B*c-A*d)*intsin9(a,b,c,d,e,f,m,Integer(-1),x)/d
                                                                                        else:
                                                                                            if eq(m,-1):
                                                                                                return intsin(c,d,a,b,e,f,A,B,n,m,x)
                                                                                            else:
                                                                                                return integrate((a+b*sin(e+f*x))**m*(A+B*sin(e+f*x))*(c+d*sin(e+f*x))**n,x)


def intsin12(a,b,c,d,e,f,A,B,C,m,n,x):
    if eq(f,0):
        return x*(a+b*sin(e))**m*(c+d*sin(e))**n*(A+B*sin(e)+C*sin(e)**2)
    else:
        if eq(C,0):
            return intsin(a,b,c,d,e,f,A,B,m,n,x)
        else:
            if eq(n,0):
                if eq(m,0):
                    return C*intsin5(e,f,Integer(1),Integer(2),x)+intsin6(A,B,e,f,Integer(1),x)
                else:
                    if eq(b,0):
                        return a**m*(C*intsin5(e,f,Integer(1),Integer(2),x)+intsin6(A,B,e,f,Integer(1),x))
                    else:
                        if eq(A*b**2-a*b*B+a**2*C,0):
                            return intsin9(a,b,b*B-a*C,b*C,e,f,Integer(1)+m,Integer(1),x)/b**2
                        else:
                            if lt(m,-1):
                                if eq(a**2-b**2,0):
                                    return intsin9(a,b,(b*B-a*C)*m+a*A*(Integer(1)+m),b*C*(Integer(1)+Integer(2)*m),e,f,Integer(1)+m,Integer(1),x)/(a**2*(Integer(1)+Integer(2)*m))+(A*b-a*B+b*C)*cos(e+f*x)*(a+b*sin(e+f*x))**m/(a*f*(Integer(1)+Integer(2)*m))
                                else:
                                    return intsin9(a,b,b*(a*A-b*B+a*C)*(Integer(1)+m),-A*b**2+a*b*B-a**2*C-b*(A*b-a*B+b*C)*(Integer(1)+m),e,f,Integer(1)+m,Integer(1),x)/(b*(a**2-b**2)*(Integer(1)+m))-(A*b**2-a*b*B+a**2*C)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b*(a**2-b**2)*f*(Integer(1)+m))
                            else:
                                return intsin9(a,b,b*C*(Integer(1)+m)+A*b*(Integer(2)+m),-a*C+b*B*(Integer(2)+m),e,f,m,Integer(1),x)/(b*(Integer(2)+m))-C*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b*f*(Integer(2)+m))
            else:
                if eq(d,0):
                    return c**n*intsin12(a,b,Integer(1),Integer(0),e,f,A,B,C,m,Integer(0),x)
                else:
                    if eq(m,0) or eq(b,0):
                        return intsin12(c,d,a,b,e,f,A,B,C,n,m,x)
                    else:
                        if eq(b*c-a*d,0):
                            if integer(m) or gt(b/d,0):
                                return (b/d)**m*intsin12(c,d,Integer(1),Integer(0),e,f,A,B,C,m+n,Integer(0),x)
                            else:
                                if integer(n) or gt(d/b,0):
                                    return (d/b)**n*intsin12(a,b,Integer(1),Integer(0),e,f,A,B,C,m+n,Integer(0),x)
                                else:
                                    return intsin12(c,d,Integer(1),Integer(0),e,f,A,B,C,m+n,Integer(0),x)*(a+b*sin(e+f*x))**m/(c+d*sin(e+f*x))**m
                        else:
                            if eq(n,1) and not eq(a**2-b**2,0):
                                if lt(m,-1):
                                    return -intsin12(a,b,c,d,e,f,b*((b*B-a*C)*(b*c-a*d)-A*b*(a*c-b*d))*(Integer(1)+m),b*B*(a**2*d+b**2*d*(Integer(1)+m)-a*b*c*(Integer(2)+m))+(b*c-a*d)*(A*b**2*(Integer(2)+m)+C*(a**2+b**2*(Integer(1)+m))),-b*(a**2-b**2)*C*d*(Integer(1)+m),Integer(1)+m,Integer(0),x)/(b**2*(a**2-b**2)*(Integer(1)+m))-(A*b**2-a*b*B+a**2*C)*(b*c-a*d)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b**2*(a**2-b**2)*f*(Integer(1)+m))
                                else:
                                    return intsin12(a,b,c,d,e,f,a*C*d+A*b*c*(Integer(3)+m),b*(B*c*(Integer(3)+m)+d*(C*(Integer(2)+m)+A*(Integer(3)+m))),-Integer(2)*a*C*d+b*(c*C+B*d)*(Integer(3)+m),m,Integer(0),x)/(b*(Integer(3)+m))-C*d*cos(e+f*x)*sin(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b*f*(Integer(3)+m))
                            else:
                                if eq(m,1) and not eq(c**2-d**2,0):
                                    return intsin12(c,d,a,b,e,f,A,B,C,n,m,x)
                                else:
                                    if eq(a**2-b**2,0):
                                        if eq(c**2-d**2,0):
                                            if lt(m,Integer(-1)/Integer(2)) or eq(Integer(2)+m+n,0) and not eq(m,Integer(-1)/Integer(2)):
                                                return intsin(a,b,c,d,e,f,B*(b*c*m+a*d*(Integer(1)+n))-C*(a*c*m+b*d*(Integer(1)+n))+A*(a*c*(Integer(1)+m)-b*d*(Integer(2)+Integer(2)*m+n)),C*(b*c*(Integer(1)+Integer(2)*m)-a*d*(Integer(-1)+m-n))+(a*A-b*B)*d*(Integer(2)+m+n),Integer(1)+m,n,x)/(b*(b*c-a*d)*(Integer(1)+Integer(2)*m))+(a*A-b*B+a*C)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(1)+n)/((b*c-a*d)*f*(Integer(1)+Integer(2)*m))
                                            else:
                                                if lt(n,Integer(-1)/Integer(2)) or eq(Integer(2)+m+n,0) and not eq(n,Integer(-1)/Integer(2)):
                                                    return intsin12(c,d,a,b,e,f,A,B,C,n,m,x)
                                                else:
                                                    if eq(n,Integer(-1)/Integer(2)):
                                                        return intsin(a,b,c,d,e,f,A+C,B,m,Integer(-1)/Integer(2),x)-Integer(2)*C*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)/(b*f*(Integer(3)+Integer(2)*m)*sqrt(c+d*sin(e+f*x)))
                                                    else:
                                                        if eq(m,Integer(-1)/Integer(2)):
                                                            return intsin12(c,d,a,b,e,f,A,B,C,n,m,x)
                                                        else:
                                                            return intsin(a,b,c,d,e,f,A*b*d*(Integer(2)+m+n)+C*(a*c*m+b*d*(Integer(1)+n)),C*(a*d*m-b*c*(Integer(1)+m))+b*B*d*(Integer(2)+m+n),m,n,x)/(b*d*(Integer(2)+m+n))-C*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*f*(Integer(2)+m+n))
                                        else:
                                            if lt(m,Integer(-1)/Integer(2)):
                                                return intsin(a,b,c,d,e,f,B*(b*c*m+a*d*(Integer(1)+n))-C*(a*c*m+b*d*(Integer(1)+n))+A*(a*c*(Integer(1)+m)-b*d*(Integer(2)+Integer(2)*m+n)),C*(b*c*(Integer(1)+Integer(2)*m)-a*d*(Integer(-1)+m-n))+(a*A-b*B)*d*(Integer(2)+m+n),Integer(1)+m,n,x)/(b*(b*c-a*d)*(Integer(1)+Integer(2)*m))+(a*A-b*B+a*C)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(1)+n)/((b*c-a*d)*f*(Integer(1)+Integer(2)*m))
                                            else:
                                                if lt(n,-1) or eq(Integer(2)+m+n,0):
                                                    return intsin(a,b,c,d,e,f,A*d*(a*d*m+b*c*(Integer(1)+n))+(c*C-B*d)*(a*c*m+b*d*(Integer(1)+n)),b*(d*(B*c-A*d)*(Integer(2)+m+n)-C*(c**2*(Integer(1)+m)+d**2*(Integer(1)+n))),m,Integer(1)+n,x)/(b*d*(c**2-d**2)*(Integer(1)+n))-(c**2*C-B*c*d+A*d**2)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*(c**2-d**2)*f*(Integer(1)+n))
                                                else:
                                                    return intsin(a,b,c,d,e,f,A*b*d*(Integer(2)+m+n)+C*(a*c*m+b*d*(Integer(1)+n)),C*(a*d*m-b*c*(Integer(1)+m))+b*B*d*(Integer(2)+m+n),m,n,x)/(b*d*(Integer(2)+m+n))-C*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*f*(Integer(2)+m+n))
                                    else:
                                        if eq(c**2-d**2,0):
                                            return intsin12(c,d,a,b,e,f,A,B,C,n,m,x)
                                        else:
                                            if gt(m,0):
                                                if lt(n,-1):
                                                    return intsin12(a,b,c,d,e,f,A*d*(b*d*m+a*c*(Integer(1)+n))+(c*C-B*d)*(b*c*m+a*d*(Integer(1)+n)),C*(b*c*d*(Integer(1)+n)-a*(c**2+d**2*(Integer(1)+n)))-d*(B*(b*d*(Integer(1)+n)-a*c*(Integer(2)+n))+A*(-b*c*(Integer(1)+n)+a*d*(Integer(2)+n))),b*(d*(B*c-A*d)*(Integer(2)+m+n)-C*(c**2*(Integer(1)+m)+d**2*(Integer(1)+n))),Integer(-1)+m,Integer(1)+n,x)/(d*(c**2-d**2)*(Integer(1)+n))-(c**2*C-B*c*d+A*d**2)*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*(c**2-d**2)*f*(Integer(1)+n))
                                                else:
                                                    if integer(n) and gt(n,0) and not integer(m) or eq(a,0) and not eq(c,0):
                                                        return intsin12(c,d,a,b,e,f,A,B,C,n,m,x)
                                                    else:
                                                        return intsin12(a,b,c,d,e,f,a*A*d*(Integer(2)+m+n)+C*(b*c*m+a*d*(Integer(1)+n)),(A*b+a*B)*d*(Integer(2)+m+n)-C*(a*c-b*d*(Integer(1)+m+n)),C*(a*d*m-b*c*(Integer(1)+m))+b*B*d*(Integer(2)+m+n),Integer(-1)+m,n,x)/(d*(Integer(2)+m+n))-C*cos(e+f*x)*(a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**(Integer(1)+n)/(d*f*(Integer(2)+m+n))
                                            else:
                                                if gt(n,0):
                                                    return intsin12(c,d,a,b,e,f,A,B,C,n,m,x)
                                                else:
                                                    if lt(m,-1):
                                                        if eq(m,Integer(-3)/Integer(2)) and eq(n,Integer(-1)/Integer(2)):
                                                            if eq(c,0):
                                                                return C*intsin9(a,b,Integer(0),d,e,f,Integer(-1)/Integer(2),Integer(1)/Integer(2),x)/(b*d)+intsin(a,b,Integer(0),d,e,f,A*b,b*B-a*C,Integer(-3)/Integer(2),Integer(-1)/Integer(2),x)/b
                                                            else:
                                                                return C*intsin9(a,b,c,d,e,f,Integer(1)/Integer(2),Integer(-1)/Integer(2),x)/b**2+intsin(a,b,c,d,e,f,A*b**2-a**2*C,b*(b*B-Integer(2)*a*C),Integer(-3)/Integer(2),Integer(-1)/Integer(2),x)/b**2
                                                        else:
                                                            if eq(a,0) and integer(m) and not integer(n) or not integer(Integer(2)*n) and lt(n,-1) and integer(n) and not integer(m) or eq(a,0):
                                                                return intsin12(a,b,c,d,e,f,(a*A-b*B+a*C)*(b*c-a*d)*(Integer(1)+m)+(A*b**2-a*b*B+a**2*C)*d*(Integer(2)+m+n),-c*(A*b**2-a*b*B+a**2*C)-(A*b-a*B+b*C)*(b*c-a*d)*(Integer(1)+m),-(A*b**2-a*b*B+a**2*C)*d*(Integer(3)+m+n),Integer(1)+m,n,x)/((a**2-b**2)*(b*c-a*d)*(Integer(1)+m))-(A*b**2-a*b*B+a**2*C)*cos(e+f*x)*(a+b*sin(e+f*x))**(Integer(1)+m)*(c+d*sin(e+f*x))**(Integer(1)+n)/((a**2-b**2)*(b*c-a*d)*f*(Integer(1)+m))
                                                            else:
                                                                return intsin12(c,d,a,b,e,f,A,B,C,n,m,x)
                                                    else:
                                                        if lt(n,-1):
                                                            return intsin12(c,d,a,b,e,f,A,B,C,n,m,x)
                                                        else:
                                                            if eq(m,-1) and eq(n,-1):
                                                                return C*x/(b*d)+(A*b**2-a*b*B+a**2*C)*intsin6(a,b,e,f,Integer(-1),x)/(b*(b*c-a*d))-(c**2*C-B*c*d+A*d**2)*intsin6(c,d,e,f,Integer(-1),x)/(d*(b*c-a*d))
                                                            else:
                                                                if eq(m,Integer(-1)/Integer(2)) and eq(n,-1):
                                                                    return C*intsin6(a,b,e,f,Integer(1)/Integer(2),x)/(b*d)-intsin(a,b,c,d,e,f,a*c*C-A*b*d,b*c*C-b*B*d+a*C*d,Integer(-1)/Integer(2),Integer(-1),x)/(b*d)
                                                                else:
                                                                    if eq(n,Integer(-1)/Integer(2)) and eq(m,-1):
                                                                        return intsin12(c,d,a,b,e,f,A,B,C,n,m,x)
                                                                    else:
                                                                        if eq(m,Integer(-1)/Integer(2)) and eq(n,Integer(-1)/Integer(2)):
                                                                            if eq(a,0) or True:
                                                                                return Integer(1)/Integer(2)*intsin12(a,b,c,d,e,f,Integer(2)*a*A*d-C*(b*c-a*d),-Integer(2)*(a*c*C-(A*b+a*B)*d),Integer(2)*b*B*d-C*(b*c+a*d),Integer(-3)/Integer(2),Integer(-1)/Integer(2),x)/d-C*cos(e+f*x)*sqrt(c+d*sin(e+f*x))/(d*f*sqrt(a+b*sin(e+f*x)))
                                                                            else:
                                                                                return intsin12(c,d,a,b,e,f,A,B,C,n,m,x)
                                                                        else:
                                                                            return integrate((a+b*sin(e+f*x))**m*(c+d*sin(e+f*x))**n*(A+B*sin(e+f*x)+C*sin(e+f*x)**2),x)

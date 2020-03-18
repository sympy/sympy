# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 23:55:57 2020

@author: Ayush Shukla
"""

from sympy import *
x,y=symbols('x y')
print('Enter the parameters a,b,g,f,h,c of equation ax^2+by^2+2gx+2fy+2hxy+c=0')
a=int(input())
b=int(input())
g=int(input())
f=int(input())
h=int(input())
c=int(input())
P=Poly((x**2)*a+(y**2)*b+2*g*x+2*h*x*y+2*f*y+c,(x,y))
S=Eq((x**2)*a+(y**2)*b+(2*g*x)+(2*h*x*y)+(2*f*y),-c)
print('Enter the point')
x1=int(input())
y1=int(input())
S1=P(x1,y1)
if S1==0:
    print('Point lies on the curve')
     T=Eq((a*x1+g+h*y1)*x+(b*y1+f+h*x1)*y,-(c+(g*x1)+(f*y1)))
    p1=plot_implicit(S,show=False)
    p2=plot_implicit(T,show=False)
    p1.extend(p2)
    p1.show()
elif S1<0:
    print('Point lies inside the curve')
else:
    print('Point lies outside the curve')

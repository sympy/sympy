from __future__ import division
import warnings

from sympy import Abs, Rational,simplify,expand, Float, S, Symbol,symbols, cos, pi, sqrt, oo
from sympy.functions.elementary.trigonometric import tan
from sympy.geometry import (Circle, Ellipse, GeometryError, Point, Point2D, Polygon, Ray, RegularPolygon, Segment, Triangle, are_similar,convex_hull, intersection, Line)
from sympy.utilities.pytest import raises, slow
from sympy.utilities.randtest import verify_numerically
from sympy.geometry.polygon import rad, deg
from sympy import integrate

def test_polygon_second_moment_of_area():
	x,y=symbols('x,y')
	x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6=symbols('x1,x2,x3,x4,x5,x6,y1,y2,y3,y4,y5,y6')
	p1,p2,p3,p4,p5,p6=Point(x1,y1),Point(x2,y2),Point(x3,y3),Point(x4,y4),Point(x5,y5),Point(x6,y6)
	p=Point(0,0)

	#Triangle (eqation of side of triangle) I_yy
	eq1_y=(((p1.y-p2.y)/(p1.x-p2.x))*(x-p1.x))+p1.y
	eq2_y=(((p3.y-p2.y)/(p3.x-p2.x))*(x-p2.x))+p2.y
	eq3_y=(((p3.y-p1.y)/(p3.x-p1.x))*(x-p3.x))+p3.y

	i1 = integrate(eq1_y*(x**2),(x,x1,x2))
	i2 = integrate(eq2_y*(x**2),(x,x2,x3))
	i3 = integrate(eq3_y*(x**2),(x,x1,x3))

	i_cal=simplify(i3-i1-i2)
	triangle=Polygon(p1,p2,p3)
	i_for=triangle.second_moment_of_area(p)[1]
	i = int(expand(i_cal-i_for))
	assert i is 0
	#I_xx
	eq1_x=(((p1.x-p2.x)/(p1.y-p2.y))*(y-p1.y))+p1.x
	eq2_x=(((p3.x-p2.x)/(p3.y-p2.y))*(y-p2.y))+p2.x
	eq3_x=(((p3.x-p1.x)/(p3.y-p1.y))*(y-p3.y))+p3.x

	i1 = integrate(eq1_x*(y**2),(y,y1,y2))
	i2 = integrate(eq2_x*(y**2),(y,y2,y3))
	i3 = integrate(eq3_x*(y**2),(y,y1,y3))

	i_cal=simplify(i1+i2-i3)  
	triangle=Polygon(p1,p2,p3)
	i_for=triangle.second_moment_of_area(p)[0]
	i = int(expand(i_cal-i_for))
	assert i is 0

	#Quadrilateral (eqation of side of quadrilateral) I_yy
	eq1_y=(((p1.y-p2.y)/(p1.x-p2.x))*(x-p1.x))+p1.y
	eq2_y=(((p2.y-p3.y)/(p2.x-p3.x))*(x-p2.x))+p2.y
	eq3_y=(((p3.y-p4.y)/(p3.x-p4.x))*(x-p3.x))+p3.y
	eq4_y=(((p4.y-p1.y)/(p4.x-p1.x))*(x-p4.x))+p4.y

	i1=integrate(eq1_y*(x**2),(x,x1,x2))
	i2=integrate(eq2_y*(x**2),(x,x2,x3))
	i3=integrate(eq3_y*(x**2),(x,x3,x4))
	i4=integrate(eq4_y*(x**2),(x,x1,x4))

	i_cal=simplify(i4-i1-i2-i3)
	triangle=Polygon(p1,p2,p3,p4)
	i_for=triangle.second_moment_of_area(p)[1]
	i = int(expand(i_cal-i_for))
	assert i is 0
	#I_xx
	eq1_x=(((p1.x-p2.x)/(p1.y-p2.y))*(y-p1.y))+p1.x
	eq2_x=(((p2.x-p3.x)/(p2.y-p3.y))*(y-p2.y))+p2.x
	eq3_x=(((p3.x-p4.x)/(p3.y-p4.y))*(y-p3.y))+p3.x
	eq4_x=(((p4.x-p1.x)/(p4.y-p1.y))*(y-p4.y))+p4.x

	i1=integrate(eq1_x*(y**2),(y,y1,y2))
	i2=integrate(eq2_x*(y**2),(y,y2,y3))
	i3=integrate(eq3_x*(y**2),(y,y3,y4))
	i4=integrate(eq4_x*(y**2),(y,y1,y4))

	i_cal=simplify(-i4+i1+i2+i3)
	triangle=Polygon(p1,p2,p3,p4)
	i_for=triangle.second_moment_of_area(p)[0]
	i = int(expand(i_cal-i_for))
	assert i is 0	

	#Pentagon (eqation of side of pentagon) I_yy
	eq1_y=(((p1.y-p2.y)/(p1.x-p2.x))*(x-p1.x))+p1.y
	eq2_y=(((p2.y-p3.y)/(p2.x-p3.x))*(x-p2.x))+p2.y
	eq3_y=(((p3.y-p4.y)/(p3.x-p4.x))*(x-p3.x))+p3.y
	eq4_y=(((p4.y-p5.y)/(p4.x-p5.x))*(x-p4.x))+p4.y
	eq5_y=(((p5.y-p1.y)/(p5.x-p1.x))*(x-p5.x))+p5.y

	i1=integrate(eq1_y*(x**2),(x,x1,x2))
	i2=integrate(eq2_y*(x**2),(x,x2,x3))
	i3=integrate(eq3_y*(x**2),(x,x3,x4))
	i4=integrate(eq4_y*(x**2),(x,x4,x5))
	i5=integrate(eq5_y*(x**2),(x,x1,x5))

	i_cal = simplify(i5-i1-i2-i3-i4)
	pentagon=Polygon(p1,p2,p3,p4,p5)
	i_for = pentagon.second_moment_of_area(p)[1]
	i = int(expand(i_cal-i_for))
	assert i is 0
	#I_xx
	eq1_x=(((p1.x-p2.x)/(p1.y-p2.y))*(y-p1.y))+p1.x
	eq2_x=(((p2.x-p3.x)/(p2.y-p3.y))*(y-p2.y))+p2.x
	eq3_x=(((p3.x-p4.x)/(p3.y-p4.y))*(y-p3.y))+p3.x
	eq4_x=(((p4.x-p5.x)/(p4.y-p5.y))*(y-p4.y))+p4.x
	eq5_x=(((p5.x-p1.x)/(p5.y-p1.y))*(y-p5.y))+p5.x

	i1=integrate(eq1_x*(y**2),(y,y1,y2))
	i2=integrate(eq2_x*(y**2),(y,y2,y3))
	i3=integrate(eq3_x*(y**2),(y,y3,y4))
	i4=integrate(eq4_x*(y**2),(y,y4,y5))
	i5=integrate(eq5_x*(y**2),(y,y1,y5))

	i_cal = simplify(-i5+i1+i2+i3+i4)
	pentagon=Polygon(p1,p2,p3,p4,p5)
	i_for = pentagon.second_moment_of_area(p)[0]
	i = int(expand(i_cal-i_for))
	assert i is 0
	#Hexagon (eqation of side of Hexagon) I_yy
	eq1_y=(((p1.y-p2.y)/(p1.x-p2.x))*(x-p1.x))+p1.y
	eq2_y=(((p2.y-p3.y)/(p2.x-p3.x))*(x-p2.x))+p2.y
	eq3_y=(((p3.y-p4.y)/(p3.x-p4.x))*(x-p3.x))+p3.y
	eq4_y=(((p4.y-p5.y)/(p4.x-p5.x))*(x-p4.x))+p4.y
	eq5_y=(((p5.y-p6.y)/(p5.x-p6.x))*(x-p5.x))+p5.y
	eq6_y=(((p6.y-p1.y)/(p6.x-p1.x))*(x-p6.x))+p6.y

	i1=integrate(eq1_y*(x**2),(x,x1,x2))	
	i2=integrate(eq2_y*(x**2),(x,x2,x3))
	i3=integrate(eq3_y*(x**2),(x,x3,x4))
	i4=integrate(eq4_y*(x**2),(x,x4,x5))
	i5=integrate(eq5_y*(x**2),(x,x5,x6))
	i6=integrate(eq6_y*(x**2),(x,x1,x6))

	i_cal=simplify(i6-i1-i2-i3-i4-i5)
	hexagon=Polygon(p1,p2,p3,p4,p5,p6)
	i_for = hexagon.second_moment_of_area(p)[1]
	i = int(expand(i_for-i_cal))
	assert i is 0
	# #I_xx
	eq1_x=(((p1.x-p2.x)/(p1.y-p2.y))*(y-p1.y))+p1.x
	eq2_x=(((p2.x-p3.x)/(p2.y-p3.y))*(y-p2.y))+p2.x
	eq3_x=(((p3.x-p4.x)/(p3.y-p4.y))*(y-p3.y))+p3.x
	eq4_x=(((p4.x-p5.x)/(p4.y-p5.y))*(y-p4.y))+p4.x
	eq5_x=(((p5.x-p6.x)/(p5.y-p6.y))*(y-p5.y))+p5.x
	eq6_x=(((p6.x-p1.x)/(p6.y-p1.y))*(y-p6.y))+p6.x

	i1=integrate(eq1_x*(y**2),(y,y1,y2))	
	i2=integrate(eq2_x*(y**2),(y,y2,y3))
	i3=integrate(eq3_x*(y**2),(y,y3,y4))
	i4=integrate(eq4_x*(y**2),(y,y4,y5))
	i5=integrate(eq5_x*(y**2),(y,y5,y6))
	i6=integrate(eq6_x*(y**2),(y,y1,y6))

	i_cal=simplify(-i6+i1+i2+i3+i4+i5)
	hexagon=Polygon(p1,p2,p3,p4,p5,p6)
	i_for = hexagon.second_moment_of_area(p)[0]
	i = int(expand(i_for-i_cal))
	assert i is 0
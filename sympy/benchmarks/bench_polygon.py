from __future__ import print_function, division

import sys
from time import time
from timeit import default_timer as clock

from sympy.geometry import (Circle, Ellipse, GeometryError, Line, Point,
                            Polygon, Ray, RegularPolygon, Segment,
                            Triangle, intersection)
from random import randint

listOfPolygons = [Polygon(Point(randint(0,10),randint(0,10)),
Point(randint(0,10),randint(0,10)),
Point(randint(0,10),randint(0,10)),
Point(randint(0,10),randint(0,10))) for w in range(10)]

def create_polygon():
    "Creating Polygon"
    p1 = Polygon(Point(0, 0), Point(3, -1),Point(6, 0), Point(4, 5))

def polygon_area(listOfPolygons):
   "Polygon.area for w in range(10)"
   [listOfPolygons[w].area for w in range(10)]

def polygon_perimeter(listOfPolygons):
    "Polygon.perimeter for w in range(10)"
    [listOfPolygons[w].perimeter for w in range(10)]

def polygon_sides(listOfPolygons):
    "Polygon.sides for w in range(10)"
    [listOfPolygons[w].sides for w in range(10)]

def polygon_centroid(listOfPolygons):
    "Polygon.centroid for w in range(10)"
    [listOfPolygons[w].centroid for w in range(10)]

def polygon_second_moment(listOfPolygons):
    "Polygon.second_moment_of_area() for w in range(10)"
    [listOfPolygons[w].second_moment_of_area() for w in range(10)]

def polygon_first_moment(listOfPolygons):
   "Polygon.first_moment_of_area for w in range(10)"
   [listOfPolygons[w].first_moment_of_area() for w in range(10)]

def polygon_polar_second_moment(listOfPolygons):
    "Polygon.polar_second_moment_of_area for w in range(10)"
    [listOfPolygons[w].polar_second_moment_of_area() for w in range(10)]

def polygon_section_modulus(listOfPolygons):
    "Polygon.section_modulus for w in range(10)"
    [listOfPolygons[w].section_modulus() for w in range(10)]

def polygon_second_moment(listOfPolygons):
    "Polygon.second_moment_of_area() for w in range(10)"
    [listOfPolygons[w].second_moment_of_area() for w in range(10)]

def polygon_is_convex(listOfPolygons):
    "Polygon.is_convex() for w in range(10)"
    [listOfPolygons[w].is_convex() for w in range(10)]

def polygon_encloses_point(listOfPolygons,Point2=Point(randint(0,10),randint(0,10))):
   "Polygon.encloses_point for w in range(10)"
   [listOfPolygons[w].encloses_point(Point2) for w in range(10)]

def polygon_arbitrary_point(listOfPolygons):
    "Polygon.arbitrary_point for w in range(10)"
    [listOfPolygons[w].arbitrary_point() for w in range(10)]

def polygon_cut_section(listOfPolygons):
    "Polygon.cut_section for w in range(10)"
    [listOfPolygons[w].cut_section(Line((randint(0,10),randint(0,10)), slope=0)) for w in range(10)]

def polygon_distance(listOfPolygons):
    "Polygon.distance() for w in range(10)"
    [listOfPolygons[w].distance(Point(randint(0,10),randint(0,10))) for w in range(10)]




if __name__ == '__main__':

    t = clock()
    create_polygon()
    t = clock() - t
    print("%s%65s: %f" % (create_polygon.__name__, create_polygon.__doc__, t))
    benchmarks = [
        polygon_area,
        polygon_perimeter,
        polygon_sides,
        polygon_centroid,
        polygon_second_moment,
        polygon_first_moment,
        polygon_polar_second_moment,
        polygon_section_modulus,
        polygon_second_moment,
        polygon_is_convex,
        polygon_encloses_point,
        polygon_arbitrary_point,
        polygon_cut_section,
        polygon_distance,
    ]
    for b in benchmarks:
        t = clock()
        try:
            b(listOfPolygons)
        except ValueError:
            t = clock()-t
            print("%s%65s: Error after %f" % (b.__name__, b.__doc__, t))
            continue
        t = clock() - t
        print("%s%65s: %f" % (b.__name__, b.__doc__, t))


from shapely import Polygon, Point, LinearRing, get_coordinates
import math
from functools import reduce
from numpy import roll, corrcoef, argmax

def signature(polygon: Polygon) -> tuple[list[float],list[float]]:
    """
    The radial signature of a polygon.
    Note1: we only consider the exterior ring of the polygon.
    Note2: we don't actually use S. Should we drop it?

    :param polygon: a polygon
    :type polygon: Polygon
    :return: Description
    :rtype: tuple[list[float], list[float]]
    """
    def distance(c1: list[float], c2: list[float]) -> float:
        return math.sqrt((c2[0] - c1[0]) ** 2 + (c2[1] - c1[1]) ** 2)
    centroid: Point = polygon.centroid
    ring: LinearRing = polygon.exterior
    # we remove the last repeated point 
    coordinate_list = get_coordinates(ring).tolist()#[:-1]
    R = [distance(coord, [centroid.x, centroid.y]) for coord in coordinate_list]
    N = max(R)
    def f(x: tuple[list[float], list[float]],y: list[float])->tuple[list[float], list[float]]:
        x[0].append(distance(x[1], y) + x[0][-1])
        return (x[0], y)
    S = list(reduce(f, coordinate_list[1:], ([0.0], coordinate_list[0])))[0]
    S = list(map(lambda x: x/S[-1], S))
    R = list(map(lambda x: x/N, R))
    return (S, R)

def correlator(signature1: list[float], signature2: list[float])-> list[float]:
    """
    Computes the correlation between all shifted versions of signature1 and signature2.
    
    :param signature1: a radial signature
    :type signature1: list[float]
    :param signature2: another radial signature
    :type signature2: list[float]
    :return: the correlation for all shifted versions of signature1
    :rtype: list[float]
    """
    return [corrcoef([roll(signature1, i).tolist(), signature2],dtype=float)[0][1] for i in range(0, len(signature1))]

def radial_distance(p1: Polygon, p2: Polygon) -> float:
    """
    The radial distance between 2 polygons.
    
    :param p1: a polygon
    :type p1: Polygon
    :param p2: another polygon
    :type p2: Polygon
    :return: radial distance between the polygons (>=0)
    :rtype: float
    """
    r1 = signature(p1)
    r2 = signature(p2)
    print(correlator(r1[1],r2[1]))
    # find the shift with the maximum correlation
    max_rho = argmax(correlator(r1[1],r2[1]))
    # get the shifted signature
    r3: list[float] = roll(r1[1], max_rho).tolist()
    # compute the distance
    return math.sqrt(sum([(x1 - x2) ** 2 for x1,x2 in zip(r2[1],r3)]) / len(r3))

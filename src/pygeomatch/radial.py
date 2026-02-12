
from shapely import Polygon, LineString, Point, LinearRing, get_coordinates
import math
from functools import reduce
from numpy import roll, corrcoef, argmax, fft, abs, conj
from itertools import chain

def interpolate(c1: list[float], c2: list[float], s: float, normalized: bool = True) -> list[float]:
    """
    Interpolation between 2 coordinates.
    
    :param c1: coordinate 1
    :type c1: list[float]
    :param c2: coordinate 2
    :type c2: list[float]
    :param s: distance at which to interpolate
    :type s: float
    :param normalized: 
        If True, the distance is a fraction of the total line length (between 0 and 1) instead of the absolute distance. 
        defaults to True.
    :type normalized: bool
    :return: Description
    :rtype: list[float]
    """
    return get_coordinates(LineString([c1, c2]).interpolate(s, normalized=True)).tolist()[0]

def distance(c1: list[float], c2: list[float]) -> float:
    """
    Distance between coordinates.
    
    :param c1: coordinate 1
    :type c1: list[float]
    :param c2: coordinate 2
    :type c2: list[float]
    :return: Euclidean distance between coordinates 1 and 2
    :rtype: float
    """
    return math.sqrt((c2[0] - c1[0]) ** 2 + (c2[1] - c1[1]) ** 2)

def signature(polygon: Polygon, fs: int = 2, nb: int = 0) -> tuple[list[float], list[float]]:
    """
    The radial signature of a polygon.
    Note1: we only consider the exterior ring of the polygon.
    Note2: we don't actually use S. Should we drop it?

    :param polygon: a polygon
    :type polygon: Polygon
    :param fs: oversampling factor (f in the paper). Default value is 2. fs=1 => no oversampling
    :type fs: int
    :param nb: nb of points used for the computation. Default value is fs*size(polygon).
    :type nb: int
    :return: Description
    :rtype: tuple[list[float], list[float]]
    """
    centroid: Point = polygon.centroid
    ring: LinearRing = polygon.exterior
    # we don't remove the last (repeated) point: it will be removed by the interpolation step
    coordinate_list = get_coordinates(ring).tolist()
    print("coordinate_list",coordinate_list)
    t = [s/fs for s in range(1, fs+1)]
    print("t",t)
    list_of_lists = [
        [interpolate(coordinate_list[i], coordinate_list[i+1], s) for s in t]
        for i in range(0, len(coordinate_list) - 1)
    ]
    XYs = list(chain(*list_of_lists))
    print("XYs",XYs)
    R = [distance(coord, [centroid.x, centroid.y]) for coord in XYs]
    N = max(R)
    def f(x: tuple[list[float], list[float]], y: list[float])->tuple[list[float], list[float]]:
        x[0].append(distance(x[1], y) + x[0][-1])
        return (x[0], y)
    S = list(reduce(f, XYs[1:], ([0.0], XYs[0])))[0]
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

def radial_distance_fft(signature1: list[float], signature2: list[float])-> list[float]:
    f1 = fft.fft(signature1)
    f2 = fft.fft(signature2)
    return abs(fft.ifft(f1 * conj(f2))).tolist()

def radial_distance(p1: Polygon, p2: Polygon, fs: int = 2, nb: int = 100, normalize: bool=False, fft: bool = False) -> float:
    """
    The radial distance between 2 polygons.
    
    :param p1: a polygon
    :type p1: Polygon
    :param p2: another polygon
    :type p2: Polygon
    :param fs: oversampling factor (f in the paper)
    :type fs: int
    :param nb: nb of points used for the computation
    :type nb: int
    :return: radial distance between the polygons (>=0)
    :rtype: float
    """
    r1 = signature(p1, fs=fs, nb=nb)
    r2 = signature(p2, fs=fs, nb=nb)
    # print(correlator(r1[1],r2[1]))
    d = radial_distance_fft(r1[1],r2[1]) if fft else correlator(r1[1],r2[1])
    # find the shift with the maximum correlation
    max_rho = argmax(d)
    # get the shifted signature
    r3: list[float] = roll(r1[1], max_rho).tolist()
    # compute the distance
    return math.sqrt(sum([(x1 - x2) ** 2 for x1,x2 in zip(r2[1],r3)]) / len(r3))


from shapely import Polygon, LineString, Point, LinearRing, get_coordinates
import math
from itertools import chain
import numpy as np

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
    return get_coordinates(LineString([c1, c2]).interpolate(s, normalized=normalized)).tolist()[0]

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

def approx(x: list[float], y: list[float], n: int) -> tuple[list[float], list[float]]:
    """
    Return a list of points which linearly interpolate given data points, or a function performing the linear (or constant) interpolation.
    Tries to mimic the R function.
    The only difference is we start at max(x)/(n) instead of min(x) to have points better spaced between 0 and max(x).
    
    :param x: numeric vectors giving the coordinates of the points to be interpolated: x values
    :type x: list[float]
    :param y: numeric vectors giving the coordinates of the points to be interpolated: y values
    :type y: list[float]
    :param n: interpolation takes place at n equally spaced points spanning the interval (max(x)/(n), max(x))
    :type n: int
    :return: numeric vectors giving the coordinates of the interpolated points (x values), (y values)
    :rtype: tuple[list[float], list[float]]
    """
    # xvals = np.linspace(min(x), max(x), n) # allows to compare directly to the R version
    xvals = np.linspace(max(x)/(n), max(x), n)
    yinterp = np.interp(xvals, x, y)
    return xvals.tolist(), yinterp.tolist()

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
    t = [s/fs for s in range(1, fs+1)]
    list_of_lists = [
        [interpolate(coordinate_list[i], coordinate_list[i+1], s) for s in t]
        for i in range(0, len(coordinate_list) - 1)
    ]
    XYs = [coordinate_list[0]]
    XYs.extend(list(chain(*list_of_lists)))
    # cumulative arc‑length
    S = np.cumsum([distance(XYs[i+1],XYs[i]) for i in range(0, len(XYs)-1)]).tolist()
    # radial distances from centroid
    R = [distance(coord, [centroid.x, centroid.y]) for coord in XYs[:-1]]
    size = len(R) if (nb == 0) else nb
	# interpolate using the nb parameter
    S, R = approx(S, R, size)
    N = max(R)
    # normalize to have all values between 0 and 1
    S = [x / S[-1] for x in S]
    # for scale‑independent distance: normalize to have all values between 0 and 1
    # R = [x / N for x in R]
    return (S, R)

def correlate(signature1: list[float], signature2: list[float])-> list[float]:
    """
    Computes the correlation between all shifted versions of signature1 and signature2.
    
    :param signature1: a radial signature
    :type signature1: list[float]
    :param signature2: another radial signature
    :type signature2: list[float]
    :return: the correlation for all shifted versions of signature1
    :rtype: list[float]
    """
    return [float(np.corrcoef([np.roll(signature1, i).tolist(), signature2],dtype=float)[0][1]) for i in range(0, len(signature1))]

def radial_distance_direct(signature1: list[float], signature2: list[float])-> list[float]:
    """
    Computes the radial_distances for all shifted versions of signature1 and signature2.
    
    :param signature1: a radial signature
    :type signature1: list[float]
    :param signature2: another radial signature
    :type signature2: list[float]
    :return: the distance for all shifted versions of signature1
    :rtype: list[float]
    """
    return [float(np.sqrt(np.mean(np.power(np.diff(np.array([signature1, np.roll(signature2, i).tolist()]), axis=0),2)))) for i in range(0, len(signature1))]

def radial_distance_fft(signature1: list[float], signature2: list[float])-> list[float]:
    """
    Computes the RMS distance between two radial signatures using the FFT (circular cross-correlation). 
    The returned list has the same length as the signatures; the smallest entry corresponds to the optimal alignment.
    
    :param signature1: signature of polygon1
    :type signature1: list[float]
    :param signature2: signature of polygon2
    :type signature2: list[float]
    :return: the distances corresponding to the different shifted versions of the signatures
    :rtype: list[float]
    """
    # Convert to NumPy arrays (float64 for numerical stability)
    a = np.asarray(signature1, dtype=np.float64)
    b = np.asarray(signature2, dtype=np.float64)
    # Compute the circular cross‑correlation via FFT
    f1 = np.fft.fft(a)
    f2 = np.fft.fft(b)
    # The inverse FFT gives the *sum* of products; divide by N to get the mean.
    N = a.shape[0] # number of samples
    corr = np.real(np.fft.ifft(f1 * np.conj(f2))) / N
    # Pre‑compute the per‑signal energy terms
    norm_a = np.mean(a ** 2)      # = (1/N) * Σ a_k²
    norm_b = np.mean(b ** 2)      # = (1/N) * Σ b_k²
    # Convert correlation to RMS distance
    rms_sq = norm_a + norm_b - 2.0 * corr
    # Numerical noise can push a few entries slightly below zero; clip to zero before sqrt.
    rms = np.sqrt(np.clip(rms_sq, a_min=0.0, a_max=None))
    return rms.tolist()

def radial_distance(p1: Polygon, p2: Polygon, fs: int = 2, nb: int = 100, fft: bool = False) -> float:
    """
    The radial distance between 2 polygons.
    If `fft=True` the fast FFT-based RMS distance is used;
    otherwise the direct O(N²) RMS distance is computed.
    
    :param p1: a polygon
    :type p1: Polygon
    :param p2: another polygon
    :type p2: Polygon
    :param fs: oversampling factor (f in the paper)
    :type fs: int
    :param nb: nb of points used for the computation
    :type nb: int
    :param fft: if True, use the fast FFT-based computation (default to False)
    :type fft: bool
    :return: radial distance between the polygons (>=0)
    :rtype: float
    """
    r1 = signature(p1, fs=fs, nb=nb)
    r2 = signature(p2, fs=fs, nb=nb)
    if fft:
        d = radial_distance_fft(r1[1],r2[1])
    else:
        d = radial_distance_direct(r1[1],r2[1])
    # find the value with the minimum distance
    rho = np.argmin(d)
    # return compute the distance
    return d[rho]

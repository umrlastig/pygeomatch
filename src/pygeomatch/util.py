import shapely
def surface_distance(geomA , geomB):
    """
    Surface distance(A,B) = 1 - inter(A,B).area/union(A,B).area
    
    :param geomA: geometry A
    :param geomB: geometry B
    """
    inter = geomA.intersection(geomB)
    union = shapely.union(geomA.buffer(0),geomB.buffer(0))
    return 1 - inter.area / union.area

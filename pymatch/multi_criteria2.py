import numpy as np
import geopandas as gpd
from itertools import groupby, chain
from collections.abc import Callable
from bitarray import frozenbitarray, util, bitarray
from typing import Iterable
from functools import partial, reduce
from operator import mul, itemgetter
import shapely
from shapely.geometry import shape
class MCMatch:
    def __init__(self, matched, not_matched, theta):
        self.matched = matched
        self.not_matched = not_matched
        self.theta = theta
    def __repr__(self):
        return "MCMatch(%s,%s,%s)" % (self.matched, self.not_matched, self.theta)

def select_candidates(ref: gpd.GeoDataFrame, comp: gpd.GeoDataFrame):
    """
    Identify candidates for ref features.
    
    :param ref: ref features
    :param comp: comp features
    """
    # The first subarray contains input geometry integer indices. The second subarray contains tree geometry integer indices.
    (refIndices, compIndices) = comp.geometry.sindex.query(ref.geometry, predicate="intersects")
    # zip both lists into tuples (refIndex, compIndex)
    zipped = zip(refIndices, compIndices)
    # group by the ref index (z[0]), map the results to dict entries with the ref index (y[0]) and keep only the second element of the grouped tuples (t[1])
    return dict(map(lambda y: (y[0], list(map(lambda t: t[1], y[1]))), groupby(zipped, lambda z: z[0])))

def get_matched_array(index, size) -> frozenbitarray:
    matched = util.zeros(size)
    matched.invert(index)
    return frozenbitarray(matched)

def get_potential_set(index: int, size: int, match: MCMatch) -> dict[bitarray,float]:
    matched = get_matched_array(index, size)
    theta = frozenbitarray(util.ones(size)) #ignorance
    return dict({matched: match.matched, frozenbitarray(~matched): match.not_matched, theta: match.theta})

# note: the threshold might not be necessary for this example but it becomes useful for very large frames
def combine(potentialSets: Iterable[dict[bitarray, float]], threshold: float = 0.00001):
    def combination2(potentialSet1: dict[bitarray, float],potentialSet2: dict[bitarray, float]):
        # print(len(potentialSet1),len(potentialSet2))
        result = dict()
        for f1, v1 in potentialSet1.items():
            def add(dict, key, value):
                if value>threshold: dict[key] = value if key not in dict else dict[key]+value
            for f2, v2 in potentialSet2.items():add(result, f1&f2, v1*v2)
        return result
    return reduce(combination2, potentialSets)

# a purely functional version
def combination_func(potentialSets, threshold: float = 0.00001):
    def combination_func2(p1,p2):
        def comb(p,x): return map(lambda y: (x[0]&y[0],x[1]*y[1]), p)
        return dict(map(lambda z: (z[0], sum(map(lambda h: h[1], z[1]))), 
                        groupby(sorted(filter(lambda z: z[1]>threshold, chain.from_iterable(map(partial(comb,p2.items()), p1.items())))), lambda t: t[0])))
    return dict(reduce(combination_func2, potentialSets))

def normalize(potentialSet: dict[bitarray, float]):
    """
    Normalize the masses according to the conflict.
    
    :param potentialSet: Description
    :type potentialSet: dict[frozenbitarray, float]
    """
    setSize = len(next(iter(potentialSet)))
    theta = frozenbitarray(util.zeros(setSize))
    conflict = potentialSet.pop(theta,None)
    # print("conflict",conflict)
    if conflict:
        return dict(map(lambda item: (item[0], item[1]/(1-conflict)), potentialSet.items()))
    return potentialSet

def pignistic_probability(potentialSet: dict[bitarray, float], threshold: float = 0.00001):
    """
    Docstring for pignisticProbability
    
    :param potentialSet: Description
    """
    result = dict[bitarray, float]()
    for f1, v1 in potentialSet.items():
        # print("prob",f1,v1)
        def add(dict: dict[bitarray, float], a_and_b: int, b: int, mass_b: float, key_b: bitarray):
            if util.subset(f1, key_b):
                # print("with",key_b,"a&b",a_and_b,"b",b,"mass",mass_b)
                v = mass_b * a_and_b / b
                if v>threshold: dict[f1] = v if f1 not in dict else dict[f1]+v
        for f2, v2 in potentialSet.items():add(result, util.count_and(f1,f2), f2.count(1), v2, f2)
    return result

def process_match(refFeature: dict, compFeatures: gpd.GeoDataFrame, criteria: list[Callable[[dict,dict],MCMatch]]) -> dict|None:
    # print("refFeature",type(refFeature),refFeature["geometry"])
    print(len(compFeatures),"candidates")
    candidates = len(compFeatures) + 1 # +1 since we'll add the 'not matched' candidate
    theta = frozenbitarray(util.ones(candidates)) #ignorance
    phi = frozenbitarray(util.zeros(candidates)) # conflict
    # we combine the criteria for all candidates
    potentialSets = [combine(list(map(lambda c: get_potential_set(i, candidates, c(refFeature, f)),criteria))) for i, f in enumerate(compFeatures.iterfeatures())]
    #flatPotentialSets = list(chain.from_iterable(potentialSets))
    # print("potentialSets",potentialSets)
    # combine and normalise the potential sets
    # fusion = normalize(combine(potentialSets))
    # print("fusion",fusion)
    # adding the not_matched hypothesis
    not_matched = get_matched_array(candidates-1, candidates)
    # the list of not_matched from other hypotheses
    #flatPotentialSets.append({not_matched:0})
    def getOrZero(index, candidates, fusion):
        if ~get_matched_array(index, candidates) in fusion:
            return fusion[~get_matched_array(index, candidates)]
        return 0.
    masses = [getOrZero(index, candidates, potentialSets) for index in range(candidates-1)]
    mass = reduce(mul, masses, 1)
    # add a source with the not_matched hypothesis
    potentialSets.append({not_matched:mass, theta: 1-mass})
    # fusion[not_matched] = mass
    final = normalize(combine(potentialSets))
    #print("final",final)
    pignistic = pignistic_probability(final)
    #print("pignistic",pignistic)
    maxPignistic = max(pignistic.items(), key=itemgetter(1))
    #print("max",maxPignistic)
    # if the max is not a simple hypothesis or is the 'not_matched' hypothesis
    if maxPignistic[0].count(1) > 1 | (maxPignistic[0] == not_matched):
        return None
    index = maxPignistic[0].find(bitarray(1))
    # FIXME this is kinda ugly
    return list(compFeatures.iterfeatures())[index-1]

def geom_criteria(a: dict, b: dict) -> MCMatch:
    # print(shape(a["geometry"]))
    # print(shape(b["geometry"]))
    geom_a = shape(a["geometry"]).buffer(0)
    geom_b = shape(b["geometry"]).buffer(0)
    inter = shapely.intersection(geom_a, geom_b)
    union = shapely.union(geom_a, geom_b)
    distance = 1 - inter.area / union.area
    T1 = 0.90
    T2 = 1.0
    E = 0.01
    S = 0.7
    K = 1 - E - S
    if distance < T1:
        app = (-(1 - E)/T2) * distance + 1 - E
        _app = E
    elif distance < T2:
        app = (-(1 - E)/T2) * distance + 1 - E
        _app = (K - E) * distance / (T2 - T1) + E - (K - E)*T1/(T2 - T1)
    else:
        app = E
        _app = K
    return MCMatch(app ,_app, 1 - app - _app)

def MCA2(ref: gpd.GeoDataFrame, comp: gpd.GeoDataFrame):
    """
    Process Multi Criteria Matching.
    
    :param ref: Ref features
    :param comp: Comp features
    """
    # get the ref features and their corresponding candidates (if they have any)
    candidateDictionary = select_candidates(ref, comp)
    results = [process_match(next(ref.iloc[[k]].iterfeatures()), comp.iloc[v], [geom_criteria]) for k,v in candidateDictionary.items()]
    return [result for result in results if result is not None]

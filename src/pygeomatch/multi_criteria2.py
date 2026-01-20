import geopandas as gpd
from itertools import groupby, chain
from collections.abc import Callable
from bitarray import frozenbitarray, util, bitarray
from typing import Iterable, Union
from functools import partial, reduce
from operator import mul, itemgetter
from shapely.geometry import shape
from pygeomatch.util import surface_distance
class MCMatch:
    def __init__(self, matched, not_matched, theta):
        self.matched = matched
        self.not_matched = not_matched
        self.theta = theta
    def __repr__(self):
        return "MCMatch(%s,%s,%s)" % (self.matched, self.not_matched, self.theta)

def select_candidates(ref: gpd.GeoDataFrame, comp: gpd.GeoDataFrame) -> dict[int, list[int]]:
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
    """
    Docstring for get_matched_array
    
    :param index: Description
    :param size: Description
    :return: Description
    :rtype: frozenbitarray
    """
    matched = util.zeros(size)
    matched.invert(index)
    return frozenbitarray(matched)

def get_potential_set(index: int, size: int, match: MCMatch) -> dict[bitarray,float]:
    """
    Get the masses for the hypotheses.
    
    :param index: Index of a feature
    :type index: int
    :param size: number of candidates
    :type size: int
    :param match: the match values
    :type match: MCMatch
    :return: a dictionay with the masses for the matched, not matched and ignorance hypotheses
    :rtype: dict[bitarray, float]
    """
    matched = get_matched_array(index, size)
    theta = frozenbitarray(util.ones(size)) #ignorance
    return dict({matched: match.matched, frozenbitarray(~matched): match.not_matched, theta: match.theta})

def combine(potentialSets: Iterable[dict[bitarray, float]], threshold: float = 0.00001):
    """
    Combine the potential sets.

    note: the threshold might not be necessary for this example but it becomes useful for very large frames
    
    :param potentialSets: Description
    :type potentialSets: Iterable[dict[bitarray, float]]
    :param threshold: Description
    :type threshold: float
    """
    def combination2(potentialSet1: dict[bitarray, float],potentialSet2: dict[bitarray, float]):
        result = dict()
        for f1, v1 in potentialSet1.items():
            def add(dict, key, value):
                if value>threshold:
                    dict[key] = value if key not in dict else dict[key]+value
            for f2, v2 in potentialSet2.items():
                add(result, f1&f2, v1*v2)
        return result
    return reduce(combination2, potentialSets)

def combination_func(potentialSets, threshold: float = 0.00001):
    """
    Combine the potential sets (a purely functional version).

    note: the threshold might not be necessary for this example but it becomes useful for very large frames
    
    :param potentialSets: Description
    :param threshold: Description
    :type threshold: float
    """
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
        def add(dict: dict[bitarray, float], a_and_b: int, b: int, mass_b: float, key_b: bitarray):
            if util.subset(f1, key_b):
                v = mass_b * a_and_b / b
                if v>threshold:
                    dict[f1] = v if f1 not in dict else dict[f1]+v
        for f2, v2 in potentialSet.items():
            add(result, util.count_and(f1,f2), f2.count(1), v2, f2)
    return result

def process_match(refIndex: int, refFeature: dict, compFeatures: gpd.GeoDataFrame, criteria: list[Callable[[dict,dict],MCMatch]]) -> Union[tuple,None]:
    """
    Docstring for process_match
    
    :param refIndex: Description
    :type refIndex: int
    :param refFeature: Description
    :type refFeature: dict
    :param compFeatures: Description
    :type compFeatures: gpd.GeoDataFrame
    :param criteria: Description
    :type criteria: list[Callable[[dict, dict], MCMatch]]
    :return: Description
    :rtype: tuple[Any, ...] | None
    """
    candidates = len(compFeatures) + 1 # +1 since we'll add the 'not matched' candidate
    theta = frozenbitarray(util.ones(candidates)) #ignorance
    #phi = frozenbitarray(util.zeros(candidates)) # conflict
    # we combine the criteria for all the candidates
    potentialSets = [combine(list(map(lambda c: get_potential_set(i, candidates, c(refFeature, f)), criteria))) for i, f in enumerate(compFeatures.iterfeatures())]
    fusion_of_criteria = combine(potentialSets)
    # adding the not_matched hypothesis
    not_matched = get_matched_array(candidates-1, candidates)
    # the list of not_matched from other hypotheses
    # create all the not matched hypotheses (except for the not_matched one)
    all_not_matched_hypotheses = [~get_matched_array(index, candidates) for index in range(candidates-1)]
    if all([c in fusion_of_criteria for c in all_not_matched_hypotheses]):
        all_not_matched_hypotheses_masses = [fusion_of_criteria[c] for c in all_not_matched_hypotheses]
        # we multiply all the masses to get the mass for the not_matched hypothesis
        mass = reduce(mul, all_not_matched_hypotheses_masses, 1)
    else:
        # there was at least a candidate for which there was not mass so...
        mass = 0.0
    # TODO we might want not to add the not_matched hypothesis if its mass is zero...
    # add a source with the not_matched hypothesis
    fusion_of_candidates = normalize(combine([fusion_of_criteria, {not_matched:mass, theta: 1-mass}]))
    pignistic = pignistic_probability(fusion_of_candidates)
    maxPignistic, maxPignisticProbability = max(pignistic.items(), key=itemgetter(1))
    # if the max is not a simple hypothesis or is the 'not_matched' hypothesis
    if (maxPignistic.count(1) > 1) | (maxPignistic == not_matched):
        return None
    index = maxPignistic.find(1)
    return (refIndex, compFeatures.iloc[index].name, maxPignisticProbability)

def geom_criteria(a: dict, b: dict) -> MCMatch:
    """
    A simple geometric criteria using the surface geometry.
    
    :param a: Description
    :type a: dict
    :param b: Description
    :type b: dict
    :return: Description
    :rtype: MCMatch
    """
    distance = surface_distance(shape(a["geometry"]), shape(b["geometry"]))
    T1 = 0.9
    T2 = 1.
    E = 0.01
    S = 0.7
    K = 1 - E - S
    m = 1 - E
    if distance < T1:
        app = (-m/T2) * distance + m
        _app = E
    elif distance < T2:
        app = (-m/T2) * distance + m
        _app = (K - E) * distance / (T2 - T1) + E - (K - E)*T1/(T2 - T1)
    else:
        app = E
        _app = K
    return MCMatch(app ,_app, 1 - app - _app)

def MCA2(ref: gpd.GeoDataFrame, comp: gpd.GeoDataFrame, criteria = [geom_criteria])->list:
    """
    Process Multi Criteria Matching.
    
    :param ref: Ref features
    :param comp: Comp features
    :param criteria: a list of criteria to compare features. Default value is a list with only a surface_geometry criteria.
    """
    # get the ref features and their corresponding candidates (if they have any)
    candidateDictionary = select_candidates(ref, comp)
    results = [process_match(k, next(ref.iloc[[k]].iterfeatures()), comp.iloc[v], criteria) for k,v in candidateDictionary.items()] # type: ignore
    return [(result[0], comp.index.get_loc(result[1]),result[2]) for result in results if result is not None]

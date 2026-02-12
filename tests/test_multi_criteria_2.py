from bitarray import frozenbitarray, util
import pytest
from pygeomatch.util import surface_distance
from pygeomatch.multi_criteria2 import MCMatch, get_potential_set, combine, get_matched_array, normalize, pignistic_probability, geom_criteria, process_match, select_candidates, MCA2, combination_func
from functools import reduce
from operator import mul, itemgetter
from shapely import Polygon, geometry
import geopandas as gpd

class TestMultiCriteria2:
    @pytest.fixture(scope="class")
    def polygon1(self):
        return Polygon(((0., 0.), (0., 1.), (1., 1.), (1., 0.), (0., 0.)))

    @pytest.fixture(scope="class")
    def polygon2(self):
        return Polygon(((0., 0.), (0., 1.), (0.4, 1.), (0.4, 0.), (0., 0.)))
    
    @pytest.fixture(scope="class")
    def polygon3(self):
        return Polygon(((2., 0.), (2., 1.), (3., 1.), (3., 0.), (2., 0.)))

    @pytest.fixture(scope="class")
    def polygon4(self):
        return Polygon(((0., 0.), (0., 1.), (0.1, 1.), (0.1, 0.), (0., 0.)))

    def test_surface_distance(self, polygon1, polygon2, polygon3, polygon4):
        s = surface_distance(polygon1, polygon2)
        assert s == pytest.approx(0.6)
        s = surface_distance(polygon1, polygon3)
        assert s == pytest.approx(1.0)
        s = surface_distance(polygon1, polygon1)
        assert s == pytest.approx(0.0)
        s = surface_distance(polygon1, polygon4)
        assert s == pytest.approx(0.9)
    
    def test_geom_criteria(self, polygon1, polygon2, polygon3, polygon4):
        a = {"geometry":geometry.mapping(polygon1)}
        b = {"geometry":geometry.mapping(polygon2)}
        c = {"geometry":geometry.mapping(polygon3)}
        d = {"geometry":geometry.mapping(polygon4)}
        criteria = geom_criteria(a, b)
        print("geom_criteria",criteria)
        assert criteria.matched == pytest.approx(0.01)#0.396
        criteria = geom_criteria(a, c)
        print("geom_criteria",criteria)
        assert criteria.matched == pytest.approx(0.01)
        criteria = geom_criteria(a, d)
        print("geom_criteria",criteria)
        assert criteria.matched == pytest.approx(0.01)#0.099
        criteria = geom_criteria(a, a)
        print("geom_criteria",criteria)
        assert criteria.matched == pytest.approx(1.0)#0.99

    def test_combine(self):
        candidates = 2
        m1c1 = MCMatch(0.4,0.0,0.6)
        m2c1 = MCMatch(0.3,0.0,0.7)
        mc1 = combine([get_potential_set(0, candidates, m1c1),get_potential_set(0, candidates, m2c1)])
        c1 = get_matched_array(0, candidates)
        theta = frozenbitarray(util.ones(candidates)) #ignorance
        assert mc1[c1] == pytest.approx(0.58)
        assert mc1[theta] == pytest.approx(0.42)
        mc1 = combination_func([get_potential_set(0, candidates, m1c1),get_potential_set(0, candidates, m2c1)])
        assert mc1[c1] == pytest.approx(0.58)
        assert mc1[theta] == pytest.approx(0.42)

    def test_normalize(self):
        candidates = 2
        m1c1 = MCMatch(0.4,0.0,0.6)
        m2c1 = MCMatch(0.3,0.0,0.7)
        mc1 = combine([get_potential_set(0, candidates, m1c1),get_potential_set(0, candidates, m2c1)])
        normalized = normalize(mc1)
        print(normalized)
        c1 = get_matched_array(0, candidates)
        theta = frozenbitarray(util.ones(candidates)) #ignorance
        assert normalized[c1] == pytest.approx(0.58)
        assert normalized[theta] == pytest.approx(0.42)

    def test_criteria_fusion(self):
        """
        From Ana-Maria Olteanu's PhD Thesis, p.110-113.
        """
        candidates = 4 # including the not_matched hypothesis
        m1c1 = MCMatch(0.4,0.0,0.6)
        m2c1 = MCMatch(0.3,0.0,0.7)
        m1c2 = MCMatch(0.1,0.9,0.0)
        m2c2 = MCMatch(0.0,1.0,0.0)
        m1c3 = MCMatch(0.35,0.0,0.65)
        m2c3 = MCMatch(0.3,0.0,0.7)
        mc1 = combine([get_potential_set(0, candidates, m1c1),get_potential_set(0, candidates, m2c1)])
        mc2 = combine([get_potential_set(1, candidates, m1c2),get_potential_set(1, candidates, m2c2)])
        mc3 = combine([get_potential_set(2, candidates, m1c3),get_potential_set(2, candidates, m2c3)])
        print("c1",mc1)
        print("c2",mc2)
        print("c3",mc3)
        c1 = get_matched_array(0, candidates)
        c2 = get_matched_array(1, candidates)
        c3 = get_matched_array(2, candidates)
        theta = frozenbitarray(util.ones(candidates)) #ignorance
        phi = frozenbitarray(util.zeros(candidates)) # conflict
        print("mc1[c1]",mc1[c1])
        print("mc1[theta]",mc1[theta])
        assert mc1[c1] == pytest.approx(0.58)
        assert mc1[theta] == pytest.approx(0.42)
        print("mc2[~c2]",mc2[~c2])
        print("mc2[phi]",mc2[phi])
        assert mc2[~c2] == pytest.approx(0.9)
        assert mc2[phi] == pytest.approx(0.1)
        print("mc3[c3]",mc3[c3])
        print("mc3[theta]",mc3[theta])
        assert mc3[c3] == pytest.approx(0.545)#approximated further to 0.54 in the manuscript
        assert mc3[c3] == pytest.approx(0.54, rel=1e-2)
        assert mc3[theta] == pytest.approx(0.455)#approximated further to 0.46 in the manuscript
        assert mc3[theta] == pytest.approx(0.46, rel=2e-2)

    def test_candidate_fusion(self):
        """
        From Ana-Maria Olteanu's PhD Thesis, p.110-118.
        """
        candidates = 4
        m1c1 = MCMatch(0.4,0.0,0.6)
        m2c1 = MCMatch(0.3,0.0,0.7)
        m1c2 = MCMatch(0.1,0.9,0.0)
        m2c2 = MCMatch(0.0,1.0,0.0)
        m1c3 = MCMatch(0.35,0.0,0.65)
        m2c3 = MCMatch(0.3,0.0,0.7)
        mc1 = combine([get_potential_set(0, candidates, m1c1),get_potential_set(0, candidates, m2c1)])
        mc2 = combine([get_potential_set(1, candidates, m1c2),get_potential_set(1, candidates, m2c2)])
        mc3 = combine([get_potential_set(2, candidates, m1c3),get_potential_set(2, candidates, m2c3)])
        potentialSets = [mc1,mc2,mc3]
        intermediate = combine(potentialSets)
        # adding the not_matched hypothesis
        not_matched = get_matched_array(candidates-1, candidates)
        # the list of not_matched from other hypotheses
        def getOrZero(index, candidates, fusion):
            if ~get_matched_array(index, candidates) in fusion:
                return fusion[~get_matched_array(index, candidates)]
            return 0.
        masses = [getOrZero(index, candidates, intermediate) for index in range(candidates-1)]
        for m in masses:
            print(m)
        # we multiply all the masses to get the mass for the not_matched hypothesis
        mass = reduce(mul, masses, 1)
        # add a source with the not_matched hypothesis
        theta = frozenbitarray(util.ones(candidates)) #ignorance
        potentialSets.append({not_matched:mass, theta: 1-mass})
        final = combine(potentialSets)
        c1 = get_matched_array(0, 4)
        c2 = get_matched_array(1, 4)
        c3 = get_matched_array(2, 4)
        theta = frozenbitarray(util.ones(4)) #ignorance
        phi = frozenbitarray(util.zeros(4)) # conflict
        print("final",final)
        print("final[c1]",final[c1])
        assert final[c1] == pytest.approx(0.23751) # further approximated to 0.24 in the manuscript
        assert final[c1] == pytest.approx(0.24, rel=2e-2)
        print("final[c3]",final[c3])
        assert final[c3] == pytest.approx(0.20601) # 0.2 in the manuscript
        assert final[c3] == pytest.approx(0.20, rel=4e-2)
        print("final[~c2]",final[~c2])
        assert final[~c2] == pytest.approx(0.17199) # 0.17 in the manuscript
        assert final[~c2] == pytest.approx(0.17, rel=2e-2)
        print("final[phi]",final[phi])
        assert final[phi] == pytest.approx(0.38449) # 0.39 in the manuscript
        assert final[phi] == pytest.approx(0.39, rel=2e-2)
        norm = normalize(final)
        print("norm",norm)
        print("norm[c1]",norm[c1])
        assert norm[c1] == pytest.approx(0.385875) # further approximated to 0.4 in the manuscript
        assert norm[c1] == pytest.approx(0.40, rel=4e-2)
        print("norm[c3]",norm[c3])
        assert norm[c3] == pytest.approx(0.334698) # 0.34 in the manuscript
        assert norm[c3] == pytest.approx(0.34, rel=2e-2)
        print("norm[~c2]",norm[~c2])
        assert norm[~c2] == pytest.approx(0.279427) # 0.288 in the manuscript
        assert norm[~c2] == pytest.approx(0.288, rel=3e-2)
        pignistic = pignistic_probability(norm)
        maxPignistic, maxPignisticProbability = max(pignistic.items(), key=itemgetter(1))
        print("pignistic",pignistic)
        print("pignistic[c1]",pignistic[c1])
        assert maxPignisticProbability == pytest.approx(0.479017) # 0.5 in the manuscript
        assert maxPignisticProbability == pytest.approx(0.5, rel=5e-2)
        print("maxPignistic",maxPignistic)
        assert maxPignistic == c1
        # assert False
    
    def test_process_match(self, polygon1):
        ref = {}
        comp = gpd.GeoDataFrame({"id": [0,1,2]}, geometry=[polygon1, polygon1, polygon1])
        print(comp.head())
        def criteria1(a: dict, b: dict) -> MCMatch:
            if b["id"] == '0':
                return MCMatch(0.4,0.0,0.6)
            elif b["id"] == '1':
                return MCMatch(0.1,0.9,0.0)
            return MCMatch(0.35,0.0,0.65)
        def criteria2(a: dict, b: dict) -> MCMatch:
            if b["id"] == '0':
                return MCMatch(0.3,0.0,0.7)
            elif b["id"] == '1':
                return MCMatch(0.0,1.0,0.0)
            return MCMatch(0.3,0.0,0.7)
        res = process_match(0, ref, comp, [criteria1, criteria2])
        assert res is not None
        _, maxPignistic, maxPignisticProbability = res
        assert maxPignistic == 0
        assert maxPignisticProbability == pytest.approx(0.479017) # 0.5 in the manuscript
        assert maxPignisticProbability == pytest.approx(0.5, rel=5e-2)

    def test_select_candidates(self, polygon1, polygon2, polygon3):
        gpd1 = gpd.GeoDataFrame({"id": [0]}, geometry=[polygon1])
        gpd2 = gpd.GeoDataFrame({"id": [0,1]}, geometry=[polygon2,polygon3])
        res = select_candidates(gpd1,gpd2)
        assert len(res[0]) == 1 # only 1 match
        assert res[0][0] == 0 # the match is 0-0

    def test_MCA2(self, polygon1, polygon2):
        gpd1 = gpd.GeoDataFrame({"id": [0]}, geometry=[polygon1])
        gpd2 = gpd.GeoDataFrame({"id": [0,1,2]}, geometry=[polygon1,polygon2,polygon2])
        def criteria1(a: dict, b: dict) -> MCMatch:
            if b["id"] == '0':
                return MCMatch(0.4,0.0,0.6)
            elif b["id"] == '1':
                return MCMatch(0.1,0.9,0.0)
            return MCMatch(0.35,0.0,0.65)
        def criteria2(a: dict, b: dict) -> MCMatch:
            if b["id"] == '0':
                return MCMatch(0.3,0.0,0.7)
            elif b["id"] == '1':
                return MCMatch(0.0,1.0,0.0)
            return MCMatch(0.3,0.0,0.7)
        res = MCA2(gpd1,gpd2, [criteria1,criteria2])
        assert len(res) == 1 # only 1 match
        assert res[0][1] == 0 # the match is 0-0
        assert res[0][2] == pytest.approx(0.5, rel=5e-2) # pignistic probability is 0.5

    def test_MCA2_not_matched(self, polygon1, polygon2):
        gpd1 = gpd.GeoDataFrame({"id": [0]}, geometry=[polygon1])
        gpd2 = gpd.GeoDataFrame({"id": [0]}, geometry=[polygon2])
        def criteria1(a: dict, b: dict) -> MCMatch:
            return MCMatch(0.1,0.6,0.3)
        def criteria2(a: dict, b: dict) -> MCMatch:
            return MCMatch(0.0,0.5,0.5)
        res = MCA2(gpd1,gpd2, [criteria1,criteria2])
        assert len(res) == 0 # no match
    
    # def test_MCA2_GT(self):
    #     gpd1 = gpd.read_file("data/popRef.shp")
    #     gpd2 = gpd.read_file("data/popComp.shp")
    #     gpd3 = gpd.read_file("data/popEvolution.shp")
    #     res = MCA2(gpd1,gpd2, [geom_criteria])

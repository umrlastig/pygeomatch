import pytest
import pygeomatch.geometric_matching_of_areas as gmoa
from shapely import Polygon
import geopandas as gpd

class TestGMA:
    @pytest.fixture(scope="class")
    def polygonA(self):
        return Polygon(((0., 0.), (0., 2.), (1., 2.), (1., 0.), (0., 0.)))

    @pytest.fixture(scope="class")
    def polygonB(self):
        return Polygon(((1., 0.), (1., 2.), (2., 2.), (2., 0.), (1., 0.)))

    @pytest.fixture(scope="class")
    def polygonC(self):
        return Polygon(((2., 0.), (2., 2.), (3., 2.), (3., 0.), (2., 0.)))

    @pytest.fixture(scope="class")
    def polygon1(self):
        return Polygon(((0., 1.), (0., 2.), (1.8, 2.), (1.8, 1.), (0., 1.)))

    @pytest.fixture(scope="class")
    def polygon2(self):
        return Polygon(((0., 0.), (0., 1.), (1.8, 1.), (1.8, 0.), (0., 0.)))

    @pytest.fixture(scope="class")
    def polygon3(self):
        return Polygon(((1.8, 0.), (1.8, 2.), (3., 2.), (3., 0.), (1.8, 0.)))

    @pytest.fixture(scope="class")
    def polygon4(self):
        return Polygon(((0., 0.), (0., 2.), (1.2, 2.), (1.2, 0.), (0., 0.)))

    @pytest.fixture(scope="class")
    def polygon5(self):
        return Polygon(((1.2, 0.), (1.2, 2.), (2., 2.), (2., 0.), (1.2, 0.)))

    @pytest.fixture(scope="class")
    def gpdA1(self, polygonA, polygonB):
        return gpd.GeoDataFrame({"id": [0, 1]}, geometry=[polygonA, polygonB])

    @pytest.fixture(scope="class")
    def gpdA2(self, polygonA, polygonB, polygonC):
        return gpd.GeoDataFrame({"id": [0, 1, 2]}, geometry=[polygonA, polygonB, polygonC])

    @pytest.fixture(scope="class")
    def gpdB1(self, polygon1, polygon2):
        return gpd.GeoDataFrame({"id": [0, 1]}, geometry=[polygon1, polygon2])

    @pytest.fixture(scope="class")
    def gpdB2(self, polygon1, polygon2, polygon3):
        return gpd.GeoDataFrame({"id": [0, 1, 2]}, geometry=[polygon1, polygon2, polygon3])

    @pytest.fixture(scope="class")
    def gpdB3(self, polygon4, polygon5):
        return gpd.GeoDataFrame({"id": [0, 1]}, geometry=[polygon4, polygon5])

    def test_accuracy(self, polygonA, polygon1, polygon4, polygon5):
        a11 = gmoa.get_accuracy(polygonA, polygonA)
        assert a11 == pytest.approx(1.0)
        a12 = gmoa.get_accuracy(polygonA, polygon1)
        assert a12 == pytest.approx(0.5)
        a13 = gmoa.get_accuracy(polygonA, polygon4)
        assert a13 == pytest.approx(1.0)
        a14 = gmoa.get_accuracy(polygonA, polygon5)
        assert a14 == pytest.approx(0.)

    def test_completeness(self, polygonA, polygon1, polygon4, polygon5):
        a11 = gmoa.get_completeness(polygonA, polygonA)
        assert a11 == pytest.approx(1.0)
        a12 = gmoa.get_completeness(polygon1, polygonA)
        assert a12 == pytest.approx(0.5)
        a13 = gmoa.get_completeness(polygon4, polygonA)
        assert a13 == pytest.approx(1.0)
        a14 = gmoa.get_completeness(polygonA, polygon5)
        assert a14 == pytest.approx(0.)

    def test_get_geom(self, polygonA, polygonB, polygon1, polygon2, gpdA1, gpdB1):
        assert gmoa.get_geom(0, gpdA1) == polygonA
        assert gmoa.get_geom(1, gpdA1) == polygonB
        assert gmoa.get_geom(0, gpdB1) == polygon1
        assert gmoa.get_geom(1, gpdB1) == polygon2

    def test_pre_match(self, gpdA1, gpdB1):
        links = gmoa.pre_match(gpdA1, gpdB1, {"minimise_surface_distance": True, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1})
        print("links",len(links))
        for l in links:
            print(l)
            if (l.ref == 0) and (l.comp == 0):
                assert l.intersection_ratio == pytest.approx(1.0/1.8)
                assert l.surface_distance == pytest.approx(1.8/2.8)
            elif (l.ref == 0) and (l.comp == 1):
                assert l.intersection_ratio == pytest.approx(1.0/1.8)
                assert l.surface_distance == pytest.approx(1.8/2.8)
            elif (l.ref == 1) and (l.comp == 0):
                assert l.intersection_ratio == pytest.approx(0.8/1.8)
                assert l.surface_distance == pytest.approx(2.2/3.0)
            elif (l.ref == 1) and (l.comp == 1):
                assert l.intersection_ratio == pytest.approx(0.8/1.8)
                assert l.surface_distance == pytest.approx(2.2/3.0)
        links = gmoa.pre_match(gpdA1, gpdB1, {"minimise_surface_distance": False, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1})
        print("links",len(links))
        for l in links:
            print(l)
            if (l.ref == 0) and (l.comp == 0):
                assert l.intersection_ratio == pytest.approx(1.0/1.8)
                assert l.accuracy == pytest.approx(0.5)
                assert l.completeness == pytest.approx(1.0/1.8)
            elif (l.ref == 0) and (l.comp == 1):
                assert l.intersection_ratio == pytest.approx(1.0/1.8)
                assert l.accuracy == pytest.approx(0.5)
                assert l.completeness == pytest.approx(1.0/1.8)
            elif (l.ref == 1) and (l.comp == 0):
                assert l.intersection_ratio == pytest.approx(0.8/1.8)
                assert l.accuracy == pytest.approx(0.8/2.0)
                assert l.completeness == pytest.approx(0.8/1.8)
            elif (l.ref == 1) and (l.comp == 1):
                assert l.intersection_ratio == pytest.approx(0.8/1.8)
                assert l.accuracy == pytest.approx(0.8/2.0)
                assert l.completeness == pytest.approx(0.8/1.8)
    
    def test_group_evaluation(self, gpdA1, gpdB1):
        params = {"minimise_surface_distance": True, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1}
        links = gmoa.pre_match(gpdA1, gpdB1, params)
        print(links)
        d = gmoa.group_evaluation(links, gpdA1, gpdB1, params)
        assert d == pytest.approx(0.1)
        print(d)
        params = {"minimise_surface_distance": False, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1}
        links = gmoa.pre_match(gpdA1, gpdB1, params)
        d = gmoa.group_evaluation(links, gpdA1, gpdB1, params)
        assert d == pytest.approx(1.9)
        print(d)

    def test_search_optimal_groups(self, gpdA1, gpdB1, gpdA2, gpdB2, gpdB3):
        print("A1-B1 minimise_surface_distance=True")
        params = {"minimise_surface_distance": True, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1, "sure_intersection_percentage": 0.9}
        links = gmoa.pre_match(gpdA1, gpdB1, params)
        assert len(links) == 4
        groups = gmoa.search_optimal_groups(links, gpdA1, gpdB1, params)
        print("search_optimal_groups",len(groups),groups)
        assert len(groups) == 4
        bools = [((group.ref == 0) and (group.ref == 0)) or ((group.ref == 0) and (group.ref == 1)) or ((group.ref == 1) and (group.ref == 0)) or ((group.ref == 1) and (group.ref == 1)) for group in groups]
        assert all(bools)
        print("A1-B3 minimise_surface_distance=True")
        params = {"minimise_surface_distance": True, "min_surface_intersection": 0.01, "min_intersection_percentage": 0.1, "sure_intersection_percentage": 0.7}
        links = gmoa.pre_match(gpdA1, gpdB3, params)
        print(links)
        assert len(links) == 3
        groups = gmoa.search_optimal_groups(links, gpdA1, gpdB3, params)
        print("search_optimal_groups",len(groups),groups)
        assert len(groups) == 3
        print("A1-B3 minimise_surface_distance=False")
        params = {"minimise_surface_distance": False, "min_surface_intersection": 0.01, "min_intersection_percentage": 0.1, "sure_intersection_percentage": 0.7}
        links = gmoa.pre_match(gpdA1, gpdB3, params)
        print(links)
        assert len(links) == 3
        groups = gmoa.search_optimal_groups(links, gpdA1, gpdB3, params)
        print("search_optimal_groups",len(groups),groups)
        assert len(groups) == 3

    def test_filter_links(self):
        link1 = gmoa.GMALink(0,0,0.0,0.0,1.0,0.0,0.0)
        link2 = gmoa.GMALink(0,0,0.0,0.0,0.0,0.0,1.0)
        link3 = gmoa.GMALink(0,0,0.0,0.0,0.0,1.0,0.0)
        link4 = gmoa.GMALink(0,0,0.0,0.0,0.0,1.0,1.0)
        params = {
            "minimise_surface_distance": True, 
            "max_surface_distance": 0.5
        }
        links = gmoa.filter_links([link1, link2, link3, link4], params, None, None)
        assert len(links) == 3
        params = {
            "minimise_surface_distance": False, 
            "min_accuracy_completeness": 0.5
        }
        links = gmoa.filter_links([link1, link2, link3, link4], params, None, None)
        assert len(links) == 1
    def test_surface_match(self, gpdA1, gpdB1, gpdA2, gpdB2, gpdB3):
        params = {
            "minimise_surface_distance": True, 
            "min_surface_intersection": 0.1, 
            "min_intersection_percentage": 0.1, 
            "sure_intersection_percentage": 0.6, 
            "use_optimal_groups": True, 
            "final_filtering": False,
            "max_surface_distance": 0.5
        }
        links = gmoa.surface_match(gpdA1, gpdB1, params)
        print("links",links)
        assert len(links) == 4
        links = gmoa.surface_match(gpdA1, gpdB2, params)
        print("links",links)
        assert len(links) == 4
        links = gmoa.surface_match(gpdA1, gpdB3, params)
        print("links",links)
        assert len(links) == 3
        links = gmoa.surface_match(gpdA2, gpdB1, params)
        print("links",links)
        assert len(links) == 4
        links = gmoa.surface_match(gpdA2, gpdB2, params)
        print("links",links)
        assert len(links) == 6
        links = gmoa.surface_match(gpdA2, gpdB3, params)
        print("links",links)
        assert len(links) == 3
        links = gmoa.surface_match(gpdA1, gpdB1, params)
        print("links without filtering",links)
        params["final_filtering"] = True
        links = gmoa.surface_match(gpdA1, gpdB1, params)
        print("links",links)
        assert len(links) == 0

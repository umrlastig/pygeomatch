import pytest
from pymatch.geometric_matching_of_areas import get_accuracy, get_completeness, pre_match, get_geom, group_evaluation, search_optimal_groups, surface_match, filter_links
from shapely import Polygon, to_geojson, geometry
import geopandas as gpd

class TestGMA:
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
    
    @pytest.fixture(scope="class")
    def gpd1(self, polygon1):
        return gpd.GeoDataFrame({"id": [0]}, geometry=[polygon1])
    
    @pytest.fixture(scope="class")
    def gpd2(self, polygon1,polygon2, polygon3):
        return gpd.GeoDataFrame({"id": [0,1,2]}, geometry=[polygon1,polygon2,polygon3])

    def test_accuracy(self, polygon1, polygon2, polygon3, polygon4):
        a11 = get_accuracy(polygon1, polygon1)
        assert a11 == pytest.approx(1.0)
        a12 = get_accuracy(polygon1, polygon2)
        assert a12 == pytest.approx(0.4)
        a13 = get_accuracy(polygon1, polygon3)
        assert a13 == pytest.approx(0.0)
        a14 = get_accuracy(polygon1, polygon4)
        assert a14 == pytest.approx(0.1)

    def test_completeness(self, polygon1, polygon2, polygon3, polygon4):
        a11 = get_completeness(polygon1, polygon1)
        assert a11 == pytest.approx(1.0)
        a12 = get_completeness(polygon1, polygon2)
        assert a12 == pytest.approx(1.0)
        a13 = get_completeness(polygon1, polygon3)
        assert a13 == pytest.approx(0.0)
        a14 = get_completeness(polygon1, polygon4)
        assert a14 == pytest.approx(1.0)

    def test_get_geom(self, polygon1, polygon2, polygon3, gpd1, gpd2):
        assert get_geom(0, gpd1) == polygon1
        assert get_geom(0, gpd2) == polygon1
        assert get_geom(1, gpd2) == polygon2
        assert get_geom(2, gpd2) == polygon3

    def test_pre_match(self, gpd1, gpd2):
        links = pre_match(gpd1, gpd2, {"minimise_surface_distance": True, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1})
        for l in links:
            print(l)
            _, id_comp, intersection_ratio, surface_distance = l
            if id_comp == 0:
                assert intersection_ratio == pytest.approx(1.0)
                assert surface_distance == pytest.approx(0.0)
            elif id_comp == 1:
                assert intersection_ratio == pytest.approx(1.0)
                assert surface_distance == pytest.approx(0.6)
        links = pre_match(gpd1, gpd2, {"minimise_surface_distance": False, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1})
        for l in links:
            print(l)
            _, id_comp, intersection_ratio, accuracy, completeness = l
            if id_comp == 0:
                assert intersection_ratio == pytest.approx(1.0)
                assert accuracy == pytest.approx(1.0)
                assert completeness == pytest.approx(1.0)
            elif id_comp == 1:
                assert intersection_ratio == pytest.approx(1.0)
                assert accuracy == pytest.approx(0.4)
                assert completeness == pytest.approx(1.0)
    
    def test_group_evaluation(self, gpd1, gpd2):
        params = {"minimise_surface_distance": True, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1}
        links = pre_match(gpd1, gpd2, params)
        d = group_evaluation(links, gpd1, gpd2, params)
        assert d == pytest.approx(0.0)
        print(d)
        params = {"minimise_surface_distance": False, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1}
        links = pre_match(gpd1, gpd2, params)
        d = group_evaluation(links, gpd1, gpd2, params)
        assert d == pytest.approx(2.0)
        print(d)

    def test_search_optimal_groups(self, gpd1, gpd2):
        params = {"minimise_surface_distance": True, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1, "sure_intersection_percentage": 0.9}
        links = pre_match(gpd1, gpd2, params)
        groups = search_optimal_groups(links, gpd1, gpd2, params)
        print("groups",groups)
        assert len(groups) == 2
        assert groups[0][0] == 0
        assert groups[0][1] == 0
        assert groups[1][0] == 0
        assert groups[1][1] == 1

    def test_surface_match(self, gpd1, gpd2):
        params = {
            "minimise_surface_distance": True, 
            "min_surface_intersection": 0.1, 
            "min_intersection_percentage": 0.1, 
            "sure_intersection_percentage": 0.9, 
            "use_optimal_groups": True, 
            "final_filtering": False,
            "min_surface_distance": 1.0
        }
        links = surface_match(gpd1, gpd2, params)
        print("links",links)
        assert len(links) == 2
        params["final_filtering"] = True
        links = surface_match(gpd1, gpd2, params)
        print("links",links)
        assert len(links) == 0

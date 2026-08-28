import pytest
import pygeomatch.geometric_matching_of_areas as gmoa
import pygeomatch.util as util
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
    def gpdA3(self, polygonA, polygonC):
        return gpd.GeoDataFrame({"id": [0, 1]}, geometry=[polygonA, polygonC])

    @pytest.fixture(scope="class")
    def gpdB1(self, polygon1, polygon2):
        return gpd.GeoDataFrame({"id": [0, 1]}, geometry=[polygon1, polygon2])

    @pytest.fixture(scope="class")
    def gpdB2(self, polygon1, polygon2, polygon3):
        return gpd.GeoDataFrame({"id": [0, 1, 2]}, geometry=[polygon1, polygon2, polygon3])

    @pytest.fixture(scope="class")
    def gpdB3(self, polygon4, polygon5):
        return gpd.GeoDataFrame({"id": [0, 1]}, geometry=[polygon4, polygon5])

    @pytest.fixture(scope="class")
    def gpdB4(self, polygon3, polygon4):
        return gpd.GeoDataFrame({"id": [0, 1]}, geometry=[polygon3, polygon4])

    def test_accuracy(self, polygonA, polygon1, polygon4, polygon5):
        a11 = util.get_accuracy(polygonA, polygonA)
        assert a11 == pytest.approx(1.0)
        a12 = util.get_accuracy(polygonA, polygon1)
        assert a12 == pytest.approx(0.5)
        a13 = util.get_accuracy(polygonA, polygon4)
        assert a13 == pytest.approx(1.0)
        a14 = util.get_accuracy(polygonA, polygon5)
        assert a14 == pytest.approx(0.)

    def test_completeness(self, polygonA, polygon1, polygon4, polygon5):
        a11 = util.get_completeness(polygonA, polygonA)
        assert a11 == pytest.approx(1.0)
        a12 = util.get_completeness(polygon1, polygonA)
        assert a12 == pytest.approx(0.5)
        a13 = util.get_completeness(polygon4, polygonA)
        assert a13 == pytest.approx(1.0)
        a14 = util.get_completeness(polygonA, polygon5)
        assert a14 == pytest.approx(0.)

    def test_get_geom(self, polygonA, polygonB, polygon1, polygon2, gpdA1, gpdB1):
        assert gmoa.get_geom(0, gpdA1) == polygonA
        assert gmoa.get_geom(1, gpdA1) == polygonB
        assert gmoa.get_geom(0, gpdB1) == polygon1
        assert gmoa.get_geom(1, gpdB1) == polygon2

    def test_pre_match(self, gpdA1, gpdB1):
        links = gmoa.pre_match(gpdA1, gpdB1, {
                               "minimise_surface_distance": True, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1})
        print("links", len(links))
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
        links = gmoa.pre_match(gpdA1, gpdB1, {
                               "minimise_surface_distance": False, "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1})
        print("links", len(links))
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
        params = {"minimise_surface_distance": True,
                  "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1}
        links = gmoa.pre_match(gpdA1, gpdB1, params)
        print("pre_match", links)
        d = gmoa.create_complex_link(links, gpdA1, gpdB1)
        print("create_complex_link", d)
        assert d.surface_distance == pytest.approx(0.1)
        params = {"minimise_surface_distance": False,
                  "min_surface_intersection": 0.1, "min_intersection_percentage": 0.1}
        links = gmoa.pre_match(gpdA1, gpdB1, params)
        d = gmoa.create_complex_link(links, gpdA1, gpdB1)
        assert d.accuracy + d.completeness == pytest.approx(1.9)
        print(d)

    def test_search_optimal_groups(self, gpdA1, gpdA2, gpdA3, gpdB1, gpdB2, gpdB3, gpdB4):
        print("A1-B1 minimise_surface_distance=True")
        params = {"minimise_surface_distance": True, "min_surface_intersection": 0.1,
                  "min_intersection_percentage": 0.1, "sure_intersection_percentage": 0.9}
        links = gmoa.pre_match(gpdA1, gpdB1, params)
        assert len(links) == 4
        groups = gmoa.search_optimal_groups(links, gpdA1, gpdB1, params)
        print("search_optimal_groups", len(groups), groups)
        assert len(groups) == 1
        complex_link = groups[0]
        if isinstance(complex_link, gmoa.ComplexGMALink):
            assert len(complex_link.simple_links) == 4
        # bools = [((group.ref == 0) and (group.ref == 0)) or ((group.ref == 0) and (group.ref == 1)) or ((group.ref == 1) and (group.ref == 0)) or ((group.ref == 1) and (group.ref == 1)) for group in groups]
        # assert all(bools)
        print("A1-B3 minimise_surface_distance=True")
        params = {"minimise_surface_distance": True, "min_surface_intersection": 0.01,
                  "min_intersection_percentage": 0.1, "sure_intersection_percentage": 0.7}
        links = gmoa.pre_match(gpdA1, gpdB3, params)
        print(links)
        assert len(links) == 3
        groups = gmoa.search_optimal_groups(links, gpdA1, gpdB3, params)
        print("search_optimal_groups", len(groups), groups)
        assert len(groups) == 1
        complex_link = groups[0]
        if isinstance(complex_link, gmoa.ComplexGMALink):
            assert len(complex_link.simple_links) == 2
        print("A1-B3 minimise_surface_distance=False")
        params = {"minimise_surface_distance": False, "min_surface_intersection": 0.01,
                  "min_intersection_percentage": 0.1, "sure_intersection_percentage": 0.7}
        links = gmoa.pre_match(gpdA1, gpdB3, params)
        print(links)
        assert len(links) == 3
        groups = gmoa.search_optimal_groups(links, gpdA1, gpdB3, params)
        print("search_optimal_groups", len(groups), groups)
        assert len(groups) == 1
        complex_link = groups[0]
        if isinstance(complex_link, gmoa.ComplexGMALink):
            assert len(complex_link.simple_links) == 2
        print("A3-B4 minimise_surface_distance=False")
        links = gmoa.pre_match(gpdA3, gpdB4, params)
        print(links)
        assert len(links) == 2
        groups = gmoa.search_optimal_groups(links, gpdA3, gpdB4, params)
        print("search_optimal_groups", len(groups), groups)
        assert len(groups) == 2

    def test_filter_links(self):
        link1 = gmoa.SimpleGMALink(0, 0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0)
        link2 = gmoa.SimpleGMALink(0, 0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0)
        link3 = gmoa.SimpleGMALink(0, 0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0)
        link4 = gmoa.SimpleGMALink(0, 0, 0.0, 0.0, 0.0, 1.0, 1.0, 0.0)
        params = {
            "minimise_surface_distance": True,
            "max_surface_distance": 0.5
        }
        links = gmoa.filter_links([link1, link2, link3, link4], params)
        assert len(links) == 3
        params = {
            "minimise_surface_distance": False,
            "min_accuracy_completeness": 0.5
        }
        links = gmoa.filter_links([link1, link2, link3, link4], params)
        assert len(links) == 1

    def test_export_link(self):
        links = gmoa.export_link(
            0, gmoa.SimpleGMALink(1, 2, 3., 4., 5., 6., 7., 8.))
        assert len(links) == 1
        ref, comp, group, measures = links[0].as_tuple()
        assert group == 0
        assert ref == 1
        assert comp == 2
        assert measures["intersection_area"] == 3.
        assert measures["intersection_ratio"] == 4.
        assert measures["surface_distance"] == 5.
        assert measures["accuracy"] == 6.
        assert measures["completeness"] == 7.
        assert measures["radial_distance"] == 8.
        simple_links = [gmoa.SimpleGMALink(
            1, 2, 3., 4., 5., 6., 7., 7.5), gmoa.SimpleGMALink(8, 9, 10., 11., 12., 13., 14., 14.5)]
        links = gmoa.export_link(42, gmoa.ComplexGMALink(
            simple_links, 0., 0., 0., 0., 0., 0.))
        assert len(links) == 2
        ref, comp, group, measures = links[0].as_tuple()
        assert group == 42
        assert ref == 1
        assert comp == 2
        assert measures["intersection_area"] == 3.
        assert measures["intersection_ratio"] == 4.
        assert measures["surface_distance"] == 5.
        assert measures["accuracy"] == 6.
        assert measures["completeness"] == 7.
        assert measures["radial_distance"] == 7.5
        ref, comp, group, measures = links[1].as_tuple()
        assert group == 42
        assert ref == 8
        assert comp == 9
        assert measures["intersection_area"] == 10.
        assert measures["intersection_ratio"] == 11.
        assert measures["surface_distance"] == 12.
        assert measures["accuracy"] == 13.
        assert measures["completeness"] == 14.
        assert measures["radial_distance"] == 14.5
        links = gmoa.export_link(42, gmoa.GMALink(0., 0., 0., 0., 0., 0.))
        assert len(links) == 0

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
        print("A1-B1", links)
        assert len(links) == 4
        links = gmoa.surface_match(gpdA1, gpdB2, params)
        print("A1-B2", links)
        assert len(links) == 4
        links = gmoa.surface_match(gpdA1, gpdB3, params)
        print("A1-B3", links)
        assert len(links) == 2
        links = gmoa.surface_match(gpdA2, gpdB1, params)
        print("A2-B1", links)
        assert len(links) == 4
        links = gmoa.surface_match(gpdA2, gpdB2, params)
        print("A2-B2", links)
        assert len(links) == 3
        links = gmoa.surface_match(gpdA2, gpdB3, params)
        print("A2-B3", links)
        assert len(links) == 2
        links = gmoa.surface_match(gpdA1, gpdB1, params)
        print("A1-B1 without filtering", links)
        params["final_filtering"] = True
        links = gmoa.surface_match(gpdA1, gpdB1, params)
        print("A1-B1 with filtering", links)
        assert len(links) == 4

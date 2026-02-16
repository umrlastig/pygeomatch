import shapely
import geopandas
from pygeomatch.util import measures, MatchingLink
from more_itertools import partition
import networkx as nx
import itertools

class GMALink:
    def __init__(self, intersection_area: float, intersection_ratio: float, surface_distance: float, accuracy: float, completeness: float):
        self.intersection_area = intersection_area
        self.intersection_ratio = intersection_ratio
        self.surface_distance = surface_distance
        self.accuracy = accuracy
        self.completeness = completeness

class SimpleGMALink(GMALink):
    def __init__(self, ref: int, comp: int, intersection_area: float, intersection_ratio: float, surface_distance: float, accuracy: float, completeness: float):
        super().__init__(intersection_area, intersection_ratio, surface_distance, accuracy, completeness)
        self.ref = ref
        self.comp = comp
    def __repr__(self):
        return "SimpleGMALink(%s,%s,%s,%s,%s,%s,%s)" % (self.ref, self.comp, self.intersection_area, self.intersection_ratio, self.surface_distance, self.accuracy, self.completeness)

class ComplexGMALink(GMALink):
    def __init__(self, simple_links: list[SimpleGMALink], intersection_area: float, intersection_ratio: float, surface_distance: float, accuracy: float, completeness: float):
        super().__init__(intersection_area, intersection_ratio, surface_distance, accuracy, completeness)
        self.simple_links = simple_links
    def __repr__(self):
        return "ComplexGMALink(%s,%s,%s,%s,%s,%s)" % (self.simple_links, self.intersection_area, self.intersection_ratio, self.surface_distance, self.accuracy, self.completeness)

def export_measures(link: GMALink) -> dict[str, float]:
    return {
        "intersection_area": link.intersection_area,
        "intersection_ratio": link.intersection_ratio,
        "surface_distance": link.surface_distance,
        "accuracy": link.accuracy,
        "completeness": link.completeness,
    }

def export_link(group:int, link: GMALink) -> list[MatchingLink]:
    if isinstance(link, SimpleGMALink):
        return [MatchingLink(link.ref, link.comp, group, export_measures(link))]
    if isinstance(link, ComplexGMALink):
        return [MatchingLink(m.ref, m.comp, group, export_measures(m)) for m in link.simple_links]
    return []

def surface_match(ref, comp, param: dict) -> list[MatchingLink]:
    """
    Match surfaces using the GMoA algorithm.
    For more information, see in particular the PhD thesis of Atef Bel Hadj (2001).
    
    :param ref: ref features
    :param comp: comp features
    :param param: algorithm parameters
    :type param: dict
    """
    links = pre_match(ref, comp , param)
    if param["use_optimal_groups"]:
        links = search_optimal_groups(links, ref, comp, param)
    # ajout petites surfaces 
    #if param["ajoutPetitesSurfaces"] : 
    #    pass
    #    liensRegroupes = ajoutPetitesSurfaces(liensRegroupes, popRef, popComp, param)
    if param["final_filtering"]:
        links = filter_links(links, param)
    return list(itertools.chain.from_iterable([export_link(idx, link) for idx, link in enumerate(links)]))

def pre_match(ref, comp, param ) -> list[SimpleGMALink]:
    """
    Pre-matches features a & b is:
    - intersection(a,b).area > min_surface_intersection
    - intersection(a,b).area / min(a.area, b.area) > min_intersection_percentage
    
    Returns 1-1 matches.
    :param ref: reference features
    :param comp: comparison features
    :param param: parameters for the algorithm
    """
    # The first subarray contains input geometry integer indices. The second subarray contains tree geometry integer indices.
    ref_indices, comp_indices = comp["geometry"].sindex.query(ref["geometry"], predicate="intersects")
    # zip both lists into tuples (refIndex, compIndex)
    zipped = zip(ref_indices, comp_indices)
    def create_link(a, b) -> SimpleGMALink:
        """
        Compute measures surface_distance, accuracy and completeness.
        TODO we might only compute the measures that will actually be used (surface_distance or accuracy and completeness).
        If we do that, we might want to create subclasses of SimpleGMALink
        :param a: index of ref feature
        :param b: index of comp feature
        """
        geom_a = ref.iloc[[a]].iloc[0]["geometry"].buffer(0)
        geom_b = comp.loc[[b]].iloc[0]["geometry"].buffer(0)
        m = measures(geom_a, geom_b)
        return SimpleGMALink(a, b, *m)
    # compure measures on all couples ref, comp
    links = list(map(lambda z: create_link(z[0], z[1]), zipped))
    def filter_links(link: SimpleGMALink):
        return not((link.intersection_area <= param["min_surface_intersection"]) | (link.intersection_ratio < param["min_intersection_percentage"]))
    return list(filter(filter_links, links))

def links_grouped_by_component(links: list[SimpleGMALink]) -> list[list[SimpleGMALink]]:
    """
    Returns a list of lists, where each inner list contains the Link
    objects that belong to the same connected component.
    """
    # ---- Build an undirected (bipartite) graph -----------------
    G = nx.Graph()
    # Add an edge for every Link (NetworkX adds the nodes automatically)
    G.add_edges_from((("ref", l.ref), ("comp", l.comp)) for l in links)
    # ---- Find connected components (as sets of node names) ---
    node_components = list(nx.connected_components(G))
    # ---- Map each node → component index for quick lookup -------
    node_to_comp = {}
    for idx, comp_nodes in enumerate(node_components):
        for n in comp_nodes:
            if (n[0] == "ref"):
                node_to_comp[n] = idx
    # ---- Bucket the original Link objects ----------------------
    groups = [[] for _ in node_components]
    for l in links:
        # Both endpoints are in the same component, so we can look up either.
        comp_idx = node_to_comp[("ref", l.ref)]
        groups[comp_idx].append(l)
    return groups

def search_optimal_groups(pre_match_links: list[SimpleGMALink], ref , comp , param) -> list[GMALink]:
    """
    On recherche les regroupements optimaux de liens de pré-traitement, pour maximiser la distance surfacique entre les groupes de référence et de comparaison.
    NB: l'appariement est symétrique

    Returns the links (liens d'appariement calculés)
    
    :param pre_match_links: links from the pre matching step
    :param ref: reference dataset
    :param comp: comparison dataset
    :param param: algorithm parameters
    """
    # create groups from connected components
    groups = links_grouped_by_component(pre_match_links)
    simple_groups, complex_groups = partition(lambda group: len(group) > 1, groups)
    final_groups:list[GMALink] = []
    for group in simple_groups:
        final_groups.append(group[0])
    for _, group in enumerate(complex_groups):
        # The first yields the items that have intersection <= param["sure_intersection_percentage"]. 
        # The second yields the items that have intersection > param["sure_intersection_percentage"].
        removable, not_removable = partition(lambda link: link.intersection_ratio > param["sure_intersection_percentage"], group)
        # converts to lists
        removable = list(removable)
        not_removable = list(not_removable)
        # print("not removable",len(not_removable),not_removable)
        # generate combinations of removable links (with a length of at most group size - 1)
        candidate_complex_links: list[ComplexGMALink] = []
        for i in range(0, min(len(removable)+1, len(group))):
            for x in itertools.combinations(removable, i):
                removable_subgroup = list(x)
                candidate_group = [link for link in group if link not in removable_subgroup]
                candidate = create_complex_link(candidate_group, ref, comp)
                candidate_complex_links.append(candidate)
                # print("candidate",len(candidate.simple_links),candidate)
        # sorted in ascending order of number of links in order to have the simplest candidates first
        candidate_complex_links = sorted(candidate_complex_links, key=lambda c: len(c.simple_links))
        if param["minimise_surface_distance"]:
            best = candidate_complex_links.index(min(candidate_complex_links, key=lambda k: k.surface_distance))
        else:
            best = candidate_complex_links.index(max(candidate_complex_links, key=lambda k: k.accuracy + k.completeness))
        best_candidate = candidate_complex_links[best]
        # print("best",best_candidate)
        # recompose
        decomposed_groups = links_grouped_by_component(candidate_complex_links[best].simple_links)
        # print("decomposed_groups",decomposed_groups)
        decomposed_links = [create_complex_link(group, ref, comp) for group in decomposed_groups]
        # do we keep the decomposed version?
        if param["minimise_surface_distance"]:
            decomposed = max(decomposed_links, key=lambda k: k.surface_distance)
            prefer_decomposed = (decomposed.surface_distance <= best_candidate.surface_distance)
        else:
            decomposed = min(decomposed_links, key=lambda k: k.accuracy + k.completeness)
            prefer_decomposed = (decomposed.accuracy + decomposed.completeness >= best_candidate.accuracy + best_candidate.completeness)
        # print("prefer_decomposed",prefer_decomposed)
        if prefer_decomposed:
            final_groups.extend(decomposed_links)
        else:
            final_groups.append(best_candidate)
    return final_groups

def get_geom(index, dataframe):
    """
    Returns the geometry for the given index in the input dataframe.
    
    :param index: index of the feature
    :param dataframe: a geodataframe
    """
    return dataframe.iloc[[index]].iloc[0]['geometry']

def list_union(list: list[int], gdf: geopandas.GeoDataFrame):
    """
    Returns the geometric union of the features with the given indices.
    
    :param list: a list of indices
    :param dataframe: a geodataframe
    """
    geom_list = [get_geom(list[k],gdf).buffer(0) for k in range(0, len(list))]
    return shapely.union_all(geom_list)

def create_complex_link(group: list[SimpleGMALink] , ref , comp) -> ComplexGMALink:
    """
    Combine the measures of the group and create a complex link.
    
    :param group: a group of links
    :param ref: reference dataset
    :param comp: comparison dataset
    :param param: algorithm parameters
    """
    # create the list of ref and comp from the input links
    l_ref, l_comp = zip(*[[g.ref, g.comp] for g in group])
    ref_union  = list_union(list(l_ref), ref)
    comp_union = list_union(list(l_comp), comp)
    return ComplexGMALink(group, *measures(ref_union, comp_union))

def filter_links(grouped_links: list[GMALink] | list[SimpleGMALink], param) -> list[GMALink]:
    """
    Filter links.

    :param grouped_links: Description
    :param param: algorithm parameters
    :param ref: reference dataset
    :param comp: comparison dataset
    """
    filtered_links: list[GMALink] = []
    for _, grouped_link in enumerate(grouped_links):
        if param["minimise_surface_distance"] : 
            if grouped_link.surface_distance <= param["max_surface_distance"]:
                filtered_links.append(grouped_link)
        else:
            if (grouped_link.accuracy >= param["min_accuracy_completeness"]) and (grouped_link.completeness >= param["min_accuracy_completeness"]):
                filtered_links.append(grouped_link)
    return filtered_links

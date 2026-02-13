import numpy as np
import shapely
import geopandas
from pygeomatch.util import surface_distance
from more_itertools import partition
import networkx as nx
import itertools

class GMALink:
    def __init__(self, ref: int, comp: int, intersection: float, intersection_ratio: float, surface_distance: float, accuracy: float, completeness: float):
        self.ref = ref
        self.comp = comp
        self.intersection = intersection
        self.intersection_ratio = intersection_ratio
        self.surface_distance = surface_distance
        self.accuracy = accuracy
        self.completeness = completeness
    def __repr__(self):
        return "GMALink(%s,%s,%s,%s,%s,%s,%s)" % (self.ref, self.comp, self.intersection, self.intersection_ratio, self.surface_distance, self.accuracy, self.completeness)

def get_accuracy(geomA, geomB):
    """
    accuracy(A,B) = surface(A inter B) / Surface (A)
    
    :param geomA: geometry A
    :param geomB: geometry B
    """
    inter = shapely.intersection(geomA , geomB)
    if inter.is_empty:
        return 0
    return inter.area / geomA.area
    
def get_completeness(geomA , geomB):
    """
    completeness(A,B) = surface(A inter B) / Surface (B)
    completeness(A,B) = accuracy(B,A)

    :param geomA: geometry A
    :param geomB: geometry B
    """
    return get_accuracy(geomB, geomA)

def surface_match(ref, comp, param: dict) -> list:
    """
    Match surfaces using the GMoA algoriithm.
    For more information, see in particular the PhD thesis of Atef Bel Hadj (2001).
    
    :param ref: ref features
    :param comp: comp features
    :param param: algorithm parameters
    :type param: dict
    """
    #TODO check consistency of all these methods: they don't return the same thing at the moment
    # pre-appariement surfaces
    links = pre_match(ref, comp , param)
    # print("pre_match",links)
    # recherche groupes optimaux 
    if param["use_optimal_groups"]:
        links = search_optimal_groups(links, ref, comp, param)
        # print("optimal_groups",links)
    # ajout petites surfaces 
    #if param["ajoutPetitesSurfaces"] : 
    #    pass
    #    liensRegroupes = ajoutPetitesSurfaces(liensRegroupes, popRef, popComp, param)
    # filtres finales  
    if param["final_filtering"]:
        links = filter_links(links, param, ref , comp)
        # print("final_filtering",links)
    return [[link.ref, link.comp, link.surface_distance] for link in links]

def pre_match(ref, comp, param ) -> list[GMALink]:
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
    def measures(a, b) -> GMALink:
        """
        Compute measures surface_distance, accuracy and completeness.
        TODO we might only compute the measures that will actually be used (surface_distance or accuracy and completeness).
        If we do that, we might want to create subclasses of GMALink
        :param a: index of ref feature
        :param b: index of comp feature
        """
        geom_a = ref.iloc[[a]].iloc[0]["geometry"].buffer(0)
        geom_b = comp.loc[[b]].iloc[0]["geometry"].buffer(0)
        intersection = shapely.intersection(geom_a , geom_b).area
        intersection_ratio = shapely.intersection(geom_a, geom_b).area / (min(geom_a.area , geom_b.area))
        _surface_distance = surface_distance(geom_a, geom_b)
        accuracy = get_accuracy(geom_a, geom_b)
        completeness = get_completeness(geom_a, geom_b)
        return GMALink(a, b, intersection, intersection_ratio, _surface_distance, accuracy, completeness)
    # compure measures on all couples ref, comp
    links = list(map(lambda z: measures(z[0], z[1]), zipped))
    def filter_links(link: GMALink):
        return not((link.intersection <= param["min_surface_intersection"]) | (link.intersection_ratio < param["min_intersection_percentage"]))
    # def map_outputs(link: GMALink):
    #     if param["minimise_surface_distance"]:
    #         return [link[0], link[1], link[3], link[4]] # keep (index_a, index_b, intersection_ratio, surface_distance)
    #     return [link[0], link[1], link[3], link[5], link[6]] # keep (index_a, index_b, intersection_ratio, accuracy, completeness)
    # return list(map(map_outputs, filter(filter_links, links)))
    return list(filter(filter_links, links))

def links_grouped_by_component(links: list[GMALink]) -> list[list[GMALink]]:
    """
    Returns a list of lists, where each inner list contains the Link
    objects that belong to the same connected component.
    """
    # ---- Build an undirected bipartite graph -----------------
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

def search_optimal_groups(pre_match_links: list[GMALink], ref , comp , param) -> list[GMALink]:
    """
    On recherche les regroupements optimaux de liens de pré-traitement, pour maximiser la distance surfacique entre les groupes de référence et de comparaison.
    NB: l'appariement est symétrique

    Returns the links (liens d'appariement calculés)
    
    :param pre_match_links: links from the pre matching step
    :param ref: reference dataset
    :param comp: comparison dataset
    :param param: algorithm parameters
    """
    groups = links_grouped_by_component(pre_match_links)
    # print("groups",len(groups),groups)
    matrix = np.zeros((len(ref), len(comp)))
    for k in range(len(pre_match_links)):
        matrix[pre_match_links[k].ref][pre_match_links[k].comp] = pre_match_links[k].intersection_ratio
    #groups_to_keep: list[list[GMALink]] = []
    #on parcrous touts les liens n-m créés
    # connected_groups = []
    # for i  in range(len(ref)):
    #     groups = [(i, j, matrix[i][j]) for j in range(len(comp)) if matrix[i][j] > 0]
    #     if len(groups) != 0:
    #         connected_groups.append(groups)
    simple_groups, complex_groups = partition(lambda group: len(group) > 1, groups)
    groups_to_keep = list(simple_groups)
    for idx, group in enumerate(complex_groups):
        # print("connected group", idx, len(group))
        # pour tous les objets isolés ou les liens 1-1, on ne fait rien de plus 
        # if len(group) == 1 : 
        #     groups_to_keep.append(group)
        #     continue 
        # pour les groupes n-m on va essayer d'enlever des arcs 
        # mais on garde à coup sûr les liens avec suffisament de recouvrement
        # The first yields the items that have intersection <= param["sure_intersection_percentage"]. 
        # The second yields the items that have intersection > param["sure_intersection_percentage"].
        removable, not_removable = partition(lambda link: link.intersection_ratio > param["sure_intersection_percentage"], group)
        # converts to lists
        removable = list(removable)
        not_removable = list(not_removable)
        # print("not removable",len(not_removable))
        if len(not_removable) == len(group) : # si on ne peut rien enlever on s'arrête la 
            groups_to_keep.append(group)
            continue
        #on cherche à enlever toutes les combinaisons possibles d'arcs virables
        # generate combinations (with a length of at most group size - 1)
        candidate_groups: list[list[GMALink]] = []
        evals: list[float] = []
        for i in range(0, min(len(removable)+1, len(group))):
            for x in itertools.combinations(removable, i):
                removable_subgroup = list(x)
                candidate_group = [link for link in group if link not in removable_subgroup]
                candidate_groups.append(candidate_group)
                evals.append(group_evaluation(candidate_group, ref, comp, param))
                # print(f"removable_subgroup({len(removable_subgroup)})",removable_subgroup)
                # print(f"candidate_group({len(candidate_group)}",candidate_group)
                # print("eval",evals[-1])
        if param["minimise_surface_distance"]:
            best = evals.index(min(evals))
        else:
            best = evals.index(max(evals))
        groups_to_keep.append(candidate_groups[best])
    L = []
    for k in range(len(groups_to_keep)):
        # print("groups_to_keep",k,groups_to_keep[k])
        if len(groups_to_keep[k]) == 2 : 
            L.append(groups_to_keep[k][0])
            L.append(groups_to_keep[k][1])
        elif len(groups_to_keep[k]) == 0 : 
            continue
        else : 
            L.append(groups_to_keep[k][0])
    return L

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
    try:
        geom_list = [get_geom(list[k],gdf).buffer(0) for k in range(0, len(list))]
        return shapely.union_all(geom_list)
    except shapely.errors.GEOSException:
        print("union list error with\n",list)
        return None

# Il me semble qu'on ajoute la zone petite à la zone plus grande pour en 
# faire une nouvelle entité 
def group_evaluation(group: list[GMALink] , ref , comp , param) -> float:
    """
    Combine the measures of the group.
    
    :param group: a group of links
    :param ref: reference dataset
    :param comp: comparison dataset
    :param param: algorithm parameters
    """
    # create the list of ref and comp from the input links
    l_ref, l_comp = zip(*[[g.ref, g.comp] for g in group])
    ref_union  = list_union(list(l_ref), ref)
    comp_union = list_union(list(l_comp), comp)
    # combine the measures of the group
    if param["minimise_surface_distance"]:
        return min(surface_distance(ref_union, comp_union), 2)
    return max(get_accuracy(ref_union , comp_union) + get_completeness(ref_union , comp_union), -1)

def filter_links(grouped_links: list[GMALink], param, ref, comp) -> list[GMALink]:
    """
    Filter links.

    :param grouped_links: Description
    :param param: algorithm parameters
    :param ref: reference dataset
    :param comp: comparison dataset
    """
    filtered_links: list[GMALink] = []
    for _, grouped_link in enumerate(grouped_links):
        # print("link",i,grouped_links[i])
        if param["minimise_surface_distance"] : 
            # distSurf = grouped_links[i][2]
            # distSurf = surface_distance(get_geom(grouped_link.ref,ref), get_geom(grouped_link.comp,comp))
            # print("distSurf",grouped_link.surface_distance)
            if surface_distance <= param["max_surface_distance"]:
                # link.append(grouped_link.ref)
                # link.append(grouped_link.comp)
                # link.append(distSurf)
                filtered_links.append(grouped_link)
        else:
            # accuracy = get_accuracy(get_geom(grouped_link.ref,ref).buffer(0), get_geom(grouped_link.comp,comp).buffer(0) )
            # completeness = get_completeness(get_geom(grouped_link.ref,ref).buffer(0), get_geom(grouped_link.comp,comp).buffer(0) )
            if (grouped_link.accuracy >= param["min_accuracy_completeness"]) and (grouped_link.completeness >= param["min_accuracy_completeness"]):
                # link.append(grouped_links[i].ref)
                # link.append(grouped_links[i].comp)
                # link.append(accuracy + completeness)
                filtered_links.append(grouped_link)
    return filtered_links

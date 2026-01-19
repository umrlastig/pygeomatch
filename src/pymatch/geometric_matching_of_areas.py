import numpy as np
import shapely
import geopandas
from pymatch.util import surface_distance
from more_itertools import partition

"""""""""""""" """ GMA """ """"""""""""""""""
#############################################

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

def surface_match(ref, comp, param: dict):
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
    # recherche groupes optimaux 
    if param["use_optimal_groups"] : 
        links = search_optimal_groups(links, ref, comp, param)
    # ajout petites surfaces 
    #if param["ajoutPetitesSurfaces"] : 
    #    pass
    #    liensRegroupes = ajoutPetitesSurfaces(liensRegroupes, popRef, popComp, param)
    # filtres finales  
    if param["final_filtering"] : 
        links = filter_links(links, param, ref , comp)
    return links

def pre_match(ref, comp, param ):
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
    refIndices, compIndices = comp["geometry"].sindex.query(ref["geometry"], predicate="intersects")
    # zip both lists into tuples (refIndex, compIndex)
    zipped = zip(refIndices, compIndices)
    def measures(a, b):
        """
        Compute measures surface_distance, accuracy and completeness.
        TODO we might only compute the measures that will actually be used (surface_distance or accuracy and completeness).
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
        return (a, b, intersection, intersection_ratio, _surface_distance, accuracy, completeness)
    # compure measures on all couples ref, comp
    links = list(map(lambda z: measures(z[0], z[1]), zipped))
    def filter_links(link):
        return not((link[2] <= param["min_surface_intersection"]) | (link[3] < param["min_intersection_percentage"]))
    def map_outputs(link):
        if param["minimise_surface_distance"]:
            return [link[0], link[1], link[3], link[4]] # keep (index_a, index_b, intersection_ratio, surface_distance)
        return [link[0], link[1], link[3], link[5], link[6]] # keep (index_a, index_b, intersection_ratio, accuracy, completeness)

    return list(map(map_outputs, filter(filter_links, links)))


def search_optimal_groups(pre_match_links, ref , comp , param):
    """
    On recherche les regroupements optimaux de liens de pré-traitement, pour maximiser la distance surfacique entre les groupes de référence et de comparaison .
    NB l'appariement est symétrique 
    Returns the links (liens d'appariement calculés)
    
    :param pre_match_links: links from the pre matching step
    :param ref: reference dataset
    :param comp: comparison dataset
    :param param: algorithm parameters
    """
    matrix = np.zeros((len(ref), len(comp)))
    for k in range(len(pre_match_links)):
        matrix[pre_match_links[k][0]][pre_match_links[k][1]] = pre_match_links[k][2]
    groups_to_keep = []
    #on parcrous touts les liens n-m créés
    connected_groups = []
    for i  in range(len(ref)):
        groups = [(i,j,matrix[i][j]) for j in range(len(comp)) if matrix[i][j] > 0]
        if len(groups) != 0:
            connected_groups.append(groups)
    for i in range(len(connected_groups)):
        # pour tous les objets isolés ou les liens 1-1, on ne fait rien de plus 
        if len(connected_groups[i]) == 1 : 
            groups_to_keep.append(connected_groups[i])
            continue 
        # pour les groupes n-m on va essayer d'enlever des arcs 
        # mais on garde à coup sûr les liens avec suffisament de recouvrement
        # The first yields the items that have intersection <= param["sure_intersection_percentage"]. 
        # The second yields the items that have intersection > param["sure_intersection_percentage"].
        removable, not_removable = partition(lambda x: x[2] > param["sure_intersection_percentage"], connected_groups[i])
        # converts to lists
        removable = list(removable)
        not_removable = list(not_removable)
        if len(not_removable) == len(connected_groups[i]) : # si on ne peut rien enlever on s'arrête la 
            groups_to_keep.append(connected_groups[i])
            continue 
        #on cherche à enlever toutes les combinaisons possibles d'arcs virables
        distSurfMin = 2
        distExacMax = 0
        to_remove = []                
        dist = group_evaluation(removable, ref, comp, param)
        if param["minimise_surface_distance"]:
            if dist < distSurfMin :
                distSurfMin = dist  # gros problemes de logique 
                to_remove.append(removable)
        else :
            if dist > distExacMax : 
                distExacMax = dist 
                to_remove.append(removable)
        if (len(to_remove)) == 0 : 
            continue 
        else : 
            for k in range(len(to_remove[0])):
                connected_groups[i].remove(to_remove[0][k])                    
            groups_to_keep.append(connected_groups[i])
    L = []
    for k in range(len(groups_to_keep)):
        print(k,groups_to_keep[k])
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
def group_evaluation(group , ref , comp , param):
    """
    Combine the measures of the group.
    
    :param group: a group of links
    :param ref: reference dataset
    :param comp: comparison dataset
    :param param: algorithm parameters
    """
    # create the list of ref and comp from the input links
    l_ref, l_comp = zip(*[[g[0],g[1]] for g in group])
    unionRef  = list_union(list(l_ref), ref)
    unionComp = list_union(list(l_comp), comp)
    # combine the measures of the group
    if param["minimise_surface_distance"]:
        return min(surface_distance(unionRef, unionComp), 2)
    return max(get_accuracy(unionRef , unionComp) + get_completeness(unionRef , unionComp), -1)

def filter_links(grouped_links, param, ref, comp):
    """
    Filter links.
    TODO check the modification of the "else" indentation level.

    :param grouped_links: Description
    :param param: algorithm parameters
    :param ref: reference dataset
    :param comp: comparison dataset
    """
    filtered_links = []
    for i in range(len(grouped_links)):
        link = []
        if param["minimise_surface_distance"] : 
            distSurf = grouped_links[i][2]
            print("distSurf",distSurf)
            if distSurf < param["min_surface_distance"]:
                link.append(grouped_links[i][0])
                link.append(grouped_links[i][1])
                link.append(distSurf)
        else:
            accuracy = get_accuracy(get_geom(grouped_links[i][0],ref).buffer(0), get_geom(grouped_links[i][1],comp).buffer(0) )
            completeness = get_completeness(get_geom(grouped_links[i][0],ref).buffer(0), get_geom(grouped_links[i][1],comp).buffer(0) )
            if (accuracy > param["min_accuracy_completeness"]) and (completeness > param["min_accuracy_completeness"]):
                link.append(grouped_links[i][0])
                link.append(grouped_links[i][1])
                link.append(accuracy + completeness)
        if len(link) != 0 :
            filtered_links.append(link)
    return filtered_links 

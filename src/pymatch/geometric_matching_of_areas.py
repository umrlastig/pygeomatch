import numpy as np
import shapely

"""""""""""""" """ GMA """ """"""""""""""""""
#############################################

# -> Distances
def surface_distance(geomA , geomB):
    """
    Surface distance(A,B) = 1 - inter(A,B).area/union(A,B).area
    
    :param geomA: geometry A
    :param geomB: geometry B
    """
    inter = geomA.intersection(geomB)
    union = shapely.union(geomA.buffer(0),geomB.buffer(0))
    return 1 - inter.area / union.area

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


def search_optimal_groups(preAppLiens, ref , comp , param):
    """
    On recherche les regroupements optimaux de liens de pré-traitement, pour maximiser la distance surfacique entre les groupes de référence et de comparaison .
    NB l'appariement est symétrique 
    Returns the links (liens d'appariement calculés)
    
    :param preAppLiens: liens issus du pré-appariement 
    :param ref: Description
    :param comp: Description
    :param param: paramètres de l'appariement 
    """
    
    matrice = np.zeros((len(ref), len(comp)))
    for k in range(len(preAppLiens)):
        matrice[preAppLiens[k][0]][preAppLiens[k][1]] = preAppLiens[k][2]
    
    groupesGardes = []
    
    #on parcrous touts les liens n-m créés
    groupesConnexes = []
    for i  in range(len(ref)):
        groupes = []
        for j in range(len(comp)): 
            if matrice[i][j] > 0 : 
                groupes.append((i,j,matrice[i][j]))
        
        if len(groupes) != 0 : 
            groupesConnexes.append(groupes)
            
    
    for i in range(len(groupesConnexes)):
        # pour tous les objets isolés ou les liens 1-1, on ne fait rien de plus 
        if len(groupesConnexes[i]) == 1 : 
            groupesGardes.append(groupesConnexes[i])
            continue 
        # pour les groupes n-m on va essayer d'enlever des arcs 
        # mais on garde à coup sûr les liens avec suffisament de recouvrement
    
        arcNonEnlevables = []
        arcEnlevables    = []
        
        for j in range(len(groupesConnexes[i])):
            if groupesConnexes[i][j][2] > param["sure_intersection_percentage"]:
                arcNonEnlevables.append(groupesConnexes[i][j])
            else :
                arcEnlevables.append(groupesConnexes[i][j])
                
        if len(arcNonEnlevables) == len(groupesConnexes[i]) : # si on ne peut rien enlever on s'arrête la 
            groupesGardes.append(groupesConnexes[i])
            continue 
        
        
        #on cherche à enlever toutes les combinaisons possibles d'arcs virables
        distSurfMin = 2
        distExacMax = 0
        # combinaisons = arcEnlevables
        arcDuGroupeEnlevesFinal = []
        # comb = 0
        
        #for j in range(len(arcEnlevables)):
        #    dist = mesureEvaluationGroupe(arcEnlevables[j], popRef, popComp, param)
        #    if param["minimise_surface_distance"]:
        #        if dist < distSurfMin :
        #            distSurfMin = dist 
        #            arcDuGroupeEnlevesFinal.append(arcEnlevables[j])
        #    else :
        #        if dist > distExacMax : 
        #            distExacMax = dist 
        #            arcDuGroupeEnlevesFinal.append(arcEnlevables[j])
        
        dist = group_evaluation(arcEnlevables, ref, comp, param)
        if param["minimise_surface_distance"]:
            if dist < distSurfMin :
                distSurfMin = dist  # gros problemes de logique 
                arcDuGroupeEnlevesFinal.append(arcEnlevables)
        else :
            if dist > distExacMax : 
                distExacMax = dist 
                arcDuGroupeEnlevesFinal.append(arcEnlevables)
        
        #groupesPreGardes = []
        #for j in range(len(groupesConnexes[i])):
        #    compteur = 0 
        #    if (len(arcDuGroupeEnlevesFinal)) == 0 : 
        #        continue 
        #    for k in range(len(arcDuGroupeEnlevesFinal[0])):
        #        if groupesConnexes[i][j] == arcDuGroupeEnlevesFinal[0][k]:
        #            compteur = 1
        #    if compteur == 0 :
         #       groupesPreGardes.append(groupesConnexes[i][j])
         #groupesGardes.append(GroupesPRegardes)
        
        if (len(arcDuGroupeEnlevesFinal)) == 0 : 
            continue 
        else : 
            for k in range(len(arcDuGroupeEnlevesFinal[0])):
                groupesConnexes[i].remove(arcDuGroupeEnlevesFinal[0][k])
                    
            groupesGardes.append(groupesConnexes[i])
    
    L = []
    for k in range(len(groupesGardes)):
        
        if len(groupesGardes[k]) == 2 : 
            L.append(groupesGardes[k][0])
            L.append(groupesGardes[k][1])
        elif len(groupesGardes[k]) == 0 : 
            continue
        else : 
            L.append(groupesGardes[k][0])
            
    return L
    
def list_union(liste,pop):
    try:
        list = [get_geom(liste[k],pop).buffer(0) for k in range(0, len(liste))]
        return shapely.union_all(list)
    except shapely.errors.GEOSException:
        print("union list error with\n",list)
        return None
    # union = getGeom(liste[0],pop)
    # for k in range(1, len(liste)):
    #     try:
    #         union = shapely.union(union , getGeom(liste[k],pop))
    #     except shapely.errors.GEOSException:
    #         print("union list error with\n", union, "\n",getGeom(liste[k],pop))

    # return union

def get_geom(indice,pop):
    # return pop.loc[indice,'geometry']
    return pop.iloc[[indice]].iloc[0]['geometry']

# Il me semble qu'on ajoute la zone petite à la zone plus grande pour en 
# faire une nouvelle entité 
def group_evaluation(groupe , popRef , popComp , param):
    if param["minimise_surface_distance"]:
        result = 2 
    else:
        result = -1
    
    listRef = []
    listComp = []

    for j in range(len(groupe)):
        listRef.append(groupe[j][0])
        listComp.append(groupe[j][1])
    unionRef  = list_union(listRef,popRef)  
    unionComp = list_union(listComp,popComp)
        
    #if len(groupe) == 3 : # A changer si on change le nombre dans groupe 
    #    unionRef  = popRef[groupe[0]]["geometry"]
    #    unionComp = popComp[groupe[1]]["geometry"]

    #else :
    #    for j in range(len(groupe)):
    #        print("groupe =", groupe)
    #        listRef.append(groupe[j][0])
    #        listComp.append(groupe[j][1])
    #    unionRef  = unionListe(listRef,popRef)  
    #    unionComp = unionListe(listComp,popComp) 
        
        
    geomRef = unionRef # peut etre que il y a un pb de type 
    geomComp = unionComp

    #on combine les mesures des parties connexes 
    if param["minimise_surface_distance"]:
        value= surface_distance(geomRef, geomComp)
        return min(value, result)
    value = get_accuracy(geomRef , geomComp) + get_completeness(geomRef , geomComp)
    return max(value, result)       

def filter_links(grouped_links, param, ref, comp):
    """
    Filter links.
    TODO check the modification of the "else" indentation level.

    :param grouped_links: Description
    :param param: Description
    :param ref: Description
    :param comp: Description
    """
    filtered_links = []
    for i in range(len(grouped_links)):
        lien = []
        if param["minimise_surface_distance"] : 
            distSurf = grouped_links[i][2]
            if distSurf < param["min_surface_distance"]:
                lien.append(grouped_links[i][0])
                lien.append(grouped_links[i][1])
                lien.append(distSurf)
        else:
            accuracy = get_accuracy(get_geom(grouped_links[i][0],ref).buffer(0), get_geom(grouped_links[i][1],comp).buffer(0) )
            completeness = get_completeness(get_geom(grouped_links[i][0],ref).buffer(0), get_geom(grouped_links[i][1],comp).buffer(0) )
            if accuracy > param["min_accuracy_completeness"] and completeness > param["min_accuracy_completeness"]: 
                lien.append(grouped_links[i][0])
                lien.append(grouped_links[i][1])
                lien.append(accuracy + completeness)
        if len(lien) != 0 :
            filtered_links.append(lien)
    return filtered_links 

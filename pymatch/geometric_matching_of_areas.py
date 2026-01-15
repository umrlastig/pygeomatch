import geopandas as gpd
import numpy as np
import shapely
from shapely import Geometry, area
from itertools import groupby, chain


"""""""""""""" """ GMA """ """"""""""""""""""
#############################################

# -> Distances
def surface_distance(geomA , geomB) :
    inter = geomA.intersection(geomB)
    # if no intersection, return 2
    if inter == None : return 2
    try:
        union = shapely.union(geomA.buffer(0),geomB.buffer(0))
    except shapely.errors.GEOSException:
        return 1
    # if no union, return 1
    if union == None : return 1
    return 1 - inter.area / union.area

# exactitude = surface( A inter B) / Surface (A)
def getExactitude(geomA, geomB):
    inter = shapely.intersection(geomA , geomB)
    if inter.is_empty is True : return 0
    return inter.area / geomA.area
    
# completude = surface( A inter B) / Surface (B) 
def getCompletude(geomA , geomB): 
    return getExactitude(geomB, geomA)

# Appariement entre deux ensemble de surfaces. Processus inspiré de celui
# defini dans la thèse de Atef Bel Hadj (2001)

def appariementSurfaces(ref, comp, param: dict):
    """
    Docstring for appariementSurfaces
    
    :param ref: Description
    :param comp: Description
    :param param: Description
    :type param: dict
    """
    #TODO check consistency of all these methods: they don't return the same thing at the moment
    # pre-appariement surfaces
    liensPreApp = preAppariementsSurfaces_avec_index_spatiale(ref, comp , param)
    
    # recherche groupes optimaux 
    liensRegroupes = liensPreApp    
    if param["regroupementOptimal"] : 
        liensRegroupes = rechercheRegroupementsOptimaux(liensPreApp, ref, comp, param)
    
    # ajout petites surfaces 
    #if param["ajoutPetitesSurfaces"] : 
    #    pass
    #    liensRegroupes = ajoutPetitesSurfaces(liensRegroupes, popRef, popComp, param)
    
    # a = time.time()
    # filtres finales  
    liensFiltres = liensRegroupes    
    if param["filtrageFinal"] : 
        liensFiltres = filtresLiens(liensRegroupes, param, ref , comp)

    # deja considerer dans writeshapefile 
    #liensFiltres = creerGeometriesDesLiens(liensFiltres, param["persistant"])
    
    return liensFiltres

# 2 surfaces sont pré-appariées si :
# 1) l'intersection des surfaces a une taille supèrieure au seuil " surf_min" ET
# 2) l'intersection fait au moins la taille d'une des surfaces multipliée par le paramètre 
# poucentage min 
# NB 1 par construction : chaque lien pointe vers UN SEUL objet de la population
# de référenes et vers un SEUL objet de la population de comparaison 
# popRef  : population des objets de référence 
# popComp : population des objets de comparaison
# param   : paramètres de l'appariement 
# return lien pré-appariement 

def preAppariementsSurfaces_avec_index_spatiale(ref, comp, param ):
    # The first subarray contains input geometry integer indices. The second subarray contains tree geometry integer indices.
    refIndices, compIndices = comp["geometry"].sindex.query(ref["geometry"], predicate="intersects")
    # zip both lists into tuples (refIndex, compIndex)
    zipped = zip(refIndices, compIndices)
    # group by the ref index (z[0]), map the results to dict entries with the ref index (y[0]) and keep only the second element of the grouped tuples (t[1]) to compure measures
    def measures(a, b):
        # geom_a = ref.loc[a,"geometry"].buffer(0)
        # geom_b = comp.loc[b,"geometry"].buffer(0)
        geom_a = ref.iloc[[a]].iloc[0]["geometry"].buffer(0)
        geom_b = comp.loc[[b]].iloc[0]["geometry"].buffer(0)
        intersection = shapely.intersection(geom_a , geom_b).area
        intersection_ratio = shapely.intersection(geom_a, geom_b).area / (min(geom_a.area , geom_b.area))
        _surface_distance = surface_distance(geom_a, geom_b)
        exactitude = getExactitude(geom_a, geom_b)
        completude = getCompletude(geom_a, geom_b)
        return (a, b, intersection, intersection_ratio, _surface_distance, exactitude, completude)

    links = list(map(lambda z: measures(z[0], z[1]), zipped))
    def filter_links(link):
        return not((link[2] <= param["surface_min_intersection"]) | (link[3] < param["pourcentage_min_intersection"]))
    def map_outputs(link):
        if param["minimiseDistanceSurfacique"]:
            return [link[0], link[1], link[3], link[4]] # keep (index_a, index_b, intersection_ratio, surface_distance)
        return [link[0], link[1], link[3], link[5], link[6]] # keep (index_a, index_b, intersection_ratio, exactitude, completude)

    return list(map(map_outputs, filter(filter_links, links)))


def rechercheRegroupementsOptimaux (preAppLiens, ref , comp , param):
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
            if groupesConnexes[i][j][2] > param["pourcentage_intersection_sur"]:
                arcNonEnlevables.append(groupesConnexes[i][j])
            else :
                arcEnlevables.append(groupesConnexes[i][j])
                
        if len(arcNonEnlevables) == len(groupesConnexes[i]) : # si on ne peut rien enlever on s'arrête la 
            groupesGardes.append(groupesConnexes[i])
            continue 
        
        
        #on cherche à enlever toutes les combinaisons possibles d'arcs virables
        distSurfMin = 2
        distExacMax = 0
        combinaisons = arcEnlevables
        arcDuGroupeEnlevesFinal = []
        comb = 0
        
        #for j in range(len(arcEnlevables)):
        #    dist = mesureEvaluationGroupe(arcEnlevables[j], popRef, popComp, param)
        #    if param["minimiseDistanceSurfacique"]:
        #        if dist < distSurfMin :
        #            distSurfMin = dist 
        #            arcDuGroupeEnlevesFinal.append(arcEnlevables[j])
        #    else :
        #        if dist > distExacMax : 
        #            distExacMax = dist 
        #            arcDuGroupeEnlevesFinal.append(arcEnlevables[j])
        
        dist = mesureEvaluationGroupe(arcEnlevables, ref, comp, param)
        if param["minimiseDistanceSurfacique"]:
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
    
def unionListe(liste,pop):
    try:
        list = [getGeom(liste[k],pop).buffer(0) for k in range(0, len(liste))]
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

def getGeom(indice,pop):
    # return pop.loc[indice,'geometry']
    return pop.iloc[[indice]].iloc[0]['geometry']

# Il me semble qu'on ajoute la zone petite à la zone plus grande pour en 
# faire une nouvelle entité 
def mesureEvaluationGroupe(groupe , popRef , popComp , param):
    
    if param["minimiseDistanceSurfacique"]:
        result = 2 
    else : result = -1
    
    # groupe = [(1, 10, 0.9558426479762911), (1, 26, 0.9999474003823016)]
    
    listRef = []
    listComp = []

    for j in range(len(groupe)):
        listRef.append(groupe[j][0])
        listComp.append(groupe[j][1])
    unionRef  = unionListe(listRef,popRef)  
    unionComp = unionListe(listComp,popComp)
        
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
    if param["minimiseDistanceSurfacique"]:
        value= surface_distance(geomRef, geomComp)
        result = min(value, result)
    else : 
        value = getExactitude(geomRef , geomComp) + getCompletude(geomRef , geomComp)
        result = max(value, result)       
    
    return result 
        
def filtresLiens(liensRegroupes, param, popRef , popComp):    
    liensFiltres = []
    for i in range(len(liensRegroupes)):
        lien = []
        if param["minimiseDistanceSurfacique"] : 
            distSurf = liensRegroupes[i][2]
            if distSurf < param["distSurfMaxFinal"]:
                lien.append(liensRegroupes[i][0])
                lien.append(liensRegroupes[i][1])
                lien.append(distSurf)
            else : 
                exactitude = getExactitude(getGeom(liensRegroupes[i][0],popRef).buffer(0), getGeom(liensRegroupes[i][1],popComp).buffer(0) )
                completude = getCompletude(getGeom(liensRegroupes[i][0],popRef).buffer(0), getGeom(liensRegroupes[i][1],popComp).buffer(0) )
                if exactitude > param["completudeExactitudeMinFinal"] and completude > param["completudeExactitudeMinFinal"]: 
                    lien.append(liensRegroupes[i][0])
                    lien.append(liensRegroupes[i][1])
                    lien.append(exactitude + completude)
        if len(lien) != 0 :
            liensFiltres.append(lien)
    return liensFiltres 

import geopandas as gpd
import numpy as np
import shapely

"""""""""""""" """ GMA """ """"""""""""""""""
#############################################
"""""""""""""" """ GMA """ """"""""""""""""""
#############################################
"""""""""""""" """ GMA """ """"""""""""""""""
#############################################




# -> Operateurs
# to remove?
def intersectionRobuste(geomA, geomB, minRes , maxRes):
    inter = shapely.intersection(geomA.buffer(0), geomB.buffer(0))
    return inter
    # if not inter.is_empty:
    #     return inter
    # for i in range(0,10):
    #     seuilDouglas = minRes + i * (maxRes - minRes)/10
    #     Amodif = DouglasPeucker(geomA , seuilDouglas)
    #     Bmodif = DouglasPeucker(geomB , seuilDouglas)
    #     inter = shapely.intersection( Amodif.buffer(0) , Bmodif.buffer(0))
    #     if not inter.is_empty:
    #         return inter
    # return None 
            
# methode Douglas-Peucker sur un polygon cf B Xiong et alt 2016
# already implemented in shapely so...
def DouglasPeucker(geom, seuil):
    return geom.simplify(seuil, preserve_topology=True)
    
    # geom_mapped = shapely.geometry.mapping(geom)
    # geom = geom_mapped['coordinates'][0]
    
    
    # #recherche de la diagonale la plus grande par rapport au point initiale
    # dmax = 0 
    # A = geom[0]
    # noeudFin = 1
    # for j in range(1, len(geom)): 
    #     B = geom[j]
    #     dist = ((A[0] - B[0])**2 + (A[1] - B[1])**2)**.5
    #     if dist > dmax:
    #         noeudFin = j
    #         dmax = dist 
    
    # #creation des listes
    # listA, listB = [], []
    # for i in range(len(geom)):
    #     if i<=noeudFin : listA.append(geom[i])
    #     else :           listB.append(geom[i])
            
    # #DouglasPeucker sur une ligne 
    # ligneA = DouglasPeuckerLigne(listA, seuil)
    # ligneB = DouglasPeuckerLigne(listB, seuil)
    
    # ligne = ligneA + ligneB

    # if len(ligne) <3 : 
    #     return shapely.geometry.Polygon(geom)
    # geom = shapely.geometry.Polygon(ligne)
    # return geom
    
def DouglasPeuckerLigne(geom, seuil): 
    print(geom)
    listFinal = [geom[0]]
    A = geom[0]
    C = geom[len(geom)-1]
    for j in range(1, len(geom)): 
        B = geom[j]
        if A[0] == B[0] : 
            continue
        dist = projete(A, B, C)
        if dist > seuil :
            listFinal.append(B)
            A = B
        
    return listFinal

def projete(A, B , C) : 
    a = ( B[1] - A[1] )/ ( B[0] - A[0] )
    #b = A[1] - a*A[0]

    # calcul du point othoganl de objet_interre sur la droite y :-> ax + b
    vect = ( B[0] - A[0] , B[1] - A[1])
    AH = ( ( C[0] - A[0] )*vect[0] + ( C[1] - A[1] )*vect[1] ) / ((vect[0]**2 + vect[1]**2 )**.5)
    xH = A[0] + (AH*vect[0])/((vect[0]**2 + vect[1]**2 )**.5)
    yH = A[1] + (AH*vect[1])/((vect[0]**2 + vect[1]**2 )**.5)
    
    return ((C[0] - xH)**2 + (C[1] - yH)**2)**.5
        

# -> Distances
def distanceSurfaciqueRobuste(geomA , geomB, minRes, maxRes ) :
    inter = intersectionRobuste(geomA, geomB, minRes, maxRes)
    # en cas de problème d'intersection avec JTS, la methode retourne 2 
    if inter == None : return 2
    try:
        union = shapely.union(geomA.buffer(0),geomB.buffer(0))
    except shapely.errors.GEOSException:
        print("union error with\n",geomA,"\n",geomB)
        return 1
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

def appariementSurfaces(popRef , popComp, param):
    # pre-appariement surfaces
    liensPreApp = preAppariementsSurfaces_avec_index_spatiale(popRef , popComp , param )
    
    # recherche groupes optimaux 
    liensRegroupes = liensPreApp    
    if param["regroupementOptimal"] : 
        liensRegroupes = rechercheRegroupementsOptimaux(liensPreApp, popRef, popComp, param)
    
    # ajout petites surfaces 
    #if param["ajoutPetitesSurfaces"] : 
    #    pass
    #    liensRegroupes = ajoutPetitesSurfaces(liensRegroupes, popRef, popComp, param)
    
    # a = time.time()
    # filtres finales  
    liensFiltres = liensRegroupes    
    if param["filtrageFinal"] : 
        liensFiltres = filtresLiens(liensRegroupes, param , popRef , popComp)

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

def preAppariementsSurfaces (popRef , popComp , param ):
    preAppLiens = [] 
    
    for i in range(len(popRef)):
        #geomRef = popRef[i].getGeom() 
        geomRef = popRef[i]['geometry']
        
        # test d'association sur tous les objets comp intersectant l'objet ref 
        for j in range(len(popComp)):
            #geomComp = popComp[j].getGeom()
            geomComp = popComp[j]['geometry']
            
            
            # creation eventuelle d'un nouveau lien de pré-appariement
            inter = intersectionRobuste(geomRef , geomComp , param["resolutionMin"],
                                        param["resolutionMax"])
            if (inter == None):
         
                continue 
            surfaceIntersection = inter.area
            
            if surfaceIntersection <= param["surface_min_intersection"]:
           
                continue 
            
            pourcentageRecouvrement = max(surfaceIntersection/ geomRef.area,
                                               surfaceIntersection/ geomComp.area)
            if pourcentageRecouvrement < param["pourcentage_min_intersection"]:
         
                continue #intersection pas suffisante 
                
            Lien = []
            #Lien.append(popRef[i]) #peut etre a revoir 
            #Lien.append(popComp[j])
            """ passe par l'indice spatial """
            
            
            Lien.append(popRef [i]["id_spatial"])
            Lien.append(popComp[j]["id_spatial"])
            Lien.append(pourcentageRecouvrement)
       
            if param["minimiseDistanceSurfacique"]:
                Lien.append(distanceSurfaciqueRobuste(geomRef , geomComp , param["resolutionMin"],
                                            param["resolutionMax"]) ) # getDistanceSurfacique
                
            else : 
                Lien.append(getExactitude(geomRef , geomComp)) #getExactitude()
                Lien.append(getCompletude(geomRef , geomComp)) #getCompletude()
            
            preAppLiens.append(Lien)
    return preAppLiens

def preAppariementsSurfaces_avec_index_spatiale (popRef , popComp , param ):

    preAppLiens = [] 
    
    geomRef = [ popRef[i]["geometry"] for i in range(len(popRef)) ]
    geomCom = [ popComp[i]["geometry"] for i in range(len(popComp)) ]
    ref = gpd.GeoSeries( geomRef )
    com = gpd.GeoSeries( geomCom )
    
    # inter = intersectionRobuste_index(geomRef , geomComp , param["resolutionMin"],
    #                            param["resolutionMax"])
    inter = com.sindex.query(ref, predicate="intersects")
    lien = [(inter[0][i] , 
             inter[1][i] , 
             shapely.intersection(geomRef[inter[0][i]].buffer(0) , geomCom[inter[1][i]].buffer(0)).area,
             shapely.intersection(geomRef[inter[0][i]].buffer(0) , geomCom[inter[1][i]].buffer(0)).area / (min(geomRef[inter[0][i]].buffer(0).area , geomCom[inter[1][i]].buffer(0).area)),
             distanceSurfaciqueRobuste(geomRef[inter[0][i]].buffer(0) , geomCom[inter[1][i]].buffer(0) , param["resolutionMin"],
                                         param["resolutionMax"]) ,
             getExactitude(geomRef[inter[0][i]].buffer(0) , geomCom[inter[1][i]].buffer(0)),
             getCompletude(geomRef[inter[0][i]].buffer(0) , geomCom[inter[1][i]].buffer(0))
             ) for i in range(len(inter[0]))]
    
    # lien = [ idRef , idCom , inter.area , pourcentageRecouvrement , distSurfaciqueRobuste,  exactitude , completude]
    # verifier si popRef[inter[0][i]]["geometry"] = geomRef[inter[0][i]]
    
    
    for i in range(len(lien)):
        lienFiltre = []
        if lien[i][2] <= param["surface_min_intersection"]:
            continue 
        if lien[i][3] < param["pourcentage_min_intersection"]:
            continue
        if param["minimiseDistanceSurfacique"]:
            lienFiltre = ( lien[i][0] , lien[i][1] , lien[i][3] , lien[i][4] )
        else : 
            lienFiltre = ( lien[i][0] , lien[i][1] , lien[i][3] ,lien[i][5] , lien[i][6] ) 
        preAppLiens.append(lienFiltre)
        
    return preAppLiens

# On recherche les regroupements optimaux de liens de pré-traitement, pour
# maximiser la distance surfacique entre les groupes de référence et de 
# comparaison 
# NB l'appariement est symétrique 
# param   : paramètres de l'appariement 
# preAppLiens : liens issus du pré-appariement 
# return liens d'appariement calculés 


def rechercheRegroupementsOptimaux (preAppLiens, popRef , popComp , param):
    
    matrice = np.zeros((len(popRef), len(popComp)))
    for k in range(len(preAppLiens)):
        matrice[preAppLiens[k][0]][preAppLiens[k][1]] = preAppLiens[k][2]
    
    groupesGardes = []
    
    #on parcrous touts les liens n-m créés
    groupesConnexes = []
    for i  in range(len(popRef)):
        groupes = []
        for j in range(len(popComp)): 
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
        
        dist = mesureEvaluationGroupe(arcEnlevables, popRef, popComp, param)
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
    return pop[indice]['geometry']

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
        value= distanceSurfaciqueRobuste(geomRef, geomComp, param["resolutionMin"], param["resolutionMax"])
        result = min(value, result)
    else : 
        value = getExactitude(geomRef , geomComp) + getCompletude(geomRef , geomComp)
        result = max(value, result)       
    
    return result 
        
def filtresLiens(liensRegroupes, param , popRef , popComp):
    # liensRegroupes = [ [ idRef, idCom, dsi]]
    
    liensFiltres = []
    
    for i in range(len(liensRegroupes)):
        lien = []
        if param["minimiseDistanceSurfacique"] : 
            distSurf = liensRegroupes[i][2]
            if distSurf < param["distSurfMaxFinal"]:
                lien.append(liensRegroupes[i][0])
                lien.append(liensRegroupes[i][1])
                lien.append(distSurf)
                #liens.append(getArcs)
            else : 
                exactitude = getExactitude(popRef[liensRegroupes[i][0]]["geometry"].buffer(0) , popComp[liensRegroupes[i][1]]["geometry"].buffer(0) )
                completude = getCompletude(popRef[liensRegroupes[i][0]]["geometry"].buffer(0) , popComp[liensRegroupes[i][1]]["geometry"].buffer(0) )
                #exactitude = liensRegroupes[2]
                #completude = liensRegroupes[3]
                if exactitude > param["completudeExactitudeMinFinal"] and completude > param["completudeExactitudeMinFinal"]: 
                    lien.append(liensRegroupes[i][0])
                    lien.append(liensRegroupes[i][1])
                    lien.append(exactitude + completude)
                    #liens.append(getArcs)
        if len(lien) != 0 :
            liensFiltres.append(lien)
    return liensFiltres 

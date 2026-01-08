#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import shapely
from shapely.geometry import Point, LineString
import geopandas as gpd
import pandas as pd
# remove to prevent confusion with the function here
# from Operateurs import intersectionRobuste  
import numpy as np
from shapely.geometry import Point, LineString
import time
import datetime

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
    
    a = time.time()
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


# early attempt
def readShapefile (url, id_column, generatedIds):
    data  = gpd.read_file(url)
    columns = [data.columns[i] for i in range(len(data.columns)) ]
    popRef = []
    for i in range(len(data)):
        L = {}
        for j in range(len(columns)):
            L[columns[j]] = data[columns[j]][i]
        L['id_spatial'] = i
        if generatedIds:
            L[id_column] = i
        print(L)
        popRef.append(L)
    return popRef
 

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
    
def returnID(popRef, popComp, liste):
    
    L = []
    for i in range(len(liste)):
        L.append((popRef[liste[i][0]]['ID'] , popComp[liste[i][1]]['ID']))
    
    return L



"""""""""""""" """ MCA """ """"""""""""""""""
#############################################
"""""""""""""" """ MCA """ """"""""""""""""""
#############################################
"""""""""""""" """ MCA """ """"""""""""""""""
#############################################

#-> Decision
def dempster(liste_critere):
    
    # ajout candidat non app
    
    masseCan = []
    
    
    for i in range(len(liste_critere)):
        masseCan.append(fusionCritereConj(liste_critere[i]))     
        
    # fusion candidat 
    ( resultList , fusCandidat ) = fusionCandidat(masseCan)
    
    # decision 
    deci = decision(resultList , fusCandidat)
    
    return deci 
    

def fusionCritereConj(Candidat) : 
    
    dic = {}
    result = {}
    resultList = []
    transit = {}
    
    
    if len(Candidat) == 1  : 
        Candidat = Candidat[0]
        result["app"] = Candidat[0]
        result["-app"] = Candidat[1]
        result["theta"] = Candidat[2]
        result["phi"] = 0
        return result
        
        
    
    for i in range(len(Candidat)) : 
        dic["app"+str(i+1)] = Candidat[i][0]
        dic["-app"+str(i+1)] = Candidat[i][1]
        dic["theta"+ str(i+1)] = Candidat[i][2]
    result["app"] = 0
    result["-app"] = 0
    result["theta"] = 0
    result["phi"] = 0
    transit["app"] = 0
    transit["-app"] = 0
    transit["theta"] = 0
    transit["phi"] = 0
    
    matrice = [["app" , "phi" , "app"],
               ["phi" , "-app" , "-app"],
               ["app" , "-app" , "theta"]]
    
    for k in range(len(Candidat)-1):
        
        result["app"] = 0
        result["-app"] = 0
        result["theta"] = 0
        result["phi"] = 0
        
        if k == 0 : 
            
            colone = ["app1" , "-app1" , "theta1" ]
            
            ligne = ["app2" , "-app2" , "theta2" ]
            
            matrice = [["app" , "phi" , "app"],
                       ["phi" , "-app" , "-app"],
                       ["app" , "-app" , "theta"]]
        
            for i in range(np.shape(matrice)[0]):
                for j in range(np.shape(matrice)[1]):
                    jeu = matrice[i][j]
                    if dic[colone[j]] * dic[ligne[i]] != 0 :
                        result[jeu] += dic[colone[j]] * dic[ligne[i]]
        
        else : 
            
            ligne = ["app"+str(k+2) , "-app"+str(k+2) , "theta"+str(k+2) ]
            
            colone = ["app" , "-app" , "theta" , "phi" ]
            
            matrice = [["app" , "phi" , "app"],
                       ["phi" , "-app" , "-app"],
                       ["app" , "-app" , "theta"],
                       ["phi" , "phi" , "phi"]]
            
            
            for i in range(np.shape(matrice)[0]):
                for j in range(np.shape(matrice)[1]):
                    
                    jeu = matrice[i][j]
                    
                    
                    if colone[i][0] == 't':
                        result[jeu] += transit['theta'] * dic[ligne[j]]
                    elif colone[i][0] == 'p':
                        result[jeu] += transit['phi'] * dic[ligne[j]]
                    elif colone[i][0] == 'a':
                        result[jeu] += transit['app'] * dic[ligne[j]]
                    elif colone[i][0] == '-':
                        result[jeu] += transit['-app'] * dic[ligne[j]]
                    
                        
        
        transit["app"] = result["app"] 
        transit["-app"] = result["-app"] 
        transit["theta"] = result["theta"]
        transit["phi"] = result["phi"]
        
    return result
        


def fusionCandidat(candidat):
    
    # Changer pa rapport au fait que masseCan est un dictionnaire 
    
    # intialisation des dictionnaires
    
    
    if len(candidat) > 9 : 
        result = {}
        resultList = []
        result["C1"] = 0
        result["-C1"] = 0
        result["theta"] = 0
        result["phi"] = 0
        resultList.append("C1")
        resultList.append("-C1")
        resultList.append("theta")
        resultList.append("phi")
        return (resultList , result)
    

    
    
    dic = {}
    transit = {}
    result = {}
    resultList = []
    for i in range(len(candidat)) : 
        dic["C"+str(i+1)] = candidat[i]["app"]
        dic["-C"+str(i+1)] = candidat[i]["-app"]
        dic["theta"+ str(i+1)] = candidat[i]["theta"]
        dic["phi"+ str(i+1)] = candidat[i]["phi"]
        transit["C"+str(i+1)] = 0
        transit["-C"+str(i+1)] = 0
        result["C"+str(i+1)] = 0
        result["-C"+str(i+1)] = 0
        resultList.append("C"+str(i+1))
        resultList.append("-C"+str(i+1))
    transit["theta"] = 0
    transit["phi"] = 0
    result["theta"] = 0
    result["phi"] = 0
    #result["NA"] = 0
    resultList.append("theta")
    resultList.append("phi")
    #resultList.append("NA")
    
    
    if len(candidat) == 1 : 
        
        result["C1"] = candidat[i]["app"]
        result["-C1"] = candidat[i]["-app"]
        result["theta"] = candidat[i]["theta"]
        result["phi"] = candidat[i]["phi"]
    
    for k in range(len(candidat)-1):
        
        (ligne , colone , matrice) =  mat( k + 2 , len(candidat))
        
        listmodif = []
        
        
        for i in range(len(resultList)):
            result[resultList[i]] = 0
        
        
        for i in range(np.shape(matrice)[0]):
            for j in range(np.shape(matrice)[1]):
                jeu = str(matrice[i][j])[2:len(str(matrice[i][j]))-1]
                listmodif.append(jeu)
                #if dic[colone[j]] * dic[ligne[i]] != 0 :
                if k == 0 : 
                    coeff = dic[colone[j]] * dic[ligne[i]]
                else : 
                    if ligne[i][0] == 't':
                        coeff = transit['theta'] * dic[colone[j]]
                    elif ligne[i][0] == 'p':
                        coeff = transit['phi'] * dic[colone[j]]
                    else :
                        coeff = transit[ligne[i]] * dic[colone[j]]
                             
                # ajout NA
                if jeu == "NA":
                    lis = [k+1 for k in range(len(candidat))]
                    a = int(colone[j][::-1][0])
                    b = int(ligne[i][::-1][0])
                    reste= []
                    for n in range(len(lis)):
                        if lis[n]  in (a,b) : 
                            pass
                        else :
                            reste.append(lis[n])
                    
                    
                    if reste == [] : 
                        result[jeu] = coeff
                    else : 
                        string = 'NA'
                        string += str(reste) 
                        result[string] = coeff
                    
                else : 
                   result[jeu] += coeff


        for i in range(len(listmodif)):
            if listmodif[i] != "phi" and listmodif[i] != "theta" and listmodif[i] != "NA":
                transit[listmodif[i]] = result[listmodif[i]] 
        transit["theta"] = result["theta"] 
        transit["phi"] = result["phi"] 
        
    
    
    a = result['phi']
    if a != 1:
        for i in range(len(result)):
            result[list(result.keys())[i]] = result[list(result.keys())[i]] / ( 1 - a )
      
    result['phi'] = 0
    
    """
    sum_ = 0
    for i in range(len(result)):
        string = str(list(result.keys())[i])
        sum_ += result[string]
    """
    
    return ( resultList , result )
        

def mat( index , nbCandidat):
    
    colone = ("C"+ str(index) , "-C"+ str(index)  , "theta"+ str(index) , "phi"+ str(index))
    
    ligne = []
    
    for i in range(index-1):
        
        ligne.append("C"+str(i+1))
        ligne.append("-C"+str(i+1))
    
            
    ligne.append("theta"+ str(1))
    ligne.append("phi"+ str(1))
        
    matrice = np.chararray((len(ligne),len(colone)), 5)
    for i in range(len(colone)):
        for j in range(len(ligne)):
            matrice[j][i] = loi(ligne[j] , colone[i])
            
    return (ligne , colone , matrice)




def loi(x,y):
    
    # x = Cx , -Cx , theta , phi 
    
    if x == y : 
        return x
    
    if len(x) > 2 :
        if x[0:3] == "phi" :
                return "phi"
        
        if x[0:5] == "theta":
            if len(y) >2 and y[0:3] == "phi" : 
                return "phi"
            else : 
                return y
        
        if x[0] == 'C' :
            if y[0] == 'C' : 
                return "phi"
            if y[0] == '-' :
                return x
        
    
    if len(y) > 2:
        if y[0:3] == "phi" : 
            return "phi"
        if y[0:5] == "theta":
            return x
        
        if y[0] == 'C' :
            if x[0] == 'C' : 
                return "phi"
            if x[0] == '-' :
                return y
                

    
    if x[0] == "C" : 
        
        if y[0] == "C" :    
            if x[1] != y[1] : 
                return "phi"
            
        if y[0] == "-" : 
            
            if x[1] != y[2] : 
                return x
            elif len(x) == len(y) - 2 :
                return x
            else : 
                return "NA"
        
        
    if x[0] == "-" : 
        
        if y[0] == "C" :    
            if x[2] != y[1] : 
                return y
            else : 
                return "NA"
            
        if y[0] == "-" : 
            return "NA"
    

import numpy
def decision(resultList , fusion):
    
    proba = {}
    
    
        
    a = int((len(resultList) - 2 )/ 2)
    
    if a == 1 :
        if fusion["C1"] > fusion["-C1"] or fusion["C1"] > fusion["theta"] :
            return "C1"
        else : return "NA"
    
    app  = []
    _app = []
    
    for i in range(a):
        app.append(fusion['C' + str(i+1)])
        _app.append(fusion['-C' + str(i+1)])
        
    result  = np.array(app) @ numpy.eye(a)
    result += np.array(_app) @ (1/(a) * ( numpy.ones(a) - numpy.eye(a) ))
    
    
    
    ligneNA = []
    matrice = []
    
    
    for i in range(len(fusion)):
        if list(fusion.keys())[i][0:2] == 'NA': 
            ligneNA.append(fusion[list(fusion.keys())[i]])
            ligne = [0 for i in range(a)]
            for j in range(int( (len(list(fusion.keys())[i]) ) /3  ) ):
                """
                if 3*j + 3 == 24 : 
                    chiffre = 10
                    ligne[int(chiffre) - 1] = 1/(a - 1)
                """
                chiffre = list(fusion.keys())[i][3*j +3]
                ligne[int(chiffre) - 1] = 1/(a - 1)
                
            
            matrice.append(ligne)
    
    
                
    result +=   np.array(ligneNA) @ matrice
    
    
    resultNA = np.array(ligneNA) @ ( 1/(a - 1) * np.ones(len(ligneNA)))
    resultNA += np.array(_app) @ (1/(a) * np.ones(a))

    
    
    for i in range(a):
        proba['C' + str(i+1)] = result[i]
        proba['C' + str(i+1)] += fusion['theta']/(a+1)
    proba['NA'] = resultNA
    proba['NA'] += fusion['theta']/(a+1)
    
    
    sum_ = 0
    for i in range(len(proba)):
        string = str(list(proba.keys())[i])
        sum_ += proba[string]
    
    max_ = 0
    for i in range(len(list(proba.keys()))):
        if proba[str(list(proba.keys())[i])] > max_ : 
            max_ = proba[str(list(proba.keys())[i])]
            result = str(list(proba.keys())[i])
        
    
    if result == 'NA' : 
        return 'NA'
    
    epsilon = 0.01
    if proba['NA'] - epsilon < max_ and max_ < proba['NA'] + epsilon :
        return "NA"
    
    
    return result

# added id column names
def MCA(popRef, popComp, id_ref, id_comp):
    
    listeCandidat = SelectionCandidatInter(popRef, popComp, id_ref, id_comp)
    
    listPopRef = listeCandidat[0]
    listPopComp = listeCandidat[1]
    
    App = []
    
    for i in range(len(listPopRef)):
        
        ( listCritere_i, listCrit , liste ) = doAppariement(listPopRef[i], listPopComp[i])
        
        
        if liste == "NA":
            continue
            
        elif liste == "theta":
            continue
        
        else : 
            decision = int(liste[1]) - 1
            App.append((listPopRef[i][0] , listPopComp[i][decision][0] , listCritere_i[decision]))
        
    return App

def doAppariement(listPopRef, listPopComp):

    listCritere = []
    geomRef = listPopRef[1]
    
    listCritere_i = []
    
    for i in range(len(listPopComp)):
        
        list_i = []
        
        
        geomComp = listPopComp[i][1]
        
        # Critere surfacique
        geomRef = geomRef.buffer(0)
        geomComp = geomComp.buffer(0)
        inter = shapely.intersection(geomRef , geomComp)
        union = shapely.union(geomRef,geomComp)
        
        ds =  1 - inter.area /union.area 

        distance = ds
        
        tableau = []
        
        T1 = 0.90
        T2 = 1.0
        E = 0.01
        S = 0.7
        K = 1 - E - S
        if distance < T1 : 
            app = (-(1 - E)/T2) * distance + 1 - E
            _app = E
            tableau.append(app)
            tableau.append(_app)
            tableau.append(1 - app - _app)
        elif distance < T2 : 
            app = (-(1 - E)/T2) * distance + 1 - E
            _app = ( K - E )* distance/ (T2 - T1)  + E - ( K - E )*T1/(T2 - T1) 
            tableau.append(app)
            tableau.append(_app)
            tableau.append(1 - app - _app)
        else :
            app = E
            _app = K
            tableau.append(E)
            tableau.append(_app)
            tableau.append(1 - app - _app)
            
            
        cs = tableau
        listCritere_i.append(ds)
        list_i.append(cs)
        list_i.append(cs)
        
        listCritere . append(list_i)
    
    
    lres = dempster(listCritere)
    
    return (listCritere_i, listCritere , lres )

# added id column names
def SelectionCandidatInter(popRef, popComp, id_ref, id_comp):
    
    idRef   = [ popRef[i][id_ref] for i in range(len(popRef)) ]
    idComp   = [ popComp[i][id_comp] for i in range(len(popComp)) ]
    geomCom = [ popComp[i]["geometry"] for i in range(len(popComp)) ]
    geomRef = [ popRef[i]["geometry"] for i in range(len(popRef)) ]
    
    ref = gpd.GeoSeries( geomRef )
    com = gpd.GeoSeries( geomCom )
    
    inter = com.sindex.query(ref, predicate="intersects")
    
    # creation liste popRef 
    listeRef = [(idRef[inter[0][0]] , geomRef[inter[0][0]])]
    for i in range(len(inter[0])-1):
        if inter[0][i] == inter[0][i+1]:
            continue 
        else :
            listeRef.append((idRef[inter[0][i+1]] , geomRef[inter[0][i+1]]))
    
    # creation liste popComp 
    listeComp   = []
    listeComp_i = [(idComp[inter[1][0]] , geomCom[inter[1][0]])]
    
    for i in range(len(inter[0])-1):
        if inter[0][i] == inter[0][i+1]:
            listeComp_i . append((idComp[inter[1][i+1]] , geomCom[inter[1][i+1]]))
        else :
            listeComp . append(listeComp_i)
            listeComp_i = [(idComp[inter[1][i+1]] , geomCom[inter[1][i+1]])]
    listeComp . append(listeComp_i)
    
    
    return (listeRef , listeComp)
    
    
"""""""""""""" """ COMMUN """ """"""""""""""""""
#############################################
"""""""""""""" """ COMMUN """ """"""""""""""""""
#############################################
"""""""""""""" """ COMMUN """ """"""""""""""""""
#############################################

def writeShapefile(popRef, popComp , data , url ) :
    
    geomRef = []
    for i in range(len(popRef)):
        for j in range(len(data)):
            if data[j][0] == popRef[i]['ID']: 
                if data[j][1] != 0 :  
                    if data[j][1] != 1 :  
                        geomRef.append(popRef[i]["geometry"])
    
    geomCom = [ popComp[i]["geometry"] for i in range(len(popComp)) ]
    #geomRef = [ popRef[i]["geometry"] for i in range(len(popRef)) ]
    
    ref = gpd.GeoSeries( geomRef )
    com = gpd.GeoSeries( geomCom )
    
    shape = []
    for i in range(len(data)):
        lien = {}
        
        if data[i][1] == 1 : 
            #construction
            
            for j in range(len(popComp)):
                if data[i][0] == popComp[j]['ID']: 
                    sourcePosition = popComp[j]["geometry"]
            sourceID = None
            targetID = data[i][0]
            
            compseul = sourcePosition
            inter = ref.sindex.query(compseul, predicate="intersects")
            
            if len(inter) == 0 : 
                nature = 'destruction'
            
            else : 
                nature = 'rien'
            """
            geomRef = compseul.buffer(0)
            geomComp = inter[0].buffer(0)
            inter = shapely.intersection(geomRef , geomComp)
            union = shapely.union(geomRef,geomComp)
            
            ds =  1 - inter.area /union.area 
            
            if ds < 0.5 : 
            
                nature = 'construction'
            """
            
            #if len(inter) == 0 : 
            #    nature2 = 'construction'
            #else : 
            #    nature2 = 'expansion'
            
        else : 
            #appariement ou deconstruction
            for j in range(len(popRef)):
                if data[i][0] == popRef[j]['ID']: 
                    sourcePosition = popRef[j]["geometry"]
            if data[i][1] == 0 : 
                targetID = None
                nature = 'construction'
                
                #refseul = sourcePosition
                #inter = ref.sindex.query(refseul, predicate="intersects")
                
                #if len(inter) == 0 : 
                #    nature2 = 'destruction'
                #else : 
                #    nature2 = 'reduction'
                    
            else : 
                targetID = data[i][1]
                nature = 'appariement'
                #nature2 = 'appariement'
            sourceID = data[i][0]
            
        ds =  data[i][2]
        lien["sourceID"] = sourceID
        lien["targetID"] = targetID
        lien["nature"]   = nature
        #lien["nature2"]   = nature2
        lien["geometry"] = sourcePosition
        lien["ds"] = ds
        shape.append(lien)
    
    listesourceID=[]
    listetargetID=[]
    listegeometry=[]
    listInfo = []
    #listInfo2 = []
    listds = []
    
    for k in range(len(shape)):
        listesourceID.append(shape[k]["sourceID"])
        listetargetID.append(shape[k]["targetID"])
        listegeometry.append(shape[k]["geometry"])
        listInfo     .append(shape[k]["nature"])
        #listInfo2    .append(shape[k]["nature2"])
        listds .append(shape[k]["ds"])

    df = pd.DataFrame(
    {
        "sourceID": listesourceID,
        "targetID": listetargetID,
        "geometry": listegeometry,
        "nature"  : listInfo,
        #"nature2"  : listInfo2,
        "ds" : listds
    })
    
    gdf = gpd.GeoDataFrame(
    df, geometry="geometry", crs="EPSG:2154"
    )

    gdf.to_file(url)   

# added id column names
def separer(popRef , popComp, id_ref, id_comp): 
    
    popRef4GMA = []
    popRef4MCA = []
    
    listeCandidat = SelectionCandidatInter(popRef, popComp, id_ref, id_comp)

    listPopRef = listeCandidat[0]
    
    listPopComp = listeCandidat[1]
    
    for i in range(len(listPopRef)):
        
        popRef_i = {}
        
        popRef_i[id_ref] = listPopRef[i][0]
        popRef_i['geometry'] = listPopRef[i][1]
        
        if len(listPopComp[i]) == 0 : 
            popRef4MCA.append(popRef_i)
        
        else : 
            
            a = 0 
            
            for j in range(len(listPopComp[i])) : 
                
                geomRef = listPopRef[i][1].buffer(0)
                geomComp = listPopComp[i][j][1].buffer(0)
                inter = shapely.intersection(geomRef , geomComp)
                union = shapely.union(geomRef,geomComp)
                
                ds =  1 - inter.area /union.area 
                
                if ds < 0.7 : a = 1
                
            if a == 1 : popRef4MCA.append(popRef_i)
            else : popRef4GMA.append(popRef_i)
        
    return ( popRef4GMA , popRef4MCA )

# added id column names
def complete(popRef , popComp, liste, id_ref, id_comp):
    
    for i in range(len(popRef)):
        a = 0
        for j in range(len(liste)):
            
            if popRef[i][id_ref] == liste[j][0] : 
                a = 1
        
        if a != 1 : 
            liste.append((popRef[i][id_ref] , 0 , 1))
            a = 0
            
    for i in range(len(popComp)):
        a = 0
        for j in range(len(liste)):
            
            if popComp[i][id_comp] == liste[j][0] : 
                a = 1
        
        if a != 1 : 
            liste.append(( popComp[i][id_comp] , 1 ,1))
            a = 0
    
    return  liste 
    
   
# remove because duplicated
# def readShapefile (url) :
#     data  = gpd.read_file(url)
#     columns = [data.columns[i] for i in range(len(data.columns)) ]
#     popRef = []
#     for i in range(len(data)):
#         L = {}
#         for j in range(len(columns)):
#             L[columns[j]] = data[columns[j]][i]
#         L['id_spatial'] = i
#         popRef.append(L) 
#     return popRef

# converts the dataframe to the expected array
def toArray(df):
    columns = [df.columns[i] for i in range(len(df.columns)) ]
    popRef = []
    for i in range(len(df)):
        L = {}
        for j in range(len(columns)):
            L[columns[j]] = df[columns[j]][i]
        L['id_spatial'] = i
        popRef.append(L) 
    return popRef

# added id column names, output_file and generatedIds
def main__(workdirectory, ref, comp, id_ref, id_comp, consigne, output_file, generatedIds=False):
    print(datetime.datetime.now(),"Running",consigne)
    url1 = workdirectory + str(ref)
    url2 = workdirectory + str(comp)
    
    param = {}
    param["surface_min_intersection"] = 1;
    param["pourcentage_min_intersection"] = 0.1;
    param["pourcentage_intersection_sur"] = 0.8;
    param["minimiseDistanceSurfacique"] = True;
    param["distSurfMaxFinal"] = 0.25;
    param["completudeExactitudeMinFinal"] = 0.8;
    param["regroupementOptimal"] = True;
    param["filtrageFinal"] = True;
    param["ajoutPetitesSurfaces"] = True;
    param["seuilPourcentageTaillePetitesSurfaces"] = 0.1;
    param["persistant"] = False;
    param["resolutionMin"] = 1;
    param["resolutionMax"] = 11;

    # read using geopandas to reuse the dataframes in the export
    gpd1  = gpd.read_file(url1)
    gpd2  = gpd.read_file(url2)
    # generated ids based on row number => same as 'id_spatial'...
    if generatedIds:
        gpd1[id_ref] = range(0, len(gpd1))
        gpd2[id_comp] = range(0, len(gpd2))
    popRef = toArray(gpd1)
    popComp = toArray(gpd2)

    # popRef = readShapefile (url1, id_ref, generatedIds)
    # popComp = readShapefile (url2, id_comp, generatedIds) 
    # url = "data/merged01.shp"
    
    if consigne == 'GMA' : 
        Appariement = appariementSurfaces(popRef, popComp, param)
        Appariement2 = appariementSurfaces(popComp,popRef, param)
        Appariement.extend(list(map(lambda m:[m[1],m[0],m[2]],Appariement2)))
        #TODO remove duplicates
        
    elif consigne == 'MCA' : 
        Appariement  = MCA(popRef, popComp, id_ref, id_comp)
        Appariement2 = MCA(popComp, popRef, id_comp, id_ref)
        Appariement.extend(list(map(lambda m:[m[1],m[0],m[2]],Appariement2)))
        #TODO remove duplicates
    
    elif consigne == 'Multi' : 
        ( popRef4GMA , popRef4MCA ) = separer(popRef , popComp, id_ref, id_comp)
        Appariement    = MCA(popRef4MCA, popComp, id_ref, id_comp)
        AppariementGMA = appariementSurfaces(popRef4GMA, popComp, param)
        for i in range(len(AppariementGMA)):
            # convert to a similar link type
            # careful: GMA and MCA do not handle the same ids: GMA returns row indices and MCA return actual ids (from the given column)
            match = AppariementGMA[i]
            Appariement.append([popRef4GMA[match[0]][id_ref],popComp[match[1]][id_comp],match[2]])
        # TODO make the matching the other way too?
    else : 
        return 'la consigne n est pas clair. Vous devez renseigner GMA, MCA ou multi'

    # create the link for the output
    def createLink(match):
        # print("LINK=",match[0],match[1],match[2])
        id1 = match[0]
        id2 = match[1]
        if consigne == 'GMA' : 
            # careful: GMA and MCA do not handle the same ids: GMA returns row indices and MCA return actual ids (from the given column)
            id1 = popRef[match[0]][id_ref]
            geom1 = popRef[match[0]]['geometry']
            id2 = popComp[match[1]][id_comp]
            geom2 = popComp[match[1]]['geometry']
        else :
            geom1 = gpd1[gpd1[id_ref]==match[0]]['geometry'].iloc[0]
            geom2 = gpd2[gpd2[id_comp]==match[1]]['geometry'].iloc[0]
        return (id1,id2,match[2],LineString([geom1.centroid,geom2.centroid]))
    links = list(map(createLink,Appariement))
    id1,id2,value,geom = list(zip(*links))
    id1List=list(id1)
    id2List=list(id2)
    valueList=list(value)
    geomList=list(geom)
    df = pd.DataFrame(
        {
            "ID1": id1List,
            "ID2": id2List,
            "prob": valueList,
        }
    )
    # remove duplicates (ignoring prob just in case)
    # df = df.drop_duplicates(subset=['ID1', 'ID2'])
    gdf = gpd.GeoDataFrame( df, geometry=geomList, crs="EPSG:2154" )
    gdf.to_file(output_file, layer=consigne, driver="GPKG")
    # for i, line in enumerate(Appariement):
    #     print(i,line)
    # for a in Appariement:
    #     print(a)

    # liste = complete(popRef, popComp , Appariement)
    # print(len(liste))
    #writeShapefile(popRef, popComp , Appariement , url)
    print(datetime.datetime.now(),"Done with",consigne)

# TODO ajouter paramètres cli (paramètre booléen pour appariement dans les 2 sens notamment)
if __name__ == "__main__" : 
    workdirectory = "data/"
    # test data
    # ref = 'popRef.shp'
    # comp = 'popComp.shp'
    # id_ref = 'ID'
    # id_comp = 'ID'
    # generatedIds = False
    ref = "Buildings_2011_02mp_buildings_with_fields_DP.shp"
    comp = "Buildings_2024_02mp_buildings_with_fields_DP.shp"
    id_ref = 'ID'
    id_comp = 'ID'
    # output_file = 'links.gpkg'
    output_file = 'links_iasi.gpkg'
    # generatedIds is a simple way to create the id columns based on the row number (no unique id on iasi data)
    generatedIds = True
    print(datetime.datetime.now(),"Let's go")
    main__(workdirectory, ref, comp, id_ref, id_comp, 'MCA', output_file, generatedIds)
    main__(workdirectory, ref, comp, id_ref, id_comp, 'GMA', output_file, generatedIds)
    main__(workdirectory, ref, comp, id_ref, id_comp, 'Multi', output_file, generatedIds)
    # saving ref & comp layers
    url1 = workdirectory + str(ref)
    url2 = workdirectory + str(comp)
    gpd1  = gpd.read_file(url1)
    gpd2  = gpd.read_file(url2)
    if generatedIds:
        # remove the id column because of conflict with ID. Ok, I could have changed the id_ref & id_comp to something else...
        gpd1 = gpd1.drop(columns=['id'])
        gpd2 = gpd2.drop(columns=['id'])
        # add the generated ids (I know, its being done 4 times now)
        gpd1[id_ref] = range(0, len(gpd1))
        gpd2[id_comp] = range(0, len(gpd2))
    # save the ref & comp layers
    gpd1.to_file(output_file, layer="ref", driver="GPKG")
    gpd2.to_file(output_file, layer="comp", driver="GPKG")
    print(datetime.datetime.now(),"All done")

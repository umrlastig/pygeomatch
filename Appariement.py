#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import shapely
from shapely.geometry import LineString
import geopandas as gpd
import pandas as pd
import datetime
import pymatch

"""""""""""""" """ COMMUN """ """"""""""""""""""
#############################################

# added id column names
def separate(popRef , popComp): 
    popRef4GMA = []
    popRef4MCA = []
    listPopRef, listPopComp = pymatch.selectCandidates(popRef, popComp)
    for i in range(len(listPopRef)):
        refIndex, refGeometry = listPopRef[i]
        if len(listPopComp[i]) == 0:
            # there is no candidate: can it really happen (I mean: does the function actually returns such a thing)?
            popRef4MCA.append(refIndex)
        else:
            def condition(comp):
                geomRef = refGeometry.buffer(0)
                geomComp = comp[1].buffer(0)
                inter = shapely.intersection(geomRef , geomComp)
                union = shapely.union(geomRef,geomComp)
                ds =  1 - inter.area /union.area 
                return ds < 0.7
            # a = 0
            # for j in range(len(listPopComp[i])) : 
            #     geomRef = listPopRef[i][1].buffer(0)
            #     geomComp = listPopComp[i][j][1].buffer(0)
            #     inter = shapely.intersection(geomRef , geomComp)
            #     union = shapely.union(geomRef,geomComp)
            #     ds =  1 - inter.area /union.area 
            #     if ds < 0.7 : a = 1
            # if a == 1 : popRef4MCA.append(listPopRef[i][0])
            # else : popRef4GMA.append(listPopRef[i][0])
            if any(condition(c) for c in listPopComp[i]):
                popRef4MCA.append(refIndex)
            else:
                popRef4GMA.append(refIndex)
    return popRef.iloc[popRef4GMA], popRef.iloc[popRef4MCA]

# added id column names
def complete(popRef , popComp, liste, id_ref, id_comp):    
    # count links from popRef
    for i in range(len(popRef)):
        a = 0
        for j in range(len(liste)):
            if popRef[i][id_ref] == liste[j][0] : 
                a = 1
        # if no link, add a 1, 0 link
        if a != 1 : 
            liste.append((popRef[i][id_ref] , 0 , 1))
            a = 0
            
    # count links from popComp
    for i in range(len(popComp)):
        a = 0
        for j in range(len(liste)):
            if popComp[i][id_comp] == liste[j][0] : 
                a = 1
        # if no link, add a 0, 1 link
        if a != 1 : 
            liste.append(( popComp[i][id_comp] , 1 ,1))
            a = 0
    
    return  liste 
    
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
def main__(workdirectory, ref, comp, consigne, output_file):
    print(datetime.datetime.now(),"Running",consigne)
    url1 = workdirectory + str(ref)
    url2 = workdirectory + str(comp)
    
    param = {}
    param["surface_min_intersection"] = 1
    param["pourcentage_min_intersection"] = 0.1
    param["pourcentage_intersection_sur"] = 0.8
    param["minimiseDistanceSurfacique"] = True
    param["distSurfMaxFinal"] = 0.25
    param["completudeExactitudeMinFinal"] = 0.8
    param["regroupementOptimal"] = True
    param["filtrageFinal"] = True
    param["ajoutPetitesSurfaces"] = True
    param["seuilPourcentageTaillePetitesSurfaces"] = 0.1
    param["persistant"] = False
    param["resolutionMin"] = 1
    param["resolutionMax"] = 11

    # read using geopandas to reuse the dataframes in the export
    gpd1  = gpd.read_file(url1)
    gpd2  = gpd.read_file(url2)
    def match_both_ways(f, p1, p2):
        """
        Use function f to match the features in both ways and return the merged result.
        
        :param f: matching function
        :param p1: parameters for way 1
        :param p2: parameters for way 2
        """
        matches_1 = f(*p1)
        matches_2 = f(*p2)
        # reverse the indices for the second matches
        return matches_1 + list(map(lambda m:[m[1],m[0],m[2]],matches_2))
    if consigne == 'GMA':
        Appariement = match_both_ways(pymatch.appariementSurfaces, (gpd1, gpd2, param), (gpd2, gpd1, param))        
    elif consigne == 'MCA':
        Appariement = match_both_ways(pymatch.MCA, (gpd1, gpd2), (gpd2, gpd1))
    elif consigne == 'Multi':
        def multi_match(ref, comp):
            popRef4GMA, popRef4MCA = separate(ref, comp)
            matches    = pymatch.MCA(popRef4MCA, comp)
            # this is sort of an ugly trick: we have to convert from the index (m[0]) of the separated dataframe (popRef4MCA) to the index from the global dataframe (gpd1)
            # to do that, we use the name of the index (with .name) and get the index with get_loc
            matches = list(map(lambda m: [ref.index.get_loc(popRef4MCA.iloc[m[0]].name), m[1], m[2]], matches))
            matches_GMA = pymatch.appariementSurfaces(popRef4GMA, comp, param)
            matches_GMA = list(map(lambda m: [ref.index.get_loc(popRef4GMA.iloc[m[0]].name), m[1], m[2]], matches_GMA))
            matches.extend(matches_GMA)
            return matches
        Appariement = match_both_ways(multi_match, (gpd1, gpd2), (gpd2, gpd1))
    elif consigne == 'MCA2':
        Appariement = match_both_ways(pymatch.MCA2, (gpd1, gpd2), (gpd2, gpd1))
    elif consigne == 'Multi2':
        def multi_match2(ref, comp):
            popRef4GMA, popRef4MCA = separate(ref, comp)
            matches    = pymatch.MCA2(popRef4MCA, comp)
            # this is sort of an ugly trick: we have to convert from the index (m[0]) of the separated dataframe (popRef4MCA) to the index from the global dataframe (gpd1)
            # to do that, we use the name of the index (with .name) and get the index with get_loc
            matches = list(map(lambda m: [ref.index.get_loc(popRef4MCA.iloc[m[0]].name), m[1], m[2]], matches))
            matches_GMA = pymatch.appariementSurfaces(popRef4GMA, comp, param)
            matches_GMA = list(map(lambda m: [ref.index.get_loc(popRef4GMA.iloc[m[0]].name), m[1], m[2]], matches_GMA))
            matches.extend(matches_GMA)
            return matches
        Appariement = match_both_ways(multi_match2, (gpd1, gpd2), (gpd2, gpd1))
    else :
        return 'la consigne n est pas clair. Vous devez renseigner GMA, MCA ou multi'

    df = pd.DataFrame({
        'ID1': list(map(lambda m: m[0], Appariement)),
        'ID2': list(map(lambda m: m[1], Appariement)),
        'prob': list(map(lambda m: m[2], Appariement))
    })
    # remove duplicates (ignoring prob just in case)
    df = df.drop_duplicates(subset=['ID1', 'ID2'])
    def createLinkGeometry(match):
        """
        Creates the linear geometry of the link using the indices.
        
        :param match: a link from a dataframe with columns ID1 & ID2
        """
        geom1 = gpd1.loc[match['ID1'],'geometry']
        geom2 = gpd2.loc[match['ID2'],'geometry']
        return LineString([geom1.centroid,geom2.centroid]) # type: ignore

    geomList = df.apply(createLinkGeometry, axis=1).tolist() # type: ignore
    # use crs from input file
    gdf = gpd.GeoDataFrame( df, geometry=geomList, crs=gpd1.crs )
    gdf.to_file(output_file, layer=consigne, driver="GPKG")
    # liste = complete(popRef, popComp , Appariement)
    # print(len(liste))
    #writeShapefile(popRef, popComp , Appariement , url)
    print(gpd1.index)
    print(datetime.datetime.now(),"Done with",consigne)

# TODO ajouter paramètres cli (paramètre booléen pour appariement dans les 2 sens notamment)
if __name__ == "__main__" : 
    workdirectory = "data/"
    # test data
    test = True
    # generatedIds is a simple way to create the id columns based on the row number (no unique id on iasi data)
    if test:
        ref = 'popRef.shp'
        comp = 'popComp.shp'
        output_file = 'links_gma.gpkg'
    else:
        ref = "Buildings_2011_02mp_buildings_with_fields_DP.shp"
        comp = "Buildings_2024_02mp_buildings_with_fields_DP.shp"
        output_file = 'links_iasi_2.gpkg'
    print(datetime.datetime.now(),"Let's go")
    
    # TODO we now export the integer-location based index as ID. 
    # Remains the question of the export: should we export the id selected by the used or the index selected by geopandas?

    # main__(workdirectory, ref, comp, 'GMA', output_file)
    # main__(workdirectory, ref, comp, 'MCA', output_file)
    # main__(workdirectory, ref, comp, 'MCA2', output_file)
    # main__(workdirectory, ref, comp, 'Multi', output_file)
    main__(workdirectory, ref, comp, 'Multi2', output_file)
    # saving ref & comp layers
    url1 = workdirectory + str(ref)
    url2 = workdirectory + str(comp)
    gpd1  = gpd.read_file(url1)
    gpd2  = gpd.read_file(url2)
    # add the index values used by the algorithms (for further joining or else)
    gpd1["pymatch_index"] = range(0, len(gpd1))
    gpd2["pymatch_index"] = range(0, len(gpd2))
    # save the ref & comp layers
    gpd1.to_file(output_file, layer="ref", driver="GPKG")
    gpd2.to_file(output_file, layer="comp", driver="GPKG")
    print(datetime.datetime.now(),"All done")

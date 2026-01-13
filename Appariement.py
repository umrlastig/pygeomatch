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
def separer(popRef , popComp, id_ref, id_comp): 
    popRef4GMA = []
    popRef4MCA = []
    
    listeCandidat = pymatch.selectCandidates(popRef, popComp, id_ref, id_comp)

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
        Appariement = pymatch.appariementSurfaces(gpd1, gpd2, param)
        Appariement2 = pymatch.appariementSurfaces(gpd2, gpd1, param)
        Appariement.extend(list(map(lambda m:[m[1],m[0],m[2]],Appariement2)))
        #TODO remove duplicates
        
    elif consigne == 'MCA' : 
        Appariement  = pymatch.MCA(popRef, popComp, id_ref, id_comp)
        Appariement2 = pymatch.MCA(popComp, popRef, id_comp, id_ref)
        Appariement.extend(list(map(lambda m:[m[1],m[0],m[2]],Appariement2)))
        #TODO remove duplicates
    
    elif consigne == 'Multi' : 
        ( popRef4GMA , popRef4MCA ) = separer(popRef , popComp, id_ref, id_comp)
        Appariement    = pymatch.MCA(popRef4MCA, popComp, id_ref, id_comp)
        AppariementGMA = pymatch.appariementSurfaces(popRef4GMA, popComp, param)
        for i in range(len(AppariementGMA)):
            # convert to a similar link type
            # careful: GMA and MCA do not handle the same ids: GMA returns row indices and MCA return actual ids (from the given column)
            match = AppariementGMA[i]
            Appariement.append([popRef4GMA[match[0]][id_ref],popComp[match[1]][id_comp],match[2]])
        # TODO make the matching the other way too?
    elif consigne == 'MCA2':
        matches = pymatch.MCA2(gpd1, gpd2)
        matches2 = pymatch.MCA2(gpd2, gpd1)
        import sys
        sys.exit()
    else : 
        return 'la consigne n est pas clair. Vous devez renseigner GMA, MCA ou multi'

    # create the link for the output
    def createLink(match):
        # print("LINK=",match[0],match[1],match[2])
        id1 = match[0]
        id2 = match[1]
        if consigne == 'GMA' : 
            # careful: GMA and MCA do not handle the same ids: GMA returns row indices and MCA return actual ids (from the given column)
            # id1 = popRef[match[0]][id_ref]
            # geom1 = popRef[match[0]]['geometry']
            # id2 = popComp[match[1]][id_comp]
            # geom2 = popComp[match[1]]['geometry']
            geom1 = gpd1.loc[id1,'geometry']
            geom2 = gpd2.loc[id2,'geometry']
        else :
            geom1 = gpd1[gpd1[id_ref]==match[0]]['geometry'].iloc[0]
            geom2 = gpd2[gpd2[id_comp]==match[1]]['geometry'].iloc[0]
        return (id1,id2,match[2],LineString([geom1.centroid,geom2.centroid])) # type: ignore
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
    
    # use crs from input file
    gdf = gpd.GeoDataFrame( df, geometry=geomList, crs=gpd1.crs )
    gdf.to_file(output_file, layer=consigne, driver="GPKG")

    # liste = complete(popRef, popComp , Appariement)
    # print(len(liste))
    #writeShapefile(popRef, popComp , Appariement , url)
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
        id_ref = 'ID'
        id_comp = 'ID'
        generatedIds = False
        output_file = 'links_gma.gpkg'
    else:
        ref = "Buildings_2011_02mp_buildings_with_fields_DP.shp"
        comp = "Buildings_2024_02mp_buildings_with_fields_DP.shp"
        id_ref = 'ID'
        id_comp = 'ID'
        output_file = 'links_iasi_2.gpkg'
        generatedIds = True
    print(datetime.datetime.now(),"Let's go")
    # main__(workdirectory, ref, comp, id_ref, id_comp, 'MCA2', output_file, generatedIds)
    main__(workdirectory, ref, comp, id_ref, id_comp, 'GMA', output_file, generatedIds)
    #main__(workdirectory, ref, comp, id_ref, id_comp, 'Multi', output_file, generatedIds)
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

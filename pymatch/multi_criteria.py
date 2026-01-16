import numpy
import geopandas as gpd
import numpy as np
import shapely
from itertools import groupby

"""""""""""""" """ MCA """ """"""""""""""""""
#############################################

def dempster(liste_critere):
    """
    Returns a decision.
    
    :param liste_critere: a list of candidates' criteria
    """
    # ajout candidat non app
    masseCan = [conjunctiveJoinRule(liste_critere[i]) for i in range(len(liste_critere))]
    # fusion candidat
    resultList, fusCandidat = candidateFusion(masseCan)
    # decision
    return decision(resultList, fusCandidat)

def conjunctiveJoinRule(criteria):
    """
    Smets' conjunctive join rule.
    
    :param criteria: a list of criteria
    """
    # print("criteria", criteria)
    dic = {}
    result = {}
    transit = {}
    if len(criteria) == 1:
        criteria = criteria[0]
        result["app"] = criteria[0]
        result["-app"] = criteria[1]
        result["theta"] = criteria[2]
        result["phi"] = 0
        return result
    for i in range(len(criteria)):
        dic["app"+str(i+1)] = criteria[i][0]
        dic["-app"+str(i+1)] = criteria[i][1]
        dic["theta" + str(i+1)] = criteria[i][2]
    result["app"] = 0
    result["-app"] = 0
    result["theta"] = 0
    result["phi"] = 0
    transit["app"] = 0
    transit["-app"] = 0
    transit["theta"] = 0
    transit["phi"] = 0

    matrice = [["app", "phi", "app"],
               ["phi", "-app", "-app"],
               ["app", "-app", "theta"]]

    for k in range(len(criteria)-1):

        result["app"] = 0
        result["-app"] = 0
        result["theta"] = 0
        result["phi"] = 0

        if k == 0:

            colone = ["app1", "-app1", "theta1"]

            ligne = ["app2", "-app2", "theta2"]

            matrice = [["app", "phi", "app"],
                       ["phi", "-app", "-app"],
                       ["app", "-app", "theta"]]

            for i in range(np.shape(matrice)[0]):
                for j in range(np.shape(matrice)[1]):
                    jeu = matrice[i][j]
                    if dic[colone[j]] * dic[ligne[i]] != 0:
                        result[jeu] += dic[colone[j]] * dic[ligne[i]]

        else:

            ligne = ["app"+str(k+2), "-app"+str(k+2), "theta"+str(k+2)]

            colone = ["app", "-app", "theta", "phi"]

            matrice = [["app", "phi", "app"],
                       ["phi", "-app", "-app"],
                       ["app", "-app", "theta"],
                       ["phi", "phi", "phi"]]

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
    # print("result", result)
    return result


def candidateFusion(candidates):
    """
    Fusion of candidates.
    
    :param candidates: a list of combined masses for candidates
    """
    # print("candidates",candidates)
    # Changer pa rapport au fait que masseCan est un dictionnaire
    # intialisation des dictionnaires
    if len(candidates) > 9:
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
        return (resultList, result)
    dic = {}
    transit = {}
    result = {}
    resultList = []
    for i in range(len(candidates)):
        dic["C"+str(i+1)] = candidates[i]["app"]
        dic["-C"+str(i+1)] = candidates[i]["-app"]
        dic["theta" + str(i+1)] = candidates[i]["theta"]
        dic["phi" + str(i+1)] = candidates[i]["phi"]
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
    # result["NA"] = 0
    resultList.append("theta")
    resultList.append("phi")
    # resultList.append("NA")

    if len(candidates) == 1:

        result["C1"] = candidates[i]["app"]
        result["-C1"] = candidates[i]["-app"]
        result["theta"] = candidates[i]["theta"]
        result["phi"] = candidates[i]["phi"]

    for k in range(len(candidates)-1):

        (ligne, colone, matrice) = mat(k + 2)#, len(candidat))

        listmodif = []

        for i in range(len(resultList)):
            result[resultList[i]] = 0

        for i in range(np.shape(matrice)[0]):
            for j in range(np.shape(matrice)[1]):
                jeu = str(matrice[i][j])[2:len(str(matrice[i][j]))-1]
                listmodif.append(jeu)
                # if dic[colone[j]] * dic[ligne[i]] != 0 :
                if k == 0:
                    coeff = dic[colone[j]] * dic[ligne[i]]
                else:
                    if ligne[i][0] == 't':
                        coeff = transit['theta'] * dic[colone[j]]
                    elif ligne[i][0] == 'p':
                        coeff = transit['phi'] * dic[colone[j]]
                    else:
                        coeff = transit[ligne[i]] * dic[colone[j]]

                # ajout NA
                if jeu == "NA":
                    lis = [k+1 for k in range(len(candidates))]
                    a = int(colone[j][::-1][0])
                    b = int(ligne[i][::-1][0])
                    reste = []
                    for n in range(len(lis)):
                        if lis[n] in (a, b):
                            pass
                        else:
                            reste.append(lis[n])

                    if reste == []:
                        result[jeu] = coeff
                    else:
                        string = 'NA'
                        string += str(reste)
                        result[string] = coeff

                else:
                    result[jeu] += coeff

        for i in range(len(listmodif)):
            if listmodif[i] != "phi" and listmodif[i] != "theta" and listmodif[i] != "NA":
                transit[listmodif[i]] = result[listmodif[i]]
        transit["theta"] = result["theta"]
        transit["phi"] = result["phi"]

    a = result['phi']
    if a != 1:
        for i in range(len(result)):
            result[list(result.keys())[i]] = result[list(
                result.keys())[i]] / (1 - a)

    result['phi'] = 0

    """
    sum_ = 0
    for i in range(len(result)):
        string = str(list(result.keys())[i])
        sum_ += result[string]
    """

    # print("result",result)
    return (resultList, result)

def loi(x, y):

    # x = Cx , -Cx , theta , phi

    if x == y:
        return x

    if len(x) > 2:
        if x[0:3] == "phi":
            return "phi"

        if x[0:5] == "theta":
            if len(y) > 2 and y[0:3] == "phi":
                return "phi"
            else:
                return y

        if x[0] == 'C':
            if y[0] == 'C':
                return "phi"
            if y[0] == '-':
                return x

    if len(y) > 2:
        if y[0:3] == "phi":
            return "phi"
        if y[0:5] == "theta":
            return x

        if y[0] == 'C':
            if x[0] == 'C':
                return "phi"
            if x[0] == '-':
                return y

    if x[0] == "C":

        if y[0] == "C":
            if x[1] != y[1]:
                return "phi"

        if y[0] == "-":

            if x[1] != y[2]:
                return x
            elif len(x) == len(y) - 2:
                return x
            else:
                return "NA"

    if x[0] == "-":

        if y[0] == "C":
            if x[2] != y[1]:
                return y
            else:
                return "NA"

        if y[0] == "-":
            return "NA"

def mat(index):
    colone = ("C" + str(index), "-C" + str(index), "theta" + str(index), "phi" + str(index))
    ligne = []
    for i in range(index-1):
        ligne.append("C"+str(i+1))
        ligne.append("-C"+str(i+1))
    ligne.append("theta" + str(1))
    ligne.append("phi" + str(1))
    matrice = np.chararray((len(ligne), len(colone)), 5)
    for i in range(len(colone)):
        for j in range(len(ligne)):
            matrice[j][i] = loi(ligne[j], colone[i])
    return (ligne, colone, matrice)

def decision(resultList, fusion):
    """
    Docstring for decision
    
    :param resultList: Description
    :param fusion: Description
    """
    proba = {}
    a = int((len(resultList) - 2) / 2)
    if a == 1:
        if fusion["C1"] > fusion["-C1"] or fusion["C1"] > fusion["theta"]:
            return "C1"
        else:
            return "NA"
    app = []
    _app = []
    for i in range(a):
        app.append(fusion['C' + str(i+1)])
        _app.append(fusion['-C' + str(i+1)])
    result = np.array(app) @ numpy.eye(a)
    result += np.array(_app) @ (1/(a) * (numpy.ones(a) - numpy.eye(a)))
    ligneNA = []
    matrice = []
    for i in range(len(fusion)):
        if list(fusion.keys())[i][0:2] == 'NA':
            ligneNA.append(fusion[list(fusion.keys())[i]])
            ligne = [0. for i in range(a)]
            for j in range(int((len(list(fusion.keys())[i])) / 3)):
                """
                if 3*j + 3 == 24 : 
                    chiffre = 10
                    ligne[int(chiffre) - 1] = 1/(a - 1)
                """
                chiffre = list(fusion.keys())[i][3*j + 3]
                ligne[int(chiffre) - 1] = 1/(a - 1)

            matrice.append(ligne)

    result += np.array(ligneNA) @ matrice

    resultNA = np.array(ligneNA) @ (1/(a - 1) * np.ones(len(ligneNA)))
    resultNA += np.array(_app) @ (1/(a) * np.ones(a))

    for i in range(a):
        proba['C' + str(i+1)] = result[i]
        proba['C' + str(i+1)] += fusion['theta']/(a+1)
    proba['NA'] = resultNA
    proba['NA'] += fusion['theta']/(a+1)

    # sum_ = 0
    # for i in range(len(proba)):
    #     string = str(list(proba.keys())[i])
    #     sum_ += proba[string]

    max_ = 0
    for i in range(len(list(proba.keys()))):
        if proba[str(list(proba.keys())[i])] > max_:
            max_ = proba[str(list(proba.keys())[i])]
            result = str(list(proba.keys())[i])

    epsilon = 0.01
    # if the probability for NA is close to the max, chose NA too
    if (result == 'NA') | (max_ < proba['NA'] + epsilon):
        return "NA"

    return result

def select_candidates(popRef: gpd.GeoDataFrame, popComp: gpd.GeoDataFrame):
    """
    Identify candidates for ref features.
    It returns two lists:
    - a list of (index, geometry) from the ref geodataframe
    - a list of corresponding candidates (index, geometry) from the comp geodataframe
    
    :param popRef: ref features
    :type popRef: gpd.GeoDataFrame
    :param popComp: comp features
    :type popComp: gpd.GeoDataFrame
    """
    # The first subarray contains input geometry integer indices. The second subarray contains tree geometry integer indices.
    refIndices, compIndices = popComp["geometry"].sindex.query(popRef["geometry"], predicate="intersects")
    # # creation liste popRef
    # listeRef = [(idRef[refIndices[0]], geomRef[refIndices[0]])]
    # for i in range(len(refIndices)-1):
    #     if refIndices[i] == refIndices[i+1]:
    #         continue
    #     else:
    #         listeRef.append((idRef[refIndices[i+1]], geomRef[refIndices[i+1]]))

    # # creation liste popComp
    # listeComp = []
    # listeComp_i = [(idComp[compIndices[0]], geomCom[compIndices[0]])]
    # for i in range(len(refIndices)-1):
    #     if refIndices[i] == refIndices[i+1]:
    #         listeComp_i.append(
    #             (idComp[compIndices[i+1]], geomCom[compIndices[i+1]]))
    #     else:
    #         listeComp.append(listeComp_i)
    #         listeComp_i = [(idComp[compIndices[i+1]],
    #                         geomCom[compIndices[i+1]])]
    # listeComp.append(listeComp_i)
    # return (listeRef, listeComp)
    # print(popRef.head())
    zipped = zip(refIndices, compIndices)
    # print(refIndices)
    # group by the ref index (z[0]), map the results to dict entries with the ref index (y[0]) and keep only the second element of the grouped tuples (t[1])
    # ref, comp = zip(*list(map(lambda y: ((y[0], popRef.loc[y[0],"geometry"]), list(map(lambda t: (t[1], popComp.loc[t[1],"geometry"]), y[1]))), groupby(zipped, lambda z: z[0]))))
    ref, comp = zip(*list(map(lambda y: ((y[0], popRef.iloc[[y[0]]].iloc[0]["geometry"]), list(map(lambda t: (t[1], popComp.iloc[[t[1]]].iloc[0]["geometry"]), y[1]))), groupby(zipped, lambda z: z[0]))))
    return (list(ref), list(comp))


def processMatch(featureRef, listPopComp):
    """
    Process a feature and its candidates.
    
    :param featureRef: a feature
    :param listPopComp: candidates
    """
    listCritere = []
    geomRef = featureRef[1]
    listCritere_i = []
    for i in range(len(listPopComp)):
        list_i = []
        geomComp = listPopComp[i][1]
        # Critere surfacique
        geomRef = geomRef.buffer(0)
        geomComp = geomComp.buffer(0)
        inter = shapely.intersection(geomRef, geomComp)
        union = shapely.union(geomRef, geomComp)
        ds = 1 - inter.area / union.area
        distance = ds
        tableau = []
        T1 = 0.90
        T2 = 1.0
        E = 0.01
        S = 0.7
        K = 1 - E - S
        if distance < T1:
            app = (-(1 - E)/T2) * distance + 1 - E
            _app = E
            tableau.append(app)
            tableau.append(_app)
            tableau.append(1 - app - _app)
        elif distance < T2:
            app = (-(1 - E)/T2) * distance + 1 - E
            _app = (K - E) * distance / (T2 - T1) + E - (K - E)*T1/(T2 - T1)
            tableau.append(app)
            tableau.append(_app)
            tableau.append(1 - app - _app)
        else:
            app = E
            _app = K
            tableau.append(E)
            tableau.append(_app)
            tableau.append(1 - app - _app)
        cs = tableau
        listCritere_i.append(ds)
        list_i.append(cs)
        list_i.append(cs)# added twice
        listCritere.append(list_i)
    lres = dempster(listCritere)
    return (listCritere_i, listCritere, lres)

def MCA(popRef: gpd.GeoDataFrame, popComp: gpd.GeoDataFrame):
    """
    Process Multi Criteria Matching.
    
    :param popRef: Ref features
    :param popComp: Comp features
    """
    # get the list of ref objects and their corresponding candidates
    listPopRef, listPopComp = select_candidates(popRef, popComp)
    App = []
    for i in range(len(listPopRef)):
        # match ref feature i with its candidates
        # listPopRef[i] is refIndex, refGeometry
        (listCritere_i, _, liste) = processMatch(listPopRef[i], listPopComp[i])
        if not((liste == "NA") | (liste == "theta")):
            decision = int(liste[1]) - 1
            App.append((listPopRef[i][0], listPopComp[i][decision][0], listCritere_i[decision]))
    return App

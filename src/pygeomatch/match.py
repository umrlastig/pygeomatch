# /// script
# dependencies = [
#   "typer"
# ]
# ///

import typer
from typing_extensions import Annotated
from pathlib import Path
from enum import Enum
from datetime import datetime
import pandas as pd
import geopandas as gpd
from shapely import LineString
from pygeomatch.geometric_matching_of_areas import surface_match
from pygeomatch.multi_criteria import MCA, select_candidates
from pygeomatch.multi_criteria2 import MCA2
from pygeomatch.util import surface_distance
from pygeomatch.optimal_transport import OT_matching

class MatchingAlgorithm(str, Enum):
    gmoa = "GMA"
    multi_criteria = "MCA"
    multi = "Multi"
    multi_criteria2 = "MCA2"
    multi2 = "Multi2"
    ot = "OT"

def separate(popRef , popComp) -> tuple[gpd.GeoDataFrame,gpd.GeoDataFrame]: 
    """
    Separates the data for the Multi algorithm.
    
    :param popRef: reference dataset
    :param popComp: comparison dataset
    :return: the reference dataset split into 2 datasets
    :rtype: tuple[GeoDataFrame, GeoDataFrame]
    """
    ref_gma = []
    ref_mca = []
    ref_candidates, comp_candidates = select_candidates(popRef, popComp)
    for i in range(len(ref_candidates)):
        ref_index, ref_geom = ref_candidates[i]
        if len(comp_candidates[i]) == 0:
            # there is no candidate: can it really happen (I mean: does the function actually returns such a thing)?
            ref_mca.append(ref_index)
        else:
            def condition(comp):
                return surface_distance(ref_geom, comp[1]) < 0.7
            if any(condition(c) for c in comp_candidates[i]):
                ref_mca.append(ref_index)
            else:
                ref_gma.append(ref_index)
    return popRef.iloc[ref_gma], popRef.iloc[ref_mca]

app = typer.Typer(name="pymatch", help="Pymatch matches geographical features", add_completion=True)

@app.command()
def main(
        file1: Annotated[Path, typer.Argument(help="File to use as reference")],
        file2: Annotated[Path, typer.Argument(help="File to use as comparison")],
        output_file: Annotated[Path, typer.Argument(help="File to save the results to")],
        algorithm: Annotated[MatchingAlgorithm, typer.Argument(help="The algorithm to use (GMA, MCA, Multi, MCA2, Multi2)")] = MatchingAlgorithm.multi_criteria,
        export_input: Annotated[bool, typer.Option("--export",help="If true, export the input layers in the output file.")] = True,
        do_match_both_ways: Annotated[bool, typer.Option("--both/--no-both","-b/-B", help="If true, match the features both ways (A=>B and B=>A) and merge the results.")] = True,
        min_surface_intersection:  Annotated[float, typer.Option(min=0, help="min surface intersection between fetures to consider matching")] = 1.,
        min_intersection_percentage: Annotated[float, typer.Option(min=0, help="min surface intersection between fetures to consider matching")] = 0.1,
        sure_intersection_percentage: Annotated[float, typer.Option(min=0, help="min surface intersection between fetures to consider matching")] = 0.8,
        minimise_surface_distance: Annotated[bool, typer.Option(help="min surface intersection between fetures to consider matching")] = True,
        min_surface_distance: Annotated[float, typer.Option(min=0, help="min surface distance between fetures to consider matching in the final evaluation (only used if minimise_surface_distance is true)")] = 0.25,
        min_accuracy_completeness: Annotated[float, typer.Option(min=0, help="min accuracy and completeness between fetures to consider matching in the final evaluation (only used if minimise_surface_distance is false)")] = 0.8,
        use_optimal_groups: Annotated[bool, typer.Option(help="if true, searches for optimal groups")] = True,
        final_filtering: Annotated[bool, typer.Option(help="if true, filter the final links using min_surface_distance or min_accuracy_completeness depending on minimise_surface_distance")] = True
        # ajoutPetitesSurfaces: Annotated[bool, typer.Option(help="min surface intersection between fetures to consider matching")] = True,
        # seuilPourcentageTaillePetitesSurfaces: Annotated[float, typer.Option(min=0, help="min surface intersection between fetures to consider matching")] = 0.1
):
    print(datetime.now(),"Running",algorithm)
    # read using geopandas to reuse the dataframes in the export
    gpd1  = gpd.read_file(file1)
    gpd2  = gpd.read_file(file2)
    param = {
        "min_surface_intersection": min_surface_intersection,
        "min_intersection_percentage": min_intersection_percentage,
        "sure_intersection_percentage": sure_intersection_percentage,
        "minimise_surface_distance": minimise_surface_distance,
        "min_surface_distance": min_surface_distance,
        "min_accuracy_completeness": min_accuracy_completeness,
        "use_optimal_groups": use_optimal_groups,
        "final_filtering": final_filtering,
        # "ajoutPetitesSurfaces": True,
        # "seuilPourcentageTaillePetitesSurfaces": 0.1
    }
    def match_both_ways_or_not(f, p1, p2):
        """
        Use function f to match the features in both ways and return the merged result.
        
        :param f: matching function
        :param p1: parameters for way 1
        :param p2: parameters for way 2
        """
        matches_1 = f(*p1)
        if not do_match_both_ways:
            return matches_1
        matches_2 = f(*p2)
        # reverse the indices for the second matches
        return matches_1 + list(map(lambda m:[m[1],m[0],m[2]],matches_2))
    def do_match():
        if algorithm == MatchingAlgorithm.gmoa:
            return match_both_ways_or_not(surface_match, (gpd1, gpd2, param), (gpd2, gpd1, param))
        elif algorithm == MatchingAlgorithm.multi_criteria:
            return match_both_ways_or_not(MCA, (gpd1, gpd2), (gpd2, gpd1))
        elif algorithm == MatchingAlgorithm.multi:
            def multi_match(ref, comp):
                popRef4GMA, popRef4MCA = separate(ref, comp)
                matches    = MCA(popRef4MCA, comp)
                # this is sort of an ugly trick: we have to convert from the index (m[0]) of the separated dataframe (popRef4MCA) to the index from the global dataframe (gpd1)
                # to do that, we use the name of the index (with .name) and get the index with get_loc
                matches = list(map(lambda m: [ref.index.get_loc(popRef4MCA.iloc[m[0]].name), m[1], m[2]], matches))
                matches_GMA = surface_match(popRef4GMA, comp, param)
                matches_GMA = list(map(lambda m: [ref.index.get_loc(popRef4GMA.iloc[m[0]].name), m[1], m[2]], matches_GMA))
                matches.extend(matches_GMA)
                return matches
            return match_both_ways_or_not(multi_match, (gpd1, gpd2), (gpd2, gpd1))
        elif algorithm == MatchingAlgorithm.multi_criteria2:
            return match_both_ways_or_not(MCA2, (gpd1, gpd2), (gpd2, gpd1))
        elif algorithm == MatchingAlgorithm.multi2:
            def multi_match2(ref, comp):
                ref_gma, ref_mca = separate(ref, comp)
                matches    = MCA2(ref_mca, comp)
                # this is sort of an ugly trick: we have to convert from the index (m[0]) of the separated dataframe (popRef4MCA) to the index from the global dataframe (gpd1)
                # to do that, we use the name of the index (with .name) and get the index with get_loc
                matches = list(map(lambda m: [ref.index.get_loc(ref_mca.iloc[m[0]].name), m[1], m[2]], matches))
                matches_GMA = surface_match(ref_gma, comp, param)
                matches_GMA = list(map(lambda m: [ref.index.get_loc(ref_gma.iloc[m[0]].name), m[1], m[2]], matches_GMA))
                matches.extend(matches_GMA)
                return matches
            return match_both_ways_or_not(multi_match2, (gpd1, gpd2), (gpd2, gpd1))
        elif algorithm == MatchingAlgorithm.ot:
            # TODO: consider a symmetric match? Not sure if the OT transport we use is perfectly symmetric.
            return OT_matching(gpd1, gpd2)
        return None
    # use unpacking to tranform a list((a,b,c)) into a list(a), list(b), list(c)
    res_match = do_match()
    if res_match is None:
        print("Unknown algorithm")
        return
    id1, id2, prob = zip(*res_match)
    df = pd.DataFrame({'ID1': id1,'ID2': id2,'prob': prob})
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
    geom_list = df.apply(createLinkGeometry, axis=1).tolist() # type: ignore
    # use crs from input file
    gdf = gpd.GeoDataFrame( df, geometry=geom_list, crs=gpd1.crs )
    gdf.to_file(output_file, layer=str(algorithm), driver="GPKG")
    print(datetime.now(),"Done with",algorithm)
    if export_input:
        # add the index values used by the algorithms (for further joining or else)
        gpd1["pymatch_index"] = range(0, len(gpd1))
        gpd2["pymatch_index"] = range(0, len(gpd2))
        # save the ref & comp layers
        gpd1.to_file(output_file, layer="ref", driver="GPKG")
        gpd2.to_file(output_file, layer="comp", driver="GPKG")
        print(datetime.now(),"All done")

if __name__ == "__main__":
    # typer.run(main)
    app()

import logging

import geopandas as gpd
import numpy as np

# POT
import ot
import ot.partial

INV_POINT_DENSITY = 100 # inverse of point density, i.e. number of points per squared meter


def OT_matching(popRef: gpd.GeoDataFrame, popComp: gpd.GeoDataFrame, *params):
    """
    Perform polygon matching based on (partial) optimal transport

    :param popRef: Ref features
    :param popComp: Comp features
    """

    # Sample points proportional to surface area
    # More points should theoretically mean more accurate distance (i.e. how much we
    # need to "warp" a polygon to match another), however it's also more computationally
    # expensive, especially in RAM (cost matrix is O(n²)).
    # We put at least 1 point for each polygon
    n_points_ref = 1 + popRef.area.astype("int32") // INV_POINT_DENSITY
    n_points_comp = 1 + popComp.area.astype("int32") // INV_POINT_DENSITY
    
    # TODO: use non-uniform sampling?
    # TODO: pretty sure the correct thing to do would be to sample points in a regular grid
    # But also it might be more difficult to ensure that every polygon has at least one point.
    # It could also lead to unexpected issues with OT because the regular grid would put all
    # points to the same distances, and discretize the Euclidean distance.
    # So all in all, I don't know.
    points_ref = popRef.sample_points(size=n_points_ref).explode()
    points_comp = popComp.sample_points(size=n_points_comp).explode()
    
    building_idx_ref = np.concatenate([np.ones(n, dtype="int32") * idx for idx, n in enumerate(n_points_ref)])
    building_idx_comp = np.concatenate([np.ones(n, dtype="int32") * idx for idx, n in enumerate(n_points_comp)])

    # Convert GeoSeries into numpy array
    arr_ref = np.stack([points_ref.x, points_ref.y], axis=1)
    arr_comp = np.stack([points_comp.x, points_comp.y], axis=1)

    # Compute cost matrix (based on Euclidean distance, so hopefully you use a meter CRS)
    cost_matrix = ot.dist(arr_ref, arr_comp)
    if max(cost_matrix.shape) > 2000:
        logging.warn(f"Cost matrix is {cost_matrix.shape}, this might take a while")


    #solved_transport = ot.solve(cost_matrix, unbalanced=1e1)#, weights_comp, weights_ref)
    # FIXME: set mass based on surface ratio instead of points ratio
    mass_ratio = min(len(arr_ref)/len(arr_comp), len(arr_comp)/len(arr_ref))
    # Compute the transport plan
    partial_transport = ot.partial.partial_wasserstein(ot.utils.unif(len(arr_ref)), ot.utils.unif(len(arr_comp)), cost_matrix, m=mass_ratio)

    # Create a dict that matches "index of polygon from popRef" to a list of tuples (index of polygon in popComp, mass to move)
    pairings = {}

    # Loop over all polygons in the source DataFrame
    for source_idx in popRef.index:
        # Get a source polygon
        building = popRef.iloc[source_idx]
        # Get all points that have been sampled for this polygon
        points_indices = np.nonzero(building_idx_ref == source_idx)[0]
    
        targets_indices = np.zeros_like(points_indices)

        for i, p in enumerate(points_indices):
            # Because we use partial OT, the point can be split so we need to find for this (source) point
            # to which target point most of its mass has been moved to by the transport plan
            max_target_point = np.argmax(partial_transport[p])
            # If the maximum transported mass is zero, then this point has no match
            if partial_transport[p][max_target_point] == 0.:
                target = -99999 # FIXME
            # Otherwise find its match
            else:
                target = max_target_point
            targets_indices[i] = target
        # Find all points with a valid match
        valid_targets = targets_indices != -99999 # FIXME
        #points_affected_ratio = np.count_nonzero(valid_targets) / len(targets_indices)
        # Get the polygon indices in the target DataFrame
        target_buildings_idx = building_idx_comp[targets_indices[valid_targets]]
    
        # Find out how many points of the source polygon has been matched to the target polygons
        target_buildings, weights = np.unique(target_buildings_idx, return_counts=True)
        # Normalize the mass moved to obtain a "weight" in [0,1] of the match at the polygon level
        normalized_weights = weights / len(target_buildings_idx)
        pairings[source_idx] = list(zip(target_buildings, normalized_weights))
        # TODO: at this stage, we could do some kind of post-processing and first stage decision making
        # e.g.:
        # - very weak matches (weight < 0.05) are unlikely, so we could remove them
        # - very strong matches (weight < 0.9) are very likely perfect matches, so we could round them to 1 and remove the rest
        # - (source) polygons with barely any matches are likely to have been removed
        # - target polygons with barely any matches are likely to be new ones

    # TODO: refactor this and the for loop above in a single loop
    matchings = []
    for source_idx, targets in pairings.items():
        for target_idx, pairing_weight in targets:
            matchings.append((source_idx, target_idx, pairing_weight))
    
    return matchings

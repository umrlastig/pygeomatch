import shapely
import networkx as nx
from itertools import groupby

def surface_distance(geomA, geomB) -> float:
    """
    Surface distance(A,B) = 1 - inter(A,B).area/union(A,B).area
    
    :param geomA: geometry A
    :param geomB: geometry B
    """
    inter = geomA.intersection(geomB)
    union = shapely.union(geomA.buffer(0),geomB.buffer(0))
    return 1 - inter.area / union.area

def get_accuracy(geomA, geomB) -> float:
    """
    accuracy(A,B) = surface(A inter B) / Surface (A)
    
    :param geomA: geometry A
    :param geomB: geometry B
    """
    inter = shapely.intersection(geomA , geomB)
    if inter.is_empty:
        return 0
    return inter.area / geomA.area
    
def get_completeness(geomA , geomB) -> float:
    """
    completeness(A,B) = surface(A inter B) / Surface (B)
    completeness(A,B) = accuracy(B,A)

    :param geomA: geometry A
    :param geomB: geometry B
    """
    return get_accuracy(geomB, geomA)

def measures(geomA , geomB) -> list[float]:
    """
    Returns all measures: (intersection, 
    :param geomA: geometry A
    :param geomB: geometry B
    """
    inter = geomA.intersection(geomB)
    intersectionArea = inter.area
    union = shapely.union(geomA.buffer(0),geomB.buffer(0))
    unionArea = union.area
    geomA_area = geomA.area
    geomB_area = geomB.area
    minArea = min(geomA_area, geomB_area)
    intersection_ratio = intersectionArea / minArea
    surf_distance = 1 - intersectionArea / unionArea
    accuracy = intersectionArea / geomA_area
    completeness = intersectionArea / geomB_area
    return [intersectionArea, intersection_ratio, surf_distance, accuracy, completeness]

class MatchingLink:
    def __init__(self, ref: int, comp: int, group: int, measures: dict[str, float]):
        self.ref = ref
        self.comp = comp
        self.group = group
        self.measures = measures
    def __repr__(self):
        return "MatchingLink(%s,%s,%s,%s)" % (self.ref, self.comp, self.group, self.measures)
    def as_tuple(self)-> tuple[int, int, int, dict]:
        return (self.ref, self.comp, self.group, self.measures)
    
def links_grouped_by_component(links: list[MatchingLink]) -> list[list[MatchingLink]]:
    """
    Returns a list of lists, where each inner list contains the Link
    objects that belong to the same connected component.
    """
    # ---- Build an undirected (bipartite) graph -----------------
    G = nx.Graph()
    # Add an edge for every Link (NetworkX adds the nodes automatically)
    G.add_edges_from((("ref", l.ref), ("comp", l.comp)) for l in links)
    # add edges between nodes of the same group
    data = sorted(links, key=lambda l: l.group)
    for _, g in groupby(data, key=lambda l: l.group):
        group = list(g)
        first = group[0]
        for other in group[1:]:
            G.add_edge(("ref", first), ("ref", other))
    # ---- Find connected components (as sets of node names) ---
    node_components = list(nx.connected_components(G))
    # ---- Map each node → component index for quick lookup -------
    node_to_comp = {}
    for idx, comp_nodes in enumerate(node_components):
        for n in comp_nodes:
            if (n[0] == "ref"):
                node_to_comp[n] = idx
    # ---- Bucket the original Link objects ----------------------
    groups = [[] for _ in node_components]
    for l in links:
        # Both endpoints are in the same component, so we can look up either.
        comp_idx = node_to_comp[("ref", l.ref)]
        groups[comp_idx].append(l)
    return groups

from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopExp import TopExp_Explorer, topexp
from OCC.Core.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCC.Core.TopoDS import topods
from OCC.Core.BRepAdaptor import BRepAdaptor_Surface
from OCC.Core.TopTools import TopTools_IndexedDataMapOfShapeListOfShape
from OCC.Core.STEPControl import STEPControl_Reader
import networkx as nx
from OCC.Core.BRepGProp import brepgprop
from OCC.Core.GProp import GProp_GProps
from OCC.Core.Bnd import Bnd_Box
from OCC.Core.BRepBndLib import brepbndlib
from OCC.Core.BRepAdaptor import BRepAdaptor_Curve
from OCC.Core.GCPnts import GCPnts_AbscissaPoint
import pymetis
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np


def get_edge_length(edge):
    try:
        curve_adaptor = BRepAdaptor_Curve(edge)
        length = GCPnts_AbscissaPoint.Length(curve_adaptor)
        return length if length > 0 else 1e-6
    except:
        return 1e-6


def load_step(filename):
    reader = STEPControl_Reader()
    reader.ReadFile(filename)
    reader.TransferRoots()
    return reader.OneShape()


def brep_to_face_adjacency_graph(shape):
    G = nx.Graph()

    # --- Build edge-to-face map first ---
    # this maps every edge -> list of faces that contain it
    edge_face_map = TopTools_IndexedDataMapOfShapeListOfShape()
    # map every topological edge to its face, makes adjacency lookup direct during edge pass
    # topological edge to ancestor faces
    topexp.MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edge_face_map)

    # --- Node pass: one node per face ---
    # computes per-face attributes (surface type, area, bounding box midpoint as centroid proxy
    # and stores them alongside graph nodes)
    face_explorer = TopExp_Explorer(shape, TopAbs_FACE)
    face_index = {}
    i = 0

    # iterates through all faces
    while face_explorer.More():
        # returns a generic TopoDS shape and downcasts it so face-specific APIs can be used
        face = topods.Face(face_explorer.Current())
        # creates a geometry adaptor over underlying surface
        adaptor = BRepAdaptor_Surface(face)
        # identifying geometric type from a given shape - used for graph attributes or filtering later on
        surf_type = adaptor.GetType()

        # container that holds geometric properties
        props = GProp_GProps()
        # traverses face and fills the container - integrates over 2D surface (surface integral)
        brepgprop.SurfaceProperties(face, props)
        # retrieves the area
        area = props.Mass()

        # bounding box approximation falls apart when geometry isn't symmetric
        bbox = Bnd_Box()
        brepbndlib.Add(face, bbox)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
        centroid = ((xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2)

        # i may not be deterministic across different runs or machines, global ID must be intrinsic
        # to CAD model itself i.e a hash or a deterministic ordering based on geometric properties (like centroid)
        # graph with generic properties, can add as much info as required
        G.add_node(i, shape=face, surface_type=surf_type, area=area,
                   bbox=bbox, centroid=centroid, global_id=i,
                   is_boundary=False, partition_id=None)

        face_index[face.TShape()] = i
        i += 1
        face_explorer.Next()

    # --- Edge pass: connect faces sharing a topological edge ---
    # creates an explorer on the shape and the type of shapes to search
    # edges are visited to establish which faces are adjacent
    # edges are topological connectors, and aren't being stored in the graph
    # the only thing extracted from each edge is its length, used as the edge weight
    edge_explorer = TopExp_Explorer(shape, TopAbs_EDGE)
    visited_edges = set()

    while edge_explorer.More():
        # current returns the current shape in exploration, casts it to the specific type that I need
        # grab the current object and make it useable as an edge
        edge = topods.Edge(edge_explorer.Current())
        # TShape returns the underlying topological entity which is a pointer to the shared geometry
        # two TopoDS objects that look different but represent the same edge will have the same TShape
        # list of pointers to uniquely identify each edge
        edge_hash = edge.TShape()

        # check that you have not traversed this edge before
        # avoid duplication of graph edges
        if edge_hash not in visited_edges:
            visited_edges.add(edge_hash)
            adjacent_faces = []

            # looks up the edge in the high-level map built by topexp.MapShapesAndAncestors
            # takes each edge, retrieves all the faces that contain it, casts each from a generic
            # shape to a face, and collects them in a list
            # this is how the code knows which faces share which edge
            if edge_face_map.Contains(edge):
                adjacent_face_list = edge_face_map.FindFromKey(edge)
                for f in adjacent_face_list:
                    # cast it back to type 'face' so I can use it as a face
                    adjacent_faces.append(topods.Face(f))

            # Case 1: interior edge, two adjacent faces
            # add graph edge between their nodes, weighted by the physical length of the shared curve
            # this is what tells METIS later how expensive it would be to cut this connection
            # longer shared boundaries = more inter-rank communication if cut
            if len(adjacent_faces) == 2:
                f0_hash = adjacent_faces[0].TShape()
                f1_hash = adjacent_faces[1].TShape()
                if f0_hash in face_index and f1_hash in face_index:
                    # add a connection between node f0 and f1 with this edge length as the weight
                    G.add_edge(face_index[f0_hash],
                               face_index[f1_hash],
                               weight=get_edge_length(edge))

            # Case 2: if only one adjacent face, it sits on the outer boundary of the model
            # no graph edge is added because there is nothing to connect to
            # face is flagged as a boundary face (during meshing: exposed edges with no geometric neighbour)
            elif len(adjacent_faces) == 1:
                f_hash = adjacent_faces[0].TShape()
                if f_hash in face_index:
                    G.nodes[face_index[f_hash]]['is_boundary'] = True

            # Case 3: non-manifold edge - more than two faces share this edge
            # should not happen in a valid solid CAD model — connect all pairs defensively
            elif len(adjacent_faces) > 2:
                print(f"Warning: non-manifold edge with {len(adjacent_faces)} faces")
                for idx1 in range(len(adjacent_faces)):
                    for idx2 in range(idx1+1, len(adjacent_faces)):
                        f0_hash = adjacent_faces[idx1].TShape()
                        f1_hash = adjacent_faces[idx2].TShape()
                        if f0_hash in face_index and f1_hash in face_index:
                            G.add_edge(face_index[f0_hash],
                                       face_index[f1_hash],
                                       weight=get_edge_length(edge))

        edge_explorer.Next()

    return G


# convert NetworkX graph into CSR array for pymetis
def build_csr(G):
    xadj = [0]
    adjacency_list = []
    weight_list = []
    for i in range(len(G.nodes())):
        for y in G.neighbors(i):
            adjacency_list.append(y)
            weight_list.append(G[i][y]['weight'])
        xadj.append(len(adjacency_list))
    return xadj, adjacency_list, weight_list


# scales float edge weights into integers and stores partition IDs on nodes
def partition_graph(G, n_partitions=4):
    """
    Partition graph using METIS.
    Edge weights: shared curve lengths (minimise communication cost).
    Node weights: face areas (balance computational load).
    Note: previously only edge weights were used. Adding vweights aligns
    the partitioning objective with the load balance metric.
    """
    print(f"Partitioning into {n_partitions} partitions")

    xadj, adjacency_list, weight_list = build_csr(G)

    # edge weights — shared curve lengths scaled to integers
    # METIS minimises total cut edge weight (communication cost)
    int_eweights = [int(w * 1000) for w in weight_list]

    # node weights — face areas scaled to integers
    # METIS balances total node weight per partition (computational load)
    # scaling: divide by min area so smallest face has weight 1
    # this prevents integer overflow for large area differences
    areas = [G.nodes[i]['area'] for i in range(len(G.nodes()))]
    min_area = max(min(areas), 1e-6)   # guard against zero area
    int_vweights = [max(1, int(a / min_area)) for a in areas]

    result = pymetis.part_graph(n_partitions,
                                xadj=xadj,
                                adjncy=adjacency_list,
                                eweights=int_eweights,
                                vweights=int_vweights)

    for node in G.nodes():
        # index is node ID, value is the partition it was assigned to,
        # this is stored as an attribute of the graph
        G.nodes[node]['partition_id'] = int(result.vertex_part[node])

    return G, result


def visualize_graph(G):
    # use actual face centroids as 2D positions (drop z for now)
    pos = {node: (data['centroid'][0], data['centroid'][2])
           for node, data in G.nodes(data=True)}
    coords = np.array(list(pos.values()))
    coords -= coords.min(axis=0)
    coords /= coords.max(axis=0)
    pos = {node: tuple(coords[i]) for i, node in enumerate(G.nodes())}
    nx.draw(G, pos, with_labels=True, node_color='lightblue',
            node_size=800, edge_color='gray', font_size=12, font_weight='bold')
    plt.savefig("wing_connectivity.png")


def visualize_partition(G, membership):
    # changed this as most of the wing lives in z
    pos = {node: (data['centroid'][0], data['centroid'][2])
           for node, data in G.nodes(data=True)}
    coords = np.array(list(pos.values()))
    coords -= coords.min(axis=0)
    coords /= coords.max(axis=0)
    pos = {node: tuple(coords[i]) for i, node in enumerate(G.nodes())}
    n_parts = len(set(membership))
    cmap = cm.get_cmap('tab20', n_parts)
    colors = [cmap(membership[node]) for node in G.nodes()]
    nx.draw(G, pos, with_labels=False, node_color=colors,
            node_size=50, edge_color='gray')
    plt.savefig('Complex_partitioned_adjacency.png')


def assess_load_balance(G, n_partitions):
    total_area = sum(data['area'] for _, data in G.nodes(data=True))
    ideal_load = total_area / n_partitions

    # initialise from actual partition IDs — METIS may not use all IDs
    # if range(n_partitions) is used, a KeyError occurs when a partition
    # number is skipped by METIS
    all_partition_ids = set(data['partition_id']
                            for _, data in G.nodes(data=True))
    partition_loads = {pid: 0.0 for pid in all_partition_ids}

    for _, data in G.nodes(data=True):
        partition_loads[data['partition_id']] += data['area']

    # count cut edges
    cut_edges = 0
    total_cut_weight = 0.0
    for u, v, data in G.edges(data=True):
        if G.nodes[u]['partition_id'] != G.nodes[v]['partition_id']:
            cut_edges += 1
            total_cut_weight += data['weight']

    print(f"--- Partition Quality Report ---")
    print(f"Total area:               {total_area:.4f}")
    print(f"Ideal load per partition: {ideal_load:.4f}")
    for pid, load in sorted(partition_loads.items()):
        imbalance = load / ideal_load
        print(f"  Partition {pid}: area={load:.4f}, imbalance={imbalance:.4f}")
    print(f"Cut edges: {cut_edges} / {G.number_of_edges()} total")
    print(f"Cut ratio: {cut_edges/G.number_of_edges():.4f}")
    print(f"Total cut weight: {total_cut_weight:.4f}")
    print(f"--------------------------------")


def detect_dominant_faces(G, n_partitions, threshold=2.0):
    """
    Flags faces whose area exceeds threshold * ideal_load.
    These faces will cause load imbalance and need subdivision
    via NURBS knot vector partitioning.
    """
    total_area = sum(data['area'] for _, data in G.nodes(data=True))
    ideal_load = total_area / n_partitions

    dominant_faces = []
    for node, data in G.nodes(data=True):
        if data['area'] > threshold * ideal_load:
            dominant_faces.append({
                'node': node,
                'area': data['area'],
                'centroid': data['centroid'],
                'partition_id': data['partition_id'],
                'ratio': data['area'] / ideal_load
            })

    print(f"\n--- Dominant Face Report ---")
    print(f"Ideal load per partition: {ideal_load:.2f}")
    print(f"Threshold: {threshold}x ideal = {threshold*ideal_load:.2f}")
    print(f"Dominant faces found: {len(dominant_faces)}")
    for f in dominant_faces:
        print(f"  Node {f['node']}: area={f['area']:.2f}, "
              f"ratio={f['ratio']:.2f}x ideal, "
              f"partition={f['partition_id']}, "
              f"centroid={f['centroid']}")
    print(f"----------------------------\n")
    return dominant_faces


if __name__ == "__main__":
    N_PARTITIONS = 3

    shape = load_step("/Users/jude/Documents/REMODEL/WP1_NURBS_Partition/data/NACA_complex_file.stp")
    G = brep_to_face_adjacency_graph(shape)
    G, result = partition_graph(G, n_partitions=N_PARTITIONS)

    assess_load_balance(G, n_partitions=N_PARTITIONS)
    detect_dominant_faces(G, n_partitions=N_PARTITIONS)
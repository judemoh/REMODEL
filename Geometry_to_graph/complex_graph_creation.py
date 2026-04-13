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
# from to_vtk import 

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
     # map eveery topological egde to its face, makes adjacency lookup direct direct during edge pass
    #  topological edge to ancestor faces
    topexp.MapShapesAndAncestors(shape, TopAbs_EDGE, TopAbs_FACE, edge_face_map)
    
    # --- Node pass: one node per face ---
    # computes per-face attributes (surface type, area, bounding boz midpoint as centroid prox and stores them alongside graph nodes)
    face_explorer = TopExp_Explorer(shape, TopAbs_FACE)
    face_index = {}
    i = 0
    # iterates through all faces 
    while face_explorer.More():
        # returns a generic Topods shape and downcasts it so face-sepcific APIs can be use d
        face = topods.Face(face_explorer.Current())
        # creates a geometry adaptor over undelrying surface
        adaptor = BRepAdaptor_Surface(face)
        # identifying geometric type from a given shape - used for graph attributes or filtering later on  
        surf_type = adaptor.GetType() 
        # compute compound properties of a global geometric system
       
    #    container that holds geometric properties 
        props = GProp_GProps() 
        # traverses face and fills the container - integrates over 2D surface (surface integral)
        brepgprop.SurfaceProperties(face, props)
        # retrieves the areas 
        area = props.Mass()
        # bounding box approxomation falls apart when geoemtry isnt symmetric 
        bbox = Bnd_Box()
        brepbndlib.Add(face, bbox)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
        centroid = ((xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2)
        print(area)
        print(centroid)
# i may not be deterministic across different runs or machines, global ID must be intrinsic to CAD model itself i.e a hash or a deterministic ordering based on geometric properties (like centroid)
        # graph with generic properties, can add as much info as required
        G.add_node(i, shape=face, surface_type=surf_type, area=area, 
           bbox=bbox, centroid=centroid, global_id=i, 
           is_boundary=False, partition_id=None)  
        
        face_index[face.TShape()] = i
        i += 1

        face_explorer.Next()
    # print(f"Extracted {G.number_of_nodes()} faces, {G.number_of_edges()} edges")
    # print(f"Average degree: {sum(d for _,d in G.degree())/G.number_of_nodes():.2f}")

  
    # --- Edge pass: connect faces sharing a topological edge ---
    # creates an explorer on the shape and the type of shapes to search
    edge_explorer = TopExp_Explorer(shape, TopAbs_EDGE)
    visited_edges = set()
    
    while edge_explorer.More():
        # current returns the current shape in exploration, casts it to the specific type that I need 
        # grab the current object and make it useable as an egde 
        edge = topods.Edge(edge_explorer.Current())
        # thsape reutnrs the underlying topologicalentity which is a pointr to the shared geometry - two topods objects that look different but represent the same edge will have the same tshape
        edge_hash = edge.TShape()
        print(f"Edge hash is {edge_hash}")
        
        # check that you have not traversed this edge before
        if edge_hash not in visited_edges:
            visited_edges.add(edge_hash)
            adjacent_faces = []
            
            if edge_face_map.Contains(edge):
                adjacent_face_list = edge_face_map.FindFromKey(edge)
                for f in adjacent_face_list:
                    # cast it back to type 'face' so I can use it as a face
                    adjacent_faces.append(topods.Face(f))
            
            if len(adjacent_faces) == 2:
                # normal interior edge - connect the two faces
                f0_hash = adjacent_faces[0].TShape()
                f1_hash = adjacent_faces[1].TShape()
                if f0_hash in face_index and f1_hash in face_index:
                    # add a connection between niode f0 and f1 with this edge length as the weight 
                    G.add_edge(face_index[f0_hash],
                            face_index[f1_hash],
                            weight=get_edge_length(edge))
                            
            elif len(adjacent_faces) == 1:
                # boundary edge - mark face as boundary
                f_hash = adjacent_faces[0].TShape()
                if f_hash in face_index:
                    G.nodes[face_index[f_hash]]['is_boundary'] = True
                    
            elif len(adjacent_faces) > 2:
                # non-manifold edge - connect all pairs and flag
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
        print(adjacent_face_list)
        print(adjacent_faces)
    return G

# convert network graph into CSR array for pymetis 
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

# scales flaot edge weight into integrers and stores partition IDs on nodes 
def partition_graph(G, n_partitions=None):
    if n_partitions is None:
        n_partitions = 4
    print(f"Partitioning into {n_partitions} partitions")
    
    xadj, adjacency_list, weight_list = build_csr(G)
    # metis only accepts integers
    int_eweights = [int(w * 1000) for w in weight_list]
    result = pymetis.part_graph(n_partitions,
                                xadj=xadj,
                                adjncy=adjacency_list,
                                eweights=int_eweights)
    for node in G.nodes():
        # index is node ID, value is the partition it was assigned to , this is stored as an attribute of the graph
        G.nodes[node]['partition_id'] = int(result.vertex_part[node])
    
    # print(f"Partitioning complete")
    return G, result


def visualize_graph(G):
    # # use actual face centroids as 2D positions (drop z for now)
    # pos = {node: (data['centroid'][0], data['centroid'][1]) 
    #    for node, data in G.nodes(data=True)}
    pos = {node: (data['centroid'][0], data['centroid'][2]) 
       for node, data in G.nodes(data=True)}
    coords = np.array(list(pos.values()))
    coords -= coords.min(axis=0)
    coords /= coords.max(axis=0)
    pos = {node: tuple(coords[i]) for i, node in enumerate(G.nodes())}

    nx.draw(G, pos,
            with_labels=True,
            node_color='lightblue',
            node_size=800,
            edge_color='gray',
            font_size=12,
            font_weight='bold')
    
    plt.savefig("wing_connectivity.png")

import matplotlib.cm as cm
import numpy as np
def visualize_partition(G, membership):
    # changed this as most of the wing lives in z 
    # pos = {node: (data['centroid'][0], data['centroid'][1]) 
    #    for node, data in G.nodes(data=True)}
    pos = {node: (data['centroid'][0], data['centroid'][2]) 
       for node, data in G.nodes(data=True)}
    coords = np.array(list(pos.values()))
    coords -= coords.min(axis=0)
    coords /= coords.max(axis=0)
    pos = {node: tuple(coords[i]) for i, node in enumerate(G.nodes())}

    n_parts = len(set(membership))
    cmap = cm.get_cmap('tab20', n_parts)
    colors = [cmap(membership[node]) for node in G.nodes()]
    
    nx.draw(G, pos,
            with_labels=False,
            node_color=colors,
            node_size=50,
            edge_color='gray')
    
    plt.savefig('Complex_partitioned_adjacency.png')

def assess_load_balance(G, n_partitions):
    total_area = sum(data['area'] for _, data in G.nodes(data=True))
    ideal_load = total_area / n_partitions
    
    partition_loads = {i: 0.0 for i in range(n_partitions)}
    for _, data in G.nodes(data=True):
        partition_loads[data['partition_id']] += data['area']
    
    # count cut edges
    cut_edges = 0
    total_cut_weight = 0.0
    for u, v, data in G.edges(data=True):
        if G.nodes[u]['partition_id'] != G.nodes[v]['partition_id']:
            cut_edges += 1
            total_cut_weight += data['weight']
    
    # print(f"--- Partition Quality Report ---")
    # print(f"Total area: {total_area:.4f}")
    # print(f"Ideal load per partition: {ideal_load:.4f}")
    # for pid, load in partition_loads.items():
    #     imbalance = load / ideal_load
    #     print(f"Partition {pid}: load = {load:.4f}, imbalance = {imbalance:.4f}")
    # print(f"Cut edges: {cut_edges} / {G.number_of_edges()} total")
    # print(f"Total cut weight: {total_cut_weight:.4f}")
    # print(f"Cut ratio: {cut_edges/G.number_of_edges():.4f}")
    # print(f"--------------------------------")

def detect_dominant_faces(G, n_partitions, threshold=2.0):
    """
    Flags faces whose area exceeds threshold * ideal_load.
    These faces will cause load imbalance and need subdivision.
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
    
    # print(f"\n--- Dominant Face Report ---")
    # print(f"Ideal load per partition: {ideal_load:.2f}")
    # print(f"Threshold: {threshold}x ideal = {threshold*ideal_load:.2f}")
    # print(f"Dominant faces found: {len(dominant_faces)}")
    # for f in dominant_faces:
    #     print(f"  Node {f['node']}: area={f['area']:.2f}, "
    #           f"ratio={f['ratio']:.2f}x ideal, "
    #           f"partition={f['partition_id']}, "
    #           f"centroid={f['centroid']}")
    # print(f"----------------------------\n")
    # return dominant_faces

import pymetis
import matplotlib.pyplot as plt

# if __name__ == "__main__":
shape = load_step("NACA_complex_file.step")
G = brep_to_face_adjacency_graph(shape)
G, result = partition_graph(G)
# export_partitioned_vtk(G, "partitioned_wing.vtk") 
# visualize_graph(G)
# visualize_partition(G, result.vertex_part)
# print(f"Nodes (faces): {G.number_of_nodes()}")
# print(f"Edges (adjacencies): {G.number_of_edges()}")
# print(f"Degree of each node (expected 4 for all):")
# for node, deg in G.degree():
#     print(f"  Face {node}: degree {deg}")
# for node, data in G.nodes(data=True):
#     print(f"Face {node}: area = {data['area']:.4f}")
#     print(f"Face {node}: area = {data['centroid']}")
# for u, v, data in G.edges(data=True):
#     print(f"Edge ({u},{v}): length = {data['weight']:.4f}")
n_parts = len(set(result.vertex_part))
assess_load_balance(G, n_partitions=n_parts)

# call it after building the graph

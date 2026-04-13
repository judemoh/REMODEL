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


#this code creates a graph from a cube geometry
def get_edge_length(edge):
    # uniform interface to the underlying geoometry, edge is a topological citizen, and adaptor translates it to a geometric entity
    curve_adaptor = BRepAdaptor_Curve(edge)
    # integratie from edge to edge - arch length method, samples points along curve and approximates the integral
    # needs the adaptror to evaluate the curve at arebitrary parameter values
    length = GCPnts_AbscissaPoint.Length(curve_adaptor)
    return length


def load_step(filename):
    reader = STEPControl_Reader()
    reader.ReadFile(filename)
    # transfers step to BREP
    reader.TransferRoots()
    # one shape containing the whole geometry
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
    # itrates through all faces 
    while face_explorer.More():
        face = topods.Face(face_explorer.Current())
        adaptor = BRepAdaptor_Surface(face)
        # identifying geometric type from a given shape 
        surf_type = adaptor.GetType() 
        # compute compound properties of a global geometric system
        props = GProp_GProps()
        brepgprop.SurfaceProperties(face, props)
        area = props.Mass()

        bbox = Bnd_Box()
        brepbndlib.Add(face, bbox)
        xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()
        centroid = ((xmin+xmax)/2, (ymin+ymax)/2, (zmin+zmax)/2)
# i may not be deterministic across different runs or machines, global ID must be intrinsic to CAD model itself i.e a hash or a deterministic ordering based on geometric properties (like centroid)
        # graph with generic properties, can add as much info as required
        G.add_node(i, shape=face, surface_type=surf_type, area=area, bbox =bbox , centroid=centroid, global_id=i)
        face_index[face.TShape()] = i
        i += 1

        face_explorer.Next()

  
    # --- Edge pass: connect faces sharing a topological edge ---
    edge_explorer = TopExp_Explorer(shape, TopAbs_EDGE)
    visited_edges = set()
    
    while edge_explorer.More():
        edge = topods.Edge(edge_explorer.Current())
        edge_hash = edge.TShape()
        
        # avoid processing the same edge twice
        if edge_hash not in visited_edges:
            visited_edges.add(edge_hash)
            adjacent_faces= []
            # look up which faces contain this edge
            if edge_face_map.Contains(edge):
                adjacent_face_list = edge_face_map.FindFromKey(edge)
                
            for f in adjacent_face_list:
                adjacent_faces.append(topods.Face(f))
                
                # an interior edge borders exactly 2 faces
                if len(adjacent_faces) == 2:
                    f0_hash = adjacent_faces[0].TShape()
                    f1_hash = adjacent_faces[1].TShape()
                    
                    if f0_hash in face_index and f1_hash in face_index:
                        G.add_edge(face_index[f0_hash], 
                                   face_index[f1_hash])
        G.add_edge(face_index[f0_hash], 
        face_index[f1_hash],
        weight=get_edge_length(edge))
        edge_explorer.Next()
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
def partition_graph(G, n_partitions=2):
    xadj, adjacency_list, weight_list = build_csr(G)
    int_eweights = [int(w * 1000) for w in weight_list]
    result = pymetis.part_graph(n_partitions,
                                xadj=xadj,
                                adjncy=adjacency_list,
                                eweights=int_eweights)
    for node in G.nodes():
        # assign partition label to every graph node after partitioning
        G.nodes[node]['partition_id'] = int(result.vertex_part[node])
        print(G)
    return G, result


def visualize_graph(G):
    # use actual face centroids as 2D positions (drop z for now)
    pos = {node: (data['centroid'][0], data['centroid'][1]) 
           for node, data in G.nodes(data=True)}
    
    nx.draw(G, pos,
            with_labels=True,
            node_color='lightblue',
            node_size=800,
            edge_color='gray',
            font_size=12,
            font_weight='bold')
    
    plt.title("Face Adjacency Graph - Cube")
    plt.savefig("cube_connectivity.png")
def visualize_partition(G, membership):
    pos = {node: (data['centroid'][0], data['centroid'][1]) 
           for node, data in G.nodes(data=True)}
    
    colors = ['lightblue' if membership[node] == 0 else 'lightcoral' 
              for node in G.nodes()]
    
    nx.draw(G, pos,
            with_labels=True,
            node_color=colors,
            node_size=800,
            edge_color='gray',
            font_size=12,
            font_weight='bold')
    
    plt.title(f"Partitioned Face Adjacency Graph - 2 partitions")
    plt.savefig('Partitioned_adjacency.png')

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
    
    print(f"--- Partition Quality Report ---")
    print(f"Total area: {total_area:.4f}")
    print(f"Ideal load per partition: {ideal_load:.4f}")
    for pid, load in partition_loads.items():
        imbalance = load / ideal_load
        print(f"Partition {pid}: load = {load:.4f}, imbalance = {imbalance:.4f}")
    print(f"Cut edges: {cut_edges} / {G.number_of_edges()} total")
    print(f"Total cut weight: {total_cut_weight:.4f}")
    print(f"Cut ratio: {cut_edges/G.number_of_edges():.4f}")
    print(f"--------------------------------")

import pymetis
import matplotlib.pyplot as plt

if __name__ == "__main__":
    shape = load_step("cube.step")
    G = brep_to_face_adjacency_graph(shape)
    G, result = partition_graph(G)
    visualize_graph(G)
    visualize_partition(G, result.vertex_part)
    print(f"Nodes (faces): {G.number_of_nodes()}  -- expected 6")
    print(f"Edges (adjacencies): {G.number_of_edges()}  -- expected 12")
    print(f"Degree of each node (expected 4 for all):")
    for node, deg in G.degree():
        print(f"  Face {node}: degree {deg}")
    for node, data in G.nodes(data=True):
        print(f"Face {node}: area = {data['area']:.4f}")
        print(f"Face {node}: area = {data['centroid']}")
    for u, v, data in G.edges(data=True):
        print(f"Edge ({u},{v}): length = {data['weight']:.4f}")
    G, result = partition_graph(G)
    assess_load_balance(G, n_partitions=2)

# call it after building the graph

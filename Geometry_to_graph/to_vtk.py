from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Core.BRep import BRep_Tool
from OCC.Core.TopAbs import TopAbs_FACE
from OCC.Core.TopExp import TopExp_Explorer
from OCC.Core.TopoDS import topods
import numpy as np

def export_partitioned_vtk(G, filename="partitioned_wing.vtk"):
    all_points = []
    all_triangles = []
    all_partition_ids = []
    all_areas = []
    offset = 0

    for node, data in G.nodes(data=True):
        face = data['shape']
        partition_id = data['partition_id']
        area = data['area']

        # tessellate the face
        mesh = BRepMesh_IncrementalMesh(face, 0.1)
        mesh.Perform()

        location = face.Location()
        triangulation = BRep_Tool.Triangulation(face, location)

        if triangulation is None:
            continue

        # extract vertices
        n_nodes = triangulation.NbNodes()
        for i in range(1, n_nodes + 1):
            p = triangulation.Node(i)
            all_points.append([p.X(), p.Y(), p.Z()])

        # extract triangles
        n_triangles = triangulation.NbTriangles()
        for i in range(1, n_triangles + 1):
            tri = triangulation.Triangle(i)
            n1, n2, n3 = tri.Get()
            all_triangles.append([offset + n1 - 1, 
                                   offset + n2 - 1, 
                                   offset + n3 - 1])
            all_partition_ids.append(partition_id)
            all_areas.append(area)

        offset += n_nodes

    # write VTK
    points = np.array(all_points)
    triangles = np.array(all_triangles)

    with open(filename, 'w') as f:
        f.write("# vtk DataFile Version 3.0\n")
        f.write("Partitioned Wing Geometry\n")
        f.write("ASCII\n")
        f.write("DATASET POLYDATA\n")

        # points
        f.write(f"POINTS {len(points)} float\n")
        for p in points:
            f.write(f"{p[0]:.6f} {p[1]:.6f} {p[2]:.6f}\n")

        # triangles
        f.write(f"POLYGONS {len(triangles)} {len(triangles) * 4}\n")
        for t in triangles:
            f.write(f"3 {t[0]} {t[1]} {t[2]}\n")

        # cell data - partition ID and area
        f.write(f"CELL_DATA {len(triangles)}\n")

        f.write("SCALARS partition_id int 1\n")
        f.write("LOOKUP_TABLE default\n")
        for pid in all_partition_ids:
            f.write(f"{pid}\n")

        f.write("SCALARS face_area float 1\n")
        f.write("LOOKUP_TABLE default\n")
        for a in all_areas:
            f.write(f"{a:.4f}\n")

    print(f"Written to {filename}")
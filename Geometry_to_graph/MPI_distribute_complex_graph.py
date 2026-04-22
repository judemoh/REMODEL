from mpi4py import MPI
import numpy as np




# This code is a parallel assignment system
# Tells each rank which faces it is responsible
# No rank recieves any geometry or evaluates anything 
# the next step would be to upload actual NURBS geometry on rank 0, which would carry actual geometry rather than metadata





# launches 4 identical copies of the script simiultaneously - each copy executres these three lines
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print(f"Rank {rank}: started", flush=True)

#rank 0: load STEP, build graph, partition, package arrays
if rank == 0:
    # read the step file via opencascade
    # do all the operations from the previous file: load step, build the rgaph, partition graph, assess load balance, and detect dominant faces with flagging
    from complex_graph_creation import (load_step, brep_to_face_adjacency_graph, 
                                 partition_graph, assess_load_balance, detect_dominant_faces)
    from to_vtk import export_partitioned_vtk
    print(f"Rank {rank}: imports done", flush=True)
    
    # load and partition with n_partitions matching number of MPI ranks
    shape = load_step("/Users/jude/Documents/REMODEL/WP1_NURBS_Partition/data/NACA_complex_file.stp")
    G, result = partition_graph(brep_to_face_adjacency_graph(shape), 
                                n_partitions=size)
    
    n_parts = len(set(result.vertex_part))
    assess_load_balance(G, n_partitions=n_parts)
    detect_dominant_faces(G, n_partitions=n_parts)
    # Then packages the results into five flat numpy arrays — one  fparvalue per face
    # package metadata into numpy arrays
    n_faces = G.number_of_nodes()
    global_ids    = np.array([G.nodes[i]['global_id']      for i in range(n_faces)])
    partition_ids = np.array([G.nodes[i]['partition_id']    for i in range(n_faces)])
    centroids     = np.array([G.nodes[i]['centroid']        for i in range(n_faces)])
    areas         = np.array([G.nodes[i]['area']            for i in range(n_faces)])
    is_boundary   = np.array([G.nodes[i]['is_boundary']     for i in range(n_faces)])
    
    print(f"Rank {rank}: graph built, {n_faces} faces, "
          f"{n_parts} partitions", flush=True)
else:
    global_ids    = None
    partition_ids = None
    centroids     = None
    areas         = None
    is_boundary   = None
    
    print(f"Rank {rank}: waiting for broadcast", flush=True)

# synchronize before broadcast - all four ranks much reach this face before the metadata can be dsitributed - barrier makes them finish
comm.Barrier()

# broadcast all arrays from rank 0 to all ranks - all ranks have copies for all five arrays containing metadata
global_ids    = comm.bcast(global_ids,    root=0)
partition_ids = comm.bcast(partition_ids, root=0)
centroids     = comm.bcast(centroids,     root=0)
areas         = comm.bcast(areas,         root=0)
is_boundary   = comm.bcast(is_boundary,   root=0)

print(f"Rank {rank}: broadcast received", flush=True)

# each rank filters its own face via a boolean mask, whcih selects only the faces assigned to it
# No communication happens here — each rank just filters the same array it already has locally using its own rank number
owned_mask       = partition_ids == rank
owned_global_ids = global_ids[owned_mask]
owned_centroids  = centroids[owned_mask]
owned_areas      = areas[owned_mask]
owned_boundary   = is_boundary[owned_mask]

# The second barrier synchronises before the final print so the output isn't mixed periodically 
comm.Barrier()
print(f"Rank {rank} owns {len(owned_global_ids)} faces, "
      f"total area = {owned_areas.sum():.2f}, "
      f"boundary faces = {owned_boundary.sum()}", flush=True)
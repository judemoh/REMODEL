from mpi4py import MPI
import numpy as np

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

print(f"Rank {rank}: started", flush=True)

if rank == 0:
    from complex_graph_creation import load_step, brep_to_face_adjacency_graph, partition_graph, assess_load_balance, get_edge_length
    
    print(f"Rank {rank}: imports done", flush=True)
    shape = load_step("cube.step")
    G, result = partition_graph(brep_to_face_adjacency_graph(shape))
    assess_load_balance(G, n_partitions=2)    
    n_faces = G.number_of_nodes()
    global_ids    = np.array([G.nodes[i]['global_id']    for i in range(n_faces)])
    # print(global_ids)
    partition_ids = np.array([G.nodes[i]['partition_id'] for i in range(n_faces)])
    # print(partition_ids)
    centroids     = np.array([G.nodes[i]['centroid']     for i in range(n_faces)])
    # print(centroids)
    areas         = np.array([G.nodes[i]['area']         for i in range(n_faces)])
    # print(areas)
    
    print(f"Rank {rank}: graph built and partitioned", flush=True)
else:
    global_ids    = None
    partition_ids = None
    centroids     = None
    areas         = None
    
    print(f"Rank {rank}: waiting for broadcast", flush=True)

comm.Barrier()
print(f"Rank {rank}: reached barrier", flush=True)
# distribute data from one source process to all other processes
global_ids    = comm.bcast(global_ids,    root=0)
partition_ids = comm.bcast(partition_ids, root=0)
centroids     = comm.bcast(centroids,     root=0)
areas         = comm.bcast(areas,         root=0)

print(f"Rank {rank}: broadcast received", flush=True)

owned_mask       = partition_ids == rank
owned_global_ids = global_ids[owned_mask]
owned_centroids  = centroids[owned_mask]
owned_areas      = areas[owned_mask]

comm.Barrier()
print(f"Rank {rank} owns {len(owned_global_ids)} faces: {owned_global_ids}", flush=True)
print(f"Rank {rank} centroids: {owned_centroids}", flush=True)
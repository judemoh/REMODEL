#include "mpi.h"
#include <stdio.h>
// int main(int argc, char** argv) {
//     int rank, size;
//     MPI_Init(&argc, &argv);
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//     MPI_Comm_size(MPI_COMM_WORLD, &size);
//     printf("Hello from process %d of %d\n", rank, size);
//     MPI_Finalize();
//     return 0;
// }
//no data exchange
// no coordiantion between processes

//The following sums contributions from each process 
// each process determines its rank then knitializes the value that will contribute to the sum

// int main(int argc, char** argv) {
//     int rank, size, valIn, valOut;
//     MPI_Init(&argc, &argv);
//     MPI_Status status;
//     //determine rank and size of all cooperating processes
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank); //rank of calling process, value between 0 and size-1
//     MPI_Comm_size(MPI_COMM_WORLD, &size);
//     // pick a simple value to contribute to the sum
//     valIn = rank;
//     //receive partial sum from the right processes (this is the sum from i = rank+1 to size-1)
//     if (rank < size - 1) {
//         // address of data, number of items, type of data, source rank, tag (a way to use a single integer with the data, this is not needed here), communicator, status
//         MPI_Recv(&valOut, 1, MPI_INT, rank + 1, 0, MPI_COMM_WORLD, &status);
//         valIn += valOut; // add the received value to the local value
//     }
//     // send the partial sum to the left (rank-1) process 
//     if (rank > 0) {
//         MPI_Send(&valIn, 1, MPI_INT, rank - 1, 0, MPI_COMM_WORLD);
//     } else { // rank 0 will have the total sum
//         printf("Total sum of ranks is: %d\n", valIn);
//     }
//     MPI_Finalize();
//     return 0;

// }

//time taken is proportion to number of processes as each process must wait for its neighbor to gfinish before it can send data
//eg most operations on PDe meshes are best done via point-to-point comm  as data exhcnagers are between pairs of neighboring processes 
//and this matches the point-to-point communication model
// int main(int argc, char** argv) {
//     int rank, size, valIn, valOut;
//     MPI_Init(&argc, &argv);
//     MPI_Status status;
//     //determine rank and size of all cooperating processes
//     MPI_Comm_rank(MPI_COMM_WORLD, &rank); //rank of calling process, value between 0 and size-1
//     valIn = rank;
//     //receive partial sum from the right processes (this is the sum from i = rank+1 to size-1), reduce the process to 0 by summing the contribtions from all processes
//     MPI_Reduce(&valIn, &valOut, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
//     // send the partial sum to the left (rank-1) process 
//     if (rank ==0) {
//         printf("Total sum of ranks is: %d\n", valOut);
//     }
//     MPI_Finalize();
//     return 0;
// }

//parallel I/O
//performance requriement: creation of file access in parallel 
//reading files is like sending message to the system
//parallel IO needs specify group of processes, also need a communicator, and sometimes a way to describe where the data is on the disk


// //Declarations
// MPI_File fh; //file handle
// int val;
// //Start MPI
// MPI_Init (&argc, &argv);

// //open file: communicator, input file, access style, and a parameter to pass additional data
// MPI_File_open(MPI_COMM_WORLD, "input.dat", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);

// //use all proccesses to reaad from the file
// MPI_File_read_all(fh, &val, 1, MPI_INT, MPI_STATUS_IGNORE);

// //close the file when finished
// MPI_File_Close(&fh);

//use of collective IO to combine with file views to write data from many processors to a single file that provides natural order for the data
//each processor writes ARRAY_Size double precision values to the file, ordered by the MPI rank of the process


#define ARRAY_SIZE 1000

//Declarations 
MPI_File fh;
int rank;
int solution[ARRAY_SIZE];
//Start MPI
MPI_Init (&argc, &argv);
//Open the file for reading only
MPI_File_open(MPI_COMM_WORLD, "output.dat", MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
//Define where each process writes in the file
MPI_Comm_rank(MPI_COMM_WORLD, &rank);
MPI_File_set_view(fh, rank*ARRAY_SIZE*sizeof(double), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

//perform the file write
MPI_File_write_all(fh, solution, ARRAY_SIZE, MPI_DOUBLE, MPI_STATUS_IGNORE);
//close the file when finished
MPI_File_close(&fh);
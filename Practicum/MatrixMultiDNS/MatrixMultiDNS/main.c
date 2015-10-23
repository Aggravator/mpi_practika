#include <mpi.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define iD 0
#define jD 1
#define kD 2

MPI_Comm mesh3d, mesh_ik, ring_j, ring_i,mesh_ij, ring_k;
int rank,size,mesh_size,matrix_size,block_size;
int* counts,*displs;
int mesh_coords[3];
MPI_Datatype resizedtype;

void createTopology(){
	int dims[3], i, periods[3];
	int keep_dims[3];

  
	int dimension = 3;
	for( i = 0; i < dimension; i++ ) {
		dims[i] = mesh_size;
		periods[i] = 1;
	}

	MPI_Cart_create(MPI_COMM_WORLD, dimension, dims, periods, 0, &mesh3d);
	MPI_Cart_coords(mesh3d, rank, dimension, mesh_coords);


	keep_dims[iD] = 1; // i
	keep_dims[jD] = 1; // j
	keep_dims[kD] = 0; // k
	MPI_Cart_sub(mesh3d, keep_dims, &mesh_ij);

	keep_dims[iD] = 1;
	keep_dims[jD] = 0;
	keep_dims[kD] = 1;
	MPI_Cart_sub(mesh3d, keep_dims, &mesh_ik);

	keep_dims[iD] = 0;
	keep_dims[jD] = 0;
	keep_dims[kD] = 1;
	MPI_Cart_sub(mesh3d, keep_dims, &ring_k);

	keep_dims[iD] = 0;
	keep_dims[jD] = 1;
	keep_dims[kD] = 0;
	MPI_Cart_sub(mesh3d, keep_dims, &ring_j);

	keep_dims[iD] = 1;
	keep_dims[jD] = 0;
	keep_dims[kD] = 0;
	MPI_Cart_sub(mesh3d, keep_dims, &ring_i);
}
void distribute(float* tmp_array, float *A, float *B) {
  int i, j;
  MPI_Status status;

  counts = (int *)malloc(sizeof(int) * size);
  displs = (int *)malloc(sizeof(int) * size);
  
  for( i = 0; i <size; i++ ) {
    counts[i] = 1;
  }
  
  for( i = 0; i < mesh_size; i++ ) {
    for(j = 0; j < mesh_size; j++ ) {
      displs[i * mesh_size + j] = j + i * matrix_size; 
    }
  }

  MPI_Datatype type;
  int sizes[2]    = {matrix_size,matrix_size};  
  int subsizes[2] = {block_size,block_size};  
  int starts[2]   = {0,0}; 
  
  MPI_Type_create_subarray(2, sizes, subsizes, starts, MPI_ORDER_C, MPI_FLOAT, &type);  

  MPI_Type_create_resized(type, 0, block_size * sizeof(float), &resizedtype);
  MPI_Type_commit(&resizedtype);


  if( mesh_coords[2] == 0 ) {
    MPI_Scatterv(tmp_array, counts, displs, resizedtype, A, block_size * block_size, MPI_FLOAT, 0, mesh_ij); 
    MPI_Scatterv(tmp_array, counts, displs, resizedtype, B, block_size * block_size, MPI_FLOAT, 0, mesh_ij); 
  }
  
  MPI_Type_free(&type);
}
void block_mult(float *A, float *B, float *C, int _q) {
  int i, j, k;
  for( i = 0; i < _q; i++ ) {
    for( j = 0; j < _q; j++ ) { 
      for( k = 0; k < _q; k++ ) {
        C[i * _q + j] += A[i * _q + k] * B[k * _q + j];
      }
    }
  }
}

void  doo(float *A, float *B, float *C) {
  int i;
  float *AB;
  AB = (float *) malloc(sizeof(float) * block_size * block_size);
  memset(AB,0,sizeof(float) * block_size * block_size);
  block_mult(A, B, AB, block_size);
  fflush(stdout);
  MPI_Reduce(AB, C, block_size * block_size, MPI_FLOAT, MPI_SUM, 0, ring_k);
  free(AB);
}
float fRand(float fMin, float fMax){
    float f = (float)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}
void generateArry(float *arr,int size){
	int i,j;
	for( i = 0; i < size; i++ ) {
		for( j = 0; j < size; j++ ) {
			arr[i * size + j] = fRand(-200,200);
		}
	}
}
int main(int argc, char **argv){
	double t_start, t_end;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	for(mesh_size=1;mesh_size*mesh_size*mesh_size<size;++mesh_size);

	if(mesh_size*mesh_size*mesh_size!=size){
		MPI_Finalize();
		return 1;
	}

	createTopology();

	matrix_size=atoi(argv[1]);
	block_size=matrix_size/mesh_size;

	float *A, *B, *C, *result;
    
	A = (float *) malloc(sizeof(float) * block_size * block_size);
	B = (float *) malloc(sizeof(float) * block_size * block_size);
	C = (float *) malloc(sizeof(float) * block_size * block_size);
	result = (float *) malloc(sizeof(float) * matrix_size * matrix_size);

	float *temp;
	if(rank==0){
		temp=(float *)malloc(sizeof(float) * (matrix_size * matrix_size));
	}
	t_start=MPI_Wtime();
	distribute(temp, A, B);
	if(rank==0)free(temp);
	
	MPI_Bcast(B, block_size * block_size, MPI_FLOAT, 0, ring_k);
	MPI_Bcast(B, block_size * block_size, MPI_FLOAT, mesh_coords[2], ring_i);
	MPI_Bcast(A, block_size * block_size, MPI_FLOAT, 0, ring_k);
	MPI_Bcast(A,block_size * block_size, MPI_FLOAT, mesh_coords[2], ring_j);
	
	fflush(stdout);

	doo(A, B, C);

	MPI_Gatherv(C, block_size * block_size, MPI_FLOAT, result, counts, displs, resizedtype, 0, mesh_ij); 
	
	t_end=MPI_Wtime();
	if(rank==0) printf("Time=%lf",t_end-t_start);
    free(A);
    free(B);
    free(C);
    free(result);

	MPI_Comm_free(&mesh_ij);
    MPI_Comm_free(&mesh_ik);
    MPI_Comm_free(&ring_k);
    MPI_Comm_free(&ring_i);
    MPI_Comm_free(&ring_j);
	MPI_Comm_free(&mesh3d);
	MPI_Finalize();
	return 0;
}
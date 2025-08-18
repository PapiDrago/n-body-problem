#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <mpi.h>

typedef struct{
	double x,y,z;
}vector;

MPI_Datatype MPI_VECTOR;

double GravConstant = 39.47;
double dt = 0.01;
int N = 100; // number of iterations
double *masses;
vector *positions, *velocities, *accelerations;

double rand_uniform(unsigned int *seed, double min, double max) {
    return min + (max - min) * ((double)rand_r(seed) / RAND_MAX);
}

unsigned int global_start_index(int r, int bodies, int rank){
    unsigned int start_index;
    if (rank < r) {
    		return rank*bodies;
    } else {
    		return r*(bodies+1) + (rank - r)*bodies;
      }

}

void initiateSystem(vector *local_positions, double *local_masses, vector *velocities, vector *accelerations, int n_body, unsigned int start_index){	
	for (int i = 0; i < n_body; i++) {
          unsigned int seed = start_index + i;

          local_positions[i].x = rand_uniform(&seed, -1.0, 1.0);
          local_positions[i].y = rand_uniform(&seed, -1.0, 1.0);
          local_positions[i].z = rand_uniform(&seed, -1.0, 1.0);

          velocities[i].x = rand_uniform(&seed, -0.5, 0.5);
          velocities[i].y = rand_uniform(&seed, -0.5, 0.5);
          velocities[i].z = rand_uniform(&seed, -0.5, 0.5);

          local_masses[i] = 1.0;

          accelerations[i].x = 0;
          accelerations[i].y = 0;
          accelerations[i].z = 0;
        }
}

void computeAccelerations(int bodies, int own_bodies, int start_global_index, vector *global_positions, double *global_masses, vector *accelerations){
	int i,j;
        const double epsilon = 1e-5;
	for(i=0;i<own_bodies;i++){
		accelerations[i].x = 0;
		accelerations[i].y = 0;
		accelerations[i].z = 0;
		for(j=0;j<bodies;j++){
		        unsigned int global_index = start_global_index + i;
			if(global_index!=j){
				
				double dx = global_positions[j].x - global_positions[global_index].x;
				double dy = global_positions[j].y - global_positions[global_index].y;
				double dz = global_positions[j].z - global_positions[global_index].z;

				double dist = sqrt(dx * dx + dy * dy + dz * dz) + epsilon;

				double invDistCubed = 1.0 / (dist * dist * dist);

				double scalar = GravConstant * global_masses[j] * invDistCubed;

				accelerations[i].x += scalar * dx;
				accelerations[i].y += scalar * dy;
				accelerations[i].z += scalar * dz;
			}
		}
	}
}

void computeVelocities(int own_bodies, vector *accelerations, vector *velocities, double dt) {
	for (int i = 0; i < own_bodies; i++) {
		velocities[i].x = velocities[i].x + dt * accelerations[i].x;
                velocities[i].y = velocities[i].y + dt * accelerations[i].y;
                velocities[i].z = velocities[i].z + dt * accelerations[i].z;
                }
}

void computePositions(vector *local_positions, vector *velocities, int own_bodies, double dt) {
	for (int i = 0; i < own_bodies; i++) {
		local_positions[i].x = local_positions[i].x + dt * velocities[i].x;
		local_positions[i].y = local_positions[i].y + dt * velocities[i].y;
		local_positions[i].z = local_positions[i].z + dt * velocities[i].z;
		}
}


void simulate(int global_bodies, int own_bodies, unsigned int start_global_index, double dt, vector *global_positions, double *global_masses, vector *accelerations, vector *velocities, vector *local_positions, double* t4, double* t5, double* t6, double* cumul_time_accel, double* cumul_time_velpos){
	*t4 = MPI_Wtime();
	computeAccelerations(global_bodies, own_bodies, start_global_index, global_positions, global_masses, accelerations);
	*t5 = MPI_Wtime();
	computeVelocities(own_bodies, accelerations, velocities, dt);
	computePositions(local_positions, velocities, own_bodies, dt);
	*t6 = MPI_Wtime();
	*cumul_time_accel += *t5-*t4; 
	*cumul_time_velpos += *t6-*t5;
}

void logPositions(FILE *fp, int global_bodies, vector *global_positions) {
    for (int i = 0; i < global_bodies; i++) {
        fprintf(fp, "%lf,%lf,%lf", global_positions[i].x, global_positions[i].y, global_positions[i].z);
        if (i < global_bodies - 1)
            fprintf(fp, ",");
    }
    fprintf(fp, "\n");
}


void createVectorType() {
    int count = 3; // My custom MPI_Datatype will have 3 fields (x, y, z)
    int blocklengths[3] = {1, 1, 1}; // Each field is just one element
    MPI_Aint offsets[3]; // To assure full architecture portability we have to tell MPI how elements in our struct are stored.
    offsets[0] = offsetof(vector, x); // the first element in my MPI_Datatype object has to be distant from the base an amount of byte equal to the base of vector and x, i.e 0 byte
    offsets[1] = offsetof(vector, y); // 8 byte
    offsets[2] = offsetof(vector, z); // 16 byte
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // Each element of my custom MPI_Datatype is an MPI_DOUBLE

    MPI_Type_create_struct(count, blocklengths, offsets, types, &MPI_VECTOR);
    MPI_Type_commit(&MPI_VECTOR); // MPI finalizes the definition and use it
}

void computeAllGathervParams(int* recvcounts, int* displs, int size, int q, int r) {
  for (int i = 0; i < size; i++) {
      if (i < r) {
         recvcounts[i] = q + 1;	// if i < r, process i will receive from MPI_AllGatherv q+1 items
         displs[i] = i * (q + 1);	// // if i < r, incoming data from process i will have a displacement i*(q+1) relative to the starting adress of 'recvbuf', i.e. global_positions or global_masses in this case
     } else {
        recvcounts[i] = q;	// if i >= r, process i will receive from MPI_AllGatherv q items
        displs[i] = r * (q + 1) + (i - r) * q;
    }
  }
}


int main (int argc, char *argv[]) {
  int myrank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  createVectorType();
  
  
  
  FILE *output = NULL;
  
  int bodies;
  double t0 = 0;
  double t1 = 0;
  double t2 = 0;
  double t3 = 0;
  double t4 = 0;
  double t5 = 0;
  double t6 = 0;
  double t7 = 0;
  double t8 = 0;
  double t9 = 0;
  double cumul_time_allgatherv3 = 0;
  double cumul_time_accel = 0;
  double cumul_time_velpos = 0;
  
  if (argc > 1) {
        bodies = atoi(argv[1]);
        if (bodies <= 0) {
        	if (myrank == 0) {
            		fprintf(stderr, "Error: number of bodies must be a positive integer.\n");
            		}
            return 1;
        }
  } else {
        bodies = 16;  // Default value
        	/**if (myrank == 0) {
        		printf("No argument provided. Using default: %d bodies.\n", bodies);
        		}**/
    }
  
  
  
  
  
  if(myrank == 0) {
    output = fopen("trajectory_parallel.csv", "w");
  }
  
  int q = bodies / size;
  int r = bodies % size;
  
  vector *global_positions = (vector*)malloc(bodies*sizeof(vector));
  double *global_masses = (double*)malloc(bodies*sizeof(double));
  
  int recvcounts[size];
  int displs[size];
  computeAllGathervParams(recvcounts, displs, size, q, r);
  int b;
  
  if (myrank < r) {
    b = q + 1; // First r processes get an extra body
  } else {
    b = q;
  }


  vector *local_positions = (vector*)malloc(b*sizeof(vector));
  double *local_masses = (double*)malloc(b*sizeof(double));
  vector *velocities = (vector*)malloc(b*sizeof(vector));
  vector *accelerations = (vector*)malloc(b*sizeof(vector));
    
  unsigned int start_index = global_start_index(r, b, myrank);
  
  t0 = MPI_Wtime();
  
  initiateSystem(local_positions, local_masses, velocities, accelerations, b, start_index);
  t1 = MPI_Wtime();
  MPI_Allgatherv(local_positions, b, MPI_VECTOR, global_positions, recvcounts, displs, MPI_VECTOR, MPI_COMM_WORLD);
  t2 = MPI_Wtime();
  MPI_Allgatherv(local_masses, b, MPI_DOUBLE, global_masses, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
  t3 = MPI_Wtime();
  
  if (myrank == 0) {
      logPositions(output, bodies, global_positions);
    }
  
  for (int i = 0; i<N; i++) {
    simulate(bodies, b, start_index, dt, global_positions, global_masses, accelerations, velocities, local_positions, &t4, &t5, &t6, &cumul_time_accel, &cumul_time_velpos);
    t7 = MPI_Wtime();
    MPI_Allgatherv(local_positions, b, MPI_VECTOR, global_positions, recvcounts, displs, MPI_VECTOR, MPI_COMM_WORLD);
    t8 = MPI_Wtime();
    cumul_time_allgatherv3 += t8-t7;
    if (myrank == 0) {
      logPositions(output, bodies, global_positions);
    }
  }
  t9 = MPI_Wtime();
  
  free(global_positions);
  free(global_masses);
  free(local_positions);
  free(local_masses);
  free(velocities);
  free(accelerations);
  
  if (myrank == 0) {
      fclose(output);
      printf("Allgatherv1: %f s\n", t2 - t1);
      printf("Allgatherv2:  %f s\n", t3 - t2);
      printf("ComputeAccelerations:  %f s\n", t5 - t4);
      printf("ComputeAccelerations total time:  %f s\n", cumul_time_accel);
      printf("Vel+pos: %f s\n", t6 - t5);
      printf("Vel+pos total time: %f s\n", cumul_time_velpos);
      printf("Allgatherv3:   %f s\n", t8 - t7);
      printf("Allgatherv3 total time:   %f s\n", cumul_time_allgatherv3);
      printf("Total:   %f s\n", t9 - t0);
    }
  
  
  MPI_Finalize();
}


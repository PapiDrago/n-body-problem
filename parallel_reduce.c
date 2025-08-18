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

void initiateSystem(vector *positions, double *masses, vector *velocities, vector *accelerations, int n_body, int rank){	
	for (int i = 0; i < n_body; i++) {
		unsigned int seed = i;

		  positions[i].x = rand_uniform(&seed, -1.0, 1.0);
		  positions[i].y = rand_uniform(&seed, -1.0, 1.0);
		  positions[i].z = rand_uniform(&seed, -1.0, 1.0);
		if (rank == 0) {
			velocities[i].x = rand_uniform(&seed, -0.5, 0.5);
			velocities[i].y = rand_uniform(&seed, -0.5, 0.5);
			velocities[i].z = rand_uniform(&seed, -0.5, 0.5);
			}
		  masses[i] = 1.0;

		  accelerations[i].x = 0;
		  accelerations[i].y = 0;
		  accelerations[i].z = 0;
        }
        
}

/**void computeAccelerations(int bodies, int own_bodies, int start_global_index, vector *global_positions, double *global_masses, vector *accelerations){
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
}**/

void computeAccelerations(int bodies, int rank, int size, vector *partial_accelerations, vector *positions, double *masses) {
    const double epsilon = 1e-5;
    
    int q = bodies / size;
    int r = bodies % size;
    int iterations;
    int start_index;
    if (rank < r) {
    		iterations = q + 1; // First r processes get an extra iteration
    		start_index = rank*iterations;
    } else {
    		iterations = q;
    		start_index = r*(iterations+1) + (rank - r)*iterations;
    	}
    
    
    for (int i = 0; i < bodies; i++) {
        partial_accelerations[i].x = 0.0;
        partial_accelerations[i].y = 0.0;
        partial_accelerations[i].z = 0.0;
  	
  	
	
        for (int j = start_index; j < start_index+iterations; j++) {
            if (i != j) {
                double dx = positions[j].x - positions[i].x;
                double dy = positions[j].y - positions[i].y;
                double dz = positions[j].z - positions[i].z;

                double dist = sqrt(dx * dx + dy * dy + dz * dz) + epsilon;

                double invDistCubed = 1.0 / (dist * dist * dist);

                double scalar = GravConstant * masses[j] * invDistCubed;

                partial_accelerations[i].x += scalar * dx;
                partial_accelerations[i].y += scalar * dy;
                partial_accelerations[i].z += scalar * dz;
            }
        }
    }
}

void computeVelocities(int bodies, vector *accelerations, vector *velocities, double dt) {
	for (int i = 0; i < bodies; i++) {
		velocities[i].x = velocities[i].x + dt * accelerations[i].x;
                velocities[i].y = velocities[i].y + dt * accelerations[i].y;
                velocities[i].z = velocities[i].z + dt * accelerations[i].z;
                }
}

void computePositions(vector *positions, vector *velocities, int bodies, double dt) {
	for (int i = 0; i < bodies; i++) {
		positions[i].x = positions[i].x + dt * velocities[i].x;
		positions[i].y = positions[i].y + dt * velocities[i].y;
		positions[i].z = positions[i].z + dt * velocities[i].z;
		}
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
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // Each element of my custom MPI_Datatype is an MPI_double

    MPI_Type_create_struct(count, blocklengths, offsets, types, &MPI_VECTOR);
    MPI_Type_commit(&MPI_VECTOR); // MPI finalizes the definition and use it
}

// Custom reduction function
void vector_sum(void *in, void *inout, int *len, MPI_Datatype *dtype) {
    vector *vin = (vector*) in;
    vector *vout = (vector*) inout;
    for (int i = 0; i < *len; i++) {
        vout[i].x += vin[i].x;
        vout[i].y += vin[i].y;
        vout[i].z += vin[i].z;
    }
}


int main (int argc, char *argv[]) {
  int myrank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Request request = MPI_REQUEST_NULL;
  MPI_Status status;
  MPI_Op MPI_VECTOR_SUM;
  MPI_Op_create(&vector_sum, 1, &MPI_VECTOR_SUM);

  
  
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
  double cumul_time_reduce = 0;
  double cumul_time_accel = 0;
  double cumul_time_velpos = 0;
  double cumul_time_bcast = 0;
  
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
    output = fopen("trajectory_parallel_reduce.csv", "w");
  }
  
  
  int q = bodies / size;
  int r = bodies % size;
  
  
  vector *positions = (vector*)malloc(bodies*sizeof(vector));
  double *masses = (double*)malloc(bodies*sizeof(double));
  
  //int recvcounts[size];
  //int displs[size];
  //computeAllGathervParams(recvcounts, displs, size, q, r);

  vector *accelerations = (vector*)malloc(bodies*sizeof(vector));
  vector *partial_accelerations = (vector*)malloc(bodies*sizeof(vector));
  vector *velocities = NULL;
  if (myrank == 0) {
  	velocities = (vector*)malloc(bodies*sizeof(vector));
  }
  
  t0 = MPI_Wtime();  
  initiateSystem(positions, masses, velocities, accelerations, bodies, myrank);
  if (myrank == 0) {
      logPositions(output, bodies, positions);
    }
    
  for (int i = 0; i<N; i++) {
  	t1 = MPI_Wtime();
  	computeAccelerations(bodies, myrank, size, partial_accelerations, positions, masses);
  	t2 = MPI_Wtime();
  	cumul_time_accel += t2 - t1;
  	MPI_Reduce(partial_accelerations, accelerations, bodies, MPI_VECTOR, MPI_VECTOR_SUM, 0, MPI_COMM_WORLD);
  	if (myrank == 0) {
  		t3 = MPI_Wtime();
  		computeVelocities(bodies, accelerations, velocities, dt);
  		computePositions(positions, velocities, bodies, dt);
  		t4 = MPI_Wtime();
  		cumul_time_velpos += t4-t3;
  		cumul_time_reduce += t3-t2;
  	}
  	t5 = MPI_Wtime();
  	MPI_Ibcast(positions, bodies, MPI_VECTOR, 0, MPI_COMM_WORLD, &request);
  	if (myrank == 0) {
  		logPositions(output, bodies, positions);
  	}
  	MPI_Wait(&request, &status);
  	t6 = MPI_Wtime();
  	cumul_time_bcast += t6 - t5;
  }
  t7 = MPI_Wtime();
  free(positions);
  free(masses);
  free(partial_accelerations);
  if (velocities != NULL) {
    free(velocities);
    }

  free(accelerations);
  
  if (myrank == 0) {
      fclose(output);
      printf("Compute: %f s\n", t2 - t1);
      printf("Compute total time:  %f s\n", cumul_time_accel);
      printf("Reduce:  %f s\n", t3 - t2);
      printf("Reduce total time:  %f s\n", cumul_time_reduce);
      printf("Vel+pos: %f s\n", t4 - t3);
      printf("Vel+pos total time:  %f s\n", cumul_time_velpos);
      printf("Ibcast:   %f s\n", t6 - t5);
      printf("Ibcast total time:  %f s\n", cumul_time_bcast);
      printf("Total:   %f s\n", t7 - t0);
    }
  
  MPI_Op_free(&MPI_VECTOR_SUM);

  
  MPI_Finalize();
}

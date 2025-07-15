#include<stdlib.h>
#include<math.h>
#include<time.h>
#include <stdio.h>
#include <mpi.h>

typedef struct{
	double x,y,z;
}vector;

MPI_Datatype MPI_VECTOR;

int bodies = 10;
double GravConstant = 39.47;
double dt = 0.01;
int N = 5000; // number of iterations
double *masses;
vector *positions, *velocities, *accelerations;

vector addVectors(vector a,vector b){
	vector c = {a.x+b.x,a.y+b.y,a.z+b.z};

	return c;
}

vector scaleVector(double b,vector a){
	vector c = {b*a.x,b*a.y,b*a.z};

	return c;
}

vector subtractVectors(vector a,vector b){
	vector c = {a.x-b.x,a.y-b.y,a.z-b.z};

	return c;
}

double mod(vector a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

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

void initiateSystem(vector *local_positions, vector *velocities, vector *accelerations, double *local_masses, int n_body, unsigned int start_index){	
	for (int i = 0; i < n_body; i++) {
          unsigned int seed = start_index + i;

          local_positions[i].x = rand_uniform(&seed, -1.0, 1.0);
          local_positions[i].y = rand_uniform(&seed, -1.0, 1.0);
          local_positions[i].z = rand_uniform(&seed, -1.0, 1.0);

          velocities[i].x = rand_uniform(&seed, -0.5, 0.5);
          velocities[i].y = rand_uniform(&seed, -0.5, 0.5);
          velocities[i].z = rand_uniform(&seed, -0.5, 0.5);

          local_masses[i] = 1.0; // For simplicity

          // Initialize accelerations to zero
          accelerations[i].x = 0;
          accelerations[i].y = 0;
          accelerations[i].z = 0;
        }
}

void computeAccelerations(int bodies, int own_bodies, int start_global_index, vector *global_positions, vector *global_masses, vector *accelerations){
	int i,j;
        const double epsilon = 1e-5;
	for(i=0;i<own_bodies;i++){
		accelerations[i].x = 0;
		accelerations[i].y = 0;
		accelerations[i].z = 0;
		for(j=0;j<bodies;j++){
		        unsigned int global_index = start_global_index + i;
			if(global_index!=j){
			        double dist = mod(subtractVectors(global_positions[global_index], global_positions[j])) + epsilon;
				accelerations[i] = addVectors(accelerations[i],scaleVector(GravConstant*global_masses[j]/pow(dist,3), subtractVectors(global_positions[j],global_positions[global_index])));
			}
		}
	}
}

void computeVelocities(int own_bodies, vector *accelerations, vector *velocities, double dt) {
	for (int i = 0; i < own_bodies; i++)
		velocities[i] = addVectors(velocities[i], scaleVector(dt, accelerations[i])); //semi-implicit euler
}

void computePositions(vector *local_positions, vector *velocities, int own_bodies, double dt) {
	for (int i = 0; i < own_bodies; i++)
		local_positions[i] = addVectors(local_positions[i], scaleVector(dt, velocities[i])); //semi-implicit euler
}


void simulate(int global_bodies, int own_bodies, unsigned int start_global_index, double dt, vector *global_positions, vector *global_masses, vector *accelerations, vector *velocities, vector *local_positions){
	computeAccelerations(global_bodies, own_bodies, start_global_index, global_positions, global_masses, accelerations);
	computeVelocities(own_bodies, accelerations, velocities, dt);
	computePositions(local_positions, velocities, own_bodies, dt);
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
    offsets[0] = offsetof(vector, x); // the first element in my MPI_Datatype object has to be distant from the base an amount of byte equal to the base of vector and x, i.e. 0                                        byte
    offsets[1] = offsetof(vector, y); // 8 byte
    offsets[2] = offsetof(vector, z); // 16 byte
    MPI_Datatype types[3] = {MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE}; // Each element of my custom MPI_Datatype is an MPI_double

    MPI_Type_create_struct(count, blocklengths, offsets, types, &MPI_VECTOR);
    MPI_Type_commit(&MPI_VECTOR); // MPI finalizes the definition and use it
}

void computeAllGathervParams(int* recvcounts, int* displs, int size, int q, int r) {
  for (int i = 0; i < size; i++) {
      if (i < r) {
         recvcounts[i] = q + 1;
         displs[i] = i * (q + 1);
     } else {
        recvcounts[i] = q;
        displs[i] = r * (q + 1) + (i - r) * q;
    }
  }
}


int main (int argc, char *argv[]) {
  int myrank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  
  FILE *output = NULL;
  
  if(myrank == 0) {
    output = fopen("trajectory.csv", "w");
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
  
  initiateSystem(local_positions, local_masses, velocities, accelerations, b, start_index);
  MPI_Allgatherv(local_positions, b, MPI_VECTOR, global_positions, recvcounts, displs, MPI_VECTOR, MPI_COMM_WORLD);
  MPI_Allgatherv(local_masses, b, MPI_DOUBLE, global_masses, recvcounts, displs, MPI_DOUBLE, MPI_COMM_WORLD);
  for (int i = 0; i<N; i++) {
    simulate(bodies, b, start_index, dt, global_positions, global_masses, accelerations, velocities, local_positions);
    MPI_Allgatherv(local_positions, b, MPI_VECTOR, global_positions, recvcounts, displs, MPI_VECTOR, MPI_COMM_WORLD);
    if (myrank == 0) {
      logPositions(output, bodies, global_positions);
    }
  }
  
  
  free(global_positions);
  free(global_masses);
  free(local_positions);
  free(local_masses);
  free(velocities);
  free(accelerations);
  
  if (myrank == 0) {
      fclose(output);
    }
  
  
  MPI_Finalize();
}

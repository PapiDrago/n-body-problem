#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>

typedef struct{
	double x,y,z;
} vector;

//int bodies = 10;
double GravConstant = 39.47;
double dt = 0.01;
int N = 100; // number of iterations
double *masses;
vector *positions, *velocities, *accelerations;

double rand_uniform(unsigned int *seed, double min, double max) {
    return min + (max - min) * ((double)rand_r(seed) / RAND_MAX); //5+16
}

void initiateSystem(int bodies){
	int i;

	masses = (double*)malloc(bodies*sizeof(double));
	positions = (vector*)malloc(bodies*sizeof(vector));
	velocities = (vector*)malloc(bodies*sizeof(vector));
	accelerations = (vector*)malloc(bodies*sizeof(vector));
	
	for (int i = 0; i < bodies; i++) { //3*bodies
          unsigned int seed = i; // 1*bodies

          positions[i].x = rand_uniform(&seed, -1.0, 1.0); // (16+21+3+3)*bodies
          positions[i].y = rand_uniform(&seed, -1.0, 1.0); // (16+21+3+3)*bodies
          positions[i].z = rand_uniform(&seed, -1.0, 1.0); // (16+21+3+3)*bodies

          velocities[i].x = rand_uniform(&seed, -0.5, 0.5);// (16+21+3+3)*bodies
          velocities[i].y = rand_uniform(&seed, -0.5, 0.5);// (16+21+3+3)*bodies
          velocities[i].z = rand_uniform(&seed, -0.5, 0.5);// (16+21+3+3)*bodies

          masses[i] = 1.0; // For simplicity //(3+3)*bodies

          // Initialize accelerations to zero
          accelerations[i].x = 0; //(3+3)*bodies
          accelerations[i].y = 0; //(3+3)*bodies
          accelerations[i].z = 0; //(3+3)*bodies
        }
}

void computeAccelerations(int bodies) {
    const double epsilon = 1e-5;
    
    for (int i = 0; i < bodies; i++) {
        accelerations[i].x = 0.0;
        accelerations[i].y = 0.0;
        accelerations[i].z = 0.0;

        for (int j = 0; j < bodies; j++) {
            if (i != j) {
                double dx = positions[j].x - positions[i].x;
                double dy = positions[j].y - positions[i].y;
                double dz = positions[j].z - positions[i].z;

                double dist = sqrt(dx * dx + dy * dy + dz * dz) + epsilon;

                double invDistCubed = 1.0 / (dist * dist * dist);

                double scalar = GravConstant * masses[j] * invDistCubed;

                accelerations[i].x += scalar * dx;
                accelerations[i].y += scalar * dy;
                accelerations[i].z += scalar * dz;
            }
        }
    }
}


void computeVelocities(int bodies) {
	for (int i = 0; i < bodies; i++) {
		velocities[i].x = velocities[i].x + dt * accelerations[i].x;
                velocities[i].y = velocities[i].y + dt * accelerations[i].y;
                velocities[i].z = velocities[i].z + dt * accelerations[i].z;
                }
}

void computePositions(int bodies) {
	for (int i = 0; i < bodies; i++) {
		positions[i].x = positions[i].x + dt * velocities[i].x;
		positions[i].y = positions[i].y + dt * velocities[i].y;
		positions[i].z = positions[i].z + dt * velocities[i].z;
		}
}


void simulate(int bodies){ //195*bodies^2+164*bodies
	computeAccelerations(bodies);//(179+16)*bodies^2
	computeVelocities(bodies); //(66+16)*bodies
	computePositions(bodies); //(66+16)*bodies
}

void logPositions(FILE *fp, int bodies) {
    for (int i = 0; i < bodies; i++) {
        fprintf(fp, "%lf,%lf,%lf", positions[i].x, positions[i].y, positions[i].z);
        if (i < bodies - 1)
            fprintf(fp, ",");
    }
    fprintf(fp, "\n");
}


int main(int argc,char* argv[])
{
  int bodies;
  if (argc > 1) {
        bodies = atoi(argv[1]);  // Convert the first argument to an int
        if (bodies <= 0) {
            fprintf(stderr, "Error: number of bodies must be a positive integer.\n");
            return 1;
        }
  } else {
        bodies = 10;  // Default value
        printf("No argument provided. Using default: %d bodies.\n", bodies);
    }
	int i,j;
	FILE *output = fopen("trajectory_serial.csv", "w");
	//clock_t start = clock();  // Start timing here
    	initiateSystem(bodies); //286*bodies
    	logPositions(output, bodies);
      	printf("Body   :     x              y               z           |       vx              vy              vz   ");
		for(i=0;i<N;i++){//195*T*bodies^2+164*T*bodies
			printf("\nCycle %d\n",i+1);
			simulate(bodies); //195*bodies^2+164*bodies
			logPositions(output, bodies);
			/**for(j=0;j<bodies;j++)
				printf("Body %d : %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n",j+1,positions[j].x,positions[j].y,positions[j].z,velocities[j].x,velocities[j].y,velocities[j].z);**/
	}
	//clock_t end = clock();    // Stop timing here
	//double elapsed_time = (double)(end - start) / CLOCKS_PER_SEC;
	//printf("Total time: %.6f seconds, Bodies: %d\n", elapsed_time, bodies);
	//fclose(output);
	
	free(masses);
	free(positions);
	free(velocities);
	free(accelerations);

	return 0;
}

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<time.h>

typedef struct{
	double x,y,z;
} vector;

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

void initiateSystem(){
	int i;

	masses = (double*)malloc(bodies*sizeof(double));
	positions = (vector*)malloc(bodies*sizeof(vector));
	velocities = (vector*)malloc(bodies*sizeof(vector));
	accelerations = (vector*)malloc(bodies*sizeof(vector));
	
	for (int i = 0; i < bodies; i++) {
          unsigned int seed = i; // Deterministic seed per body

          positions[i].x = rand_uniform(&seed, -1.0, 1.0);
          positions[i].y = rand_uniform(&seed, -1.0, 1.0);
          positions[i].z = rand_uniform(&seed, -1.0, 1.0);

          velocities[i].x = rand_uniform(&seed, -0.5, 0.5);
          velocities[i].y = rand_uniform(&seed, -0.5, 0.5);
          velocities[i].z = rand_uniform(&seed, -0.5, 0.5);

          masses[i] = 1.0; // For simplicity

          // Initialize accelerations to zero
          accelerations[i].x = 0;
          accelerations[i].y = 0;
          accelerations[i].z = 0;
        }
}

void computeAccelerations(){
	int i,j;
        const double epsilon = 1e-5;
	for(i=0;i<bodies;i++){
		accelerations[i].x = 0;
		accelerations[i].y = 0;
		accelerations[i].z = 0;
		for(j=0;j<bodies;j++){
			if(i!=j){
			        double dist = mod(subtractVectors(positions[i], positions[j])) + epsilon;
				accelerations[i] = addVectors(accelerations[i],scaleVector(GravConstant*masses[j]/pow(dist,3), subtractVectors(positions[j],positions[i])));
			}
		}
	}
}

void computeVelocities() {
	for (int i = 0; i < bodies; i++)
		velocities[i] = addVectors(velocities[i], scaleVector(dt, accelerations[i])); //semi-implicit euler
}

void computePositions() {
	for (int i = 0; i < bodies; i++)
		positions[i] = addVectors(positions[i], scaleVector(dt, velocities[i])); //semi-implicit euler
}


void simulate(){
	computeAccelerations();
	computeVelocities(dt);
	computePositions(dt);
}

void logPositions(FILE *fp) {
    for (int i = 0; i < bodies; i++) {
        fprintf(fp, "%lf,%lf,%lf", positions[i].x, positions[i].y, positions[i].z);
        if (i < bodies - 1)
            fprintf(fp, ",");
    }
    fprintf(fp, "\n");
}


int main(int argC,char* argV[])
{
	int i,j;
	FILE *output = fopen("trajectory_serial.csv", "w");
    	initiateSystem(bodies);
    	logPositions(output);
      	printf("Body   :     x              y               z           |       vx              vy              vz   ");
		for(i=0;i<N;i++){
			printf("\nCycle %d\n",i+1);
			simulate(dt);
			logPositions(output);
			/**for(j=0;j<bodies;j++)
				printf("Body %d : %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n",j+1,positions[j].x,positions[j].y,positions[j].z,velocities[j].x,velocities[j].y,velocities[j].z);**/
	}
	
	fclose(output);
	
	free(masses);
	free(positions);
	free(velocities);
	free(accelerations);

	return 0;
}

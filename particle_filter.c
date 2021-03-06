/*-------------------------------------INCLUDES-------------------------------------*/
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <stdlib.h>
#include <string.h>
#include "error.h"
/*-------------------------------------DEFINES-------------------------------------*/
#define REFERENCE_POINTS 3
#define NUM_PARTICLE 100
#define PI 3.141592654
void sense(float x, float y, float direction, float* measurement_ref);

double gaussrand();

struct particle{
	float x;
	float y;
	float vx;
	float vy;
	float direction;
	float weight;
};


typedef struct particle particle;
particle** particles;
float calculate_weight(particle* p, float* measurement_ref);

float state_vector[4] = {0,0,0,0};
float reference_points[3][2] = {{0,0}, {10, 10}, {10, 0}};

int create_particles(){
	int i;
	particles = (particle**)malloc(NUM_PARTICLE*sizeof(particle*));
	//particles = malloc(NUM_PARTICLE * sizeof(particle));
	for(i=0; i<NUM_PARTICLE; i++){
		particles[i] = (particle*)malloc(sizeof(particle));
		particles[i]->x = rand()%10;
		particles[i]->y = rand()%10;
		particles[i]->direction = (float)(rand()/RAND_MAX)*2*PI;
		particles[i]->weight = 1/(float)NUM_PARTICLE;
		//printf("create_particles: particles[%d]->x,y=%f, %f\n", i, particles[i]->x, particles[i]->y);
	}
}

int main(){
	int i;
	float total_weight = 0;
	float measurement_ref[REFERENCE_POINTS];
	create_particles();
	memset(measurement_ref, 0.0, sizeof(measurement_ref));
	sense(1, 1, 1, measurement_ref);
	printf("sensed value is %f, %f, %f\n", measurement_ref[0], measurement_ref[1], measurement_ref[2]);
	for(i=0; i<NUM_PARTICLE; i++){
		float weight = calculate_weight(particles[i], measurement_ref);
		particles[i]->weight = weight;
		total_weight += weight;
		if(weight > 0.01)
			printf("main: particles[%d]->weight=%f\n", i, particles[i]->weight);	
	}
	printf("main: total_weight=%f\n", total_weight);
}
float calculate_gaussian_probability(float mean, float sigma, float x){
	float value = exp(- (pow((mean - x), 2) / (pow(sigma, 2)) / 2.0)) / sqrt(2.0 * PI * (pow(sigma, 2)));
	//value = exp(value);
	//printf("value is %f\n", value);
	return value;

}



void normalize_weight(particle** particles, float total_weight){
	int i;
	for(i=0; i<NUM_PARTICLE; i++){
		particles[i]->weight = particles[i]->weight/total_weight;
	}

}



int resample_particles(particle** particles, double max_weight){
	/*particle** particles_resampled = (particle**)malloc(NUM_PARTICLE*sizeof(particle*));
	for(i=0; i<NUM_PARTICLE; i++){
		float rand = (float)rand()/RAND_MAX;
		particles_resampled[i] = (particle*) malloc(sizeof(particle));
		
		if(particles[i]->weight > rand){
			particles_resampled[i] = particles[i];
		}
	}*/

	//Resample the particles
	int i;
	int index = (int)rand()%NUM_PARTICLE;
	double beta = 0;
	printf("index=%d\n", index);
	//int num0 = 0;
	//int num1 = 0;
	particle** resampled = (particle**)malloc(NUM_PARTICLE*sizeof(particle*));
	for(i=0; i<NUM_PARTICLE; i++)
	{

		beta += ((float)rand()/RAND_MAX)*2.0*max_weight;
		//printf("beta=%f, particles[%d]->weight=%f\n", beta, index, particles[index]->weight);
		while(beta>particles[index]->weight)
		{
			beta = beta-particles[index]->weight;
			//printf("beta=%f\n", beta);
			index = (index+1)%NUM_PARTICLE;
		}
		resampled[i] = (particle*)malloc(sizeof(particle));
		copy_particle(resampled[i], particles[index]);
		//particles[i] = particles[index];

		/*if(index==0)
			num0++;
		else if(index == 1)
			num1++;
		printf("selected index=%d, num0=%d, num1=%d\n", index, num0, num1);*/
	}
	free_particles(particles);
	particles = resampled;
	//print_particles(particles);
}

void sense(float x, float y, float direction, float* measurement_ref){
	int i;
	
	for(i=0; i<REFERENCE_POINTS; i++){
		float dist_from_ref = sqrt(pow((x-reference_points[i][0]), 2) + pow((y-reference_points[i][1]), 2));
		measurement_ref[i] = dist_from_ref; //+ gaussrand();
	}
	return;
	
}


float calculate_weight(particle* p, float* measurement_ref){
	int i;
	float weight = 1;
	//printf("calculate_weight entered\n");
	for(i=0; i<REFERENCE_POINTS; i++){
		float dist_from_ref = sqrt(pow((p->x - reference_points[i][0]), 2) + pow((p->y - reference_points[i][1]), 2));
		//printf("calculate_weight calculated dist\n");
		//printf("dist_from_ref=%f\n", dist_from_ref);
		//printf("measurement_ref=%f\n", measurement_ref[i]);
		weight = weight * calculate_gaussian_probability(dist_from_ref, 0.5, measurement_ref[i]);
		//printf("calculate_weight=%f\n", weight);
	}
	//printf("calculate_weight going out\n");
	return weight;
}

double gaussrand()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = (double)rand() / RAND_MAX;
			double U2 = (double)rand() / RAND_MAX;

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;
	printf("returning X from gaussrand=%f\n", X);
	return X;
}

int copy_particle(particle* dst, particle* src)
{
	if(dst == NULL || src == NULL)
	{
		printf("copy_particle: either src or dst pointer is null.\n");
		return NULL_POINTER_ENCOUNTERED;
	}
	dst->x = src->x;
	dst->y = src->y;
	dst->vx = src->vx;
	dst->vy = src->vy;
	dst->direction = src->direction;
	dst->weight = src->weight;
	return SUCCESS;
}


int free_particles(particle** particles)
{
	int i;
	if(particles == NULL)
	{
		printf("free_particles: particles pointer is null.\n");
		return NULL_POINTER_ENCOUNTERED;
	}
	for(i=0; i<NUM_PARTICLE; i++)
	{
		if(particles[i])
		{
			free(particles[i]);
		}
	}
	free(particles);
	return SUCCESS;
}

int print_particles(particle** particles)
{
	int i;
	if(particles == NULL)
	{
		printf("print_particle: particles pointer is null.\n");
		return NULL_POINTER_ENCOUNTERED;
	}
	printf("index          x          y          vx          vy          direction          weight\n");
	for(i=0; i<NUM_PARTICLE; i++)
	{
		if(particles[i])
		{
			printf("%d    %f    %f    %f    %f    %f    %f\n", i, particles[i]->x, particles[i]->y, particles[i]->vx, particles[i]->vy, particles[i]->direction, particles[i]->weight);
		}
	
	}

}





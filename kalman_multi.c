#include <stdio.h>
#include <math.h>
#include <omp.h>
#include <sys/time.h>
#include "matrix.h"
#define SAMPLE_TIME 1
#define NUM_MEASUREMENTS 18
#define PI 3.14159265358979323846
double* update(double old_mean, double new_mean, double old_variance, double new_variance);
double* predict(double old_mean, double motion_mean, double old_variance, double motion_variance);

double updated_values[2] = {1, 100000};
double initial_position[2] = {0, 0};
double initial_velocity[2] = {0, 0};
double x_matrix[4][1] = {{-1}, {0}, {0}, {0}};
double** updated_state_vector;
double F_matrix[4][4] = {{1, 0, SAMPLE_TIME, 0}, {0, 1, 0, SAMPLE_TIME}, {0, 0, 1, 0}, {0, 0, 0, 1}};
double H_matrix[2][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}};
double P_matrix[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1000, 0}, {0, 0, 0, 1000}};
double R_matrix[2][2] = {{0.1, 0}, {0, 0.1}};
double I_matrix[4][4] = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

double measurements[2][NUM_MEASUREMENTS] = { {10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180}, {1*PI/4, 1*PI/4, 2*PI/4, 2*PI/4, 3*PI/4, 3*PI/4, 4*PI/4, 4*PI/4, 5*PI/4, 5*PI/4, 6*PI/4, 6*PI/4, 7*PI/4, 7*PI/4, 0, 0, 1*PI/4, 1*PI/4}};
double* values = updated_values;
int main(){
        int i,j,n;
		printf("sin30=%f\n", sin(3.14/6));
       /* double measurement[10] = {10, 19, 29, 38, 51, 59, 69, 80, 91, 101};
        double measurement_variance = 1;
        double motion[10] = {10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
        double motion_variance = 1;
		double transpose[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
		calculate_transpose(*P_matrix, *transpose, 4, 4);
		printf("print the transpose: %f, %f, %f, %f\n", transpose[0][0], transpose[0][1], transpose[0][2], transpose[0][3]);
		double matrix_mul[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
		calculate_multiplication(*P_matrix, *transpose, 4, 4, 4, 4, *matrix_mul);
		double matrix_add[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
		calculate_addition(*P_matrix, *transpose, 4, 4, 4, 4, *matrix_add);
		printf(" matrix_add: %f, %f, %f, %f\n",  matrix_add[0][0],  matrix_add[0][1],  matrix_add[0][2],  matrix_add[0][3]);
		double matrix_tobeinversed[2][2] = {{1, 0}, {0, 1}};
		double matrix_inverse[2][2] = {{0, 0}, {0, 0}};
		calculate_square_inverse(*matrix_tobeinversed, *matrix_inverse);
		printf(" matrix_inverse: %f, %f, %f, %f\n",  matrix_inverse[0][0],  matrix_inverse[0][1],  matrix_inverse[1][0],  matrix_inverse[1][1]);

        for(i=0; i<sizeof(measurement)/sizeof(double); i++){
                        values = update(values[0], values[1], measurement[i], measurement_variance);
                        values = predict(values[0], values[1], motion[i], motion_variance );
        }*/
		double F_transpose[4][4];
		calculate_transpose(*F_matrix, *F_transpose, 4, 4);
		//print_matrix("F_transpose", *F_transpose, 4, 4);

		double H_transpose[4][2];
		calculate_transpose(*H_matrix, *H_transpose, 2, 4);
		//print_matrix("H_transpose", *H_transpose, 4, 2);
		//double start_time = omp_get_wtime ();

		struct timeval tv;
		gettimeofday(&tv, 0);
		double start_mill = (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ;
		//#pragma omp parallel 
		{

		//#pragma omp for num_threads(1)
		for(n=0; n<NUM_MEASUREMENTS; n++){ 
			//int tid = omp_get_thread_num();
			//int num_threads = omp_get_num_threads ();
			//printf("tid=%d, num_threads=%d\n", tid, num_threads);
			//printf("this is the iteration number %d\n", n);
		    // prediction
		    //x = (F * x) + u
			double matrix_mul[4][1] = {{0}, {0}, {0}, {0}};
			calculate_multiplication(*F_matrix, *x_matrix, 4, 4, 4, 1, *matrix_mul);
			//print_matrix("matrix_mul", *matrix_mul, 4, 1);
			copy_matrix(*matrix_mul, *x_matrix, 4, 1);
			////print_matrix("x_matrix", *x_matrix, 4, 1);

		    //P = F * P * F.transpose()
			double matrix_mult_PFt[4][4] = {{0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}, {0, 0, 0, 0}};
			calculate_multiplication(*P_matrix, *F_transpose, 4, 4, 4, 4, *matrix_mult_PFt);
			//print_matrix("matrix_mult_PFt", *matrix_mult_PFt, 4, 4);
			calculate_multiplication(*F_matrix, *matrix_mult_PFt, 4, 4, 4, 4, *P_matrix);
			//print_matrix("P_matrix", *P_matrix, 4, 4);

		    
		    //# measurement update
		    //Z = matrix([measurements[n]])

			double Z_matrix[2][1];
			/*for(i=0; i<2; i++){
					Z_matrix[i][0] = measurements[i][n] ;
			}*/
			if(measurements[1][n] >= 0  && measurements[1][n] < PI/2)
			{
				Z_matrix[0][0] = measurements[0][n]/sqrt((pow(tan(measurements[1][n]), 2) +1));
				Z_matrix[1][0] = (measurements[0][n]*tan(measurements[1][n]))/sqrt((pow(tan(measurements[1][n]), 2) +1));
			}
			else if(measurements[1][n] >= PI/2 && measurements[1][n] < PI)
			{
				Z_matrix[0][0] = -1*measurements[0][n]/sqrt((pow(tan(measurements[1][n]), 2) +1));
				Z_matrix[1][0] = -1*(measurements[0][n]*tan(measurements[1][n]))/sqrt((pow(tan(measurements[1][n]), 2) +1));
			}
			else if(measurements[1][n] >= PI && measurements[1][n] < 3*PI/2)
			{
				Z_matrix[0][0] = -1*measurements[0][n]/sqrt((pow(tan(measurements[1][n]), 2) +1));
				Z_matrix[1][0] = -1*(measurements[0][n]*tan(measurements[1][n]))/sqrt((pow(tan(measurements[1][n]), 2) +1));
			}
			else if(measurements[0][n] >= 3*PI/2 && measurements[0][n] < 2*PI)
			{
				Z_matrix[0][0] = measurements[0][n]/sqrt((pow(tan(measurements[1][n]), 2) +1));
				Z_matrix[1][0] = (measurements[0][n]*tan(measurements[1][n]))/sqrt((pow(tan(measurements[1][n]), 2) +1));
			}
			print_matrix("Z_matrix", *Z_matrix, 2, 1);
			//y = Z.transpose() - (H * x)
			/*H_matrix[0][0] = atan2(x_matrix[1][0], x_matrix[0][0])/x_matrix[0][0];
			H_matrix[1][1] = sqrt(pow(x_matrix[0][0], 2) + pow(x_matrix[1][0], 2))/x_matrix[1][0];
			printf("H00=%lf, H11=%lf\n",H_matrix[0][0], H_matrix[1][1] );*/
			double matrix_mult_Hx[2][1];
			calculate_multiplication(*H_matrix, *x_matrix, 2, 4, 4, 1, *matrix_mult_Hx);
			//print_matrix("matrix_mult_Hx", *matrix_mult_Hx, 2, 1);
			double y_matrix[2][1];
			calculate_subtraction(*Z_matrix, *matrix_mult_Hx, 2, 1, 2, 1, *y_matrix);
			//print_matrix("y_matrix", *y_matrix, 2, 1);

		    //S = H * P * H.transpose() + R
			double matrix_mult_PHt[4][2];
			calculate_multiplication(*P_matrix, *H_transpose, 4, 4, 4, 2, *matrix_mult_PHt);
			//print_matrix("matrix_mult_PHt", *matrix_mult_PHt, 4, 2);
			double matrix_mult_HPHt[2][2];
			calculate_multiplication(*H_matrix, *matrix_mult_PHt, 2, 4, 4, 2, *matrix_mult_HPHt);
			//print_matrix("matrix_mult_HPHt", *matrix_mult_HPHt, 2, 2);
			double S_matrix[2][2];
			calculate_addition(*matrix_mult_HPHt, *R_matrix, 2, 2, 2, 2, *S_matrix);
			//print_matrix("S_matrix", *S_matrix, 2, 2);

		    //K = P * H.transpose() * S.inverse()
			double S_inverse[2][2];
			calculate_square_inverse(*S_matrix, *S_inverse);
			//print_matrix("S_inverse", *S_inverse, 2, 2);
			double matrix_mult_HtSi[4][2];
			calculate_multiplication(*H_transpose, *S_inverse, 4, 2, 2, 2, *matrix_mult_HtSi);
			//print_matrix("matrix_mult_HtSi", *matrix_mult_HtSi, 4, 2);
			double K_matrix[4][2];
			calculate_multiplication(*P_matrix, *matrix_mult_HtSi, 4, 4, 4, 2, *K_matrix);
			//print_matrix("K_matrix", *K_matrix, 4, 2);

		    //x = x + (K * y)
			double matrix_mult_Ky[4][1];
			calculate_multiplication(*K_matrix, *y_matrix, 4, 2, 2, 1, *matrix_mult_Ky);
			//print_matrix("matrix_mult_Ky", *matrix_mult_Ky, 4, 1);
			double matrix_add_x_Ky[4][1];
			calculate_addition(*x_matrix, *matrix_mult_Ky, 4, 1, 4, 1, *matrix_add_x_Ky);
			//print_matrix("matrix_add_x_Ky", *matrix_add_x_Ky, 4, 1);
			copy_matrix(*matrix_add_x_Ky, *x_matrix, 4, 1);
			print_matrix("x_matrix", *x_matrix, 4, 1);
		
		    //P = (I - (K * H)) * P
			double matrix_mult_KH[4][4];
			calculate_multiplication(*K_matrix, *H_matrix, 4, 2, 2, 4, *matrix_mult_KH);
			//print_matrix("matrix_mult_KH", *matrix_mult_KH, 4, 4);	
			double matrix_sub_I_KH[4][4];
			calculate_subtraction(*I_matrix, *matrix_mult_KH, 4, 4, 4, 4, *matrix_sub_I_KH);
			//print_matrix("matrix_sub_I_KH", *matrix_sub_I_KH, 4, 4);
			double matrix_mult_IKHP[4][4];
			calculate_multiplication(*matrix_sub_I_KH, *P_matrix, 4, 4, 4, 4, *matrix_mult_IKHP);
			//print_matrix("matrix_mult_IKHP", *matrix_mult_IKHP, 4, 4);
			copy_matrix(*matrix_mult_IKHP, *P_matrix, 4, 4);
			////print_matrix("P_matrix", *P_matrix, 4, 2);
		}

		}
		gettimeofday(&tv, 0);
		double finish_mill = (tv.tv_sec) * 1000 + (tv.tv_usec) / 1000 ;
		printf("start_mill = %lf\n", start_mill);
		printf("finish_mill = %lf\n", finish_mill);
		printf("total time for the  kalman loop execution = %lf\n", finish_mill-start_mill);

		/*double finish_time = omp_get_wtime();
		double total_time = finish_time - start_time;
		printf("total time for the  kalman loop execution = %lf\n", total_time);*/
}


double* update(double old_mean, double old_variance, double new_mean, double new_variance){
        //double values[2];
        updated_values[0] = (old_mean*new_variance + new_mean*old_variance)/(old_variance + new_variance);
        updated_values[1] = 1/(1/old_variance + 1/new_variance);
        printf("updated position=%f and updated variance=%f\n", updated_values[0], updated_values[1]);
        return updated_values;
}

double* predict(double old_mean, double old_variance, double motion_mean, double motion_variance){
        //double predicted_values[2];
        updated_values[0] = old_mean + motion_mean;
        updated_values[1] = old_variance + motion_variance;
        printf("predicted position=%f and predicted variance=%f\n", updated_values[0], updated_values[1]);
        return updated_values;
}


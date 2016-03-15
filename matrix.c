#include <stdio.h>
#include <math.h>
#include "matrix.h"


#define SUCCESS 0
#define FAILURE 1

void calculate_transpose(double* matrix, double* transpose, int dimx, int dimy){
	int i,j;
	//int dimx = 4;//sizeof(matrix)/sizeof(matrix[0]);
	//int dimy = 4;//sizeof(matrix)/sizeof(matrix[0][0])/dimx;
	//printf("source matrix: %f, %f, %f, %f\n", *(matrix+3*dimy+0), *(matrix+3*dimy+1), *(matrix+3*dimy+2), *(matrix+3*dimy+3));
	for(i=0; i<dimy; i++)
		for(j=0; j<dimx; j++)
			*(transpose+i*dimx+j) = *(matrix+j*dimy+i);
	//printf("print the matrix: %f, %f, %f, %f\n", matrix[0][0], matrix[1][0], matrix[2][0], matrix[3][0]);
	//printf("print the transpose: %f, %f, %f, %f\n", transpose[0][0], transpose[0][1], transpose[0][2], transpose[0][3]);
	//printf("print the transpose: %f, %f, %f, %f\n", transpose[1][0], transpose[1][1], transpose[1][2], transpose[1][3]);
	//printf("print the transpose: %f, %f, %f, %f\n", transpose[2][0], transpose[2][1], transpose[2][2], transpose[2][3]);
	//printf("print the transpose: %f, %f, %f, %f\n", transpose[3][0], transpose[3][1], transpose[3][2], transpose[3][3]);
}

int calculate_multiplication(double* matrix_a, double* matrix_b, int dimx_a, int dimy_a, int dimx_b, int dimy_b, double *matrix_mul){
	int i,j,k;
	if(dimy_a != dimx_b){
		printf("Dimensions of the matrices not suitable for multiplication.");
		return FAILURE;
	}
	double value = 0;
	for(i=0; i<dimx_a; i++){
		for(j=0; j<dimy_b; j++){
			for(k=0; k<dimy_a; k++){
				value += (*(matrix_a+i*dimy_a+k)) *  (*(matrix_b+k*dimy_b+j));
			}
			*(matrix_mul+i*dimy_b+j) = value;
			value = 0;
		}
	}
	//printf("multiplied matrix: %f, %f, %f, %f\n", *(matrix_mul+0*dimy_b+0), *(matrix_mul+0*dimy_b+1), *(matrix_mul+0*dimy_b+2), *(matrix_mul+0*dimy_b+3));
	return SUCCESS;
}

int calculate_addition(double* matrix_a, double* matrix_b, int dimx_a, int dimy_a, int dimx_b, int dimy_b, double *matrix_add){
	int i,j;
	if(dimx_a != dimx_b && dimy_a != dimy_b){
		printf("Dimensions of the matrices do not match.");
		return FAILURE;
	}
	for(i=0; i<dimx_a; i++){
		for(j=0; j<dimy_b; j++){
			*(matrix_add+i*dimy_a+j) = (*(matrix_a+i*dimy_a+j)) +  (*(matrix_b+i*dimy_b+j));
		}
	}
	//printf("added matrix: %f, %f, %f, %f\n", *(matrix_add+2*dimy_a+0), *(matrix_add+2*dimy_a+1), *(matrix_add+2*dimy_a+2), *(matrix_add+2*dimy_a+3));
	return SUCCESS;
}

int calculate_subtraction(double* matrix_a, double* matrix_b, int dimx_a, int dimy_a, int dimx_b, int dimy_b, double *matrix_sub){
	int i,j;
	if(dimx_a != dimx_b && dimy_a != dimy_b){
		printf("Dimensions of the matrices do not match.");
		return FAILURE;
	}
	for(i=0; i<dimx_a; i++){
		for(j=0; j<dimy_b; j++){
			*(matrix_sub+i*dimy_a+j) = (*(matrix_a+i*dimy_a+j)) -  (*(matrix_b+i*dimy_b+j));
		}
	}
	//printf("subtracted matrix: %f, %f, %f, %f\n", *(matrix_sub+2*dimy_a+0), *(matrix_sub+2*dimy_a+1), *(matrix_sub+2*dimy_a+2), *(matrix_sub+2*dimy_a+3));
	return SUCCESS;
}

int calculate_square_inverse(double* matrix_a, double *matrix_inverse){
	double determinant = (*(matrix_a+0))*(*(matrix_a+3)) - (*(matrix_a+1))*(*(matrix_a+2));
	printf("determinant is %f\n", determinant);
	*(matrix_inverse+0) = *(matrix_a+3)/determinant;
	*(matrix_inverse+1) = (-1)*(*(matrix_a+1))/determinant;
	*(matrix_inverse+2) = (-1)*(*(matrix_a+2))/determinant;
	*(matrix_inverse+3) = *(matrix_a+0)/determinant;
	//printf("inversed matrix: %f, %f, %f, %f\n", *(matrix_inverse+0), *(matrix_inverse+1), *(matrix_inverse+2), *(matrix_inverse+3));
	return SUCCESS;
}

int copy_matrix(double* matrix_src, double *matrix_dst, int dimx, int dimy){
	int i,j;
	for(i=0; i<dimx; i++){
		for(j=0; j<dimy; j++){
			*(matrix_dst+i*dimy+j) = *(matrix_src+i*dimy+j);
		}
	}
	return SUCCESS;
}

int print_matrix(char* name, double* matrix, int dimx, int dimy){
	int i,j;
	printf("printing matrix %s\n", name);
	for(i=0; i<dimx; i++){
		for(j=0; j<dimy; j++){
			if(j==0 && j==(dimy-1)){
				printf("[%f]\n ", *(matrix+i*dimy+j));
			}
			else if(j==0){
				printf("[%f ", *(matrix+i*dimy+j));
			}
			else if(j==(dimy-1)){
				printf("%f]\n", *(matrix+i*dimy+j));
			}
			else{
				printf("%f ", *(matrix+i*dimy+j));
			}
		}
	}
}

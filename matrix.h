#include <stdio.h>

typedef struct matrix{
	int dimx;
	int dimy;
	double** matrix_value;
} matrix;

void calculate_transpose(double* matrix, double* transpose, int dimx, int dimy);
int calculate_multiplication(double* matrix_a, double* matrix_b, int dimx_a, int dimy_a, int dimx_b, int dimy_b, double *matrix_mul);
int calculate_addition(double* matrix_a, double* matrix_b, int dimx_a, int dimy_a, int dimx_b, int dimy_b, double *matrix_add);
int calculate_subtraction(double* matrix_a, double* matrix_b, int dimx_a, int dimy_a, int dimx_b, int dimy_b, double *matrix_sub);
int calculate_square_inverse(double* matrix_a, double *matrix_inverse);
int copy_matrix(double* matrix_src, double *matrix_dst, int dimx, int dimy);
int print_matrix(char* name, double* matrix, int dimx, int dimy);

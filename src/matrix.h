/* matrix.h : header file for matrix routines.*/

#ifndef MATRIX_HEADER_IN
#define MATRIX_HEADER_IN
/* The basic matrix structure */
#define TOL 1e-10
#define VEC M[0]

typedef struct
{ int vec, r,c,original_r,original_c;long mem;double **M,*V;} matrix;

extern matrix null_mat;

extern long matrallocd;
#endif
/* The user routines */

void mcopy(matrix *A,matrix *B);
matrix initmat(int rows,int cols);
void freemat(matrix A);
void vmult(matrix *A,matrix *b,matrix *c,int t);
void matmult(matrix C,matrix A,matrix B,int tA,int tB);
void invert(matrix *a);
double dot(matrix a,matrix b);
double enorm(matrix d);
void householder(matrix *u,matrix a,matrix b,int t1);
void Hmult(matrix C,matrix u);
void HQmult(matrix C,matrix U,int p,int t);
void QT(matrix Q,matrix A,int Qfull);
void Rsolv(matrix *R,matrix *p,matrix *y, int transpose);
int QR(matrix *Q,matrix *R);
void OrthoMult(matrix *Q,matrix *A,int off,int rows,int t,int pre,int o_pre);
void matrixintegritycheck(void);
void msort(matrix a);
void RArrayFromMatrix(double *a,int r,matrix *M);
matrix Rmatrix(double *A,int r,int c);
matrix initvec(int rows);



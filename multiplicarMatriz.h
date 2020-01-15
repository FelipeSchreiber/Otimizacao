/*
Universidade Federal do Rio de Janeiro
Escola Politecnica
Autor:Felipe Schreiber Fernandes 

$Author$
$Date$
$Log$
*/

#include <math.h>
#ifndef _MultiplicarMatriz_
#define _MultiplicarMatriz_ "@(#)multiplicarMatriz.h $Revision$"
/*#define NUMERO_MAXIMO_LINHAS1 2 //caso for fazer minimos quadrados esse valor deve ser 2, senao a quantidade de linhas da matriz de entrada
#define NUMERO_MAXIMO_COLUNAS1 4//quantidade de colunas da matriz de entrada
#define INCOMPATIVEL 5*/

//void leastSquaresTest();

float* leastSquares(float *xValues,float *yValues,int n);

void eigenvectorsJacobiMethodTest(float **matrix, int n);

float**eigenvectorsJacobiMethod(float **matrix, int n,float tolerance);

//void powerMethodTest(float **matrix,int n);

float* powerMethod(float **matriz, float *chute,int n,float residuoAceitavel);

//void gaussSeidelMethodTest(float **matrix,float *termosIndependentes,int n);

void gaussSeidelMethod(float **matriz, float *vector, float *termosIndependentes,int n,int normaDoVetor,float residuoAceitavel);

void gaussJordanEliminationTests(float **matrix,float *termosIndependentes,int n);

float** eliminacaoGaussJordan(float **matriz1, float *termosIndependentes, int n);

void gaussEliminationTests(float **matrix,float *termosIndependentes,int n);

float** eliminacaoGauss(float **matriz1,float **history ,float *termosIndependentes, int n);

//void jacobiMethodTest(float **matrix,float *termosIndependentes,int n);

void jacobiMethod(float **matriz, float *vector, float *termosIndependentes,int n,int normaDoVetor,float residuoAceitavel);

void choleskyTest(float**matrix,float* termosIndependentes,int n);

float** choleskyDecomposition(float **matriz, int n);

float** getRegressorsMatrix(float *xValues,int n);

void deallocateArray(float ** ptr, int row);

float** criarMatrizIdentidade(int n);

float calcularResiduo(float *X1, float* X0, int n, int norma);

float** eliminacaoJordan(float **matriz1,float **history, float *termosIndependentes, int n);

int swap(float **matriz, float *termosIndependentes,int column, int n);

int swap2(float **matriz, float *termosIndependentes,int column, int n);

void imprimeMatriz(float** MatrizVertice,int n,int m);

void MultiplicarMatrizVetor(float **matriz1, float *matriz2, float *matrizResultado, int n, int m);

void MultiplicarMatrizes(float **matriz1, float **matriz2, float **matrizResultado,int m, int n, int o);

//float** matrizInit(int n,float matriz[NUMERO_MAXIMO_LINHAS1][NUMERO_MAXIMO_COLUNAS1]);

//float* vectorInit(int n,float matriz[NUMERO_MAXIMO_LINHAS1]);

void copyMatrix(float **matrix1,float **matrix2,int n);

float** calculateRotationMatrix(float **matrix,int n);

int* findGreatestValueIndex(float **matrix,int n);

float calculatePhi(float **matrix, int i, int j);

float** getTranspose(float **matrix,int m,int n);

float det(float **A, int n);

int isDiagonalDominant(float **matrix, int n);

int checkSimmetry(float**matrix, int n);

int checkEigenvalues(float**matrix, int n);

int isPositiveDefinite(float**matrix, int n);

void multiplicarMatrizPorEscalar(float**matrix,float alpha,int n);

float* retroSubstituicao(float**matrix,float*termosIndependentes,int n);

float* forwardSubstituicao(float**matrix,float*termosIndependentes,int n);

void multiplicarVetorPorEscalar(float*vector,float beta,int n);

#endif 
/*$RCSfile$*/

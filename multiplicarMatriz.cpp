/*
Universidade Federal do Rio de Janeiro
Escola Politecnica
Autor: Felipe Schreiber Fernandes

$Author$
$Date$
$Log$

*/

#include <stdio.h>
#include <stdlib.h>
#include "multiplicarMatriz.h"

/*
void leastSquaresTest()
{
 float beta = 0.58;
 float yValues[NUMERO_MAXIMO_COLUNAS1] = {0,0.693,2.197,2.996};
 float *YValues = vectorInit(NUMERO_MAXIMO_COLUNAS1, yValues);
 multiplicarVetorPorEscalar(YValues,beta,NUMERO_MAXIMO_COLUNAS1);
 float xValues[NUMERO_MAXIMO_COLUNAS1] = {0,0.693,1.099,1.386};
 float *XValues = vectorInit(NUMERO_MAXIMO_COLUNAS1, xValues);
 float*B = leastSquares(XValues,YValues,NUMERO_MAXIMO_COLUNAS1);
 free(B);
}

void powerMethodTest(float **matrix,int n)
{
  float chuteInicial[NUMERO_MAXIMO_LINHAS1] = {3.0,3.0};
  float* chute= vectorInit(n, chuteInicial);    
  float residuoMin = 0.001;
  powerMethod(matrix,chute,n,residuoMin);
  printf("\npower method solution eigenvector:\n\n");
    for (int i = 0; i<n; i++)
     printf("%.3f\n", chute[i]);
    free(chute);      
}

void gaussSeidelMethodTest(float **matrix,float *termosIndependentes,int n)
{    
    float chuteInicial[NUMERO_MAXIMO_LINHAS1] = {1.0,1.0};
    float* chute= vectorInit(n, chuteInicial);    
    float residuoMin = 0.001;
    int norma = 2;
    gaussSeidelMethod(matrix,chute,termosIndependentes,n,norma,residuoMin);
    printf("\nGauss-Seidel solution vector:\n\n");
    for (int i = 0; i<n; i++)
     printf("%.3f\n", chute[i]);
    free(chute);    
}

void jacobiMethodTest(float **matrix,float *termosIndependentes,int n)
{    
    float chuteInicial[NUMERO_MAXIMO_LINHAS1] = {1.0,1.0};
    float* chute= vectorInit(n, chuteInicial);    
    float residuoMin = 0.001;
    int norma = 2;
    jacobiMethod(matrix,chute,termosIndependentes,n,norma,residuoMin);
    printf("\nJacobi solution vector:\n\n");
    for (int i = 0; i<n; i++)
     printf("%.3f\n", chute[i]);
    free(chute);    
}


//Funcao que converte a matriz do tipo float(*)[NUMERO_MAXIMO_LINHAS1] pro tipo float** 
float** matrizInit(int n,float matriz[NUMERO_MAXIMO_LINHAS1][NUMERO_MAXIMO_COLUNAS1])
{
  float**Matrix = (float**) malloc(n*sizeof(float*));
  for(int i = 0; i < n; i++)
  {
   Matrix[i] = (float*) malloc(n*sizeof(float));
   for(int j = 0; j < n; j++)
   {
    Matrix[i][j] = matriz[i][j];
   }
  }
  return Matrix;
}

float* vectorInit(int n,float matriz[NUMERO_MAXIMO_LINHAS1])
{
  float* vector = (float*) malloc(n*sizeof(float));
   for(int j = 0; j < n; j++)
   {
    vector[j] = matriz[j];
   }
  
  return vector;
}*/

float* leastSquares(float *xValues,float *yValues,int n)
{
 float determinante =0;
 float**P = getRegressorsMatrix(xValues,n);
 float**result = criarMatrizIdentidade(2);
 float**A = criarMatrizIdentidade(2);
 float**inverseA = criarMatrizIdentidade(2);
 float*C = (float*)malloc(2*sizeof(float));
 float*B = (float*)malloc(2*sizeof(float));
 for(int i = 0; i<n;i++)
 {
  for(int j = 0; j<2; j++)
   printf("%.3f ",P[i][j]);
  printf("\n");
 } 
 printf("\n\n");
 float**PTranspose = getTranspose(P,n,2);
 for(int i = 0; i<2;i++)
 {
  //for(int j = 0; j<n; j++)
   //printf("%.3f ",PTranspose[i][j]);
  printf("\n");
 }
 printf("\n\n\n"); 
 MultiplicarMatrizes(PTranspose,P,result,2, n, 2); 
 copyMatrix(result,A,2);
 determinante = det(A,2);
 //NUMERO_MAXIMO_LINHAS1
 //printf("\nDet: %.3f\n",determinante);
 inverseA[0][0] = A[1][1]/determinante;
 inverseA[1][1] = A[0][0]/determinante;
 inverseA[0][1] = (-1)*A[0][1]/determinante;
 inverseA[1][0] = (-1)*A[1][0]/determinante;
 //imprimeMatriz(A, 2, 2);
 //imprimeMatriz(inverseA, 2, 2);
 MultiplicarMatrizes(A,inverseA,result,2, 2, 2);
 //imprimeMatriz(result, 2, 2);
 MultiplicarMatrizVetor(PTranspose,yValues,C,2,n); 
 printf("C[0] = %.3f \n",C[0]);
 printf("C[1] = %.3f \n\n",C[1]);
 B[0] = inverseA[0][0]*C[0] + inverseA[0][1]*C[1];
 B[1] = inverseA[1][0]*C[0] + inverseA[1][1]*C[1];
 printf("B[0] = %.3f \n",B[0]);
 printf("B[1] = %.3f \n",B[1]);
 deallocateArray(P,n);
 deallocateArray(PTranspose,2);
 deallocateArray(result,2);
 deallocateArray(inverseA,2);
 deallocateArray(A,2);
 free(C);
 return B;
}

float** getRegressorsMatrix(float *xValues,int n)
{
 int i;
 float**P = (float**) malloc(n*sizeof(float*));
 for(i = 0; i < n; i++)
  P[i] = (float*) malloc(2*sizeof(float));
 for(i = 0; i<n; i++)
 {
  P[i][0] = 1;
  P[i][1] = xValues[i];
 }
 return P;
}

void eigenvectorsJacobiMethodTest(float **matrix, int n)
{
 float tolerance = 0.001; 
 float**eigenvectors = eigenvectorsJacobiMethod(matrix,n,tolerance);  
 if(eigenvectors)    
 { 
	 printf("\nEigenvectors:\n");
	 imprimeMatriz(eigenvectors,n,n);
	 printf("\n\n");
	 printf("\nEigenvalues:\n");
	 imprimeMatriz(matrix,n,n);
	 printf("\n\n");
	 deallocateArray(eigenvectors,n);
 }
}

float**eigenvectorsJacobiMethod(float **matrix, int n,float tolerance)
{
 if(checkSimmetry(matrix, n))
 {       int count = 0; 
	 int *i_j;
	 int i,j;
	 float**eigenvectors = criarMatrizIdentidade(n);
	 float**P;
	 float**PTranspose;
	 float**result = criarMatrizIdentidade(n);
	 do
	 {
          count++;
	  P = calculateRotationMatrix(matrix,n);
	  PTranspose = getTranspose(P,n,n);
	  MultiplicarMatrizes(matrix,P,result,n,n,n);
	  copyMatrix(result,matrix,n);
	  MultiplicarMatrizes(PTranspose,matrix,result,n,n,n); 
	  copyMatrix(result,matrix,n);
	  /*printf("\nPTranspose x A x P:\n");
	  imprimeMatriz(matrix,n,n);
	  printf("\n\n");*/
	  MultiplicarMatrizes(eigenvectors,P,result,n,n,n);
	  copyMatrix(result,eigenvectors,n);
	  i_j = findGreatestValueIndex(matrix,n);
	  i = i_j[0];
	  j = i_j[1];
	 }while(fabs(matrix[i][j]) > tolerance);
	 free(i_j);
	 deallocateArray(P, n);
	 deallocateArray(PTranspose, n);
	 deallocateArray(result, n);
         //printf("\nTotal iteracoes:%i\n",count); 
       	return eigenvectors;
 }
 else
 {
  //printf("\n\nMATRIZ NAO SIMETRICA, METODO DE JACOBI INDETERMINADO \n\n");
  return NULL;
 }
}

void copyMatrix(float **matrix1,float **matrix2, int n)
{
 int i,j;
 for(i = 0; i<n; i++)
  for(j = 0; j<n; j++)
   matrix2[i][j] = matrix1[i][j];
}

//original matrix has m lines and n columns
float** getTranspose(float **matrix,int m,int n)
{
 int i,j;
 float**Transpose = (float**) malloc(n*sizeof(float*));
 for(i = 0; i<n ; i++)
  Transpose[i] = (float*) malloc(m*sizeof(float));
 for(i = 0; i<n; i++)
  for(j = 0; j<m; j++)
   Transpose[i][j] = matrix[j][i];
 return Transpose;
}
void multiplicarMatrizPorEscalar(float**matrix,float alpha,int n)
{
 int i = 0;
 int j = 0;
 for(i = 0; i<n; i++)
 {
  for(j = 0; j<n;j++)
  {
   matrix[i][j] *= alpha;
  }
 }
}

float** calculateRotationMatrix(float **matrix,int n)
{
 float** P = criarMatrizIdentidade(n);
 int *i_j = findGreatestValueIndex(matrix,n);
 int i = i_j[0];
 int j = i_j[1];
 //printf("\nTermo a ser zerado: A(%i,%i)\n",i,j);
 free(i_j);
 float phi = calculatePhi(matrix,i,j);
 //printf("phi: %.3f",phi);
 P[j][j] = P[i][i] = (float)cos((float)phi);
 P[i][j] = (-1)*(float)sin((float)phi);
 P[j][i] = (-1)*P[i][j];
 return P;
}

int* findGreatestValueIndex(float **matrix,int n)
{
 int i,j;
 float greatestValue = 0;
 int *greatestValueIndex = (int*)malloc(2*sizeof(int));
 for(i = 0;i<n;i++)
  for(j = 0;j<n;j++)
  {
   if((i != j) && (fabs(matrix[i][j]) > greatestValue))
   {
    greatestValue = fabs(matrix[i][j]);
    greatestValueIndex[0] = i;
    greatestValueIndex[1] = j;
   }
  }
 return greatestValueIndex;
}

float calculatePhi(float **matrix, int i, int j)
{
 float phi = 0;
 if(matrix[i][i] != matrix[j][j])
 {
  //phi = atan2((double)2*matrix[i][j], (matrix[i][i] - matrix[j][j])  )/2;
  phi = atan((double)2*matrix[i][j]/(matrix[i][i] - matrix[j][j])  )/2;
 }
 else
 {
  phi = (float)M_PI/4;
 }
 return phi;
}

float* retroSubstituicao(float**matrix,float*termosIndependentes,int n)
{
 int i = 0;
 int j = 0;
 float*solution = (float*)malloc(n*sizeof(float));
 for(i = n-1; i>=0; i--)
 {
  solution[i] = termosIndependentes[i];
  for(j = n-1; j>i; j--)
  {
   solution[i] -= matrix[i][j]*solution[j];
   //printf("\nmatrix[%i][%i] = %.3f\n",i,j,matrix[i][j]);
  }
  solution[i] /= matrix[i][i];
 }
 return solution;
}

float* forwardSubstituicao(float**matrix,float*termosIndependentes,int n)
{
 int i = 0;
 int j = 0;
 float*solution = (float*)malloc(n*sizeof(float));
 for(i = 0; i<n; i++)
 {
  solution[i] = termosIndependentes[i];
  for(j = 0; j<i; j++)
  {
   //printf("\nmatrix[%i][%i] = %.3f\n",i,j,matrix[i][j]);
   solution[i] -= matrix[i][j]*solution[j];
  }
  solution[i] /= matrix[i][i];
 }
 //printf("\n\nFim forward\n\n");
 return solution;
}

float * powerMethod(float **matriz, float *chute,int n,float residuoAceitavel)
{
 int i,j;
 int count = 0;
 float *eigenvector = (float *)malloc(n*sizeof(float));
 float *err = (float *)malloc(n*sizeof(float));
 float *result = (float *)malloc(n*sizeof(float));
 float errMax = 0;
 float lambda = 0;
 do
 {
  count++;
  for(i=0;i<n;i++)
  {
   eigenvector[i]=0;
   for(j=0;j<n;j++)
    eigenvector[i]+=matriz[i][j]*chute[j];
  }
  lambda = fabs(eigenvector[0]);
  for(i=0;i<n;i++)
  {
   if(fabs(eigenvector[i])>fabs(lambda))
    lambda=fabs(eigenvector[i]);
  }
  errMax = 0;
  for(i=0;i<n;i++)
  {
   eigenvector[i]/=lambda;
   err[i] = 0;
   err[i] = fabs(fabs(eigenvector[i]) - fabs(chute[i]))/fabs(eigenvector[i]); 
   if(err[i]>errMax)
    errMax = err[i];  
   chute[i] = eigenvector[i];
  }
 }while(errMax>residuoAceitavel);
 printf("\nLambda: %.3f\n",lambda);
 MultiplicarMatrizVetor(matriz,eigenvector,result,n,n);
 printf("\nA x Eigenvector:\n");
 for(i = 0;i< n;i++)
  printf("%.3f ",result[i]);
 printf("\nEigenvalue x Eigenvector:\n");
 for(i = 0;i< n;i++)
  printf("%.3f ",lambda*eigenvector[i]);
 printf("\n\n"); 
 printf("\nTotal iteracoes power method: %i\n",count);
 return eigenvector;
}

void gaussJordanEliminationTests(float **matrix,float *termosIndependentes,int n)
{
 printf("\nOriginal:\n");
 imprimeMatriz((float **)matrix, n, n);
 float** inversa = eliminacaoGaussJordan(matrix, termosIndependentes, n); 
 printf("\nMatrix D:\n");
 imprimeMatriz((float **)matrix, n, n);
 printf("\nInversa de A:\n");
 imprimeMatriz((float **)inversa, n, n);	
 printf("\nTermos Independentes:\n");
 for (int i = 0;i <n; i++)
  printf("%.4f\n",termosIndependentes[i]);
 printf("\nSOLUCAO:\n");
 for (int i = 0;i <n; i++)
 {
  termosIndependentes[i] /= matrix[i][i];
  printf("%.4f\n",termosIndependentes[i]);
 } 
}

void gaussEliminationTests(float **matrix,float *termosIndependentes,int n)
{
 float **result = criarMatrizIdentidade(n);
 float **HISTORY = criarMatrizIdentidade(n);
 float *Y = (float*)malloc(n*sizeof(float));
 float *X = (float*)malloc(n*sizeof(float));
 printf("\nOriginal:\n");
 imprimeMatriz((float **)matrix, n, n);
 float **L = eliminacaoGauss(matrix, HISTORY,termosIndependentes, n); 
 Y = forwardSubstituicao(L,termosIndependentes,n);
 X = retroSubstituicao(matrix,Y,n);
 printf("\nX:\n");
 for (int i = 0;i <n; i++)
  printf("%.4f\n",X[i]);
 printf("\nMatrix U:\n");
 imprimeMatriz((float **)matrix, n, n);	
 printf("\nMatriz P x L:\n");
 imprimeMatriz(L, n, n);
 printf("\nP x L x U\n");
 MultiplicarMatrizes(L, matrix, result,n,n,n);//retorna o produto matricial LU na matriz result
 imprimeMatriz(result, n, n);
 deallocateArray(result, n);
 deallocateArray(L, n);
}

void choleskyTest(float**matrix,float*termosIndependentes,int n)
{
 int i = 0;
 float *Y = (float*)malloc(n*sizeof(float));
 float *X = (float*)malloc(n*sizeof(float));
 /*
 for(i = 0;i<n;i++)
  Y[i] = termosIndependentes[i];*/
 float **L = choleskyDecomposition(matrix,n);
 if(L)
 {
  float **LTransposta = criarMatrizIdentidade(n);
  float **result = criarMatrizIdentidade(n);
  for(i = 0; i<n; i++)
  {
   for(int j=0; j<n; j++)
   {
    LTransposta[i][j] = L[j][i];
   }
  }
  printf("\nMatriz Original:\n\n");
  imprimeMatriz((float **)matrix, n, n);	
  printf("\n\n");
  printf("\nMatriz L:\n\n");       
  imprimeMatriz((float **)L, n, n);	
  printf("\n\n");
  printf("\nMatriz L Transposta:\n\n");  
  imprimeMatriz((float **)LTransposta, n, n);	
  printf("\n\n"); 
  printf("\nMatriz L x L*\n\n");      
  //MultiplicarMatrizes(L,LTransposta, result);
  MultiplicarMatrizes(L,LTransposta, result,n,n,n);
  imprimeMatriz(result, n, n);
  Y = forwardSubstituicao(L,termosIndependentes,n);
  
  //gaussJordanEliminationTests(L,Y,n);
  //gaussJordanEliminationTests(LTransposta,Y,n);
 
  X = retroSubstituicao(LTransposta,Y,n);
  printf("\nX:\n"); 
  for(i = 0; i<n;i++)
   printf("\n%.3f\n",X[i]);
  printf("\nSolucao\n"); 
  for(i = 0; i<n;i++)
   printf("\n%.3f\n",X[i]);
  deallocateArray(result, n);
  deallocateArray(LTransposta, n);
  deallocateArray(L, n);
 }
}

void multiplicarVetorPorEscalar(float*vector,float beta,int n)
{
 for(int i = 0;i<n;i++)
  vector[i] *= beta;
}

float** choleskyDecomposition(float **matriz, int n)
{
 int i = 0;
 int j = 0;
 int k = 0;
 if( isPositiveDefinite(matriz, n) )
 {
  float **L = criarMatrizIdentidade(n);
  for(i = 0; i<n; i++)
  {
   L[i][i] = matriz[i][i];
   for(k = 0; k < i; k++)
   {
    L[i][i] -= L[i][k]*L[i][k];
    //printf("\nL[%i][%i] <= %.3f\n",i,i,L[i][i]);
   }
   L[i][i] = sqrt(L[i][i]);
   //printf("\nL[i][i] = %.3f\n",L[i][i]);
   for(j = i+1; j < n; j++)
   {
    L[j][i] = matriz[i][j];
    for(k = 0; k<i; k++)
    {
     L[j][i] -= L[i][k]*L[j][k]; 
    }
    L[j][i] /= L[i][i]; 
   }
  }
  return L;
 }
 else
  return NULL;
}

void deallocateArray(float ** ptr, int row)
{
	for(int i = 0; i < row; i++)
	{
		free(ptr[i]);
	}
	free(ptr);
}

float** criarMatrizIdentidade(int n)
{
  float**identityMatrix = (float**) malloc(n*sizeof(float*));
  for(int i = 0; i < n; i++)
  {
   identityMatrix[i] = (float*) malloc(n*sizeof(float));
   for(int j = 0; j < n; j++)
   {
    if(i == j)
     identityMatrix[i][j] = 1;
    else
     identityMatrix[i][j] = 0;
   }
  }
  return identityMatrix;
}

/*Faz a permutacao entre as linhas quando necessario*/
int swap(float **matriz, float *termosIndependentes,int column, int n)
{
 int i;
 printf("Column:%i\n",column);
 for(i = column+1; i<n; i++) //percorre as linhas da coluna cujo pivot eh zero
 {
  if( matriz[i][column] != 0)
  {
   float aux = 0;
   for(int j= 0; j<n; j++)
   {
    aux = matriz[column][j];
    matriz[column][j] = matriz[i][j];
    matriz[i][j] = aux;
   }
   
   if(termosIndependentes)
   {
    aux = termosIndependentes[column];
    termosIndependentes[column] = termosIndependentes[i];
    termosIndependentes[i] = aux;
   }
   break;
  }  
 }
 printf("\nSWAP OK\n");
 if(i>=n)
  i=column;
 return i;
}

/*De maneira analoga, mas percorre no sentido inverso pro caso da elimicao Gauss Jordan*/
int swap2(float **matriz, float *termosIndependentes,int column, int n)
{
 int i;
 for(i = column-1; i>0; i--) //percorre as linhas da coluna cujo pivot eh zero
 {
  if( matriz[i][column] != 0)
  {
   float aux = 0;
   for(int j= 0; j<n; j++)
   {
    aux = matriz[column][j];
    matriz[column][j] = matriz[i][j];
    matriz[i][j] = aux;
   }
   
   if(termosIndependentes)
   {
    aux = termosIndependentes[column];
    termosIndependentes[column] = termosIndependentes[i];
    termosIndependentes[i] = aux;
   }
   break;
  }  
 }
 if(i<0)
  i=column;
 return i;
}

float** eliminacaoGauss(float **matriz1,float **history,float *termosIndependentes, int n)
{ 
  int j = 0;
  int k = 0;
  float aux = 0;
  int lineChanged;
  float **resultado = criarMatrizIdentidade(n);
  float **permutation = criarMatrizIdentidade(n);
  float **L = criarMatrizIdentidade(n);
  float *termosInd = (float *)malloc(n*sizeof(float));
  for(int i = 0; i<n-1; i++)
  {
   //printf("\nloop\n");
   if(matriz1[i][i] == 0)
   {
     printf("Swap\n");
     //imprimeMatriz(matriz1, n, n);
     //printf("-------------->>\n");
     if(termosIndependentes)
     {
      lineChanged = swap(matriz1, termosIndependentes, i, n);
      //printf("-------------->>line:%i\n",lineChanged);
     }
     else
     { 
      lineChanged = swap(matriz1, NULL, i, n);
     }
     //printf("Ok ate aqui\n");
     for(j= 0; j<n; j++)
     {
      aux = permutation[lineChanged][j];
      permutation[lineChanged][j] = permutation[i][j];
      permutation[i][j] = aux;
      //L[lineChanged][j] = matriz1[i][j];
      //matriz1[i][j] = aux;
      aux = history[lineChanged][j];
      history[lineChanged][j] = history[i][j];
      history[i][j] = aux;
     }
     //imprimeMatriz(matriz1, n, n);
     //printf("-------------->>>>>\n");
   }
    float **identidade = criarMatrizIdentidade(n);
    for(j = i+1; j<n; j++)
    { 
     if(matriz1[i][i] != 0)
     {    
      //printf("Matriz1[%i][%i] = %0.3f\n",j,i,matriz1[j][i]);
      identidade[j][i] = -matriz1[j][i]/matriz1[i][i];
     }
     else
     {
      identidade[j][i] = 0; 
     }
     L[j][i] = -identidade[j][i];
    }
   //printf("M:\n");
   //imprimeMatriz(identidade, n, n);
   // printf("\n\n");
   //imprimeMatriz(matriz1, n, n);
    MultiplicarMatrizes(identidade, matriz1, resultado,n,n,n);
    //printf("\nOK\n");
    //imprimeMatriz(resultado, n, n);
    if(termosIndependentes)
     MultiplicarMatrizVetor(identidade, termosIndependentes, termosInd,n,n); 
    //printf("TermosInd:\n");
    for(k = 0; k < n; k++)
    {
     if(termosIndependentes)
      termosIndependentes[k] = termosInd[k];
     //printf("%.2f\n",termosIndependentes[k]);
    //printf("\nFoi\n");
     for(j = 0; j < n; j++)
      matriz1[j][k] = resultado[j][k];    
    }
   MultiplicarMatrizes(identidade, history, resultado,n,n,n);
   for(k = 0; k < n; k++)
   {
    for(j = 0; j < n; j++)
     history[j][k] = resultado[j][k];    
   }
   deallocateArray(identidade,n);
   //printf("\nHistory:\n");
   //imprimeMatriz(history, n, n);
  }
 MultiplicarMatrizes(permutation, L, resultado,n,n,n);
 for(k = 0; k < n; k++)
 {
  for(j = 0; j < n; j++)
   L[j][k] = resultado[j][k];    
 }
 deallocateArray(resultado,n);
 deallocateArray(permutation,n);
 free(termosInd);
 //printf("L:\n");
 //imprimeMatriz(L,n,n);
 return L;
 //printf("\nSegFault?\n");
}

float** eliminacaoJordan(float **matriz1,float **history, float *termosIndependentes, int n)
{ 
  int j = 0;
  int k = 0;
  float aux = 0;
  int lineChanged;
  float **U = criarMatrizIdentidade(n);
  float **resultado = criarMatrizIdentidade(n);
  float **permutation = criarMatrizIdentidade(n);
  float *termosInd = (float *)malloc(n*sizeof(float));
  
  for(int i = n-1; i>0; i--)
  {
   //printf("\nloop2\n");
   if(matriz1[i][i] == 0)
   {
    // imprimeMatriz(matriz1, n, n);
     //printf("Swap2\n");
     lineChanged = swap2(matriz1, termosIndependentes, i, n);
     for(j= 0; j<n; j++)
     {
      aux = permutation[lineChanged][j];
      permutation[lineChanged][j] = permutation[i][j];
      permutation[i][j] = aux;
      //U[lineChanged][j] = matriz1[i][j];
      //matriz1[i][j] = aux;
      aux = history[lineChanged][j];
      history[lineChanged][j] = history[i][j];
      history[i][j] = aux;
     }
     // imprimeMatriz(matriz1, n, n);
   }
   
    float **identidade = criarMatrizIdentidade(n);
    for(j = i-1; j>=0; j--)
    { 
     if(matriz1[i][i] != 0)    
      identidade[j][i] = -matriz1[j][i]/matriz1[i][i];
     else
      identidade[j][i] = 0;
      U[j][i] = -identidade[j][i];
    }
    //imprimeMatriz(identidade, n, n);
   // printf("\n\n");
    //imprimeMatriz(matriz1, n, n);
    MultiplicarMatrizes(identidade, matriz1, resultado,n,n,n);
    //printf("\nOK\n");
    //imprimeMatriz(resultado, n, n);
    
   if(termosIndependentes)
    MultiplicarMatrizVetor(identidade, termosIndependentes, termosInd,n,n);
    //printf("TermosInd:\n");
    for(k = 0; k < n; k++)
    {
     if(termosIndependentes)
      termosIndependentes[k] = termosInd[k];
    // printf("%.2f",termosIndependentes[k]);
   //  printf("\nFoi\n");
     for(j = 0; j < n; j++)
      matriz1[j][k] = resultado[j][k];    
    }
    MultiplicarMatrizes(identidade, history, resultado,n,n,n);
   for(k = 0; k < n; k++)
   {
    for(j = 0; j < n; j++)
     history[j][k] = resultado[j][k];    
   }
   deallocateArray(identidade,n);
  //printf("\nHistory:\n"); 
  //imprimeMatriz(history, n, n);
 }
 MultiplicarMatrizes(permutation, U, resultado,n,n,n);
 for(k = 0; k < n; k++)
 {
  for(j = 0; j < n; j++)
   U[j][k] = resultado[j][k];    
 }
 deallocateArray(resultado,n);
 deallocateArray(permutation,n);
 free(termosInd);
 //printf("U:\n");
 //imprimeMatriz(U, n, n);
 return U;
 //printf("\nSegFault?\n");
}

float calcularResiduo(float *X1, float* X0, int n, int norma)
{
 float numerador = 0;
 float denominador = 0;
 for(int i = 0; i<n; i++){ 
  numerador += (float)pow( (double)fabs(X1[i] - X0[i]) , (double)norma);
  denominador += (float)pow( (double)fabs(X1[i]), (double)norma ); 
 }
 //printf("Numerador: %.6f\n",numerador);
 //printf("Denominador: %.6f\n",denominador);
 return (float)pow( (double)numerador/denominador, (double)(1/(float)norma) );;
}

/*Vector eh o chute inicial e tambem o X0 para toda iteracao (solucao antiga)*/
void jacobiMethod(float **matriz, float *vector, float *termosIndependentes,int n,int normaDoVetor,float residuoAceitavel)
{
  int i,j,p;
  int ok = p = 0;
  float a = 0;
  if ( isDiagonalDominant(matriz, n) )
  {
   float *solution = (float *)malloc(n*sizeof(float));//Solution eh o X1 para toda iteracao (solucao mais recente)
   while(ok != 1)
   {
    p++;
    for(i=0; i<n; i++)
    {
     solution[i] = termosIndependentes[i];   
     for(j=0; j<n; j++)
     {
      if(i != j)
       solution[i] -= matriz[i][j]*vector[j];
     }
     solution[i] /= matriz[i][i];
     //printf("%.2f\n", solution[i]);
    }
    a = calcularResiduo(solution, vector, n, normaDoVetor);
    if(a < residuoAceitavel)
     ok = 1;
    /*else*/
     //printf("%.3f\n",a);
    for(i = 0; i<n; i++)
    {
     //printf("%.3f\n", solution[i]);
     vector[i] = solution[i];
    }
   }
   free(solution); 
   printf("Numero de iteracoes necessarias: %i\n",p); 
   printf("Residuo Final: %.6f\n",a);
  }
 else
 {
  printf("\nJACOBI NAO CONVERGE\n");
 }
}

void gaussSeidelMethod(float **matriz, float *vector, float *termosIndependentes,int n,int normaDoVetor,float residuoAceitavel)
{
  int i,j,p;
  int ok = p = 0;
  float a = 0; 
  if( isDiagonalDominant(matriz, n) || isPositiveDefinite(matriz, n) )
  {
   float *solution = (float *)malloc(n*sizeof(float));
   while(ok != 1)
   {
    p++;
    for(i=0; i<n; i++)
    {
     solution[i] = termosIndependentes[i];   
     for(j=0; j<n; j++)
     {
      if((i != j) && (j < i))
      {
       //printf("solution[%i](%.3f) -= matriz[%i][%i](%.3f)*solution[j](%.3f)\n\n",i,solution[i],i,j,matriz[i][j],solution[j]);
       solution[i] -= matriz[i][j]*solution[j];
      }
      else if (i != j)
      {
       //printf("solution[%i](%.3f) -= matriz[%i][%i](%.3f)*solution[j](%.3f)\n\n",i,solution[i],i,j,matriz[i][j],vector[j]);
       solution[i] -= matriz[i][j]*vector[j];
      }
     }
     solution[i] /= matriz[i][i];
     //printf("Solucao da iteracao %i \n%.2f\n", p,solution[i]);
    }
    a = calcularResiduo(solution, vector, n, normaDoVetor);
    if(a < residuoAceitavel)
     ok = 1;
    else
     //printf("Residuo: %.6f\n",a);
    for(i = 0; i<n; i++)
    {
     //printf("%.3f\n", solution[i]);
     vector[i] = solution[i];
    }
   }
   free(solution);
   printf("Numero de iteracoes necessarias: %i\n",p); 
   printf("Residuo Final: %.6f\n",a);
  }
 else
  printf("\nGauss-Seidel NAO CONVERGE\n");
}

float** eliminacaoGaussJordan(float **matriz1, float *termosIndependentes, int n)
{ 
  int i,j;
  float **HISTORY = criarMatrizIdentidade(n);
  float **copy = criarMatrizIdentidade(n);
  float **L;
  float **U;
  for(i = 0;i<n;i++)
   for(j = 0;j<n;j++)
    copy[i][j] = matriz1[i][j];
  if(termosIndependentes)
  {
   L = eliminacaoGauss(matriz1, HISTORY,termosIndependentes, n);
   U = eliminacaoJordan(matriz1, HISTORY,termosIndependentes, n);
  }
  else
  {
   printf("Term == NULL\n");
   L = eliminacaoGauss(matriz1, HISTORY,NULL, n);
   U = eliminacaoJordan(matriz1, HISTORY,NULL, n);
  }
  float **result = criarMatrizIdentidade(n);
  float**M = criarMatrizIdentidade(n);
  for(i= 0; i<n; i++)
  {
   if(matriz1[i][i] != 0) 
    M[i][i] /= matriz1[i][i];
   else
    M[i][i] = 1;
  }
  /*printf("\nU\n");
  imprimeMatriz(U, n, n);
  printf("\nL\n");
  imprimeMatriz(L, n, n);
  printf("\nM\n");
  imprimeMatriz(M, n, n);
  printf("\nIDENTITY\n");*/
  //imprimeMatriz(M, n, n);
  MultiplicarMatrizes(M,HISTORY, result,n,n,n);
  MultiplicarMatrizes(result,copy,M,n,n,n);
  //imprimeMatriz(M, n, n);
  //imprimeMatriz(result, n, n);
  deallocateArray(L, n);
  deallocateArray(U, n);
  deallocateArray(M, n);
  deallocateArray(copy, n);
  deallocateArray(HISTORY, n);
  return result; //retorna a inversa de A
  //imprimeMatriz(M, n, n);
}

//Matriz 1 possui m linhas e n colunas e Matriz2 possui o colunas
void MultiplicarMatrizes(float **matriz1, float **matriz2, float **matrizResultado,int m, int n, int o)
{
 int linha1,coluna1,linha2,coluna2;
 for(linha1 = 0;linha1 < m;linha1++)
 {
  for(coluna2 = 0;coluna2 < o;coluna2++)
  {
   matrizResultado[linha1][coluna2] = 0;
   linha2 = 0;
   for(coluna1 = 0;coluna1 < n;coluna1++)
   {
    matrizResultado[linha1][coluna2] += matriz1[linha1][coluna1]*matriz2[linha2][coluna2];
    linha2++;   
   }
  }
 }
}

//Multiplica o vetor coluna por uma matriz
void
MultiplicarMatrizVetor (float **matriz1, float *matriz2, float *matrizResultado, int n, int m)
{
 int linha1,coluna1;
 for(linha1 = 0;linha1 < n;linha1++)
 {
   matrizResultado[linha1]= 0;
   for(coluna1 = 0;coluna1 < m;coluna1++)
   {
    matrizResultado[linha1] += matriz1[linha1][coluna1]*matriz2[coluna1];
   }
 }
}

void imprimeMatriz(float** MatrizVertice, int n, int m){
    int j,i;
    for(i = 0; i < n; i++){
        for(j = 0; j< m; j++){
            printf("%.9f ",MatrizVertice[i][j]);
        }
        printf("\n");
    }
}

int isPositiveDefinite(float**matrix, int n)
{
 if(checkSimmetry(matrix,n) && checkEigenvalues(matrix,n))
  return 1;
 else
 {
  printf("\n\nMATRIZ NAO EH POSITIVA DEFINIDA\n\n");
 } 
  return 0;
}

int checkEigenvalues(float**matrix, int n)
{
 float**copy = criarMatrizIdentidade(n);
 copyMatrix(matrix,copy,n);
 eigenvectorsJacobiMethodTest(copy, n);
 for(int i = 0; i<n; i++)
  if(copy[i][i] <= 0)
   return 0;
 deallocateArray(copy,n);
 return 1;
}

int checkSimmetry(float**matrix, int n)
{
 for(int i = 0; i<n; i++)
  for(int j = 0; j<n; j++)
   if(matrix[i][j] != matrix[j][i])
   {
    printf("\n\nmatrix[%i][%i] = %.9f e matrix[%i][%i] = %.9f\n\n",i,j,matrix[i][j],j,i,matrix[j][i]);
    return 0;
   }
 return 1;
}

int isDiagonalDominant(float **matrix, int n)
{
 float *sumColumnsPerLine = (float*)malloc(n*sizeof(float));
 float *sumLinesPerColumn = (float*)malloc(n*sizeof(float));
 for(int i = 0; i<n; i++)
 {
  sumColumnsPerLine[i] = 0;
  sumLinesPerColumn[i] = 0;
  for(int j = 0; j<n; j++)
  {
   if(j != i)
   {
    sumColumnsPerLine[i] += fabs(matrix[i][j]);
    sumLinesPerColumn[i] += fabs(matrix[j][i]);
   }
  }
  if( (fabs(matrix[i][i]) < sumColumnsPerLine[i]) || (matrix[i][i] < sumLinesPerColumn[i]) )
  {
   printf("\n\nNAO EH DIAGONAL DOMINANTE\n\n");
   return 0;
  }
 }
 return 1;  
}

float det(float **A, int n)
{
 float**copy = criarMatrizIdentidade(n);
 copyMatrix(A,copy,n);
 float tol = 0.00001;
 float result = 1;
 eigenvectorsJacobiMethod(copy, n,tol);
 for(int i = 0; i<n;i++)
  result *= copy[i][i];
 return fabs(result);
}

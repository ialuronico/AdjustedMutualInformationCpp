#include "mex.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <time.h>
#include <limits>
#include <vector>
#include <string>
#include <cstdlib>
#include <time.h>
#include <string.h>
#include <list>


#include <stdlib.h>
#include <string.h>
#include <stdio.h>


using namespace std;

#define MAX(a, b) (((a) > (b) ? (a) : (b)))
#define MIN(a, b) (((a) < (b) ? (a) : (b)))

int** Contingency(double* A, double* B,int n_stateA, int n_stateB, int N);
double Mutu_Info(int** T, int n_stateA,int n_stateB);
double Adjusted_Mutu_Info(int** T, int n_stateA,int n_stateB);
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{

#define AMI_OUT plhs[0]

#define A_IN prhs[0]
#define B_IN prhs[1]


int N;

N = mxGetN(A_IN);

double *A = mxGetPr(A_IN);
double *B = mxGetPr(B_IN);
int n_stateA = 0;
int n_stateB = 0;

// Get number of clusters
for(int i=0;i<N;i++)
   if (n_stateA < A[i]) n_stateA = A[i];

// Get number of clusters
for(int i=0;i<N;i++)
  if (n_stateB < B[i]) n_stateB = B[i];

int** T = Contingency(A, B, n_stateA, n_stateB, N);
/*
for(int i=0;i<n_stateA;i++){
  for(int j=0;j<n_stateB;j++){
    printf("%d ",T[i][j]);
  }
  printf("\n");
}
*/

double AMI = Adjusted_Mutu_Info(T, n_stateA, n_stateB);

AMI_OUT = mxCreateDoubleScalar(AMI);
  
for(int i=0;i<n_stateA;i++)
  free(T[i]);
free(T);

return;

}


/* Adjusted Mutual Information */


int** Contingency(double* A, double* B,int n_stateA,int n_stateB, int N){

  int i,j;

  int** T = (int**)malloc(n_stateA*sizeof(int*));

  for(i=0;i<n_stateA;i++)
    T[i] = (int*)malloc(n_stateB*sizeof(int));

  for(i=0;i<n_stateA;i++)
    for(j=0;j<n_stateB;j++)
       T[i][j]=0;

    for(i =0;i<N;i++){
       T[(int)A[i]-1][(int)B[i]-1]++;
    }
   return T;
}

//get the mutual information from a contingency table
double Mutu_Info(int **T, int n_stateA,int n_stateB){  
    //n_state: #rows n_stateC:#cols
    double MI=0;

    int *a = (int*)malloc(n_stateA*sizeof(int));
    int *b = (int*)malloc(n_stateB*sizeof(int));
    int N = 0;
    int i,j;

    for(i=0;i<n_stateA;i++){ //row sum
        a[i] = 0;
        for(j=0;j<n_stateB;j++)
        {a[i]+=T[i][j];}
    }

    for(i=0;i<n_stateB;i++){ //col sum
        b[i] = 0;
        for(j=0;j<n_stateA;j++)
        {b[i]+=T[j][i];}
    }

    for(i=0;i<n_stateA;i++) {N+=a[i];}

    for(i=0;i<n_stateA;i++){
        for(j=0;j<n_stateB;j++){
            if(T[i][j]>0){
                MI+= T[i][j]*log((double)T[i][j]*N/a[i]/b[j]);
            }
        }
    }

    free(a);
    free(b);
    
    if(N>0) return MI/N;
    else return 0;
}
//get the AMI for a contigency table
double Adjusted_Mutu_Info(int **T, int n_stateA, int n_stateB){ 
    //n_stateA: #rows n_stateB:#cols
    
    double MI=0;
    double Ha=0;
    double Hb=0;

    double *a = (double*)malloc(n_stateA*sizeof(double));
    double *b = (double*)malloc(n_stateB*sizeof(double));
    double N = 0;
    int i,j,i_;

    // Compute marginals and entropy for rows
    for(i=0;i<n_stateA;i++){ //row sum
        a[i]=0;
        for(j=0;j<n_stateB;j++)
          {a[i]+=T[i][j];}
        if (a[i] > 0)
          Ha += a[i]*log(a[i]);
    }
    
    for(i=0;i<n_stateA;i++) {N+=a[i];}
    
    Ha = log(N)-Ha/N;
    
    // Compute marginals and entropy for columns
    for(j=0;j<n_stateB;j++){ //col sum
        b[j]=0;
        for(i=0;i<n_stateA;i++)
          {b[j]+=T[i][j];}
        if (b[j] > 0)
          Hb += b[j]*log(b[j]);
    }
    Hb = log(N) - Hb/N;

    for(i=0;i<n_stateA;i++){
        for(j=0;j<n_stateB;j++){
            if(T[i][j]>0){
                MI+= T[i][j]*log((double)T[i][j]*N/a[i]/b[j]);
            }
        }
    }

    if(N>0) MI = MI/N;
    else return 0;


    // Compute the expected value
    double EMI = 0;
    for(i=0;i<n_stateA;i++){
        for(j=0;j<n_stateB;j++){
          double nij = MAX(1.0,a[i]+b[j]-N);

          double x_1, x_2, v;
          double p = 1;

          // Two linked lists in order to avoid precision errors
          
          // sizetmp is O(max(ai,bj))
          // I can use N as maximum size
          
          int sizetmp = N;                    
          
          double *num = (double*)malloc(sizetmp*sizeof(double));
          double *den = (double*)malloc(sizetmp*sizeof(double));

          // Some simplifications on numerator and denominator
          x_1 = MIN(nij, N - a[i] - b[j] + nij);
          x_2 = MAX(nij, N - a[i] - b[j] + nij);


          int i_num = 0;
          int i_den = 0;
          // numerator
          v = a[i];

          while (v > a[i] - nij)
              num[i_num++] = v--;

          v = b[j];
          while (v > b[j] - nij)
            num[i_num++] = v--;

          // denominator
          v = N;
          while (v > N - a[i])
            den[i_den++] = v--;
          v = x_1;
          while (v > 0)
            den[i_den++] = v--;

          // Due to the simplifications
          if (N - b[j] > x_2)   { // add a term to numerator
              v = N - b[j];
              while (v > x_2)
                num[i_num++] = v--;
          } else {            // add a term to denominator
              v = x_2;
              while (v > N - b[j])
                den[i_den++] = v--;
          }

          for(i_ = 0; i_ < i_num; i_++)
            p*= num[i_]/den[i_];
                    
          free(num);
          free(den);
          num = NULL;
          den = NULL;

          while (nij <= MIN(a[i],b[j])){
            EMI += nij*log((double)nij*N/a[i]/b[j])*p;
            p *= (a[i]-nij)*(b[j]-nij)/(nij+1)/(N-a[i]-b[j]+nij+1);
            nij++;
          }
        }
    }
    EMI = EMI/N;
        
    free(a);
    free(b);
    a = NULL;
    b = NULL;
        
    return (MI - EMI)/(MAX(Ha,Hb) - EMI);    
}


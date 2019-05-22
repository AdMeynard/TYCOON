#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
  double *InData;
  int NRow,NCol;
  double lambda;

  int *Freq;
  double *FreqOut;
  double *Energy, *FVal;
  double minval, val;
  
  int i,j,k,j1,j2;
 
  /* Portal to matlab */
  InData = mxGetPr(prhs[0]);
  NRow = mxGetM(prhs[0]);
  NCol = mxGetN(prhs[0]);

  lambda = mxGetScalar(prhs[1]);

  plhs[0] = mxCreateNumericArray(1,&NRow,mxDOUBLE_CLASS,mxREAL);
  FreqOut = mxGetPr(plhs[0]);

  printf("lambda = %lf\n", lambda);
  printf("NRow = %d\n", NRow);
  printf("NCol = %d\n", NCol);

  /* Main operations start here */
  const double eps = 1e-8;
  double sum = 0;

  /* Initialization */
  Energy = (double *)mxMalloc(sizeof(double)*NRow*NCol);
  FVal = (double *)mxMalloc(sizeof(double)*NRow*NCol*NCol);
  Freq = (int *)mxMalloc(sizeof(int)*NRow);

  for (i=0;i<NRow;++i)
    for (j=0;j<NCol;++j){
      Energy[i*NCol+j] = InData[i+j*NRow];
      sum += Energy[i*NCol+j];
    }

  /* Dynamic programming */
  for (i=0;i<NRow;++i)
    for (j=0;j<NCol;++j)
      Energy[i*NCol+j] = -log(Energy[i*NCol+j]/sum+eps);
  for (j1=0;j1<NCol;++j1)
    for (j2=0;j2<NCol;++j2)
      FVal[NCol*NCol+j1*NCol+j2] = Energy[j1]+Energy[NCol+j2];

  for (i=2;i<NRow;++i){
    for (j1=0;j1<NCol;++j1)
      for (j2=0;j2<NCol;++j2){
        FVal[i*NCol*NCol+j1*NCol+j2] = 1e13;
        for (k=0;k<NCol;++k)
          if (FVal[i*NCol*NCol+j1*NCol+j2] > \
              FVal[(i-1)*NCol*NCol+k*NCol+j1]+lambda*(k+j2-2*j1)*(k+j2-2*j1))
            FVal[i*NCol*NCol+j1*NCol+j2] =  \
              FVal[(i-1)*NCol*NCol+k*NCol+j1]+lambda*(k+j2-2*j1)*(k+j2-2*j1);
        FVal[i*NCol*NCol+j1*NCol+j2] += Energy[i*NCol+j2];
      }
  }

  /* Traceback */
  minval = FVal[(NRow-1)*NCol*NCol];
  Freq[NRow-1] = 0; Freq[NRow-2] = 0;
  for (j1=1;j1<NCol;++j1)
    for (j2=1;j2<NCol;++j2)
      if (FVal[(NRow-1)*NCol*NCol+j1*NCol+j2]<minval){
        minval = FVal[(NRow-1)*NCol*NCol+j1*NCol+j2];
        Freq[NRow-1] = j2; Freq[NRow-2] = j1;
      }
  for (i=NRow-3;i>=0;--i){
    val = FVal[(i+2)*NCol*NCol+Freq[i+1]*NCol+Freq[i+2]] - \
      Energy[(i+2)*NCol+Freq[i+2]];
    for (j=0;j<NCol;++j)
      if (fabs(val - FVal[(i+1)*NCol*NCol+j*NCol+Freq[i+1]] - \
               lambda*(j+Freq[i+2]-2*Freq[i+1])*(j+Freq[i+2]-2*Freq[i+1]))<eps){
        Freq[i] = j;
        break;
      }
  }
          
  for (i=0;i<NRow;++i){
    FreqOut[i] = (double)Freq[i]+1.;
  }


  mxFree(Freq);
  mxFree(FVal);
  mxFree(Energy); 
  return;
}

 

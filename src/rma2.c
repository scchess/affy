/************************************************************************
 **
 ** rma.c
 **
 ** created by: B. M. Bolstad
 ** created on: June 26, 2002
 ** 
 ** last modified: Dec 26, 2002 (Laurent)
 **
 ** License: GPL V2 or later (same as the rest of the Affy package)
 **
 ** version 0.2 - Initial release
 **
 ** Version History
 ** 0.1 - Initial version released on July 14, 2002. Implements median
 **       polish RMA method with 
 ** 0.2 - background implemented in c with the density estimation still carried
 **       out by the R function density()
 ** 0.25 - correct background implementation, version 0.2 is broken.
 **        background is implemented in rma_background.c
 ** 0.30 - Have a copy and none copy path. ie we can either work inplace or on
 **        duplicates. the purpose of this is to reduce memory overhea. For 
 **        someone with an interest only in expression estimates this should not be a problem
 **
 ** a c language implementation of the RMA method as given in the RMA.R file I 
 ** received from Rafael at an earlier point, but assume already had background 
 ** correction to PM's at somepoint (perhaps in the c code) bg will be written in later.
 ** Possibly another background method will be inserted in at this stage.
 **
 ** Note that the normalization code that is used in this algorithm is updated
 ** from that in the affy version 1.1.1 (there will be slight differences in the 
 ** expression values because of this), there is also slight differences in the 
 ** ordering of the results.
 **
 ** Ideally and at some later point a more modular approach that can be called 
 ** in a better manner from R than this will be written. This is a quick and 
 ** dirty approach to get something that will run acceptably.in terms of memory 
 ** and speed. From a software engineering viewpoint expect the code to be poorly 
 ** constructed.
 **
 ** Input to the function should be process from within R
 **
 ** Note that the qnorm code here will not be the development tree
 ** LG: what do you mean ?
 **
 ** Nov 2, 2002 - modify so that it will work efficently with affy2
 ** Nov 3, 2002 - More modifications, remove cruft from old version
 ** Nov 4, 2002 - testing, check docs etc
 ** Nov 10,2002 - remove pesky debug printf()
 **
 ** Dec 26, 2002 - '//' is not a valid way to comment out (and some C compilers complain about it)
 **                (Laurent)
 ************************************************************************/


/* #include "rma_structures.h" */
#include "rma_common.h"
#include "rma_background2.h"
/* #include "qnorm.h" */

#include <R.h> 
#include <Rdefines.h>
#include <Rmath.h>
#include <Rinternals.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void do_RMA(double *PM, char **ProbeNames, int *rows, int * cols,double *results,char **outNames,int nps);


/**************************************************************************
 **
 ** double median(double *x, int length)
 **
 ** double *x - vector
 ** int length - length of *x
 **
 ** returns the median of *x
 **
 *************************************************************************/

double  median(double *x, int length){
  int i;
  int half;
  double med;
  double *buffer = malloc(length*sizeof(double));
  
  for (i = 0; i < length; i++)
    buffer[i] = x[i];
  
  qsort(buffer,length,sizeof(double), (int(*)(const void*, const void*))sort_double);
  half = (length + 1)/2;
  if (length % 2 == 1){
    med = buffer[half - 1];
  } else {
    med = (buffer[half] + buffer[half-1])/2.0;
  }
  
  free(buffer);
  return med;
}

/*******************************************************************************
 **
 ** double sum_abs(double *z, int rows, int cols)
 **
 ** double *z - matrix of doubles
 ** int rows - dimension of matrix
 ** int cols - dimension of matrix
 **
 ** returns the sum of the absolute values of elements of the matrix *z
 **
 ******************************************************************************/

double sum_abs(double *z, int rows, int cols){
 
  int i, j;
  double sum = 0.0;

  for (i=0; i < rows; i++)
    for (j=0; j < cols; j++)
      sum+=fabs(z[j*rows+i]);

  return sum;
}

/********************************************************************************
 **
 ** void get_row_median(double *z, double *rdelta, int rows, int cols)
 **
 ** double *z - matrix of dimension  rows*cols
 ** double *rdelta - on output will contain row medians (vector of length rows)
 ** int rows, cols - dimesion of matrix
 **
 ** get the row medians of a matrix 
 **
 ********************************************************************************/

void get_row_median(double *z, double *rdelta, int rows, int cols){
  int i,j;
  double *buffer = malloc(cols*sizeof(double));

  for (i = 0; i < rows; i++){ 
    for (j = 0; j < cols; j++){
      buffer[j] = z[j*rows + i];
    }
    rdelta[i] = median(buffer,cols);
  }
  
  free(buffer);
}

/********************************************************************************
 **
 ** void get_col_median(double *z, double *cdelta, int rows, int cols)
 **
 ** double *z - matrix of dimension  rows*cols
 ** double *cdelta - on output will contain col medians (vector of length cols)
 ** int rows, cols - dimesion of matrix
 **
 ** get the col medians of a matrix 
 **
 ********************************************************************************/

void get_col_median(double *z, double *cdelta, int rows, int cols){
  
  int i, j;
  
  double *buffer = malloc(rows*sizeof(double));
  for (j = 0; j < cols; j++){
    for (i = 0; i < rows; i++){  
      buffer[i] = z[j*rows + i];
    }
    cdelta[j] = median(buffer,rows);
  }
  
  free(buffer);

}

/***********************************************************************************
 **
 ** void subtract_by_row(double *z, double *rdelta, int rows, int cols)
 ** 
 ** double *z - matrix of dimension rows by cols
 ** double *rdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** subtract the elements of *rdelta off each row of *z
 **
 ***********************************************************************************/

void subtract_by_row(double *z, double *rdelta, int rows, int cols){
  
  int i,j;

  for (i = 0; i < rows; i++){
    for (j = 0; j < cols; j++){
      z[j*rows +i]-= rdelta[i];
    }
  }
}


/***********************************************************************************
 **
 ** void subtract_by_col(double *z, double *cdelta, int rows, int cols)
 ** 
 ** double *z - matrix of dimension rows by cols
 ** double *cdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** subtract the elements of *cdelta off each col of *z
 **
 ***********************************************************************************/

void subtract_by_col(double *z, double *cdelta, int rows, int cols){
  
  int i,j;
  for (j = 0; j < cols; j++){
    for (i = 0; i < rows; i++){
      z[j*rows +i]-= cdelta[j];
    }
  }

}

/***********************************************************************************
 **
 ** void rmod(double *r, double *rdelta, int rows)
 ** 
 ** double *r - vector of length rows
 ** double *rdelta - vector of length rows
 ** int rows, cols dimensions of matrix
 **
 ** add elementwise *rdelta to *r
 **
 ***********************************************************************************/


void rmod(double *r, double *rdelta, int rows){
  int i;

  for (i = 0; i < rows; i++){
    r[i]= r[i] + rdelta[i];
  }
}

/***********************************************************************************
 **
 ** void cmod(double *c, double *cdelta, int cols)
 ** 
 ** double *c - vector of length rows
 ** double *cdelta - vector of length rows
 ** int cols length of vector
 **
 ** add elementwise *cdelta to *c
 **
 ***********************************************************************************/

void cmod(double *c, double *cdelta, int cols){
  int j;

  for (j = 0; j < cols; j++){
    c[j]= c[j] + cdelta[j];
  }
}


/*************************************************************************************
 **
 ** void R_median_polish(double *data, int *rs, int *cs)
 **
 ** a function to do median polish that can be called from R. Main reason is to test 
 ** that c implementation of median polish is correct
 **
 *************************************************************************************/

void R_median_polish(double *data, int *rs, int *cs){
  
  int rows = *rs, cols = *cs, nprobes=*rs;

  int i,j,iter;
  int maxiter = 10;
  double eps=0.01;
  double oldsum = 0.0,newsum = 0.0;
  double t = 0.0;
  double delta;
  double *rdelta = calloc(nprobes,sizeof(double));
  double *cdelta = calloc(cols,sizeof(double));
  
  double *r = calloc(nprobes,sizeof(double));
  double *c = calloc(cols,sizeof(double));

  double *data_matrix = malloc(nprobes*cols*sizeof(double));
  double *z = malloc(nprobes*cols*sizeof(double));
  double *results = malloc(cols*sizeof(double));

  /*   for (i=0; i < nprobes; i++){
    r[i] = 0.0;
  }

  for (j=0; j < nprobes; j++){
    c[j] = 0.0;
    } */


  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] =  log(data[j*nprobes + i])/log(2.0);
      printf("%f ",z[j*nprobes +i]);
    }
    printf("\n");
  } 
  
  
  for (iter = 1; iter <= 10; iter++){
    get_row_median(z,rdelta,nprobes,cols);
    subtract_by_row(z,rdelta,nprobes,cols);
    rmod(r,rdelta,nprobes);
    delta = median(c,cols);
    for (j = 0; j < cols; j++){
      c[j] = c[j] - delta;
    }
    t = t + delta;

    get_col_median(z,cdelta,nprobes,cols);
    subtract_by_col(z,cdelta,nprobes,cols);
    cmod(c,cdelta,cols);
    delta = median(r,nprobes);
    for (i =0; i < nprobes; i ++){
      r[i] = r[i] - delta;
    }
    t = t+delta;
    newsum = sum_abs(z,nprobes,cols);
    if (newsum == 0.0 || fabs(1.0 - oldsum/newsum) < eps)
      break;
    oldsum = newsum;
  }

  printf("%f ",t);
  for (j=0; j < cols; j++){
    results[j] = c[j]; 
    printf("%f ",results[j]);
  }
  
  printf("\n");

   for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      printf("%f ",z[j*nprobes +i]);
    }
    printf("\n");
  } 



  free(rdelta);
  free(r);
  free(cdelta);
  free(c);
  free(z);
  free(data_matrix);


}


/*************************************************************************************
 **
 ** void median_polish(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes)
 **
 ** double *data - a data matrix of dimension rows by cols (the entire PM matrix)
 ** int rows, cols - rows and columns dimensions of matrix
 ** int cur_rows - vector of length nprobes containg row indicies of *data matrix which apply for a 
 **                particular probeset
 ** double *results - a vector of length cols already allocated. on output contains expression values
 ** int nprobes - number of probes in current probeset.
 **
 ** a function to do median polish.
 **
 *************************************************************************************/

void median_polish(double *data, int rows, int cols, int *cur_rows, double *results, int nprobes){

  int i,j,iter;
  int maxiter = 10;
  double eps=0.01;
  double oldsum = 0.0,newsum = 0.0;
  double t = 0.0;
  double delta;
  double *rdelta = Calloc(nprobes,double);
  double *cdelta = Calloc(cols,double);
  
  double *r = Calloc(nprobes,double);
  double *c = Calloc(cols,double);
  double *z = Calloc(nprobes*cols,double);

  for (j = 0; j < cols; j++){
    for (i =0; i < nprobes; i++){
      z[j*nprobes + i] = log(data[j*rows + cur_rows[i]])/log(2.0);  
    }
  } 
  
  
  for (iter = 1; iter <= 10; iter++){
    get_row_median(z,rdelta,nprobes,cols);
    subtract_by_row(z,rdelta,nprobes,cols);
    rmod(r,rdelta,nprobes);
    delta = median(c,cols);
    for (j = 0; j < cols; j++){
      c[j] = c[j] - delta;
    }
    t = t + delta;
    get_col_median(z,cdelta,nprobes,cols);
    subtract_by_col(z,cdelta,nprobes,cols);
    cmod(c,cdelta,cols);
    delta = median(r,nprobes);
    for (i =0; i < nprobes; i ++){
      r[i] = r[i] - delta;
    }
    t = t+delta;
    newsum = sum_abs(z,nprobes,cols);
    if (newsum == 0.0 || fabs(1.0 - oldsum/newsum) < eps)
      break;
    oldsum = newsum;
  }
  
  for (j=0; j < cols; j++){
    results[j] =  t + c[j]; 
  }
  
  Free(rdelta);
  Free(cdelta);
  Free(r);
  Free(c);
  Free(z); 
}

/************************************************************************************
 **
 **  void do_RMA(double *PM, char **ProbeNames, int *rows, int * cols)
 **
 ** double *PM - matrix of dimension rows by cols (probes by chips) should already be 
 **              normalized and background corrected.
 ** char **ProbeNames - Probeset names, one for each probe.
 ** int *rows, *cols - dimensions of matrix
 **
 ** perform the multichip averaging. PM should be background corrected and normalized
 **
 ** assumed that Probes are sorted, by ProbeNames, so that we can just look at 
 ** consecutive rows in PM matrix when doing the median polish
 **
 ** each item is then used to create a matrix that is median polished to give
 ** expression estimates.
 **
 ************************************************************************************/

void do_RMA(double *PM, char **ProbeNames, int *rows, int *cols, double *results, char **outNames, int nps){
  int j = 0;
  int i = 0;
  int k = 0;
  int size;
  char *first;
  int first_ind;

  /* buffers of size 200 should be enough. */

  char *curname=Calloc(200,char);
  int *cur_rows=Calloc(200,int);
  int nprobes=0;

  double *cur_exprs = Calloc(*cols,double);

  double *OLDPM = NULL;

  first = ProbeNames[0];
  first_ind = 0;
  i =0;
  nprobes = 1;
  for (j = 1; j < *rows; j++){
    if ((strcmp(first,ProbeNames[j]) != 0) | (j == (*rows -1))){
      if (j == (*rows -1)){
	nprobes++;
       	for (k = 0; k < nprobes; k++){
	  cur_rows[k] = (j+1 - nprobes)+k; 
	  /* printf("%d ", (j+1 - nprobes)+k); */
	}
      } else {
	for (k = 0; k < nprobes; k++){
	  cur_rows[k] = (j - nprobes)+k; 
	  /* printf("%d ", (j - nprobes)+k); */
	}
      }
      /* printf("%d \n", nprobes); */
      median_polish(PM, *rows, *cols, cur_rows, cur_exprs, nprobes);
      for (k =0; k < *cols; k++){
	results[k*nps + i] = cur_exprs[k];
      } 
      size = strlen(first);
      outNames[i] = Calloc(size+1,char);
      strcpy(outNames[i],first);
      i++;
      first = ProbeNames[j];
      first_ind = j;
      nprobes = 0;
    }
    nprobes++;
  }

  Free(cur_exprs);
  Free(curname);
  Free(cur_rows);
}




/*************************************************************************'
 **
 ** void rma_c_call(double *PM, double *MM, char **ProbeNames, int *rows, int *cols)
 **
 ** a function to actually carry out the RMA method.
 ** this function will only call three functions, one that does background,
 ** one to quantile normalize and one to to the multichip average (median polish)
 **
 ** this is the .Call() interface. reason for this is to avoid having to using
 ** dup=FALSE in a .C() call which won't work with our example
 **
 **
 *************************************************************************/
SEXP rma_c_call(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP norm_flag){ /*, SEXP outvec, SEXP outnamesvec){ */
  
  int rows, cols;
  double *outexpr;
  double *PM,*MM;
  char **outnames;
  char **ProbeNames;
  int i,nprobesets;
  


  SEXP dim1;
  SEXP outvec,outnamesvec;
  SEXP dimnames,names;
  
  PROTECT(dim1 = getAttrib(PMmat,R_DimSymbol)); 
  rows = INTEGER(dim1)[0];
  cols = INTEGER(dim1)[1]; 

  PM = NUMERIC_POINTER(AS_NUMERIC(PMmat));
  MM = NUMERIC_POINTER(AS_NUMERIC(MMmat));
  
  nprobesets=INTEGER(N_probes)[0];
  
  /*  printf("%i\n",nprobesets); */
  /* printf("%d ",INTEGER(norm_flag)[0]); */
  if (INTEGER(norm_flag)[0]){
  /* normalize PM using quantile normalization */
    printf("Normalizing\n");
    qnorm_c(PM,&rows,&cols);
  }

  ProbeNames =malloc(rows*(sizeof(char *)));

  for (i =0; i < rows; i++)
    ProbeNames[i] = CHAR(VECTOR_ELT(ProbeNamesVec,i));
  
  outnames = malloc(nprobesets*sizeof(char *));

  /* PROTECT(outvec = NEW_NUMERIC(nprobesets*cols)); */
  
  PROTECT(outvec = allocMatrix(REALSXP, nprobesets, cols));


  outexpr = NUMERIC_POINTER(outvec);
 	    
  printf("Calculating Expression\n");
  do_RMA(PM, ProbeNames, &rows, &cols,outexpr,outnames,nprobesets);

  UNPROTECT(2);

  /* now lets put names on the matrix */

  PROTECT(dimnames = allocVector(VECSXP,2));
  PROTECT(names = allocVector(STRSXP,nprobesets));
  
  for ( i =0; i < nprobesets; i++)
    SET_VECTOR_ELT(names,i,mkChar(outnames[i]));
  
  SET_VECTOR_ELT(dimnames,0,names);
  setAttrib(outvec, R_DimNamesSymbol, dimnames);
  UNPROTECT(2);

  return outvec;
}

/*******************************************************************************************************************
 **
 ** SEXP rma_c_complete(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP densfunc, SEXP rho)
 **
 ** SEXP PMmat  - PM's
 ** SEXP MMmat - MM's
 ** SEXP ProbeNamesVec - names of probeset for each row
 ** SEXP N_probes  - number of probesets
 ** SEXP densfunc  - density function to use in computation of background
 ** SEXP rho - an r environment
 ** 
 ** Main function to be called from R. Modifies the PM matrix from the parent environment. More dangerous than the
 ** function below, but less memory intensive. This is a function that implements the complete RMA method. ie
 ** background correction, quantile normalization, then expression summarization using median polish
 **
 *******************************************************************************************************************/

SEXP rma_c_complete(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP densfunc, SEXP rho,SEXP norm_flag){
  
  PMmat = bg_correct_c(PMmat,MMmat,densfunc,rho);
  return rma_c_call(PMmat, MMmat, ProbeNamesVec,N_probes,norm_flag);
}

/********************************************************************************************************************
 **
 ** SEXP rma_c_complete_copy(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP densfunc, SEXP rho)
 **
 ** SEXP PMmat   - PM's
 ** SEXP MMmat   - MM's
 ** SEXP ProbeNamesVec - names of probeset for each row
 ** SEXP N_probes  - number of probesets
 ** SEXP densfunc - density function to use in computation of background
 ** SEXP rho - an r environment
 **
 ** Main function to be called from R. Makes a copy of the PM matrix and then works with that. Safer than the 
 ** other function above, but more memory intensive. This is the function that implements the complete RMA method.
 ** ie background correction, quantile normalization, then expression summarization using median polish
 **
 ********************************************************************************************************************/

SEXP rma_c_complete_copy(SEXP PMmat, SEXP MMmat, SEXP ProbeNamesVec,SEXP N_probes,SEXP densfunc, SEXP rho,SEXP norm_flag){
  
  PMmat = bg_correct_c_copy(PMmat,MMmat,densfunc,rho);
  return rma_c_call(PMmat, MMmat, ProbeNamesVec,N_probes,norm_flag);
}

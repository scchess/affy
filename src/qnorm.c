/**********************************************************
 **
 ** file: qnorm.c
 **
 ** aim: A c implementation of the quantile normalization method 
 **
 ** Copyright (C) 2002-2003    Ben Bolstad
 **
 ** written by: B. M. Bolstad  <bolstad@stat.berkeley.edu>
 **
 ** written: Feb 2, 2002
 ** last modified: Apr 19, 2002
 ** 
 ** This c code implements the quantile normalization method
 ** for normalizing high density oligonucleotide data as discussed
 ** in
 **
 ** Bolstad, B. M., Irizarry R. A., Astrand, M, and Speed, T. P. (2003)(2003) 
 ** A Comparison of Normalization Methods for High 
 ** Density Oligonucleotide Array Data Based on Bias and Variance.
 ** Bioinformatics 19,2,pp 185-193
 **
 ** History
 ** Feb 2, 2002 - Intial c code version from original R code
 ** Apr 19, 2002 - Update to deal more correctly with ties (equal rank)
 ** Jan 2, 2003 - Documentation/Commenting updates reformating
 ** Feb 17, 2003 - add in a free(datvec) to qnorm(). clean up freeing of dimat
 ** Feb 25, 2003 - try to reduce or eliminate compiler warnings (with gcc -Wall)
 ** Feb 28, 2003 - update reference to normalization paper in comments
 **
 ***********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*#include "rma_common.c"*/

/*************************************************************
 **
 ** the dataitem record is used to keep track of data indicies 
 ** along with data value when sorting and unsorting in the 
 ** quantile algorithm.
 **
 ************************************************************/

typedef struct{
  double data;
  int rank;
} dataitem;
  

/***********************************************************
 **  
 ** int min(int x1, int x2)							    
 **
 ** returns the minimum of x1 and x2
 **		    
 **********************************************************/

int min(int x1,int x2){
  if (x1 > x2)
    return x2;
  else
    return x1;
}

/**********************************************************
 **
 ** int sort_fn(const void *a1,const void *a2)
 **
 ** a comparison function for sorting objects of the dataitem type.
 **
 **
 **********************************************************/

int sort_fn(const void *a1,const void *a2){
  dataitem *s1, *s2;
  s1 = (dataitem *)a1;
  s2 = (dataitem *)a2;
  
  if (s1->data < s2->data)
    return (-1);
  if (s1 ->data > s2->data)
    return (1);
  return 0;
}


/************************************************************
 **
 ** dataitem **get_di_matrix(double *data, int rows, int cols)
 **
 ** given data  form a matrix of dataitems, each element of
 ** matrix holds datavalue and original index so that 
 ** normalized data values can be resorted to the original order
 **
 ***********************************************************/

dataitem **get_di_matrix(double *data, int rows, int cols){
  int i,j;
  dataitem **dimat;
  /* dataitem *xtmp; */
  
  dimat = (dataitem **)malloc((cols)*sizeof(dataitem *));
  
  if (dimat == NULL){
    printf("\nERROR - Sorry the normalization routine could not allocate adequate memory\n       You probably need more memory to work with a dataset this large\n");
  }

  /* xtmp = malloc(cols*rows*sizeof(dataitem));
     for (j=0; j < cols; j++, xtmp +=rows) dimat[j] = xtmp; */
  
  for (j=0; j < cols; j++)
    dimat[j] = malloc(rows*sizeof(dataitem));



  for (j =0; j < cols; j++)
    for (i =0; i < rows; i++){
      dimat[j][i].data = data[j*rows + i];
      dimat[j][i].rank = i;
    }

  return(dimat); 
}

/************************************************************
 **
 ** double *get_ranks(dataitem *x,int n)
 **
 ** get ranks in the same manner as R does. Assume that *x is
 ** already sorted
 **
 *************************************************************/

void get_ranks(double *rank, dataitem *x,int n){
  int i,j,k;
   
  i = 0;

  while (i < n) {
    j = i;
    while ((j < n - 1) && (x[j].data  == x[j + 1].data))
      j++;
    if (i != j) {
      for (k = i; k <= j; k++)
	rank[k] = (i + j + 2) / 2.0;
    }
    else
      rank[i] = i + 1;
    i = j + 1;
  }
  /*return rank;*/
}

/*********************************************************
 **
 ** void qnorm_c(double *data, int *rows, int *cols)
 **
 **  this is the function that actually implements the 
 ** quantile normalization algorithm. It is called from R
 **
 ********************************************************/

void qnorm_c(double *data, int *rows, int *cols){
  int i,j,ind;
  dataitem **dimat;
  double sum;
  double *row_mean = malloc((*rows)*sizeof(double));
  double *datvec = malloc(*cols*sizeof(double));
  double *ranks = malloc((*rows)*sizeof(double));

  /*# sort original columns */
  
  dimat = get_di_matrix(data, *rows, *cols);
  
  for (j=0; j < *cols; j++){
    qsort(dimat[j],*rows,sizeof(dataitem),sort_fn);
  }

  /*# calculate means */
  
  for (i =0; i < *rows; i++){
    sum = 0.0;
    for (j=0; j < *cols; j++)
      datvec[j] = dimat[j][i].data;
    /*qsort(datvec,*cols,sizeof(double),(int(*)(const void*, const void*))sort_double); */
    for (j=0; j < *cols; j++){
      sum +=datvec[j]/(double)*cols;;
    }
    row_mean[i] = sum; /*/(double)*cols;*/
  }
  
  /*# unsort mean columns */
  for (j =0; j < *cols; j++){
    get_ranks(ranks,dimat[j],*rows);
    for (i =0; i < *rows; i++){
      ind = dimat[j][i].rank;
      data[j*(*rows) + ind] = row_mean[(int)floor(ranks[i])-1];
    }
  }
  free(ranks);
  free(datvec);   

  for (j=0; j < *cols; j++){
    free(dimat[j]);
  }

  free(dimat);
  free(row_mean); 
}


/*********************************************************
 **
 ** void qnorm_robust_c(double *data, int *rows, int *cols)
 **
 ** this is the function that actually implements the 
 ** quantile normalization algorithm. It is called from R. 
 ** this function allows the user to downweight particular
 ** chips, in the calculation of the mean.
 **
 ********************************************************/

void qnorm_robust_c(double *data,double *weights, int *rows, int *cols){
  int i,j,ind;
  dataitem **dimat;
  double sum,sumweights;
  double *row_mean = malloc((*rows)*sizeof(double));
  double *datvec = malloc(*cols*sizeof(double));
  double *ranks = malloc((*rows)*sizeof(double));

  dimat = get_di_matrix(data, *rows, *cols);
  
  for (j=0; j < *cols; j++){
    qsort(dimat[j],*rows,sizeof(dataitem),sort_fn);
  }

  for (i =0; i < *rows; i++){
    sum = 0.0;
    for (j=0; j < *cols; j++)
      datvec[j] = dimat[j][i].data;
    /* qsort(datvec,*cols,sizeof(double),(int(*)(const void*, const void*))sort_double); */
    for (j=0; j < (*cols); j++){
      sum +=weights[j]*datvec[j];
    }
    sumweights = 0.0;
    for (j=0; j < (*cols); j++){
      sumweights = sumweights + weights[j];
    }
    row_mean[i] = sum/sumweights;
  }
  
  /*# unsort mean columns */
  for (j =0; j < *cols; j++){
    get_ranks(ranks,dimat[j],*rows);


    for (i =0; i < *rows; i++){
      ind = dimat[j][i].rank;
      data[j*(*rows) + ind] = row_mean[(int)floor(ranks[i])-1];
    }
  }
  free(ranks);
  free(dimat);
  free(row_mean); 

}

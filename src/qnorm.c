/* 
   qnorm.c

   A c implementation of the quantile normalization method 

   B. M. Bolstad

   written: Feb 2, 2001

*/

#include <stdio.h>
#include <stdlib.h>

typedef struct{
  double data;
  int rank;
} dataitem;
  

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


int sort_double(const double *a1,const double *a2){
  if (*a1 < *a2)
    return (-1);
  if (*a1 > *a2)
    return (1);
  return 0;
}

dataitem **get_di_matrix(double *data, int rows, int cols){
  int i,j;
  dataitem **dimat;
  dataitem *xtmp;
  
  dimat = (dataitem **)malloc((cols)*sizeof(dataitem *));
  

  xtmp = malloc(cols*rows*sizeof(dataitem));

  
  for (j=0; j < cols; j++, xtmp +=rows) dimat[j] = xtmp;
  
  for (j =0; j < cols; j++)
    for (i =0; i < rows; i++){
      dimat[j][i].data = data[j*rows + i];
      dimat[j][i].rank = i;
    }
  
  return(dimat); 
}


void qnorm_c(double *data, int *rows, int *cols){
  int i,j,ind;
  dataitem **dimat;
  double sum;
  double *row_mean = malloc((*rows)*sizeof(double));
  double *datvec = malloc(*cols*sizeof(double));


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
    qsort(datvec,*cols,sizeof(double),(int(*)(const void*, const void*))sort_double);
    for (j=0; j < *cols; j++){
      sum +=datvec[j]/(double)*cols;;
    }
    row_mean[i] = sum; 
  }
  
  /*# unsort mean columns */
  
  for (i =0; i < *rows; i++){
    for (j =0; j < *cols; j++){
      ind = dimat[j][i].rank;
      data[j*(*rows) + ind] = row_mean[i];

    }
  }
  free(dimat);
  free(row_mean); 
  free(datvec);
}


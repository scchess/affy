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
  

int min(int x1,int x2){
  if (x1 > x2)
    return x2;
  else
    return x1;
}

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



  /* printf("fred is %i\n",*rows); */
  
/*
  for (i=0; i < *rows; i++){
    for (j=0; j < *cols; j++){
      printf("%f ",data[j*(*rows) + i]);
    }
    printf("\n");
  }
  printf("\n"); */
  /*# sort original columns */
  
  dimat = get_di_matrix(data, *rows, *cols);
  
  /* qsort(data, (*rows),sizeof(double),(int(*)(const void*, const void*))comp_doub); */

/*  for (i=0; i < *rows; i++){
    for (j=0; j < *cols; j++){
      printf("%i ",dimat[j][i].rank);
    }
    printf("\n");
  }
  printf("\n"); */

  

  for (j=0; j < *cols; j++){
    qsort(dimat[j],*rows,sizeof(dataitem),sort_fn);
  }

/*  for (i=0; i < *rows; i++){
    for (j=0; j < *cols; j++){
      printf("%f ",dimat[j][i].data);
    }
    printf("\n");
  }
  printf("\n");

  for (i=0; i < *rows; i++){
    for (j=0; j < *cols; j++){
      printf("%i ",dimat[j][i].rank);
    }
    printf("\n");
  }

  printf("\n"); */
  /*# calculate means */
  
  for (i =0; i < *rows; i++){
    sum = 0.0;
    for (j=0; j < *cols; j++)
      datvec[j] = dimat[j][i].data;
    qsort(datvec,*cols,sizeof(double),(int(*)(const void*, const void*))sort_double);
    for (j=0; j < *cols; j++){
      sum +=datvec[j]/(double)*cols;;
      /*printf("%f ",datvec[j]); */
    }
    /*printf("\n"); */
    row_mean[i] = sum; /*/(double)*cols;*/
/*   printf("%f \n",row_mean[i]); */
  }
  
  /*# unsort mean columns */
  
  for (i =0; i < *rows; i++){
    for (j =0; j < *cols; j++){
      ind = dimat[j][i].rank;
/*      printf("%i \n",ind);*/
      data[j*(*rows) + ind] = row_mean[i];

    }
  }
  free(dimat);
  free(row_mean); 

/*  for (i=0; i < *rows; i++){
    for (j=0; j < *cols; j++){
      printf("%f ",data[j*(*rows) + i]);
    }
    printf("\n");
  }
  printf("\n");*/
}



void qnorm_robust_c(double *data,double *weights, int *rows, int *cols){
  int i,j,ind;
  dataitem **dimat;
  double sum,sumweights;
  double *row_mean = malloc((*rows)*sizeof(double));
  double *datvec = malloc(*cols*sizeof(double));

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
  
  for (i =0; i < *rows; i++){
    for (j =0; j < *cols; j++){
      ind = dimat[j][i].rank;
      data[j*(*rows) + ind] = row_mean[i];
    }
  }
  free(dimat);
  free(row_mean); 
}

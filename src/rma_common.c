/***********************************************************************
 **
 ** file: rma_common.c
 ** 
 ** aim: a location for commonly used utility functions
 **
 **
 ** written by: B. M. Bolstad
 **
 ** created: Oct 16, 2002
 ** last modified: Oct 16, 2002
 **
 ** history:
 ** Oct 16, 2002 - a place to put common utility code, created to help
 **                the R package build.
 **
 ***********************************************************************/

#include "rma_common.h"

/**********************************************************

 int sort_double(const void *a1,const void *a2)
 
 a sorting function for sorting doubles.

 not currently used

**********************************************************/

int sort_double(const double *a1,const double *a2){
  if (*a1 < *a2)
    return (-1);
  if (*a1 > *a2)
    return (1);
  return 0;
}



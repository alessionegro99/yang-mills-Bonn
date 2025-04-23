#ifndef ALIGNCHECK_H
#define ALIGNCHECK_H

#include"macro.h"

#include<stdlib.h>
#include<stdio.h>

inline void is_aligned(const void * p, size_t byte_align, char* file, int line)
  {
  if( (unsigned long)p % byte_align != 0)
    {
    fprintf(stderr, "Error in alignement (%s, %d)\n", file, line);
    exit(EXIT_FAILURE);
    }
  }

#endif // ALIGNCHECK_H

#include <stdio.h>
#include <stdlib.h>


void write_counts(char filename [], int* x, int* y, double* counts, int n_lines){
  FILE * fd = fopen(filename, "w");
  unsigned int i = 0;
  while(i != n_lines){
    // if counts is 0, do not write
    if(counts[i] == 0){
      i++;
      continue;
    }
    fprintf(fd, "%d\t%d\t%f\n", x[i], y[i], counts[i]);
    i++;
  }
  fclose(fd);
}


void write_counts_int(char filename [], int* x, int* y, int* counts,
                      int n_lines){
  FILE * fd = fopen(filename, "w");
  unsigned int i = 0;
  while(i != n_lines){
    // if counts is 0, do not write
    if(counts[i] == 0){
      i++;
      continue;
    }
    fprintf(fd, "%d\t%d\t%d\n", x[i], y[i], counts[i]);
    i++;
  }
  fclose(fd);
}

#include <stdio.h>
#include <stdlib.h>


int get_num_line(char filename []){
  FILE * fd = fopen(filename, "r");
  char line [1024];
  unsigned int n_lines = 0;
  while(fgets(line, 1024, fd) != NULL){
    n_lines++;
  }
  fclose(fd);
  return n_lines;
}

void read_counts(char filename [], int* array){
  FILE * fd = fopen(filename, "r");
  char line[1024];

  unsigned int i = 0;
  while(fgets(line, 1024, fd) != NULL){
    sscanf(line, "%d %d %d",
           &array[i * 3], &array[i * 3 + 1], &array[i * 3 + 2]);
      i++;
  }
  fclose(fd);
}

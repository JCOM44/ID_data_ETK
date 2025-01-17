#ifndef MEMORY_H
#define MEMORY_H

#define ARRAY_END 1
#define FREE_ARG char*

#define MDIV 301 
#define SDIV 601

void AllocationErrorSF(char error_text[]);
double **array_allocateSF(long nrl, long nrh, long ncl, long nch);
void array_freeSF(double **m, long nrl, long ncl);
void check_status_error(int status);
#endif

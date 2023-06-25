#ifndef _ARRAY_UTILS_H_
#define _ARRAY_UTILS_H_

#include <iostream>

#ifndef VIE_INF
#define VIE_INF 10000000
#endif

// Function declarations for creating arrays
int **create2DArray(int d1, int d2);
int ***create3DArray(int d1, int d2, int d3);
int ****create4DArray(int d1, int d2, int d3, int d4);
int *****create5DArray(int d1, int d2, int d3, int d4, int d5);
int ******create6DArray(int d1, int d2, int d3, int d4, int d5, int d6);

// Function declarations for deallocating arrays
void delete2DArray(int **arr, int d1);
void delete3DArray(int ***arr, int d1, int d2);
void delete4DArray(int ****arr, int d1, int d2, int d3);
void delete5DArray(int *****arr, int d1, int d2, int d3, int d4);
void delete6DArray(int ******arr, int d1, int d2, int d3, int d4, int d5);

#endif

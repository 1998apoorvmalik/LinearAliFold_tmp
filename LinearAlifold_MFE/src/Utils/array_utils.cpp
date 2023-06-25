#include "array_utils.h"

// Create arrays
int **create2DArray(int d1, int d2)
{
    int **arr = new int *[d1];
    for (int i = 0; i < d1; ++i)
    {
        arr[i] = new int[d2];
        std::fill_n(arr[i], d2, VIE_INF);
    }
    return arr;
}

int ***create3DArray(int d1, int d2, int d3)
{
    int ***arr = new int **[d1];
    for (int i = 0; i < d1; ++i)
    {
        arr[i] = create2DArray(d2, d3);
    }
    return arr;
}

int ****create4DArray(int d1, int d2, int d3, int d4)
{
    int ****arr = new int ***[d1];
    for (int i = 0; i < d1; ++i)
    {
        arr[i] = create3DArray(d2, d3, d4);
    }
    return arr;
}

int *****create5DArray(int d1, int d2, int d3, int d4, int d5)
{
    int *****arr = new int ****[d1];
    for (int i = 0; i < d1; ++i)
    {
        arr[i] = create4DArray(d2, d3, d4, d5);
    }
    return arr;
}

int ******create6DArray(int d1, int d2, int d3, int d4, int d5, int d6)
{
    int ******arr = new int *****[d1];
    for (int i = 0; i < d1; ++i)
    {
        arr[i] = create5DArray(d2, d3, d4, d5, d6);
    }
    return arr;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////// Function definitions for deallocating arrays //////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void delete2DArray(int **arr, int d1)
{
    for (int i = 0; i < d1; ++i)
    {
        delete[] arr[i];
    }
    delete[] arr;
}

void delete3DArray(int ***arr, int d1, int d2)
{
    for (int i = 0; i < d1; ++i)
    {
        delete2DArray(arr[i], d2);
    }
    delete[] arr;
}

void delete4DArray(int ****arr, int d1, int d2, int d3)
{
    for (int i = 0; i < d1; ++i)
    {
        delete3DArray(arr[i], d2, d3);
    }
    delete[] arr;
}

void delete5DArray(int *****arr, int d1, int d2, int d3, int d4)
{
    for (int i = 0; i < d1; ++i)
    {
        delete4DArray(arr[i], d2, d3, d4);
    }
    delete[] arr;
}

void delete6DArray(int ******arr, int d1, int d2, int d3, int d4, int d5)
{
    for (int i = 0; i < d1; ++i)
    {
        delete5DArray(arr[i], d2, d3, d4, d5);
    }
    delete[] arr;
}

#include "Tools.h"

int arraySize(int arr[]) {
    return (sizeof(arr) / sizeof(*arr));
}

int arraySize(char* arr[]) {
    return (sizeof(arr) / sizeof(*arr));
}

double GetCurrentTimeSec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + tv.tv_usec * 1e-6;
    //return 1.0;
}

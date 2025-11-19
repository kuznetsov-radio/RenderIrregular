#pragma once

#ifdef LINUX
#define finite isfinite
#else
#define finite _finite
#endif

inline double sqr(double x)
{
 return x*x;
}

inline int min(int a, int b)
{
 return (a<b) ? a : b;
}

inline int max(int a, int b)
{
 return (a>b) ? a : b;
}

inline double max(double a, double b)
{
 return (a>b) ? a : b;
}

inline double min(double a, double b)
{
 return (a<b) ? a : b;
}

inline void arrswap(double *a, int i, int j)
{
 double tmp=a[i];
 a[i]=a[j];
 a[j]=tmp;
}

inline void arrswap(int *a, int i, int j)
{
 int tmp=a[i];
 a[i]=a[j];
 a[j]=tmp;
}

int value_locate(double*, int, double);
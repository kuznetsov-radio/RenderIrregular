int value_locate(double *a, int n, double x)
{
 if (x<a[0]) return -1;
 if (x>=a[n-1]) return n-1;

 int j, j1, l;
 j=0; 
 j1=n-1; 
 while (j1-j>1) 
 { 
  l=(j+j1)>>1; 
  if (a[l]>x) j1=l; 
  else j=l; 
 } 
 return j;
} 
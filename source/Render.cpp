#include <stdlib.h>
#include <math.h>
#include <float.h>
#include "ExtMath.h"
#include "Render.h"
#ifndef LINUX
#include <ppl.h>
#else
#include <omp.h>
#endif

int CheckCrossPixel(double xl, double xr, double yb, double yt, double x1, double x2, double y1, double y2)
{
 if ((x1>xl) && (x1<=xr) && (y1>yb) && (y1<=yt)) return 1; //1st point inside
 
 if ((x2>xl) && (x2<=xr) && (y2>yb) && (y2<=yt)) return 1; //2nd point inside
 
 double bx=x2-x1;
 if (bx==0.0) if ((x1>xl) && (x1<=xr) && (min(y1, y2)<=yb) && (max(y1, y2)>=yt)) return 1; //vertical crossing
                                                                            else return 0;
 
 double by=y2-y1;
 if (by==0.0) if ((y1>yb) && (y1<=yt) && (min(x1, x2)<=xl) && (max(x1, x2)>=xr)) return 1; //horizontal crossing
                                                                            else return 0;
 
 double yc=y1+(xl-x1)*by/bx;
 if (finite(yc)) if ((yc>=min(y1, y2)) && (yc<=max(y1, y2)) && (yc>yb) && (yc<=yt)) return 1; //crosses the left boundary
 
 yc=y1+(xr-x1)*by/bx;
 if (finite(yc)) if ((yc>=min(y1, y2)) && (yc<=max(y1, y2)) && (yc>yb) && (yc<=yt)) return 1; //crosses the right boundary
 
 double xc=x1+(yb-y1)*bx/by;
 if (finite(xc)) if ((xc>=min(x1, x2)) && (xc<=max(x1, x2)) && (xc>xl) && (xc<=xr)) return 1; //crosses the bottom boundary
 
 xc=x1+(yt-y1)*bx/by;
 if (finite(xc)) if ((xc>=min(x1, x2)) && (xc<=max(x1, x2)) && (xc>xl) && (xc<=xr)) return 1; //crosses the top boundary
 
 return 0;
}

int RenderMain(int Nx, int Ny, int Nz, double dx, double dy, double *dz, 
               double x1, double x2, double y1, double y2, double z1, double z2, 
               int VoxOn, int OnedOn, int arrN, __int32 *VoxelID, 
               __int32 *NVoxels, __int32 *VoxList, double *ds, 
               double *x_ind, double *y_ind, double *z_ind, double *entry_point, double *exit_point)
{
 double *tmidarr=(double*)malloc(sizeof(double)*arrN);
 int *iarr=(int*)malloc(sizeof(int)*Nx*Ny);
 int *jarr=(int*)malloc(sizeof(int)*Nx*Ny);
 double *zbarr=(double*)malloc(sizeof(double)*(Nz+1));
 int *Ncrossings=(int*)malloc(sizeof(int)*Nz);
 double *tarr=(double*)malloc(sizeof(double)*Nz*2);

 int res=0;

 //----------------------------------------------

 double L=sqrt(sqr(x2-x1)+sqr(y2-y1)+sqr(z2-z1));
 
 double bx=x2-x1;
 double by=y2-y1;
 double bz=z2-z1;
 
 double tcmin= 1e100;
 double tcmax=-1e100;
 
 *NVoxels=0;
 
 int imin=max(int(floor(min(x1, x2)/dx)), 0);
 int imax=min(int(ceil(max(x1, x2)/dx)), Nx-1);
 int jmin=max(int(floor(min(y1, y2)/dy)), 0);
 int jmax=min(int(ceil(max(y1, y2)/dy)), Ny-1);
 
 if (bx==0.0) imin=max(imin-1, 0);
 if (by==0.0) jmin=max(jmin-1, 0);
 
 if ((imax>=imin) && (jmax>=jmin))
 {
  int Nxy=0;
 
  for (int i=imin; i<=imax; i++)
  {
   double xl=dx*i;
   double xr=dx*(i+1);
  
   for (int j=jmin; j<=jmax; j++)
   {
    double yb=dy*j;
    double yt=dy*(j+1);
   
    if (CheckCrossPixel(xl, xr, yb, yt, x1, x2, y1, y2))
	{
     iarr[Nxy]=i;
     jarr[Nxy]=j; 
     Nxy+=1;
	}
   }
  }
 
  if (Nxy>0)
  {
   for (int m=0; m<Nxy; m++)
   {
    int i=iarr[m];
    int j=jarr[m];

    double xl=dx*i;
    double xr=dx*(i+1);
    double yb=dy*j;
    double yt=dy*(j+1);
    
    zbarr[0]=0.0;
    for (int k=1; k<=Nz; k++) zbarr[k]=zbarr[k-1]+dz[i+(j+(k-1)*Ny)*Nx];
    
    for (int k=0; k<Nz; k++) Ncrossings[k]=0;
    
    if (bx!=0.0) 
	{
     double tc=(xl-x1)/bx;
     double yc=y1+tc*by;
     double zc=z1+tc*bz;
     
     if ((yc>yb) && (yc<=yt) && (zc>=zbarr[0]) && (zc<=zbarr[Nz]))
	 {
      int r=value_locate(zbarr, Nz+1, zc);
      if ((r>=0) && (r<Nz))
	  {
       if (Ncrossings[r]<2) tarr[r+Ncrossings[r]*Nz]=tc;
       Ncrossings[r]++;
       tcmin=min(tcmin, tc);
       tcmax=max(tcmax, tc); 
	  }
	 }
     
     tc=(xr-x1)/bx;
     yc=y1+tc*by;
     zc=z1+tc*bz;
     
     if ((yc>yb) && (yc<=yt) && (zc>=zbarr[0]) && (zc<=zbarr[Nz]))
	 {
      int r=value_locate(zbarr, Nz+1, zc);
      if ((r>=0) && (r<Nz))
	  {
       if (Ncrossings[r]<2) tarr[r+Ncrossings[r]*Nz]=tc;
       Ncrossings[r]++;
       tcmin=min(tcmin, tc);
       tcmax=max(tcmax, tc); 
	  }
	 }
	}
    
    if (by!=0.0)
	{
     double tc=(yb-y1)/by;
     double xc=x1+tc*bx;
     double zc=z1+tc*bz;
     
     if ((xc>xl) && (xc<=xr) && (zc>=zbarr[0]) && (zc<=zbarr[Nz]))
	 {
      int r=value_locate(zbarr, Nz+1, zc);
      if ((r>=0) && (r<Nz))
	  {
       if (Ncrossings[r]<2) tarr[r+Ncrossings[r]*Nz]=tc;
       Ncrossings[r]++;
       tcmin=min(tcmin, tc);
       tcmax=max(tcmax, tc); 
	  }
	 }
     
     tc=(yt-y1)/by;
     xc=x1+tc*bx;
     zc=z1+tc*bz;
     
     if ((xc>xl) && (xc<=xr) && (zc>=zbarr[0]) && (zc<=zbarr[Nz]))
	 {
      int r=value_locate(zbarr, Nz+1, zc);
      if ((r>=0) && (r<Nz))
	  {
       if (Ncrossings[r]<2) tarr[r+Ncrossings[r]*Nz]=tc;
       Ncrossings[r]++;
       tcmin=min(tcmin, tc);
       tcmax=max(tcmax, tc); 
	  }
	 }
	}
    
    if (bz!=0.0) for (int k=0; k<=Nz; k++)
	{
     double tc=(zbarr[k]-z1)/bz;
     double xc=x1+tc*bx;
     double yc=y1+tc*by ;
     
     if ((xc>xl) && (xc<=xr) && (yc>yb) && (yc<=yt)) for (int r=max(k-1, 0); r<=min(k, Nz-1); r++)
	 {
      if (Ncrossings[r]<2) tarr[r+Ncrossings[r]*Nz]=tc;
      Ncrossings[r]++;
      tcmin=min(tcmin, tc);
      tcmax=max(tcmax, tc); 
	 }
	}
    
    for (int r=0; r<Nz; r++) if ((*NVoxels<arrN) && (Ncrossings[r]==2))
	{
     VoxList[*NVoxels]=r*Nx*Ny+j*Nx+i;
     ds[*NVoxels]=fabs(tarr[r]-tarr[r+Nz])*L;
     tmidarr[*NVoxels]=(tarr[r]+tarr[r+Nz])/2;
     
     double xmid=x1+tmidarr[*NVoxels]*bx;
     double ymid=y1+tmidarr[*NVoxels]*by;
     double zmid=z1+tmidarr[*NVoxels]*bz;
     x_ind[*NVoxels]=min(max(xmid/dx-0.5, 0.0), double(Nx-1));
     y_ind[*NVoxels]=min(max(ymid/dy-0.5, 0.0), double(Ny-1));
     z_ind[*NVoxels]=min(max(r+(zmid-zbarr[r])/(zbarr[r+1]-zbarr[r])-0.5, 0.0), double(Nz-1));
     
     if (VoxOn && !OnedOn) if (VoxelID[i+(j+r*Ny)*Nx] & 1)
	 {
      x_ind[*NVoxels]=i;
      y_ind[*NVoxels]=j;
      z_ind[*NVoxels]=r;
	 } 
	 else
	 {
	  int ke=0;
	  int vox0=VoxelID[i+(j+r*Ny)*Nx];
	  for (int i1=max(i-1, 0); i1<=min(i+1, Nx-1); i1++) for (int j1=max(j-1, 0); j1<=min(j+1, Ny-1); j1++) for (int k1=max(r-1, 0); k1<=min(r+1, Nz-1); k1++)
	   if (VoxelID[i1+(j1+k1*Ny)*Nx]!=vox0) ke++;
      if (ke!=0) 
	  {
       x_ind[*NVoxels]=i;
       y_ind[*NVoxels]=j;
       z_ind[*NVoxels]=r;
	  }
	 }
	
     (*NVoxels)++; 
     if (*NVoxels>=arrN) res=1;
	}
   }
  }
 }
 
 if ((*NVoxels>0) && (*NVoxels<arrN))
 {
  for (int p=0; p<(*NVoxels-1); p++) for (int q=p+1; q<(*NVoxels); q++) if (tmidarr[p]>tmidarr[q])
  {
   arrswap(tmidarr, p, q);
   arrswap(VoxList, p, q);
   arrswap(ds, p, q);
   arrswap(x_ind, p, q);
   arrswap(y_ind, p, q);
   arrswap(z_ind, p, q);
  }
  
  if (VoxOn && OnedOn) for (int i=0; i<(*NVoxels); i++) if (VoxelID[VoxList[i]] & 1)
  {
   x_ind[i]=round(x_ind[i]);
   y_ind[i]=round(y_ind[i]);
   z_ind[i]=round(z_ind[i]);
  } 
  else 
  {
   int ke=0;
   for (int i1=max(i-1, 0); i1<=min(i+1, *NVoxels-1); i1++) if (VoxelID[VoxList[i1]]!=VoxelID[VoxList[i]]) ke++;
   if (ke!=0) 
   {
    x_ind[i]=round(x_ind[i]);
    y_ind[i]=round(y_ind[i]);
    z_ind[i]=round(z_ind[i]);
   }
  }
  
  entry_point[0]=x1+tcmin*bx; 
  entry_point[1]=y1+tcmin*by; 
  entry_point[2]=z1+tcmin*bz;
  exit_point[0]=x1+tcmax*bx;
  exit_point[1]=y1+tcmax*by;
  exit_point[2]=z1+tcmax*bz;
 }

 if (*NVoxels>=arrN) *NVoxels=0;

 //---

 free(tmidarr);
 free(iarr);
 free(jarr);
 free(zbarr);
 free(Ncrossings);
 free(tarr);

 return res;
}
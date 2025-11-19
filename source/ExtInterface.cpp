#include <stdlib.h>
#include "Render.h"
#ifndef LINUX
#include <ppl.h>
#else
#include <omp.h>
#endif

#define D2(s1, i1, i2) ((i1)+(i2)*(s1))

#ifndef LINUX
extern "C" __declspec(dllexport) int RENDER(int argc, void **argv)
#else
extern "C" int RENDER(int argc, void **argv)
#endif
{
 __int32 *Lparms=(__int32*)argv[0];
 double *pixsize=(double*)argv[1];
 double *dz=(double*)argv[2];
 double *LOS=(double*)argv[3];
 __int32 *VoxelID=(__int32*)argv[4];

 __int32 *NVoxels=(__int32*)argv[5];
 __int32 *VoxList=(__int32*)argv[6];
 double *VoxData=(double*)argv[7];
 double *epoints=(double*)argv[8];

 int Nx=Lparms[0];
 int Ny=Lparms[1];
 int Nz=Lparms[2];
 int arrN=Lparms[3];
 int VoxOn=Lparms[4];
 int OnedOn=Lparms[5];
 double dx=pixsize[0];
 double dy=pixsize[1];
 double x1=LOS[0];
 double y1=LOS[1];
 double z1=LOS[2];
 double x2=LOS[3];
 double y2=LOS[4];
 double z2=LOS[5];

 double *ds=(double*)malloc(sizeof(double)*arrN);
 double *x_ind=(double*)malloc(sizeof(double)*arrN);
 double *y_ind=(double*)malloc(sizeof(double)*arrN);
 double *z_ind=(double*)malloc(sizeof(double)*arrN);
 double entry_point[3]{};
 double exit_point[3]{};

 int res=RenderMain(Nx, Ny, Nz, dx, dy, dz, 
                    x1, x2, y1, y2, z1, z2, 
                    VoxOn, OnedOn, arrN, VoxelID, 
                    NVoxels, VoxList, ds, 
                    x_ind, y_ind, z_ind, entry_point, exit_point);

 for (int i=0; i<arrN; i++)
 {
  VoxData[D2(arrN, i, 0)]=ds[i];
  VoxData[D2(arrN, i, 1)]=x_ind[i];
  VoxData[D2(arrN, i, 2)]=y_ind[i];
  VoxData[D2(arrN, i, 3)]=z_ind[i];
 }

 for (int i=0; i<3; i++)
 {
  epoints[D2(3, i, 0)]=entry_point[i];
  epoints[D2(3, i, 1)]=exit_point[i];
 }

 free(ds);
 free(x_ind);
 free(y_ind);
 free(z_ind);

 return res;
}

#ifndef LINUX
extern "C" __declspec(dllexport) int RENDER_MULTI(int argc, void **argv)
#else
extern "C" int RENDER_MULTI(int argc, void **argv)
#endif
{
 __int32 *Lparms_M=(__int32*)argv[0];

 int arrN=Lparms_M[3];
 int NLOS=Lparms_M[6];

 int *res_arr=(int*)malloc(sizeof(int)*NLOS);

 #ifndef LINUX
 concurrency::parallel_for(0, NLOS, [&](int l)
 #else
 #pragma omp parallel for
 for (int l=0; l<NLOS; l++)
 #endif
 {
  void *ARGV[9];
  
  ARGV[0]=argv[0]; //Lparms
  ARGV[1]=argv[1]; //pixsize
  ARGV[2]=argv[2]; //dz
  ARGV[3]=(void*)(((double*)argv[3])+l*6); //LOS
  ARGV[4]=argv[4]; //VoxelID
  ARGV[5]=(void*)(((__int32*)argv[5])+l); //NVoxels
  ARGV[6]=(void*)(((__int32*)argv[6])+l*arrN); //VoxList
  ARGV[7]=(void*)(((double*)argv[7])+l*arrN*4); //VoxData
  ARGV[8]=(void*)(((double*)argv[8])+l*6); //epoints

  res_arr[l]=RENDER(23, ARGV);
 #ifndef LINUX
 });
 #else
 }
 #endif

 int res=0;
 for (int i=0; i<NLOS; i++) if (res_arr[i]) res=1;

 free(res_arr);

 return res;
}
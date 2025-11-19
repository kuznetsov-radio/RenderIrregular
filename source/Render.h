#pragma once

#ifdef LINUX
#define __int32 int32_t
#endif

int RenderMain(int Nx, int Ny, int Nz, double dx, double dy, double *dz, 
               double x1, double x2, double y1, double y2, double z1, double z2, 
               int VoxOn, int OnedOn, int arrN, __int32 *VoxelID, 
               __int32 *NVoxels, __int32 *VoxList, double *ds, 
               double *x_ind, double *y_ind, double *z_ind, double *entry_point, double *exit_point);
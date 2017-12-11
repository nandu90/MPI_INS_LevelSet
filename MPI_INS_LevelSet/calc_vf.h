/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   calc_vf.h
 * Author: nsaini3
 *
 * Created on November 23, 2016, 1:11 PM
 */

#ifndef CALC_VF_H
#define CALC_VF_H

void calc_vf(double ***phi, double *init_vf, double *vf, double *err)
{
  int i,j;
    double eps=epsilon*max(xlen/(gxelem), ylen/(gyelem));
    double ***H;
    allocator3(&H,xelem,yelem, zelem);
    
    heavy_func(H,  phi, eps);
    
    double localvf;
    (*vf)=0.0;
    for(i=2; i<xelem-2; i++)
    {
        for(j=2; j<yelem-2; j++)
        {
	  localvf += (1.0 - H[i][j][0])*area[i][j][0][0]*area[i][j][1][1];
        }
    }
    
    double totvol = xlen*ylen;
    MPI_Allreduce(&localvf,vf,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    (*vf) = (*vf)*100.0/totvol;
    
    if((*init_vf) == 0.0)
    {
      (*init_vf) = (*vf);
    }
    //if(myrank==master)printf("\n%.6f %.6f %.6f\n",(*init_vf),(*vf),eps);
    (*err) = ((*vf) - (*init_vf))*100.0/(*init_vf);
    
    deallocator3(&H,xelem,yelem, zelem);
}

#endif /* CALC_VF_H */


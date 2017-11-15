/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   bub_advect.h
 * Author: nsaini3
 *
 * Created on October 18, 2016, 7:42 PM
 */

#ifndef BUB_ADVECT_H
#define BUB_ADVECT_H




void bub_advect(struct elemsclr sclr, int iter, double deltat)
{


    ///Interpolate velocity at cell edges to cell centers
  double ***ucen;
  allocator3(&ucen,xelem,yelem,zelem);

  double ***vcen;
  allocator3(&vcen,xelem,yelem,zelem);
  
    for(int i=1; i < xelem-1; i++)
    {
        for(int j=1; j < yelem-1; j++)
        {
            ucen[i][j][0]=0.5*(sclr.u[i][j][0] + sclr.u[i-1][j][0]);
            vcen[i][j][0]=0.5*(sclr.v[i][j][0] + sclr.v[i][j-1][0]);
        }
    }
    cell_center_vel_BC(ucen,vcen);



//    Main bubble iteration loop
    if(iter == 0)
    {
        output_vtk(sclr,iter);
    }




    double **rhsx;
    double **rhsy;
    allocator(&rhsx, xelem, yelem);
    allocator(&rhsy, xelem, yelem);
    //Calculate the fluxes
    rhs_bub(rhsx, rhsy, ucen, vcen, sclr.phi);

    double ***phistar;
    allocator3(&phistar, xelem, yelem, zelem);
    
    #pragma omp parallel for schedule(dynamic)
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            phistar[i][j][0] = sclr.phi[i][j][0] + deltat * (rhsx[i][j] + rhsy[i][j]);
        }
    }

    level_setBC(phistar);

    double **rhstarx;
    double **rhstary;
    allocator(&rhstarx, xelem, yelem);
    allocator(&rhstary, xelem, yelem);
    //Calculate the star fluxes
    rhs_bub(rhstarx, rhstary, ucen, vcen, phistar);
    #pragma omp parallel for schedule(dynamic)
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            //<<phi[i][j][0];
            sclr.phi[i][j][0] = sclr.phi[i][j][0] + 0.5*deltat *(rhsx[i][j] + rhstarx[i][j] + rhsy[i][j] + rhstary[i][j]);
            //<<" "<<phi[i][j][0]<<endl;
        }
    }

    level_setBC(sclr.phi);

    deallocator(&rhstarx, xelem, yelem);
    deallocator(&rhstary, xelem, yelem);
    deallocator3(&phistar, xelem, yelem, zelem);
    deallocator(&rhsx, xelem, yelem);
    deallocator(&rhsy, xelem, yelem);
    deallocator3(&ucen,xelem,yelem,zelem);
    deallocator3(&vcen,xelem,yelem,zelem);
    
}

#endif /* BUB_ADVECT_H */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   pressure_solver.h
 * Author: nsaini3
 *
 * Created on November 3, 2016, 6:39 PM
 */

#ifndef PRESSURE_SOLVER_H
#define PRESSURE_SOLVER_H

#include "gauss_siedel.h"

void pressure(double ***ustar, double ***vstar, double ***p, double deltat)
{
    double rho=1.0;

    /****Calculate the RHS of matrix****/
    double **b;
    allocator(&b, xelem, yelem);
    #pragma omp parallel for schedule(dynamic)
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            double hx = area[i][j][1][1];
            double hy = area[i][j][0][0];
            b[i][j] = (rho*hx*hy/deltat)*(hy*(ustar[i][j][0] - ustar[i-1][j][0]) + hx*(vstar[i][j][0]-vstar[i][j-1][0]));

        }
    }

    double ***a;
    allocator3(&a, xelem, yelem, 5);
    #pragma omp parallel for schedule(dynamic)
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            double hx = area[i][j][1][1];
            double hy = area[i][j][0][0];
            a[i][j][0] = hy*hy;
            a[i][j][1] = hx*hx;
            a[i][j][2] = -2.0*(hx*hx+hy*hy);
            a[i][j][3] = hx*hx;
            a[i][j][4] = hy*hy;
        }
    }
    //<<a[1][1][0]<<" "<<a[1][1][1]<<" "<<a[1][1][2]<<" "<<a[1][1][3]<<" "<<a[1][1][4]<<" "<<endl;
    //exit(0);
    gs_solver(a,b,p);
    deallocator3(&a, xelem, yelem, 5);
    deallocator(&b, xelem, yelem);
}

#endif /* PRESSURE_SOLVER_H */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   surface_tension.h
 * Author: nsaini3
 *
 * Created on October 31, 2016, 12:08 PM
 */

#ifndef SURFACE_TENSION_H
#define SURFACE_TENSION_H

void surface(struct elemsclr sclr, double ***st_forcex, double ***st_forcey)
{
    /**Compute eps based on grid size*/
    double eps = epsilon*max(xlen/(xelem-2), ylen/(yelem-2));
    
    /**Compute Heavyside function**/
    double ***H;
    allocator3(&H, xelem, yelem, zelem);
    
    heavy_func(H, sclr.phi, eps);
    
    find_density_visc(H, sclr.rho, sclr.mu);
    //<<eps<<endl;
    
    
    
    if(sf_toggle == 1)
    {
      double ***delta;
      allocator3(&delta, xelem, yelem, zelem);
        for(int i=1; i<xelem-1; i++)
        {
            for(int j=1; j<yelem-1; j++)
            {
                if(fabs(sclr.phi[i][j][0]) > eps)
                {
                    delta[i][j][0] = 0.0;
                }
                else
                {
                    //delta[i][j][0] = 3.0*exp(-3.0*sclr.phi[i][j][0]/eps)/(1.0 * pow(1.0 + exp(-3.0*sclr.phi[i][j][0]/eps),2.0));
                    delta[i][j][0] = (1.0/2.0*eps) * (1.0 + cos(PI * sclr.phi[i][j][0]/eps));
                }
                //<<delta[i][j][0]<<" ";
            }
            //<<endl;
        }
        //exit(0);
	double ***del_scaling;
	allocator3(&del_scaling, xelem, yelem, zelem);

         //This is equivalent to marker function
         for(int i=1 ; i<xelem-1; i++)
         {
             for(int j=1; j<yelem-1; j++)
             {
                 del_scaling[i][j][0] = 2.0*H[i][j][0]*delta[i][j][0];
                 //<<del_scaling[i][j][0]<<" ";
             }
             //<<endl;
         }
         //exit(0);

         /***Now calculate gradient of level set for ultimately calculating curvature*/
	 double ***grad_phix;
	 double ***grad_phiy;
	 double ***phiRface;
	 double ***phiTface;
	 allocator3(&phiRface, xelem, yelem, zelem);
	 allocator3(&phiTface, xelem, yelem, zelem);
	 allocator3(&grad_phiy, xelem, yelem, zelem);
	 allocator3(&grad_phix, xelem, yelem, zelem);
         for(int i=0; i<xelem-1; i++)
         {
             for(int j=0; j<yelem-1; j++)
             {
                 phiRface[i][j][0] = 0.5*(sclr.phi[i+1][j][0] + sclr.phi[i][j][0]);
                 phiTface[i][j][0] = 0.5*(sclr.phi[i][j+1][0] + sclr.phi[i][j][0]);
             }
         }

         for(int j=1; j<yelem-1; j++)
         {
             for(int i=1; i<xelem-1; i++)
             {
                 grad_phix[i][j][0] = (phiRface[i][j][0] - phiRface[i-1][j][0])/area[i][j][1][1];
                 grad_phiy[i][j][0] = (phiRface[i][j][0] - phiRface[i][j-1][0])/area[i][j][0][0];
                 //if(delta[i][j][0] != 0.0){
                 //<<grad_phix[i][j][0]<<" ";
                 //}
             }
             //<<endl;
         }
         //exit(0);
         /*Need to impose BC for grad_phix and grad_phiy*/
         grad_level_setBC(grad_phix);
         grad_level_setBC(grad_phiy);

         /***Compute double and mixed derivatives***/
	 double ***grad_phixx;
	 double ***grad_phixy;
	 double ***grad_phiyy;
	 double ***phixRface;
	 double ***phiyTface;
	 double ***phixTface;
	 allocator3(&grad_phixx, xelem, yelem, zelem);
	 allocator3(&grad_phixy, xelem, yelem, zelem);
	 allocator3(&grad_phiyy, xelem, yelem, zelem);
	 allocator3(&phixRface, xelem, yelem, zelem);
	 allocator3(&phiyTface, xelem, yelem, zelem);
	 allocator3(&phixTface, xelem, yelem, zelem);

         for(int i=0; i<xelem - 1; i++)
         {
             for(int j=0; j<yelem-1; j++)
             {
                 phixRface[i][j][0] = 0.5*(grad_phix[i][j][0] + grad_phix[i+1][j][0]);
                 phiyTface[i][j][0] = 0.5*(grad_phiy[i][j][0] + grad_phiy[i][j+1][0]);
                 phixTface[i][j][0] = 0.5*(grad_phix[i][j][0] + grad_phix[i][j+1][0]);
             }
         }

         for(int i=1; i<xelem-1; i++)
         {
             for(int j=1; j<yelem-1; j++)
             {
                 grad_phixx[i][j][0] = (phixRface[i][j][0] - phixRface[i-1][j][0])/area[i][j][1][1];
                 grad_phiyy[i][j][0] = (phiyTface[i][j][0] - phiyTface[i][j-1][0])/area[i][j][0][0];
                 grad_phixy[i][j][0] = (phixTface[i][j][0] - phixTface[i][j-1][0])/area[i][j][0][0];
             }
         }
	 double ***curvature;
	 allocator3(&curvature, xelem, yelem, zelem);


         for(int j=1; j<yelem-1; j++)
         {
             for(int i=1; i<xelem-1; i++)
             {
                 //<<i<<" "<<j<<endl;
                 curvature[i][j][0] = (pow(grad_phiy[i][j][0],2.0)*grad_phixx[i][j][0] + pow(grad_phix[i][j][0],2.0)*grad_phiyy[i][j][0] - 2.0*grad_phix[i][j][0]*grad_phiy[i][j][0]*grad_phixy[i][j][0]);
                 curvature[i][j][0] = curvature[i][j][0]/pow((pow(grad_phix[i][j][0],2.0) + pow(grad_phiy[i][j][0],2.0)),1.5);
                 /*if(delta[i][j][0] != 0)
                 {
                 <<curvature[i][j][0]<<" ";
                 }*/
             }
             //<<endl;

         }
         //exit(0);

         for(int j=1; j<yelem-1; j++)
         {
             for(int i=1; i<xelem-1; i++)
             {
                 st_forcex[i][j][0] = sf_coeff*curvature[i][j][0]*del_scaling[i][j][0]*grad_phix[i][j][0];
                 st_forcey[i][j][0] = sf_coeff*curvature[i][j][0]*del_scaling[i][j][0]*grad_phiy[i][j][0];
                 //<<st_forcey[i][j][0]<<" ";
             }
             //<<endl;
         }
         //exit(0);
         grad_level_setBC(st_forcex);
         grad_level_setBC(st_forcey);

	 deallocator3(&curvature, xelem, yelem, zelem);
	 deallocator3(&grad_phixx, xelem, yelem, zelem);
	 deallocator3(&grad_phixy, xelem, yelem, zelem);
	 deallocator3(&grad_phiyy, xelem, yelem, zelem);
	 deallocator3(&phixRface, xelem, yelem, zelem);
	 deallocator3(&phiyTface, xelem, yelem, zelem);
	 deallocator3(&phixTface, xelem, yelem, zelem);
	 deallocator3(&grad_phix, xelem, yelem, zelem);
	 deallocator3(&grad_phiy, xelem, yelem, zelem);
	 deallocator3(&phiRface, xelem, yelem, zelem);
	 deallocator3(&phiTface, xelem, yelem, zelem);
	 deallocator3(&del_scaling, xelem, yelem, zelem);
	 deallocator3(&delta, xelem, yelem, zelem);
	 
    }
    deallocator3(&H, xelem, yelem, zelem);
}

#endif /* SURFACE_TENSION_H */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   re_distance.h
 * Author: nsaini3
 *
 * Created on November 7, 2016, 8:24 PM
 */

#ifndef RE_DISTANCE_H
#define RE_DISTANCE_H

#include "rhs_bub_redist.h"
void monitor_res_redist(double *ires, bool *exitflag, int iter, double ***phi,  double ***phitemp)
{
  int i,j;
    double res=0.0;
    for(i=1; i<xelem-1; i++)
    {
        for(j=1; j<yelem-1; j++)
        {
            res=res + pow(phi[i][j][0]-phitemp[i][j][0],2.0)*vol[i][j];
        }
    }
    
    double buf = res;
    MPI_Allreduce(&buf,&res,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    res=sqrt(res);
    
    
    
    if(iter == 0)
    {
      (*ires)=res;
    }
    else
    {
        //<<"Step: "<<iter<<" phi residual: "<<res/ires<<endl;//<<" V vel residual: "<<res[1]/ires[1]<<endl;
        
      if(res/(*ires) < tol)
        {
            *exitflag=true;
        }
    }
    
}



void re_distance(struct elemsclr sclr)
{
  int i,j;
    /***Store phi values in a separate matrix***/
  double ***phi2;
  allocator3(&phi2, xelem, yelem, zelem);
    for(i=0; i< xelem; i++)
    {
        for(j=0; j< yelem; j++)
        {
            phi2[i][j][0] = sclr.phi[i][j][0];
        }
    }
    
   
    /**Compute eps based on grid size*/
    double eps = epsilon*max(xlen/(xelem-2), ylen/(yelem-2));
    
     double deltat=re_time;
    
    double ires=0.0;  
    
    bool exitflag = false;
    
    int iter;
    for(iter=0; iter < re_loops; iter++)
    {
        /**Compute Heavyside function**/
      double ***H;
      allocator3(&H, xelem, yelem, zelem);
        //heavy(H,phi2,eps);
      double ***signnew;
      allocator3(&signnew, xelem, yelem, zelem);
        for(i=0; i<xelem; i++)
        {
            for(j=0; j<yelem; j++)
            {
                signnew[i][j][0] = 2.0*(H[i][j][0] -0.5);
            }
        }
        double ***temp_phi2;
	allocator3(&temp_phi2, xelem, yelem, zelem);

        for(i=1;i<xelem-1;i++)
        {
            for(j=1;j<yelem-1;j++)
            {
                temp_phi2[i][j][0]=phi2[i][j][0];
            }
        }
        
        //<<iter<<endl;
        /*Calculate convection velocites*/
	double ***grad_phix;
	double ***grad_phiy;
	double ***phiRface;
	double ***phiTface;
	allocator3(&grad_phix, xelem, yelem, zelem);
	allocator3(&grad_phiy, xelem, yelem, zelem);
	allocator3(&phiRface, xelem, yelem, zelem);
	allocator3(&phiTface, xelem, yelem, zelem);

        for(i=0; i<xelem-1; i++)
        {
            for(j=0; j<yelem-1; j++)
            {
                phiRface[i][j][0] = 0.5*(phi2[i+1][j][0] + phi2[i][j][0]);
                phiTface[i][j][0] = 0.5*(phi2[i][j+1][0] + phi2[i][j][0]);
            }
        }

        for(j=1; j<yelem-1; j++)
        {
            for(i=1; i<xelem-1; i++)
            {
                grad_phix[i][j][0] = (phiRface[i][j][0] - phiRface[i-1][j][0])/area[i][j][1][1];
                grad_phiy[i][j][0] = (phiRface[i][j][0] - phiRface[i][j-1][0])/area[i][j][0][0];
            }
        }
        periodicBC(grad_phix);
        periodicBC(grad_phiy);
        gradBC(grad_phix);
        gradBC(grad_phiy);

	double ***ucen;
	double ***vcen;
	allocator3(&ucen, xelem, yelem, zelem);
	allocator3(&vcen, xelem, yelem, zelem);
        
        for(j=1; j<yelem-1; j++)
        {
            for(i=1; i<xelem-1; i++)
            {
                double mag_phi = sqrt(pow(grad_phix[i][j][0],2.0) + pow(grad_phiy[i][j][0],2.0));
                ucen[i][j][0] = signnew[i][j][0]*grad_phix[i][j][0]/mag_phi;
                vcen[i][j][0] = signnew[i][j][0]*grad_phiy[i][j][0]/mag_phi;
                //<<grad_phiy[i][j][0]<<" ";
            }
            //<<endl;
        }
        //exit(0);
        
        periodicBC(ucen);
        periodicBC(vcen);
        zerogradBC(ucen); //Note velocity at the wall should not be zero
        zerogradBC(vcen); //Note velocity at the wall should not be zero
        
        
        
        
        /*****Now onto calculating fluxes******/
	double **rhsx;
	double **rhsy;
	allocator(&rhsx, xelem, yelem);
	allocator(&rhsy, xelem, yelem);

        
        rhs_bub(rhsx, rhsy, ucen, vcen, phi2);

	double ***phistar;
	allocator3(&phistar, xelem, yelem, zelem);

        for(i=1; i<xelem-1; i++) 
        {
            for(j=1; j<yelem-1; j++)
            {
                phistar[i][j][0] = phi2[i][j][0] + deltat * (signnew[i][j][0] + rhsx[i][j] + rhsy[i][j]);
            }
        }
        periodicBC(phistar);
        zerogradBC(phistar);
        //bothscalarBC(phistar);

	double **rhstarx;
	double **rhstary;
	allocator(&rhstarx, xelem, yelem);
	allocator(&rhstary, xelem, yelem);
        //Calculate the star fluxes
        rhs_bub(rhstarx, rhstary, ucen, vcen, phistar);

	double ***Hstar;
	allocator3(&Hstar, xelem, yelem, zelem);

//        /heavy(Hstar,phistar,eps);
	double ***signnew2;
	allocator3(&signnew2, xelem, yelem, zelem);
	

        for(i=0; i<xelem; i++)
        {
            for(j=0; j<yelem; j++)
            {
                signnew2[i][j][0] = 2.0*(Hstar[i][j][0] -0.5);
            }
        }
        
        
        for(i=1; i<xelem-1; i++)
        {
            for(j=1; j<yelem-1; j++)
            {
                phi2[i][j][0] = phistar[i][j][0] + 0.5*deltat *(signnew[i][j][0] + signnew2[i][j][0] + rhsx[i][j] + rhstarx[i][j] + rhsy[i][j] + rhstary[i][j]);
                //<<" "<<phi[i][j][0]<<endl;
            }
        }
        
        periodicBC(phi2);
        zerogradBC(phi2);
        
        if(exitflag == false)
        {
            monitor_res_redist(&ires, &exitflag, iter,  phi2,  temp_phi2);
        }
        
	deallocator3(&H, xelem, yelem, zelem);
	deallocator3(&signnew, xelem, yelem, zelem);
	deallocator3(&temp_phi2, xelem, yelem, zelem);
	deallocator3(&grad_phix, xelem, yelem, zelem);
	deallocator3(&grad_phiy, xelem, yelem, zelem);
	deallocator3(&phiRface, xelem, yelem, zelem);
	deallocator3(&phiTface, xelem, yelem, zelem);
	deallocator3(&ucen, xelem, yelem, zelem);
	deallocator3(&vcen, xelem, yelem, zelem);
	deallocator3(&phistar, xelem, yelem, zelem);
	deallocator3(&Hstar, xelem, yelem, zelem);
	deallocator3(&signnew2, xelem, yelem, zelem);
	deallocator(&rhsx, xelem, yelem);
	deallocator(&rhsy, xelem, yelem);
	deallocator(&rhstarx, xelem, yelem);
	deallocator(&rhstary, xelem, yelem);
	
    }
    
    
    /*Reassign values*/
    for(i=0; i< xelem; i++)
    {
        for(j=0; j< yelem; j++)
        {
            sclr.phi[i][j][0] = phi2[i][j][0];
        }
    }
    deallocator3(&phi2, xelem, yelem, zelem);
}

#endif /* RE_DISTANCE_H */


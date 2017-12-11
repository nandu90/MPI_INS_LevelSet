/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   re_distance2.h
 * Author: nsaini3
 *
 * Created on November 8, 2016, 12:17 AM
 */

#ifndef HYPERBOLIC_H
#define HYPERBOLIC_H


void hyperbolic(struct elemsclr sclr)
{
  //printf("%d is here1\n",myrank);
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
    double eps=epsilon*max(xlen/(gxelem), ylen/(gyelem));
    
     double deltat=re_time;
    
    double ires=0.0;  
    
    bool exitflag = false;
    
    /****Heavyside and delta functions for volume constraint***/
    double ***H, ***delta, ***grad_phi;
    allocator3(&H, xelem, yelem, zelem);
    allocator3(&delta, xelem, yelem, zelem);
    allocator3(&grad_phi, xelem, yelem, zelem);
    heavy_func(H,sclr.phi,eps);
    delta_func(delta,sclr.phi,eps);
    grad_func(grad_phi, sclr.phi);
    
    
    int iter;
    //printf("%d is here2\n",myrank);
    for(iter=0; iter < re_loops; iter++)
    {
      if(myrank == master)printf("    Re-dist iter = %d\n",iter+1);
      double ***lambda;
      allocator3(&lambda, xelem, yelem, zelem);

      double ***temp_phi2;
       allocator3(&temp_phi2, xelem, yelem, zelem);
       
        for(i=2;i<xelem-2;i++)
        {
            for(j=2;j<yelem-2;j++)
            {
                temp_phi2[i][j][0]=phi2[i][j][0];
            }
        }
        
                
        /*****Now onto calculating fluxes******/
	double **rhs;
	allocator(&rhs, xelem, yelem);
        
        rhs_redist2(rhs, phi2, sclr.phi);
        
        
        double ***phistar;
	allocator3(&phistar, xelem, yelem, zelem);
        
        for(i=2; i<xelem-2; i++) 
        {
            for(j=2; j<yelem-2; j++)
            {
                phistar[i][j][0] = phi2[i][j][0] + deltat * (rhs[i][j]);
            }
        }
        
        /*if(vf_control == 1)
        {
            vol_contraint(phistar, sclr.phi, grad_phi, delta, deltat);
	    }*/
        
	//printf("%d is here3\n",myrank);
        //bothscalarBC(phistar);
	commu2(phistar);
        level_setBC(phistar);

	//printf("%d id here4\n",myrank);
	double **rhstar;
	allocator(&rhstar, xelem, yelem);
	
        //Calculate the star fluxes
        rhs_redist2(rhstar, phistar, sclr.phi);
        
        for(i=2; i<xelem-2; i++)
        {
            for(j=2; j<yelem-2; j++)
            {
                phi2[i][j][0] = phistar[i][j][0] + 0.5*deltat *(rhs[i][j] + rhstar[i][j]);
                //<<" "<<phi[i][j][0]<<endl;
            }
        }
        //printf("%d is here5\n",myrank);
        //level_setBC(phi2);
        if(vf_control == 1)
        {
            vol_contraint(phi2, sclr.phi, grad_phi, delta, deltat);
        }
	commu2(phi2);
        level_setBC(phi2);
        //printf("%d is here6\n",myrank);
        //bothscalarBC(phistar);
        /****Apply volume constraint****/
        
        /*if(exitflag == false)
        {
            monitor_res_redist(&ires, &exitflag, iter,  phi2,  temp_phi2);
        }
        else
        {
            break;
	    }*/
        
	//printf("%d is here7\n",myrank);
        deallocator3(&lambda, xelem, yelem, zelem);
	deallocator3(&temp_phi2, xelem, yelem, zelem);
	deallocator3(&phistar, xelem, yelem, zelem);
	deallocator(&rhs, xelem, yelem);
	deallocator(&rhstar, xelem, yelem);
    }
    
    
    /*Reassign values*/
    for(i=0; i< xelem; i++)
    {
        for(j=0; j< yelem; j++)
        {
            sclr.phi[i][j][0] = phi2[i][j][0];
        }
    }

    commu2(sclr.phi);
    level_setBC(sclr.phi);
    deallocator3(&H, xelem, yelem, zelem);
    deallocator3(&delta, xelem, yelem, zelem);
    deallocator3(&grad_phi, xelem, yelem, zelem); 
    deallocator3(&phi2,xelem,yelem,zelem);
}

#endif /* RE_DISTANCE2_H */


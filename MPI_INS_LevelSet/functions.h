/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   functions.h
 * Author: nsaini3
 *
 * Created on September 27, 2016, 7:44 PM
 */

#ifndef FUNCTIONS_H
#define FUNCTIONS_H

void calcp(struct elemsclr sclr)
{
  int i,j,k;
    for(i=0; i<xelem; i++)
    {
        for(j=0; j<yelem; j++)
        {
            for(k=0; k<zelem; k++)
            {
                sclr.p[i][j][k] = (-4800*xc[i][j] + 96)*nu;
            }
        }
        //<<endl;
    }
    
    /*for(int i=1; i < xelem-1; i++)
    {
        sclr.p[i][0][0]= sclr.p[i][1][0];
        sclr.p[i][yelem-1][0] = sclr.p[i][yelem-2][0];
    }*/
    
    /*for(int j=1; j < yelem-1; j++)
    {
        sclr.p[0][j][0] = sclr.p[xelem-2][j][0];
        sclr.p[xelem-1][j][0] = sclr.p [1][j][0];
    }*/
}

void monitor_res(double *ires, bool *exitflag, int iter, struct elemsclr sclr, double ***utemp,  double ***vtemp)
{
  int i,j,k;
  double res[3];
  for (i=0; i<3; i++)
    {
      res[i] = 0.0;
    }
    for(i=1; i<xelem-1; i++)
    {
        for(j=1; j<yelem-1; j++)
        {
            res[0]=res[0] + pow(sclr.u[i][j][0]-utemp[i][j][0],2.0)*vol[i][j];
            res[1]=res[1] + pow(sclr.v[i][j][0]-vtemp[i][j][0],2.0)*vol[i][j];
        }
    }
    
    for(i=0; i<3; i++)
    {
        res[i]=sqrt(res[i]);
    }
    
    
    if(iter == 0)
    {
        for(i=0; i<3; i++)
        {
            ires[i]=res[i];
            //<<ires[i]<<endl;
        }
    }
    else
    {
      printf(" U vel residual: %.6f V vel residual: %.6f",res[0]/ires[0],res[1]/ires[1]);
      //<<" U vel residual: "<<res[0]/ires[0]<<" V vel residual: "<<res[1]/ires[1];
        
        if(res[0]/ires[0] < tol && res[1]/ires[1] <  tol)
        {
            *exitflag=true;
        }
    }
    
}


void timestep_calc(struct elemsclr sclr, double *deltat, double *cfl)
{
    /*****Find the existing maximum cfl of the domain***/
  int i,j,k;
  (*cfl) = 0.0;
    int reqi, reqj;
    for(i=1; i<xelem-1; i++)
    {
        for(j=1; j<yelem-1; j++)
        {
            double surf_int = area[i][j][0][0]* (fabs(sclr.u[i][j][0]) + fabs(sclr.u[i-1][j][0]));
            surf_int = surf_int + area[i][j][1][1]* (fabs(sclr.v[i][j][0]) + fabs(sclr.v[i][j-1][0]));
            
            /*double surf_int = area[i][j][0][0]* (fabs(sclr.u[i][j][0] - sclr.u[i-1][j][0]));
            surf_int = surf_int + area[i][j][1][1]* (fabs(sclr.v[i][j][0] - sclr.v[i][j-1][0]));*/
            
            double temp_cfl = 0.0;
            if((*cfl) == 0.0)
            {
	      (*cfl) = surf_int*(*deltat)/(area[i][j][0][0]*area[i][j][1][1]);
	      reqi = i;
	      reqj = j;
            }
            else
            {
	      double temp_cfl = surf_int*(*deltat)/(area[i][j][0][0]*area[i][j][1][1]);
	      if(temp_cfl > (*cfl))
                {
		  (*cfl) = temp_cfl;
                    reqi = i;
                    reqj = j;
                }
            }
        }
    }
    
    struct
    {
      double buf;
      int loc;
    }in,out;

    in.loc = myrank;
    in.buf = (*cfl);

    MPI_Allreduce(&in,&out,1,MPI_DOUBLE_INT,MPI_MAXLOC,MPI_COMM_WORLD);
    
    (*cfl) = out.buf;
    if((*cfl) < max_cfl)
    {
      (*cfl) = (*cfl) + 0.01*(max_cfl - (*cfl));
    }
    else
    {
      (*cfl) = max_cfl;
    }
    
    double surf_int = area[reqi][reqj][0][0]* (fabs(sclr.u[reqi][reqj][0]) + fabs(sclr.u[reqi-1][reqj][0]));
    surf_int = surf_int + area[reqi][reqj][1][1]* (fabs(sclr.v[reqi][reqj][0]) + fabs(sclr.v[reqi][reqj-1][0]));
    
    /*double surf_int = area[reqi][reqj][0][0]* (fabs(sclr.u[reqi][reqj][0] - sclr.u[reqi-1][reqj][0]));
    surf_int = surf_int + area[reqi][reqj][1][1]* (fabs(sclr.v[reqi][reqj][0] - sclr.v[reqi][reqj-1][0]));*/
    
    (*deltat) = *cfl * area[reqi][reqj][0][0]*area[reqi][reqj][1][1]/surf_int;

    MPI_Bcast(deltat,1,MPI_DOUBLE,out.loc,MPI_COMM_WORLD);
    
    
    
    //deltat = advect_deltat;
    
}
#endif /* FUNCTIONS_H */


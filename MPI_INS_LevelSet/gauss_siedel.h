/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   gauss_siedel.h
 * Author: nsaini3
 *
 * Created on November 3, 2016, 7:18 PM
 */

#ifndef GAUSS_SIEDEL_H
#define GAUSS_SIEDEL_H

void gs_solver(double ***a, double **b, double ***p)
{
  //printf("%d reached pressure solver\n",myrank);
  int i,j;
  double ***tempp;
  double **delp;
  double **d;
  allocator3(&tempp, xelem, yelem, zelem);
  allocator(&delp, xelem, yelem);
  allocator(&d, xelem, yelem);

    double ires;
	double resnorm=0.0;

    for(j=0; j<yelem; j++)
    {
        for(i=0; i<xelem; i++)
        {
            tempp[i][j][0] = p[i][j][0];
        }
    }

    int iter;
    for(iter=0; iter< 10000; iter++)
    {
        for(i=2; i<xelem-2; i++)
        {
            for(j=2; j<yelem-2; j++)
            {
                double res = b[i][j] - (a[i][j][0]*tempp[i-1][j][0] + a[i][j][1]*tempp[i][j-1][0] + a[i][j][3]*tempp[i][j+1][0] + a[i][j][4]*tempp[i+1][j][0]);
                delp[i][j] = res/a[i][j][2];
            }
        }

        for(i=2; i<xelem-2; i++)
        {
            for(j=2; j<yelem-2; j++)
            {
                tempp[i][j][0] = delp[i][j];
            }
        }
	commu2(tempp);
        pressureBC(tempp);

	double res1=0.0;
        
	double buf = 0.0;
        for(i=2; i<xelem-2; i++)
        {
            for(j=2; j<yelem-2; j++)
            {
                res1 = b[i][j] - (a[i][j][0]*tempp[i-1][j][0] + a[i][j][1]*tempp[i][j-1][0] + a[i][j][3]*tempp[i][j+1][0] + a[i][j][4]*tempp[i+1][j][0] + a[i][j][2]*tempp[i][j][0]);

                buf = buf + pow(res1,2.0)*vol[i][j];
            }
        }

	MPI_Allreduce(&buf,&resnorm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	
	
        if(iter == 0)
        {
            ires = resnorm;
            //if(myrank==master)printf("Initial residual %.6f\n",ires);
            
        }
        else
        {
	  //if(myrank == master)printf("Pressure Step: %d residual: %.6f\n",iter,resnorm/ires);
            if(resnorm / ires < ptol)
            {
	      if(myrank == master)printf("Pressure converged in %d \n",iter);
                break;
            }
        }
    }

    commu2(tempp);
    pressureBC(tempp);

    if(myrank==master)printf("Pressure iterations %d Residual %.2e\n",iter,resnorm/ires);
    for(i=0; i<xelem; i++)
    {
        for(j=0; j<yelem; j++)
        {
            p[i][j][0] = tempp[i][j][0];
        }
    }
    deallocator3(&tempp, xelem, yelem, zelem);
    deallocator(&delp, xelem, yelem);
    deallocator(&d, xelem, yelem);
}

#endif /* GAUSS_SIEDEL_H */


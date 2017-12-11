/*
Notes:
-  3d groundwater flow problem using 2d decomposition for reference - /home/nsaini3/PGREM3D/flow
-  Conjugate gradient solver - /home/nsaini3/blas_tar/solverc/source
-  BiCGstab or GMRES solver - /home/nsaini3/PGREM3D/trans/
 */
/*
 * File:   2phase.cpp
 * Author: nsaini3
 *
 * Created on November 1, 2017
 */

#define PI 3.1415926535897
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <fenv.h>
#include <stdio.h>
#include <limits.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <errno.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
#include <string.h>
#include <stdbool.h>
#include <dirent.h>


#include "common.h"


//Include user generated header files. Must come after global variable decalaration
#include "control.h"
#include "partition.h"
#include "grid.h"
#include "commu.h"
#include "output.h"
/*  #include "read_write.h"*/
#include "bound_cond.h"
/*#include "pressure_solver.h"
  #include "variable_pressure.h"*/
#include "heavy_delta.h"
#include "initial_conditions.h"
/*#include "surface_tension.h"
#include "body_force.h"
#include "rhs.h"*/
#include "functions.h"
#include "rhs_bub.h"
#include "bub_advect.h"
#include "re_distance.h"
#include "hyperbolic.h"
#include "calc_vf.h"



int main(int argc, char **argv)
{
    ///Catches mathematical exceptions
    //feenableexcept(FE_INVALID | FE_OVERFLOW |FE_DIVBYZERO);

  //MPI relevant variables/

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
  master = 0;

  double time1 = MPI_Wtime();
  int i,j,k;

  //Create necessary ouput directories//
  //Let only the master create directory
  if(myrank == master)
    {
      char* path;
      path = concat(getexepath(), "/output");
      mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      memset(path,0,strlen(path));
      path = concat(getexepath(),"/laststep");
      mkdir(path, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      free(path);
    }
/////////////////////////////////////
  
    
    //Read control file/
  //Lets not worry about openmp at the moment
  //    omp_set_num_threads(1);

  //All processors may read the control file
    control();

    
    //Routine to partition the mesh
    partition();
    
    
    xelem=xelem+4; //Include 2 ghost cells on each side
    yelem=yelem+4; //Include 2 ghost cells on each side
    xnode=xelem+1; //
    ynode=yelem+1; //
    
    ///Resize the vectors and initialize data structures//
    allocator(&x,xnode,ynode);
    allocator(&y,xnode,ynode);

    allocator(&xc,xelem,yelem);
    allocator(&yc,xelem,yelem);
    allocator(&vol,xelem,yelem);
    iallocator(&iBC,xelem,yelem);
    
    allocator4(&area,xelem,yelem,2,2);
    
    //Read grid and populate element and node vectors/
    gridread();
    genibc();
    sendptr = (double **) malloc(4 * sizeof(double *));
    recvptr = (double **) malloc(4 * sizeof(double *));
    setupcommu();
    /*for(i=0; i<yelem; i++)
    {
      if(myrank == master)printf("%d %.4f\n",i,bhai.sendrbuf[i]);
      }*/
    

    //Initialize solution vectors/
    struct elemsclr sclr;
    
    allocator3(&sclr.p,xelem,yelem,zelem);
    allocator3(&sclr.u,xelem,yelem,zelem);
    allocator3(&sclr.v,xelem,yelem,zelem);
    allocator3(&sclr.phi,xelem,yelem,zelem);
    allocator3(&sclr.rho,xelem,yelem,zelem);
    allocator3(&sclr.mu,xelem,yelem,zelem);
    
    initialize(sclr);
    
    
    
    //Read from file is startstep != 0/
    //prevfileread(sclr);
    
    if(case_tog == 1)
    {
      //Velocity field for vortex/
        for(i=0; i<xelem; i++)
        {
            for(j=0; j<yelem; j++)
            {
                sclr.u[i][j][0] = -pow(sin(PI*x[i][j]),2.0) * sin(2.0*PI*y[i][j]);
                sclr.v[i][j][0] = pow(sin(PI*y[i][j]),2.0) * sin(2.0*PI*x[i][j]);
            }
        }
    }
    else if(case_tog == 4)
    {
      //Velocity field for zalesak/
        for(i=0; i<xelem; i++)
        {
            for(j=0; j<yelem; j++)
            {
                sclr.u[i][j][0] = PI*(50.0-y[i][j])/314.0;
                sclr.v[i][j][0] = PI*(x[i][j]-50.0)/314.0;
            }
        }
    }
    
    //fast_march(sclr);
    debug = 0;
    double *ires = (double *) malloc(3 * sizeof(double));
    bool exitflag = false;
    int iter=0;
    int print_count=0;
    double deltat=advect_deltat;
    double init_vf=0.0;

    FILE *out;
    if(myrank == master)
      {
	out = fopen("sim_out.txt","w");
	if(out == NULL)
	  {
	    printf("Error opening sim_out.txt!\n");
	    exit(0);
	  }
      }

    
    for(iter=startstep; iter<itermax; iter++)
      {

	/*if(flow_solve == 1)
        {
	  double ***utemp;
	  allocator3(&utemp,xelem,yelem,zelem);
	  double ***vtemp;
	  allocator3(&vtemp,xelem,yelem,zelem);
	  
	  for(i=2;i<xelem-2;i++)
            {
	      for(j=2;j<yelem-2;j++)
                {
		  utemp[i][j][0]=sclr.u[i][j][0];
		  vtemp[i][j][0]=sclr.v[i][j][0];
                }
            }
	  double **rhsx, **rhsy;
	  allocator(&rhsx,xelem,yelem);
	  allocator(&rhsy,xelem,yelem);
	  rhscalc(sclr, rhsx, rhsy, iter, exitflag);

	  double ***ustar;
	  allocator3(&ustar,xelem,yelem,zelem);
	  double ***vstar;
	  allocator3(&vstar,xelem,yelem,zelem);
            //Predictor Step
	  for(i=2; i<xelem-2; i++)
            {
	      for(j=2; j<yelem-2; j++)
                {
		  ustar[i][j][0] = sclr.u[i][j][0]+deltat*rhsx[i][j];
		  vstar[i][j][0] = sclr.v[i][j][0]+deltat*rhsy[i][j];
                }
            }

	  commu(ustar);
	  commu(vstar);
	  vel_BC(ustar, vstar);
	  
            ///Calculate contribution from source term - Surface tension force
            //Note that surface tension force is calculated at the centre of cell at i,j (where p and phi are stored)/
	  double ***st_forcex;
	  allocator3(&st_forcex,xelem,yelem,zelem);
	  double ***st_forcey;
	  allocator3(&st_forcey,xelem,yelem,zelem);


	  surface(sclr,st_forcex, st_forcey);
	  body(sclr,st_forcex,st_forcey);
	  
	  if(p_solver == 1)
            {
	      pressure(ustar,vstar, sclr.p, deltat);
	      
            }
	  else if(p_solver == 2)
            {
	      variable_pressure(ustar, vstar, sclr.p, deltat, sclr.rho, st_forcex, st_forcey);
            }

            //Projection Step
	  for(i=2; i<xelem-2; i++)
            {
	      for(j=2; j<yelem-2; j++)
                {
		  
                    sclr.u[i][j][0] = ustar[i][j][0] - deltat*((2.0/(sclr.rho[i][j][0]+sclr.rho[i+1][j][0]))*((sclr.p[i+1][j][0]-sclr.p[i][j][0])/area[i][j][1][1] + 0.5*(st_forcex[i+1][j][0]+st_forcex[i][j][0])));
                    sclr.v[i][j][0] = vstar[i][j][0] - deltat*((2.0/(sclr.rho[i][j][0]+sclr.rho[i][j+1][0]))*((sclr.p[i][j+1][0]-sclr.p[i][j][0])/area[i][j][0][0] + 0.5*(st_forcey[i][j+1][0]+st_forcey[i][j][0])));

                }
            }

	  commu(sclr.u);
	  commu(sclr.v);
	  vel_BC(sclr.u, sclr.v);

	  if(myrank == master)
	    {
	      printf("Step: %d\n",iter+1);
	      fprintf(out,"Step: %d\n",iter+1);
	    }
	  if(exitflag == false && sol_type == 0)
            {
	      monitor_res(ires, &exitflag, iter, sclr,utemp,vtemp);
            }
	  if(exitflag == true && sol_type == 0)
            {
	      if(myrank==master) printf("Flow solution converged\n");
	      break;
            }
	  
	  deallocator3(&utemp,xelem,yelem,zelem);
	  deallocator3(&vtemp,xelem,yelem,zelem);
	  deallocator3(&ustar,xelem,yelem,zelem);
	  deallocator3(&vstar,xelem,yelem,zelem);
	  deallocator3(&st_forcex,xelem,yelem,zelem);
	  deallocator3(&st_forcey,xelem,yelem,zelem);
	  deallocator(&rhsx,xelem,yelem);
	  deallocator(&rhsy,xelem,yelem);
	  
	  }*/
	

	if(flow_solve == 0 && myrank == master)
	  {
	    printf("Step: %d\n",iter+1);
	    fprintf(out,"Step: %d\n",iter+1);
	  }
        //Bubble Advection and re-distance solvers/
        if(advect_solve == 1)
	  {
	    bub_advect(sclr, iter, deltat);
            //re_distance(sclr);
	    //printf("Value of redist_method on %d id %d\n",myrank,redist_method);
            if(redist_method == 1)
	      {
		//printf("%d calls hyperbolic\n",myrank);
		hyperbolic(sclr);
	      }
            else if(redist_method == 2)
	      {
		printf("Fast March not present. Use hyperbolic\n");
		exit(1);
		//fast_march(sclr);
	      }

            else if(redist_method == 3)
	      {
		printf("direct re-distance not present. Use hyperbolic\n");
		exit(1);
	      //direct_redist(sclr);
	      }

	  }

        print_count++;
        if(print_count == 1)
        {
            output_xml(sclr,iter);//,st_forcex, st_forcey);
        }
        if(print_count == print_gap)
        {
            print_count = 0;
        }


        //Determine time step based on CFL/

        if(time_control == 1)
        {
            double cfl;
            timestep_calc(sclr, &deltat, &cfl);
	    if(myrank == master)
	      {
		printf("CFL number: %.6f time step: %.6f\n",cfl,deltat);
		fprintf(out,"CFL number: %.6f time step: %.6f\n",cfl,deltat);
	      }
        }


        //Determine Void Fraction/
        double vf = 0.0;
        double err = 0.0;
        calc_vf(sclr.phi, &init_vf, &vf, &err);
	if(myrank == master)
	  {
	    printf("Void Fraction: %.6f%% Error in vf: %.6f%%\n\n",vf ,err);
	    fprintf(out,"Void Fraction: %.6f%% Error in vf: %.6f%%\n\n",vf ,err);

	    printf("\n");
	    fprintf(out,"\n");
	  }
    
      }



    
    
    free(ires);
    //filewrite(sclr);

    deallocator3(&sclr.p,xelem,yelem,zelem);
    deallocator3(&sclr.u,xelem,yelem,zelem);
    deallocator3(&sclr.v,xelem,yelem,zelem);
    deallocator3(&sclr.phi,xelem,yelem,zelem);
    deallocator3(&sclr.rho,xelem,yelem,zelem);
    deallocator3(&sclr.mu,xelem,yelem,zelem);
    
    deallocator(&x,xnode,ynode);
    deallocator(&y,xnode,ynode);
    ideallocator(&iBC,xelem,yelem);
    deallocator(&xc,xelem,yelem);
    deallocator(&yc,xelem,yelem);
    deallocator(&vol,xelem,yelem);

    deallocator4(&area,xelem,yelem,2,2);
    ideallocator(&io_info,nprocs,4);

    double time2 = MPI_Wtime();
    double secs = time2-time1;
    if(myrank == master)
      {
	printf("Total run time: %.6f secs\n",secs);
	fprintf(out,"Total run time: %.6f secs\n",secs);
	fclose(out);
      }
    
    
    destroycommu(); 
    free(sendptr);
    free(recvptr);
    MPI_Finalize();
 }
      

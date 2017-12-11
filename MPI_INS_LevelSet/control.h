/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   control.h
 * Author: nsaini3
 *
 * Created on September 26, 2016, 4:35 PM
 */

#ifndef CONTROL_H
#define CONTROL_H

void control()
{
  

  if(myrank == master)
    {
      FILE *controlfile;
  controlfile = fopen("control.txt","r");
  if(controlfile == NULL)
      {
	printf("Error opening control.txt!\n");
	exit(0);
      }
  
    char* line = NULL;
    ssize_t size;
    size_t len = 0;

    char delim [1] = " ";
    char* word = NULL;
    
    while((size = getline(&line, &len, controlfile)) != -1)
    {
      
      word = strtok(line, delim);
      while(word != NULL)
	{
	  //printf("%s\n",word);
	  if(strcmp(word,"xlen") == 0)
	    {
	      word = strtok(NULL,delim);
	      xlen = atof(word);
	    }
	  if(strcmp(word,"ylen") == 0)
	    {
	      word = strtok(NULL,delim);
	      ylen = atof(word);
	    }
	  if(strcmp(word,"zlen") == 0)
	    {
	      word = strtok(NULL,delim);
	      zlen = atof(word);
	    }
	  if(strcmp(word,"xelem") == 0)
	    {
	      word = strtok(NULL,delim);
	      gxelem = atoi(word);
	    }
	  if(strcmp(word,"yelem") == 0)
	    {
	      word = strtok(NULL,delim);
	      gyelem = atoi(word);
	    }
	  if(strcmp(word,"zelem") == 0)
	    {
	      word = strtok(NULL,delim);
	      gzelem = atoi(word);
	    }
	  if(strcmp(word,"Max_Iterations") == 0)
	    {
	      word = strtok(NULL,delim);
	      itermax = atoi(word);
	    }
	  if(strcmp(word,"Kinematic_viscosity") == 0)
	    {
	      word = strtok(NULL,delim);
	      nu = atof(word);
	    }
	  if(strcmp(word,"Tolerance_x") == 0)
	    {
	      word = strtok(NULL,delim);
	      tol = atof(word);
	    }
	  if(strcmp(word,"Bubble_radius") == 0)
	    {
	      word = strtok(NULL,delim);
	      rb_in = atof(word);
	    }
	  if(strcmp(word,"x_pos_of_bubble") == 0)
	    {
	      word = strtok(NULL,delim);
	      xb_in = atof(word);
	    }
	  if(strcmp(word,"y_pos_of_bubble") == 0)
	    {
	      word = strtok(NULL,delim);
	      yb_in = atof(word);
	    }
	  if(strcmp(word,"advect_steps") == 0)
	    {
	      word = strtok(NULL,delim);
	      advect_steps = atoi(word);
	    }
	  if(strcmp(word,"advect_deltat") == 0)
	    {
	      word = strtok(NULL,delim);
	      advect_deltat = atof(word);
	    }
	  if(strcmp(word,"Solution_read") == 0)
	    {
	      word = strtok(NULL,delim);
	      solnread = atoi(word);
	    }
	  if(strcmp(word,"Liquid_density") == 0)
	    {
	      word = strtok(NULL,delim);
	      rhof = atof(word);
	    }
	  if(strcmp(word,"Liquid_viscosity") == 0)
	    {
	      word = strtok(NULL,delim);
	      muf = atof(word);
	    }
	  if(strcmp(word,"Gas_density") == 0)
	    {
	      word = strtok(NULL,delim);
	      rhog = atof(word);
	    }
	  if(strcmp(word,"Gas_viscosity") == 0)
	    {
	      word = strtok(NULL,delim);
	      mug = atof(word);
	    }
	  if(strcmp(word,"Epsilon") == 0)
	    {
	      word = strtok(NULL,delim);
	      epsilon = atof(word);
	    }
	  if(strcmp(word,"Surface_tension_coefficient") == 0)
	    {
	      word = strtok(NULL,delim);
	      sf_coeff = atof(word);
	    }
	  if(strcmp(word,"GS_relaxation_factor") == 0)
	    {
	      word = strtok(NULL,delim);
	      relax = atof(word);
	    }
	  if(strcmp(word,"Tolerance_p") == 0)
	    {
	      word = strtok(NULL,delim);
	      ptol = atof(word);
	    }
	  if(strcmp(word,"Re_distance_timestep") == 0)
	    {
	      word = strtok(NULL,delim);
	      re_time = atof(word);
	    }
	  if(strcmp(word,"Re_distance_loops") == 0)
	    {
	      word = strtok(NULL,delim);
	      re_loops = atoi(word);
	    }
	  if(strcmp(word,"print_gap") == 0)
	    {
	      word = strtok(NULL,delim);
	      print_gap = atoi(word);
	    }
	  if(strcmp(word,"Start_step") == 0)
	    {
	      word = strtok(NULL,delim);
	      startstep = atoi(word);
	    }
	  if(strcmp(word,"Surface_tension") == 0)
	    {
	      word = strtok(NULL,delim);
	      sf_toggle = atoi(word);
	    }
	  if(strcmp(word,"Solve_flow") == 0)
	    {
	      word = strtok(NULL,delim);
	      flow_solve = atoi(word);
	    }
	  if(strcmp(word,"Variable_density_pressure_solver") == 0)
	    {
	      word = strtok(NULL,delim);
	      p_solver = atoi(word);
	    }
	  if(strcmp(word,"Advect_bubble") == 0)
	    {
	      word = strtok(NULL,delim);
	      advect_solve = atoi(word);
	    }
	  if(strcmp(word,"Void_fraction_control") == 0)
	    {
	      word = strtok(NULL,delim);
	      vf_control = atoi(word);
	    }
	  if(strcmp(word,"Re-distance_method") == 0)
	    {
	      word = strtok(NULL,delim);
	      redist_method = atoi(word);
	    }
	   if(strcmp(word,"max_CFL") == 0)
	    {
	      word = strtok(NULL,delim);
	      max_cfl = atof(word);
	    }
	   if(strcmp(word,"gx") == 0)
	    {
	      word = strtok(NULL,delim);
	      gx = atof(word);
	    }
	   if(strcmp(word,"gy") == 0)
	    {
	      word = strtok(NULL,delim);
	      gy = atof(word);
	    }

	   if(strcmp(word,"bub_advect_scheme") == 0)
	     {
	       word = strtok(NULL, delim);
	       if(strcmp(word,"book\n") == 0)
		 {
		   bub_conv_scheme=1;
		 }
	       else if(strcmp(word,"HJ-WENO\n") == 0)
		 {
		   bub_conv_scheme=2;
		 }
	     }

	   if(strcmp(word,"x-boundary") == 0)
	     {
	       word = strtok(NULL, delim);
	       if(strcmp(word,"no-slip\n") == 0)
		 {
		   x_bound=1;
		 }
	       else if(strcmp(word,"slip\n") == 0)
		 {
		   x_bound=2;
		 }
	       else if(strcmp(word,"periodic\n") == 0)
		 {
		   x_bound=3;
		 }
	     }

	   if(strcmp(word,"y-boundary") == 0)
	     {
	       word = strtok(NULL, delim);
	       if(strcmp(word,"no-slip\n") == 0)
		 {
		   y_bound=1;
		 }
	       else if(strcmp(word,"slip\n") == 0)
		 {
		   y_bound=2;
		 }
	       else if(strcmp(word,"periodic\n") == 0)
		 {
		   y_bound=3;
		 }
	     }
	   if(strcmp(word,"Solver_type") == 0)
	     {
	       word = strtok(NULL, delim);
	       if(strcmp(word,"steady-state\n") == 0)
		 {
		   sol_type=0;
		 }
	       else if(strcmp(word,"transient\n") == 0)
		 {
		   sol_type=1;
		 }
	     }
	   if(strcmp(word,"Time_control") == 0)
	     {
	       word = strtok(NULL, delim);
	       if(strcmp(word,"CFL-based\n") == 0)
		 {
		   time_control=1;
		 }
	       else if(strcmp(word,"constant_time\n") == 0)
		 {
		   time_control=2;
		 }
	     }
	   if(strcmp(word,"Case") == 0)
	     {
	       word = strtok(NULL, delim);
	       if(strcmp(word,"vortex\n") == 0)
		 {
		   case_tog=1;
		 }
	       else if(strcmp(word,"bubble_rise\n") == 0)
		 {
		   case_tog=2;
		 }
	       else if(strcmp(word,"bubble_break\n") == 0)
		 {
		   case_tog=3;
		 }
	       else if(strcmp(word,"zalesak\n") == 0)
		 {
		   case_tog=4;
		 }
	       else if(strcmp(word,"homework\n") == 0)
		 {
		   case_tog=5;
		 }
	     }
	   word = strtok(NULL,delim);
	}
        
     }
    fclose(controlfile);

    free(word);
    free(line);
    }

  MPI_Bcast(&xlen,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&ylen,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&zlen,1,MPI_DOUBLE,master,MPI_COMM_WORLD);

  MPI_Bcast(&gxelem,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&gyelem,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&gzelem,1,MPI_INT,master,MPI_COMM_WORLD);

  MPI_Bcast(&itermax,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&nu,1,MPI_DOUBLE,master,MPI_COMM_WORLD);

  MPI_Bcast(&rb_in,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&xb_in,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&yb_in,1,MPI_DOUBLE,master,MPI_COMM_WORLD);

  MPI_Bcast(&advect_steps,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&advect_deltat,1,MPI_DOUBLE,master,MPI_COMM_WORLD);

  MPI_Bcast(&solnread,1,MPI_INT,master,MPI_COMM_WORLD);
  
  MPI_Bcast(&rhof,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&rhog,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&muf,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&mug,1,MPI_DOUBLE,master,MPI_COMM_WORLD);

  MPI_Bcast(&epsilon,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&sf_coeff,1,MPI_DOUBLE,master,MPI_COMM_WORLD);

  MPI_Bcast(&relax,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&ptol,1,MPI_DOUBLE,master,MPI_COMM_WORLD);

  MPI_Bcast(&re_time,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&re_loops,1,MPI_INT,master,MPI_COMM_WORLD);
  
  MPI_Bcast(&print_gap,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&startstep,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&sf_toggle,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&flow_solve,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&p_solver,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&advect_solve,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&vf_control,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&redist_method,1,MPI_INT,master,MPI_COMM_WORLD);
  
  MPI_Bcast(&max_cfl,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&gx,1,MPI_DOUBLE,master,MPI_COMM_WORLD);
  MPI_Bcast(&gy,1,MPI_DOUBLE,master,MPI_COMM_WORLD);

  MPI_Bcast(&bub_conv_scheme,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&x_bound,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&y_bound,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&sol_type,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&time_control,1,MPI_INT,master,MPI_COMM_WORLD);
  MPI_Bcast(&case_tog,1,MPI_INT,master,MPI_COMM_WORLD);
  

  
    // printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",xelem, yelem, zelem, advect_steps, solnread, bub_conv_scheme, print_gap, startstep, sf_toggle, flow_solve, p_solver, x_bound, y_bound, advect_solve, sol_type, vf_control, time_control, redist_method, case_tog);

    //printf("%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f",nu, cfl, tol, xlen, ylen, zlen, rb_in, xb_in, yb_in, advect_deltat, rhof, rhog, muf, mug, epsilon, sf_coeff, relax, ptol, re_time, re_loops, gx, gy, max_cfl);
    //printf("Value of redist_method on %d id %d\n",myrank,redist_method);
    
}

#endif /* CONTROL_H */


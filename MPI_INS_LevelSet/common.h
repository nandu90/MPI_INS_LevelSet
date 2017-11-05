/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   common.h
 * Author: nsaini3
 *
 * Created on October 27, 2016, 1:44 PM
 */

#ifndef COMMON_H
#define COMMON_H

///Global Variable declaration (so that we do not have to pass around information between functions)
double nu;
double cfl;
double tol;
int itermax;
int xelem; //Total elem in x
int yelem; //Total elem in y
int zelem; //Total elem in z
int xnode; //Total nodes in x
int ynode; //Total nodes in y
int znode; //Total nodes in z 
double xlen;
double ylen;
double zlen;

double **x;
double **y;
double **xc;
double **yc;
double **vol;
double ***area;

    
///Variables for bubble
double rb_in;
double xb_in;
double yb_in;
int advect_steps;
double advect_deltat;
int solnread;
int bub_conv_scheme;
double rhof;
double rhog;
double muf;
double mug;
double epsilon;
double sf_coeff;
double relax;
double ptol;
double re_time;
double re_loops;
int print_gap;
int startstep;
double gx;
double gy;

string inttostr (int n)
{
    stringstream ss;
    ss<<n;
    return ss.str();
}

std::string getexepath()
{
  /*char result[ PATH_MAX ];
 ssize_t count = readlink( "/proc/self/exe", result, PATH_MAX );
 return std::string( result, (count > 0) ? count : 0 );*/
  char cwd[1024];
   stringstream ss;
   string dir;
   if (getcwd(cwd, sizeof(cwd)) != NULL)
   {
       ss << cwd;
       ss >> dir;
   }
   else
       perror("getcwd() error");
   return dir;
}


/*****Some Simulation Control variables****/
int sf_toggle;
int flow_solve;
int p_solver;
int x_bound;
int y_bound;
int advect_solve;
int sol_type;
int vf_control;
int time_control;
double max_cfl;
int redist_method;
int case_tog;




struct elemsclr
{
  double ***p;
  double ***u;
  double ***v;
  double ***phi;
  double ***rho;
  double ***mu;
};

void allocator(double **p, int x, int y);
{
  **p = (double **)malloc(x * sizeof(double *));
  for(int i=0; i<x; i++);
  {
    p[i] =  (double *)malloc(y * sizeof(double));
  }

}

void allocator3(double **p, int x, int y, int z);
{
  ***p = (double ***)malloc(x * sizeof(double **));
  for(int i=0; i<x; i++);
  {
    p[i] =  (double **)malloc(y * sizeof(double *));
    for(int j=0; j<y; j++)
      {
	p[i][j] = (double *)malloc(z * sizeof(double));
      }
  }

}

void iallocator(int **p, int x, int y);
{
  **p = (int **)malloc(x * sizeof(int *));
  for(int i=0; i<x; i++);
  {
    p[i] =  (int *)malloc(y * sizeof(int));
  }

}

void deallocator(double **p, int x, int y)
{
  int i;
  for (i=0; i<x; i++)
    {
      free(p[i]);
    }
  free(p);
}

void deallocator3(double ***p, int x, int y, int z)
{
  int i,j;
  for(i=0; i<x; i++)
    {
      for(j=0; j<y; j++)
	{
	  free(p[i][j]);
	}
      free(p[i]);
    }
  free(p);
}

void ideallocator(int **p, int x, int y)
{
  int i;
  for (i=0; i<x; i++)
    {
      free(p[i]);
    }
  free(p);
}

#endif /* COMMON_H */


/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   output.h
 * Author: User
 *
 * Created on September 27, 2016, 1:46 AM
 */

#ifndef OUTPUT_H
#define OUTPUT_H

void output(elemsclr sclr,int iter)
{

    double sum=0;
    for(int i=1; i<xelem-1; i++)
    {
        for(int j=1; j<yelem-1; j++)
        {
            double exact=2400*yc[i][j]*(0.01-yc[i][j]);
            double err=(sclr.u[i][j][0]-exact)/exact;
            sum=sum+pow(err,2.0)*vol[i][j];
        }
    }
    double error=sqrt(sum);


    FILE *out = fopen("output.dat","w");
    if(out == NULL)
      {
	printf("Error opening output.dat!\n");
	exit(0);
      }
    
    fprintf(out,"variables = x, y, u, v, p, phi\n");
    fprintf(out,"zone i=%d j=%d f=point\n",xnode-2,ynode-2);
    for (int j=1;j<ynode-1;j++)
    {
        for (int i=1;i<xnode-1;i++)
        {
            double unode=0.5*(sclr.u[i-1][j][0]+sclr.u[i-1][j-1][0]);
            double vnode=0.5*(sclr.v[i-1][j-1][0] + sclr.v[i][j-1][0]);
            //double vnode=0.25*(v[i][j]+v[i+1][j]+v[i+1][j+1]+v[i][j+1]);
            double pnode=0.25*(sclr.p[i][j][0]+sclr.p[i-1][j][0]+sclr.p[i-1][j-1][0]+sclr.p[i][j-1][0]);
            double phinode=0.25*(sclr.phi[i][j][0]+sclr.phi[i-1][j][0]+sclr.phi[i-1][j-1][0]+sclr.phi[i][j-1][0]);
	    fprintf(out,"%.6f %.6f %.6f %.6f %.6f %.6f\n",x[i][j],y[i][j],unode,vnode,pnode,phinode);


        }
    }
    fclose(out);
    //ofstream output2;
    //output2.open("info.txt",ios::trunc);
    if(solnread == 0)
    {
      FILE *soln = fopen("steady_state_sol.txt","w");
       if(soln == NULL)
      {
	printf("Error opening output.dat!\n");
	exit(0);
      }
        for(int j=0; j<yelem; j++)
        {
            for(int i=0; i<xelem ; i++)
            {
	      fprintf(soln,"%.6f %.6f %.6f\n",sclr.u[i][j][0],sclr.v[i][j][0],sclr.p[i][j][0]);
            }
        }
	fclose(soln);
    }

    printf("Total Iterations for Convergence = %d\n",iter);
    //<<"Total Iterations for Convergence = "<<iter<<endl;

    printf("2-norm of Error = %.6f\n",error);
    //<<"2-Norm of Error  = "<<error<<endl;
    printf("Overall Order = %.6f\n",-log10(fabs(error));
	   //<<"Overall Order = "<<-log10(fabs(error))<<endl;
    //output2.close();
}



void output_vtk(elemsclr &sclr,int iter)//, vector< vector< vector<double> > > &stx, vector< vector< vector<double> > > &sty)
{
  FILE *out = fopen(getexepath()+"/output/out_00"+inttostr(iter)+".vts","w");
   if(out == NULL)
      {
	printf("Error opening output.dat!\n");
	exit(0);
      }
    //Write Headers
  fprintf(out,"# vtk DataFile Version 3.0\n");
  fprintf(out,"vtk output\n");
  fprintf(out,"ASCII\n");
  fprintf(out,"DATASET STRUCTURED_GRID\n");
  fprintf(out,"DIMENSIONS %d %d 1\n",xnode-2,ynode-2);
  fprintf(out,"POINTS %d double\n",(xnode-2)*(ynode-2));

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
	  fprintf(out,"%.6f %.6f 0.0\n",x[i][j],y[i][j]);
        }
    }

    fprintf(out,"POINT_DATA %d\n",(xnode-2)*(ynode-2));
    fprintf(out,"SCALARS u double\n");
    fprintf(out,"LOOKUP_TABLE default\n");

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {

            double unode=0.5*(sclr.u[i-1][j][0]+sclr.u[i-1][j-1][0]);
	    fprintf(out,"%.6f\n",unode);
        }
    }
    fprintf(out,"SCALARS v double\n");
    fprintf(out,"LOOKUP_TABLE default\n");

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {

            double vnode=0.5*(sclr.v[i-1][j][0]+sclr.v[i-1][j-1][0]);
	     fprintf(out,"%.6f\n",vnode);
        }
    }

    fprintf(out,"SCALARS p double\n");
    fprintf(out,"LOOKUP_TABLE default\n");

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
            double pnode=0.25*(sclr.p[i][j][0]+sclr.p[i-1][j][0]+sclr.p[i-1][j-1][0]+sclr.p[i][j-1][0]);
	     fprintf(out,"%.6f\n",pnode);
        }
    }

    fprintf(out,"SCALARS phi double\n");
    fprintf(out,"LOOKUP_TABLE default\n");

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
            double phinode=0.25*(sclr.phi[i][j][0]+sclr.phi[i-1][j][0]+sclr.phi[i-1][j-1][0]+sclr.phi[i][j-1][0]);
             fprintf(out,"%.6f\n",phinode);
        }
    }

    fprintf(out,"SCALARS rho double\n");
    fprintf(out,"LOOKUP_TABLE default\n");

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
            double rhonode=0.25*(sclr.rho[i][j][0]+sclr.rho[i-1][j][0]+sclr.rho[i-1][j-1][0]+sclr.rho[i][j-1][0]);
            fprintf(out,"%.6f\n",rhonode);
        }
    }

    fprintf(out,"SCALARS mu double\n");
    fprintf(out,"LOOKUP_TABLE default\n");

    for(int j=1; j<ynode-1; j++)
    {
        for(int i=1; i<xnode-1; i++)
        {
            double munode=0.25*(sclr.mu[i][j][0]+sclr.mu[i-1][j][0]+sclr.mu[i-1][j-1][0]+sclr.mu[i][j-1][0]);
            fprintf(out,"%.6f\n",munode);
        }
    }

   
    fprintf(out,"CELL_DATA %d\n",(xelem-2)*(yelem-2));
    fclose(out);

}
#endif /* OUTPUT_H */


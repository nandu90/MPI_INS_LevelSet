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

void output_xml(struct elemsclr sclr,int iter)
{
  int i,j,k;
  char* dirpath;
  dirpath = concat(getexepath(), "/output/");
  char buf_dir[12];
  snprintf(buf_dir,12,"%d",iter);
  dirpath = concat(dirpath,buf_dir);
  DIR* dir = opendir("dirpath");
  if(dir)
    {
      closedir(dir);
    }
  else if(ENOENT == errno)
    {
      mkdir(dirpath, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }
  free(dirpath);


  FILE *out;
  FILE *out1;
  //if(myrank == master)
  //{
  if(myrank == master)
    {
      char* path1;
      path1 = concat(getexepath(),"/output/out_00");
      char buffer0[12];
      snprintf(buffer0,12,"%d",iter);
      path1 = concat(path1,buffer0);
      path1 = concat(path1,".pvts");
      out1 = fopen(path1,"w");
      free(path1);

      fprintf(out1,"<VTKFile type=\"PStructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
      fprintf(out1,"<PStructuredGrid WholeExtent=\"%d %d %d %d %d %d\" GhostLevel=\"1\">\n",1,io_info[nprocs-1][1]+1,1,io_info[nprocs-1][3]+1,0,0);
      fprintf(out1,"<PPoints>\n");
      fprintf(out1,"<PDataArray NumberOfComponents=\"3\" format=\"ascii\" type =\"Float32\" Name=\"mesh\"/>\n");
      fprintf(out1,"</PPoints>\n");
      fprintf(out1,"<PPointData>\n");
      fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"phi\"/>\n");
      fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"u\"/>\n");
      fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"v\"/>\n");
      fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"mu\"/>\n");
      fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"rho\"/>\n");
      fprintf(out1,"<PDataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"p\"/>\n");
      fprintf(out1,"</PPointData>\n");
      for(i=0;i<nprocs; i++)
	{
	  fprintf(out1,"<Piece Extent=\"%d %d %d %d %d %d\" Source=\"%d/out.%d.vts\"/>\n",io_info[i][0],io_info[i][1]+1,io_info[i][2],io_info[i][3]+1,0,0,iter,i);
	}
      fprintf(out1,"</PStructuredGrid>\n");
      fprintf(out1,"</VTKFile>\n");
      fclose(out1);
    }
  
      char* path;
      char buffer1[12];
      char buffer2[12];
      path = concat(getexepath(), "/output/");
      snprintf(buffer1,12,"%d",iter);
      path = concat(path,buffer1);
      path = concat(path,"/out.");
      snprintf(buffer2,12,"%d",myrank);
      path = concat(path,buffer2);
      path = concat(path,".vts");
  
  
      out = fopen(path,"w");
      free(path);
      if(out == NULL)
	{
	  printf("Error opening output.dat!\n");
	  exit(0);
	}
    //Write Headers
      fprintf(out,"<VTKFile type=\"StructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n");
      fprintf(out,"<StructuredGrid WholeExtent=\"%d %d %d %d %d %d\">\n",1,io_info[nprocs-1][1]+1,1,io_info[nprocs-1][3]+1,0,0);
      fprintf(out,"<Piece Extent=\"%d %d %d %d %d %d\">\n",io_info[myrank][0],io_info[myrank][1]+1,io_info[myrank][2],io_info[myrank][3]+1,0,0);
      fprintf(out,"<PointData></PointData>\n");
      fprintf(out,"<CellData></CellData>\n");
      fprintf(out,"<Points>\n");
      fprintf(out,"<DataArray NumberOfComponents=\"3\" format=\"ascii\" type =\"Float32\" Name=\"mesh\">\n");
      for(j=2; j<ynode-2; j++)
	{
	  for(i=2; i<xnode-2; i++)
	    {
	      fprintf(out,"%.6f %.6f 0.0\n",x[i][j],y[i][j]);
	    }
	}
      fprintf(out,"</DataArray>");
      fprintf(out,"</Points>\n");
      fprintf(out,"<PointData>\n");
      fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"phi\">\n");
      for(j=2; j<ynode-2; j++)
	{
	  for(i=2; i<xnode-2; i++)
	    {
	      double phinode=0.25*(sclr.phi[i][j][0]+sclr.phi[i-1][j][0]+sclr.phi[i-1][j-1][0]+sclr.phi[i][j-1][0]);
	      fprintf(out,"%.6f\n",phinode);
	    }
	}

      fprintf(out,"</DataArray>\n");
      fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"u\">\n");
      for(j=2; j<ynode-2; j++)
	{
	  for(i=2; i<xnode-2; i++)
	    {
	      
	      double unode=0.5*(sclr.u[i-1][j][0]+sclr.u[i-1][j-1][0]);
	      fprintf(out,"%.6f\n",unode);
	    }
	}

      fprintf(out,"</DataArray>\n");
      fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"v\">\n");
      for(j=2; j<ynode-2; j++)
	{
	  for(i=2; i<xnode-2; i++)
	    {
	      
	      double vnode=0.5*(sclr.v[i-1][j][0]+sclr.v[i-1][j-1][0]);
	      fprintf(out,"%.6f\n",vnode);
	    }
	}
      fprintf(out,"</DataArray>\n");

      fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"mu\">\n");
      for(j=2; j<ynode-2; j++)
	{
	  for(i=2; i<xnode-2; i++)
	    {
	       double munode=0.25*(sclr.mu[i][j][0]+sclr.mu[i-1][j][0]+sclr.mu[i-1][j-1][0]+sclr.mu[i][j-1][0]);
	       fprintf(out,"%.6f\n",munode);
	      
	    }
	}
      fprintf(out,"</DataArray>\n");

      fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"rho\">\n");
      for(j=2; j<ynode-2; j++)
	{
	  for(i=2; i<xnode-2; i++)
	    {
	      
	      double rhonode=0.25*(sclr.rho[i][j][0]+sclr.rho[i-1][j][0]+sclr.rho[i-1][j-1][0]+sclr.rho[i][j-1][0]);
	      fprintf(out,"%.6f\n",rhonode);
	    }
	}
      fprintf(out,"</DataArray>\n");

      fprintf(out,"<DataArray NumberOfComponents=\"1\" format=\"ascii\" type =\"Float32\" Name=\"p\">\n");
      for(j=2; j<ynode-2; j++)
	{
	  for(i=2; i<xnode-2; i++)
	    {
	      
	      double pnode=0.25*(sclr.p[i][j][0]+sclr.p[i-1][j][0]+sclr.p[i-1][j-1][0]+sclr.p[i][j-1][0]);
	      fprintf(out,"%.6f\n",pnode);
	    }
	}
      fprintf(out,"</DataArray>\n");
      fprintf(out,"</PointData>\n");
      fprintf(out,"</Piece>\n");
      fprintf(out,"</StructuredGrid>\n");
      fprintf(out,"</VTKFile>\n");
      fclose(out);
      //}
}
#endif /* OUTPUT_H */


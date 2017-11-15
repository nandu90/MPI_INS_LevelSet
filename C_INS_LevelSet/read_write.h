/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/*
 * File:   read_write.h
 * Author: nsaini3
 *
 * Created on November 11, 2016, 5:13 PM
 */

#ifndef READ_WRITE_H
#define READ_WRITE_H

void prevfileread(struct elemsclr sclr)
{
    if(startstep != 0)
    {
      
      char* path;
      
      path = concat(getexepath(),"/laststep/u.dat");
      FILE *startu = fopen(path,"r");
      if(startu == NULL)
	{
	  printf("Error openeing file laststep/u.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      path = concat(getexepath(),"/laststep/v.dat");
      FILE *startv = fopen(path,"r");
      if(startv == NULL)
	{
	  printf("Error openeing file laststep/v.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      path = concat(getexepath(),"/laststep/p.dat");
      FILE *startp = fopen(path,"r");
      if(startp == NULL)
	{
	  printf("Error openeing file laststep/p.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      path = concat(getexepath(),"/laststep/phi.dat");
      FILE *startphi = fopen(path,"r");
      if(startphi == NULL)
	{
	  printf("Error openeing file laststep/phi.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      path = concat(getexepath(),"/laststep/rho.dat");
      FILE *startrho = fopen(path,"r");
      if(startu == NULL)
	{
	  printf("Error openeing file laststep/rho.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      path = concat(getexepath(),"/laststep/mu.dat");
      FILE *startmu = fopen(path,"r");
      if(startmu == NULL)
	{
	  printf("Error openeing file laststep/mu.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      char* line = NULL;
      ssize_t size;
      size_t len = 0;

      
      int i=0;
      int j=0;
      while((size = getline(&line, &len, startu)) != -1)
	{
	  sclr.u[i][j][0] = atof(line);
	  j++;
	  if(j == yelem-1)
	    {
	      i++;
	      j=0;
	    }
	  if(i == xelem-1)break;
	}
      
      i=0;
      j=0;
      while((size = getline(&line, &len, startv)) != -1)
	{
	  sclr.v[i][j][0] = atof(line);
	  j++;
	  if(j == yelem-1)
	    {
	      i++;
	      j=0;
	    }
	  if(i == xelem-1)break;
	}

      i=0;
      j=0;
      while((size = getline(&line, &len, startp)) != -1)
	{
	  sclr.p[i][j][0] = atof(line);
	  j++;
	  if(j == yelem-1)
	    {
	      i++;
	      j=0;
	    }
	  if(i == xelem-1)break;
	}

      i=0;
      j=0;
      while((size = getline(&line, &len, startphi)) != -1)
	{
	  sclr.phi[i][j][0] = atof(line);
	  j++;
	  if(j == yelem-1)
	    {
	      i++;
	      j=0;
	    }
	  if(i == xelem-1)break;
	}

      i=0;
      j=0;
      while((size = getline(&line, &len, startrho)) != -1)
	{
	  sclr.rho[i][j][0] = atof(line);
	  j++;
	  if(j == yelem-1)
	    {
	      i++;
	      j=0;
	    }
	  if(i == xelem-1)break;
	}

      i=0;
      j=0;
      while((size = getline(&line, &len, startmu)) != -1)
	{
	  sclr.mu[i][j][0] = atof(line);
	  j++;
	  if(j == yelem-1)
	    {
	      i++;
	      j=0;
	    }
	  if(i == xelem-1)break;
	}
      

        
      free(line);
	free(path);
	fclose(startu);
	fclose(startv);
	fclose(startrho);
	fclose(startp);
	fclose(startmu);
	fclose(startphi);
    }
}

void filewrite(struct elemsclr sclr)
{
  char* path;
      
      path = concat(getexepath(),"/laststep/u.dat");
      FILE *startu = fopen(path,"w");
      if(startu == NULL)
	{
	  printf("Error openeing file laststep/u.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      path = concat(getexepath(),"/laststep/v.dat");
      FILE *startv = fopen(path,"w");
      if(startv == NULL)
	{
	  printf("Error openeing file laststep/v.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      path = concat(getexepath(),"/laststep/p.dat");
      FILE *startp = fopen(path,"w");
      if(startp == NULL)
	{
	  printf("Error openeing file laststep/p.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      path = concat(getexepath(),"/laststep/phi.dat");
      FILE *startphi = fopen(path,"w");
      if(startphi == NULL)
	{
	  printf("Error openeing file laststep/phi.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      path = concat(getexepath(),"/laststep/rho.dat");
      FILE *startrho = fopen(path,"w");
      if(startrho == NULL)
	{
	  printf("Error openeing file laststep/rho.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      path = concat(getexepath(),"/laststep/mu.dat");
      FILE *startmu = fopen(path,"w");
      if(startmu == NULL)
	{
	  printf("Error openeing file laststep/mu.dat/n");
	  exit(1);
	}
      memset(path,0,strlen(path));

      for(int j=0; j<yelem; j++)
	{
	  for(int i=0; i<xelem; i++)
	    {
	      fprintf(startu,"%.12f\n",sclr.u[i][j][0]);
	      fprintf(startv,"%.12f\n",sclr.v[i][j][0]);
	      fprintf(startrho,"%.12f\n",sclr.rho[i][j][0]);
	      fprintf(startphi,"%.12f\n",sclr.phi[i][j][0]);
	      fprintf(startmu,"%.12f\n",sclr.mu[i][j][0]);
	      fprintf(startp,"%.12f\n",sclr.p[i][j][0]);
	    }
	}
      
      free(path);
      fclose(startu);
      fclose(startv);
      fclose(startrho);
      fclose(startp);
      fclose(startmu);
      fclose(startphi);
}

#endif /* READ_WRITE_H */


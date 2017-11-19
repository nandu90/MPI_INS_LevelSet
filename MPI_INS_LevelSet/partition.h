/*Routine to partition mesh to various processors*/



#ifndef PARTITION_H
#define PARTITION_H
void partition()
{
  /*Overlay processors on the entire domain*/
  /*Determine m x n distribution of processors*/
  // Make sure number of processors are even//
  if(nprocs % 2 != 0 &&  nprocs != 1)
    {
      if(myrank == master)
	{
	  printf("Please use even number of processors.\nExiting... %d",nprocs);
	  exit(0);
	}
    }

  int eperproc;
  int m,n;
  int totale = gxelem*gyelem;
  if(nprocs > 1)
    {
       //////////Processor division
      int mp,np,temp;
      np = (int)floor(sqrt(nprocs*gyelem/gxelem));
      mp = (int)floor(nprocs/np); 
      if(myrank == master)printf("mp np = %d %d\n",mp,np);
      if(gxelem > gyelem)
	{
	  temp = mp;
	  mp = (int)max(temp,np);
	  np = (int)min(temp,np);
	}
      else
	{
	  temp = np;
	  np = (int)max(mp,np);
	  mp = (int)min(mp,temp);
	}
      if(myrank == master)printf("mp np = %d %d\n",mp,np);
      if(mp*np-nprocs > 0)
	{
	  if(gxelem > gyelem)
	    {
	      np--;
	    }
	  else
	    {
	      mp--;
	    }
	}
      if(myrank == master)printf("mp np = %d %d\n",mp,np);
      while(mp*np-nprocs != 0)
	{
	  if(gxelem > gyelem)
	    {
	      mp++;
	      np = (int)nprocs/mp;
	    }
	  else
	    {
	      np++;
	      mp = (int)nprocs/np;
	    }
	  if(mp == 1 || np == 1 || mp == nprocs || np == nprocs) break;
	}

      if(myrank == master)printf("mp np = %d %d\n",mp,np);
      if(myrank == master)
	{
	  printf("Processor matrix = %d X %d = %d\n", mp ,np,mp*np);
	  if(mp*np-nprocs !=0)
	    {
	      printf("Error dividing processors\n.Exiting...");
	      exit(0);
	    }
	}
      ///////////End of proc division////////


      /////Element division
      eperproc = (int)floor((totale)/nprocs);
      n = (int)floor(gyelem/np);
      m = (int)floor(gxelem/mp);
    
      elemm = m;
      elemn = n;

      
      if(gxelem > mp*m)
	{
	  int extra = gxelem - mp*m;
	  if((myrank+1) % mp == 0)
	    {
	      m = m+extra;
	    }
	}

      if(gyelem > np*n)
	{
	  int extra = gyelem - np*n;
	  double factor = (double)((myrank+1.0)/(mp*(np-1)));
	  if((factor-1.0) > 1.0e-12)
	    {
	      n = n+extra;
	    }
	}

      //printf(" proc m n = %d %d %d\n",(myrank+1),m,n);

      xelem = m;
      yelem = n;

      procm = mp;
      procn = np;

      int sum = xelem*yelem;
      int gsum=0;
      MPI_Reduce(&sum, &gsum,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD);
      if(myrank == master)
	{
	  if(gsum != gxelem*gyelem)
	    {
	      printf("Mismatch between distributed elements and total elements\n.Exiting...");
	      exit(0);
	    }
      }

    }
  else
    {
      xelem = gxelem;
      yelem = gyelem;
    }

}
#endif 

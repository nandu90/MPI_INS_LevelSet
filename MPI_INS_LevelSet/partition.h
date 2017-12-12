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
	      printf("Mismatch between distributed elements and total elements\n.Exiting... gsum = %d, gxelem = %d, gyelem = %d\n", gsum, gxelem, gyelem);
	      exit(0);
	    }
      }

      int i,j,k;
      int send[2];
      int recv[2];
      
      //Array to inform io routine output_xml that how many elements are on each processor
      io_info = (int **)malloc(nprocs*sizeof(int *));
      for(i=0; i<nprocs; i++)
	{
	  io_info[i] = (int *)malloc(4*sizeof(int));
	}
      
      if(myrank != master)
	{
	  send[0] = xelem;
	  send[1] = yelem;
	  MPI_Send(send, 2, MPI_INT,master,myrank, MPI_COMM_WORLD);
	}
      if(myrank == master)
	{	  
	  io_info[0][0] = 1;
	  io_info[0][1] = xelem;
	  io_info[0][2] = 1;
	  io_info[0][3] = yelem;
	  j=1;
	  k=0;
	  for(i=1; i<nprocs; i++)
	    {
	      MPI_Recv(recv,2, MPI_INT,i,i,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	      if(j<mp)
		{
		  io_info[i][0] = io_info[i-1][1]+1;
		  io_info[i][1] = io_info[i-1][1]+recv[0];
		  io_info[i][2] = io_info[i-1][2];
		  io_info[i][3] = io_info[i-1][3];
		}
	      else
		{
		  k++;
		  io_info[i][0] = io_info[i-mp][0];
		  io_info[i][1] = io_info[i-mp][1];
		  io_info[i][2] = io_info[i-mp][3]+1;
		  io_info[i][3] = io_info[i-mp][3] + recv[1];
		  j=0;
		}
	      j++;
	    }
	}
    }
  else
    {
      xelem = gxelem;
      yelem = gyelem;
    }

  int i;
  for(i=0; i<nprocs; i++)
  {
    MPI_Bcast(io_info[i],4,MPI_INT,master,MPI_COMM_WORLD);
  }
  /*printf("outside %d %d %d\n",myrank, xelem, yelem);
  if(myrank == master)
    {
      int i;
      for(i=0; i<nprocs; i++)
	{
	  printf("%d %d %d %d %d\n", i , io_info[i][0], io_info[i][1], io_info[i][2],io_info[i][3]);
	}
	}*/
}
#endif 

/*Routine to partition mesh to various processors*/



#ifndef PARTITION_H
#define PARTITION_H
void partition()
{
  /*Overlay processors on the entire domain*/
  /*Determine m x n distribution of processors*/
  // Make sure number of processors are even//
  if(nprocs % 2 != 0 || nprocs != 1)
    {
      if(myrank == master)
	{
	  printf("Please use even number of processors./nExiting...");
	  exit(1);
	}
    }

  int diff = 10000;
  int temp;
  int factor1, factor2;
  if(nprocs == 1)
    {
      xelem = gxelem;
      yelem = gyelem;
    }
  else
    {
      for(int i=2; i <= (int)nprocs/2; i++) 
	{
	  if(nprocs % i == 0)
	    {
	      factor1 = i;
	      factor2 = (int)nprocs/factor1;
	      temp = abs(factor1 - factor2);
	      if(diff >= temp)
		{
		  diff =temp;
		}
	      else
		{
		  break;
		}
	    }
	}
    }

  int m,n;
  if(gxelem > gyelem)
    {
      m = (int)max(factor1,factor2);
      n = (int)min(factor1,factor2);
    }
  else
    {
      n = (int)max(factor1,factor2);
      m = (int)min(factor1,factor2);
    }
}
#endif 

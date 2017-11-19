

#ifndef COMMU_H
#define COMMU_H

void genibc()
{
  /*This routine works only on ghost elements on each processor
    Convention: 
    if interior element iBC = 0 (IBC has been initalized to 0, so need to do this here)
    if interior boundary element iBC = 1 - Responsible for communication
    if domain boundary element iBC  = 2
  */

  //First consider all ghost elements to be 1
  int i,j;
  for(j=0; j<yelem; j++)
    {
      iBC[0][j] = 1; //left band
      iBC[xelem-1][j] = 1; //right band
    }
  for(i=0; i<xelem; i++)
    {
      iBC[i][0] = 1; //lower band
      iBC[i][yelem-1] = 1; //upper band
    }

  int left = 0;
  int right = 0;
  int up = 0;
  int down = 0;
  

  /*Now take care of ghosts of boundary processors*/
  if((int)floor(myrank/procm) == 0) //Means that this processor is aligned with lower boundary
    {
      down = 1;
      for(i=0; i<xelem; i++)
	{
	  iBC[i][0] = 2;
	}
      
    }

  if((int)floor(myrank/procm) == procn-1) //Means processor is on upper boundary
    {
      up = 1; 
      for(i=0; i<xelem; i++)
	{
	  iBC[i][yelem-1] = 2;
	}
    }

  if((myrank % procm) == 0) //Means processor is on left boundary
    {
      left = 1;
      for(j=0; j<yelem; j++)
	{
	  iBC[0][j] = 2;
	}
    }

  if((myrank % procm) == procm-1) //Means processor is on right boundary
    {
      right = 1;
      for(j=0; j<yelem; j++)
	{
	  iBC[xelem-1][j] = 2;
	}
    }

  //  printf("%d processor is on left = %d; right = %d; up = %d; down = %d\n",myrank+1,left,right,up,down);
}
#endif

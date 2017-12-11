/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   bound_cond.h
 * Author: nsaini3
 *
 * Created on September 27, 2016, 8:02 PM
 */

#ifndef BOUND_COND_H
#define BOUND_COND_H

void imposeBC(struct elemsclr sclr)
{
    //Wall BC
  int i,j;
    for( i=0; i<xelem; i++)
    {
      if(iBC[i][1] == 2)
	{
	  sclr.u[i][1][0] = -sclr.u[i][2][0];
        
	  sclr.v[i][1][0] = 0.0;
	}
    }
    for( i=0; i<xelem; i++)
    {
      if(iBC[i][yelem-2] == 2)
	{
	  sclr.u[i][yelem-2][0]= -sclr.u[i][yelem-3][0];
        
	  sclr.v[i][yelem-2][0] = 0.0;
	  sclr.v[i][yelem-3][0] = 0.0;
	}
    }
    
    //Periodic BC
    //To be handled differently
    /*
    for( j=0; j<yelem; j++)
    {
        sclr.u[0][j][0] = sclr.u[xelem-2][j][0];
        sclr.u[xelem-1][j][0] = sclr.u[1][j][0];
        
        sclr.v[0][j][0] = sclr.v[xelem-2][j][0];
        sclr.v[xelem-1][j][0] = sclr.v[1][j][0];
	}*/

    
}

void periodicBC(double ***scalar)
{
  int j;
    //Periodic BC
    for( j=0; j<yelem; j++)
    {
        scalar[0][j][0] = scalar[xelem-2][j][0];
        scalar[xelem-1][j][0] = scalar[1][j][0];
    }
}

void wallBC(double ***scalar)
{
  int i;
    //Wall BC
    for( i=0; i<xelem; i++)
    {
      if(iBC[i][1] == 2)scalar[i][1][0] = -scalar[i][2][0];
      if(iBC[i][yelem-2] == 2)scalar[i][yelem-2][0]= -scalar[i][yelem-3][0];
        
    }
}

void walluBC(double ***scalar)
{
  int i;
    //Wall BC
    for( i=0; i<xelem; i++)
    {
        if(iBC[i][1] == 2)scalar[i][1][0] = -scalar[i][2][0];
        if(iBC[i][yelem-2] == 2)scalar[i][yelem-2][0]= -scalar[i][yelem-3][0];
        
    }
}


void wallvBC(double ***scalar)
{
  int i;
    //Wall BC
    for( i=0; i<xelem; i++)
    {
        if(iBC[i][1] == 2)scalar[i][1][0] = 0.0;
        if(iBC[i][yelem-2] == 2)scalar[i][yelem-2][0]= 0.0;
        if(iBC[i][yelem-2] == 2)scalar[i][yelem-3][0]= 0.0;
        
    }
}


void zerogradBC(double ***scalar)
{
  int i;
    for( i=0; i<xelem; i++)
    {
        if(iBC[i][1] == 2)scalar[i][1][0] = scalar[i][2][0];
        if(iBC[i][yelem-2] == 2)scalar[i][yelem-2][0]= scalar[i][yelem-3][0];
        
    }
}

void gradBC(double ***scalar)
{
  int i,j;
    for( i=0; i<xelem; i++)
    {
        if(iBC[i][1] == 2)scalar[i][1][0] = 0.0;
        if(iBC[i][yelem-2] == 2)scalar[i][yelem-2][0]= 0.0;
        
    }
}

void bothscalarBC(double ***scalar)
{
  int i,j;
    //Wall BC
    for( i=0; i<xelem; i++)
    {
        if(iBC[i][1] == 2)scalar[i][1][0] = -scalar[i][2][0];
        if(iBC[i][yelem-2] == 2)scalar[i][yelem-2][0]= -scalar[i][yelem-3][0];
        
    }
    
    //Periodic BC
    for( j=0; j<yelem; j++)
    {
        scalar[0][j][0] = scalar[xelem-2][j][0];
        scalar[xelem-1][j][0] = scalar[1][j][0];
        
    }
}

void pressureBC(double ***scalar)
{
  int i,j;
    //For advection case
    /*for( i=0; i<xelem; i++)
    {
        scalar[i][0][0] = scalar[i][1][0];
        scalar[i][yelem-1][0]= scalar[i][yelem-2][0];
        
    }
    
    
    for( j=0; j<yelem; j++)
    {
        scalar[0][j][0] = (-4800*xc[0][j] + 96)*nu;//scalar[1][j][0] + 4800.0*area[1][j][1][1];
        scalar[xelem-1][j][0] = scalar[xelem-2][j][0] - 4800*nu*area[xelem-2][j][1][1];
        
    }*/
    
    //For Bubble breakup and bubble rise case
    for( i=0; i<xelem; i++)
    {
        if(iBC[i][1] == 2)scalar[i][1][0] = scalar[i][2][0];
        if(iBC[i][yelem-2] == 2)scalar[i][yelem-2][0]= 0.0;
        
    }
    
    
    for( j=0; j<yelem; j++)
    {
        if(iBC[1][j] == 2)scalar[1][j][0] = scalar[2][j][0];//scalar[1][j][0] + 4800.0*area[1][j][1][1];
        if(iBC[xelem-2][j] == 2)scalar[xelem-2][j][0] = scalar[xelem-3][j][0];
        
    }
    
    /*double line = 2.5*2.0*rb_in;
    for( i=2; i<xelem-2; i++)
    {
        for( j=2; j<yelem-2; j++)
        {
            if(yc[i][j] > line)
            {
                scalar[i][j][0] = 0.0;
            }
        }
	}*/
}

void vel_BC(double ***u, double ***v)
{
  int i,j;
    /***Taking care of left and right direction wall****/
    if(x_bound == 1)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)u[1][j][0] = 0.0;
            if(iBC[xelem-2][j] == 2)u[xelem-2][j][0] = 0.0;
            if(iBC[xelem-2][j] == 2)u[xelem-3][j][0] = 0.0;
            
            if(iBC[1][j] == 2)v[1][j][0] = -v[2][j][0];
            if(iBC[xelem-2][j] == 2)v[xelem-2][j][0] = -v[xelem-3][j][0];
        }
    }
    else if(x_bound == 2)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)u[1][j][0] = u[2][j][0];
            if(iBC[xelem-2][j] == 2)u[xelem-3][j][0] = u[xelem-4][j][0];
            if(iBC[xelem-2][j] == 2)u[xelem-2][j][0] = u[xelem-3][j][0];
            
            if(iBC[1][j] == 2)v[1][j][0] = v[2][j][0];
            if(iBC[xelem-2][j] == 2)v[xelem-2][j][0] = v[xelem-3][j][0];
        }
    }
    else if(x_bound == 3)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)u[1][j][0] = u[xelem-3][j][0];
            if(iBC[xelem-2][j] == 2)u[xelem-2][j][0] = u[2][j][0];
        
            if(iBC[1][j] == 2)v[1][j][0] = v[xelem-3][j][0];
            if(iBC[xelem-2][j] == 2)v[xelem-1][j][0] = v[1][j][0];
        }
    }
    
    /****Now take care of up and down walls****/
    if(y_bound == 1)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][1] == 2)u[i][1][0] = -u[i][2][0];
            if(iBC[i][yelem-2] == 2)u[i][yelem-2][0]= -u[i][yelem-3][0];

            if(iBC[i][1] == 2)v[i][1][0] = 0.0;
            if(iBC[i][yelem-2] == 2)v[i][yelem-2][0] = 0.0;
            if(iBC[i][yelem-2] == 2)v[i][yelem-3][0] = 0.0;
        }
    }
    else if(y_bound == 2)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][1] == 2)u[i][1][0] = u[i][2][0];
            if(iBC[i][yelem-2] == 2)u[i][yelem-2][0]= u[i][yelem-3][0];

            if(iBC[i][1] == 2)v[i][1][0] = v[i][2][0];
            if(iBC[i][yelem-2] == 2)v[i][yelem-3][0] = v[i][yelem-4][0];
            if(iBC[i][yelem-2] == 2)v[i][yelem-2][0] = v[i][yelem-3][0];
            
        }
    }
    else if(y_bound == 3)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][0] == 2)u[i][0][0] = u[i][yelem-2][0];
            if(iBC[i][yelem-1] == 2)u[i][yelem-1][0] = u[i][1][0];
            
            if(iBC[i][0] == 2)v[i][0][0] = v[i][yelem-2][0];
            if(iBC[i][yelem-1] == 2)v[i][yelem-1][0] = v[i][1][0];
        }
    }
    
}


void level_setBC(double ***scalar)
{
  int i,j;
    /****left and right walls*****/
    if(x_bound == 1 || x_bound == 2)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)scalar[1][j][0] = scalar[2][j][0];
            if(iBC[xelem-2][j] == 2)scalar[xelem-2][j][0] = scalar[xelem-3][j][0];
        }
    }
    else if(x_bound == 3)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)scalar[1][j][0] = scalar[xelem-3][j][0];
            if(iBC[xelem-2][j] == 2)scalar[xelem-2][j][0] = scalar[2][j][0];
        }
    }
    
    /****Top and bottom walls***/
    if(y_bound == 1 || y_bound == 2)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][1] == 2)scalar[i][1][0] = scalar[i][2][0];
            if(iBC[i][yelem-2] == 2)scalar[i][yelem-2][0] = scalar[i][yelem-3][0];
        }
    }
    else if(y_bound == 3)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][0] == 2)scalar[i][0][0] = scalar[i][yelem-2][0];
            if(iBC[i][yelem-1] == 2)scalar[i][yelem-1][0] = scalar[i][1][0];
        }
    }
}

void grad_level_setBC(double ***scalar)
{
  int i,j;
    /****left and right walls*****/
    if(x_bound == 1 || x_bound == 2)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)scalar[1][j][0] = 0.0;
            if(iBC[xelem-2][j] == 2)scalar[xelem-2][j][0] = 0.0;
        }
    }
    else if(x_bound == 3)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[0][j] == 2)scalar[0][j][0] = scalar[xelem-2][j][0];
            if(iBC[xelem-1][j] == 2)scalar[xelem-1][j][0] = scalar[1][j][0];
        }
    }
    
    /****Top and bottom walls***/
    if(y_bound == 1 || y_bound == 2)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][1] == 2)scalar[i][1][0] = 0.0;
            if(iBC[i][yelem-2] == 2)scalar[i][yelem-2][0] = 0.0;
        }
    }
    else if(y_bound == 3)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][0] == 2)scalar[i][0][0] = scalar[i][yelem-2][0];
            if(iBC[i][yelem-1] == 2)scalar[i][yelem-1][0] = scalar[i][1][0];
        }
    }
}

void cell_center_vel_BC(double ***u, double ***v)
{
  int i,j;
    /***Taking care of left and right direction wall****/
    if(x_bound == 1)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)u[1][j][0] = -u[2][j][0];
            if(iBC[xelem-2][j] == 2)u[xelem-2][j][0] = -u[xelem-3][j][0];
            
            if(iBC[1][j] == 2)v[1][j][0] = -v[2][j][0];
            if(iBC[xelem-2][j] == 2)v[xelem-2][j][0] = -v[xelem-3][j][0];
        }
    }
    else if(x_bound == 2)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[1][j] == 2)u[1][j][0] = u[2][j][0];
            if(iBC[xelem-2][j] == 2)u[xelem-2][j][0] = u[xelem-3][j][0];
            
            if(iBC[1][j] == 2)v[1][j][0] = v[2][j][0];
            if(iBC[xelem-2][j] == 2)v[xelem-2][j][0] = v[xelem-3][j][0];
        }
    }
    else if(x_bound == 3)
    {
        for( j=0; j<yelem; j++)
        {
            if(iBC[0][j] == 2)u[0][j][0] = u[xelem-2][j][0];
            if(iBC[xelem-1][j] == 2)u[xelem-1][j][0] = u[1][j][0];
        
            if(iBC[0][j] == 2)v[0][j][0] = v[xelem-2][j][0];
            if(iBC[xelem-1][j] == 2)v[xelem-1][j][0] = v[1][j][0];
        }
    }
    
    /****Now take care of up and down walls****/
    if(y_bound == 1)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][1] == 2)u[i][1][0] = -u[i][2][0];
            if(iBC[i][yelem-2] == 2)u[i][yelem-2][0]= -u[i][yelem-3][0];

            if(iBC[i][1] == 2)v[i][1][0] = -v[i][2][0];
            if(iBC[i][yelem-2] == 2)v[i][yelem-2][0]= -v[i][yelem-3][0];
        }
    }
    else if(y_bound == 2)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][1] == 2)u[i][1][0] = u[i][2][0];
            if(iBC[i][yelem-2] == 2)u[i][yelem-2][0]= u[i][yelem-3][0];

            if(iBC[i][1] == 2)v[i][1][0] = v[i][2][0];
            if(iBC[i][yelem-2] == 2)v[i][yelem-2][0]= v[i][yelem-3][0];
            
        }
    }
    else if(y_bound == 3)
    {
        for( i=0; i<xelem; i++)
        {
            if(iBC[i][0] == 2)u[i][0][0] = u[i][yelem-2][0];
            if(iBC[i][yelem-1] == 2)u[i][yelem-1][0] = u[i][1][0];
            
            if(iBC[i][0] == 2)v[i][0][0] = v[i][yelem-2][0];
            if(iBC[i][yelem-1] == 2)v[i][yelem-1][0] = v[i][1][0];
        }
    }
    
}
#endif /* BOUND_COND_H */


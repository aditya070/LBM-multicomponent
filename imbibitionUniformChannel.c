/*  multi component.
6	2	5
3	0	1
7	4	8
*/
// search edit to remove/add wall node
#include                <time.h>
#include                <stdio.h>
#include                <math.h>
#include                <stdlib.h>
#include                "mkl.h"
#include		<omp.h>
#include                <sys/stat.h>

#define		XLENGTH     1500
#define		YLENGTH     32
#define		WIDTH       30
#define		LENGTH	    500


#define		SUB         2
#define		Q           9
#define		w0          4.0 / 9.0
#define		w1          1.0 / 9.0
#define		w2          1.0 / 36.0
#define		RHO_R  	    1/* 0 is red */
#define		RHO_B  	    1
#define		DROP_RAD    50
#define		G           2/* interaction force duet to resultant force resulted due to absence of similar molecule in adjacent node and being replaced by a new particle at that node*/
#define		Gads0       0
#define		Gads1       -1// positive is repel for G positive is repel for Gads
#define		g           0.000000
//#define     ux0         0.1
#define		T           100001


void        tables                              ( void );
void        propagate                           ( void );
void        initialise                          ( void );
void        calc_obs                            ( void );
void        project                             ( void);
double	    y_0                                 ( void );
void        density_data                        ( int data_counter, int id );
void        readrestart                         (void);


/*double		f	      [XLENGTH][YLENGTH][Q];*/
double		r	      [XLENGTH][YLENGTH][Q];
double		b	      [XLENGTH][YLENGTH][Q];
double		back	      [XLENGTH][YLENGTH][Q];
int		next_x        [XLENGTH]                  [Q];
int		next_y	               [YLENGTH]         [Q];
int             cx                     [Q] = { 0,1,0,-1,0,1,-1,-1,1};     /* Lattice basis                            */
int             cy                     [Q] = { 0,0,1,0,-1,1,1,-1,-1};
double		    t                      [Q] = {w0,w1,w1,w1,w1,w2,w2,w2,w2};     /* Link weights (D2Q9 specific)             */
int		bounce			[Q] = {0,3,4,1,2,7,8,5,6};      //bounceback
double		psi      [XLENGTH][YLENGTH][SUB];
double		ux	     [XLENGTH][YLENGTH][SUB];
double		uy	     [XLENGTH][YLENGTH][SUB];

double		bx	     [XLENGTH][YLENGTH][SUB];/* interaction force between blue and red fluid*/
double		by	     [XLENGTH][YLENGTH][SUB];

double		rx	     [XLENGTH][YLENGTH][SUB];/* interaction force between blue or red fluid and wall*/
double		ry	     [XLENGTH][YLENGTH][SUB];

double		gx	     [XLENGTH][YLENGTH][SUB];
double		gy	     [XLENGTH][YLENGTH][SUB];

double		rho      [XLENGTH][YLENGTH][SUB];
double		rhot     [XLENGTH][YLENGTH];
double          ux_eq0       [XLENGTH][YLENGTH]; /* */
double          uy_eq0       [XLENGTH][YLENGTH];

double          ux_eq1       [XLENGTH][YLENGTH];
double          uy_eq1       [XLENGTH][YLENGTH];

double          uprime_x     [XLENGTH][YLENGTH];
double          uprime_y     [XLENGTH][YLENGTH];

double 		U		[XLENGTH][YLENGTH];
double		V		[XLENGTH][YLENGTH];
double          rhoN         [XLENGTH][YLENGTH];

double		bnode        [XLENGTH][YLENGTH];//wall
double          rnode        [XLENGTH][YLENGTH];//fluid boundary

double		OMEGA_B, OMEGA_R, U_0;
int         k;
int        xstart,xend,xfill,ystart,yend;




int	main	(void)
{
    int                      count, argc;

    tables();
    OMEGA_R = 1;
    OMEGA_B =1 ;
     
    xstart = (int)(XLENGTH/2);
    xend=xstart+LENGTH+1;
    xfill = xstart;
    ystart= (int)((YLENGTH-1)/2-WIDTH/2);
    yend= ystart+WIDTH;
	
    initialise();


    for (k = 0; k < T; k++)
    {   printf("%d\t",k);

         propagate           ();
        calc_obs( );
        project (  );

        if (k==2000 || k ==1000 || k % 5000 ==  0)                                                      /* Sample evolution                 */
        {

             density_data ( count, k );/*Presently u_prime is measured as fluid velocity which has to be changed*/

                //readrestart();
        }

    }

    scanf("%d",&k);
    return  0;
}


void    density_data ( int data_counter, int id )
{   int     x, y;
    FILE    *f1;
    FILE    *f2;
    FILE    *f3;
    char    phase[100] = "phase";
    char    velocity[100] = "velocity";
    char    density[100] = "density";
 /*  Data for Tecplot
    */

    sprintf ( phase, "phase%d.dat",id);
    f1 = fopen(phase, "w");
    fprintf( f1, "Variables = X, Y, rho0, rho1, rhoN, U, V \n");
    fprintf(f1,"Zone I=   %d,J=   %d,F=POINT\n", XLENGTH, YLENGTH);
 	for (y = YLENGTH - 1; y > -1; y--)
	{	for (x = 0; x < XLENGTH; x++)
		{
			fprintf(f1,"%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\n",x,y,rho[x][y][0],rho[x][y][1],rhoN[x][y],U[x][y],V[x][y]);
			fprintf(f1, "\n");
		}

    }
    fflush(f1);
	fclose(f1);


//   sprintf ( velocity, "velocity%d.dat",id);
  /* f2 = fopen(velocity, "a");

    fprintf(f2,"Zone I=   %d,J=   %d,F=POINT\n", XLENGTH, YLENGTH);

 	for (y = YLENGTH - 1; y > -1; y--)
	{	for (x = 0; x < XLENGTH; x++)
        fprintf(f2,"%d\t%d\t%lf\t%lf\n",x,y,uprime_x[x][y],uprime_y[x][y]);
       fprintf(f1, "\n");
   }
    fflush(f2);
	 fclose(f2);
/*
//    sprintf ( density, " density%d.dat",id);
    f3 = fopen(density, "w");
   fprintf(f3,"Zone I=   %d,J=   %d, k=   %d,F=POINT\n", XLENGTH, YLENGTH,ZLENGTH);
   for (z = ZLENGTH - 1; z > -1; z--)
 	for (y = 0 ; y < YLENGTH; y++)
	{	for (x = 0; x < XLENGTH; x++)
            fprintf(f3,"%d\t%d\t%d\t%lf\n",x,y,z,rhot[x][y][z]);
        fprintf(f3, "\n");
    }
    fflush(f3);
	 fclose(f3);*/



    return;

}
/*
void readrestart()
{
    int x,y,i;
    FILE    *f4;
    char    FileNamerestart[10]="restart";
    char    string[10]={'0'};

    f4 = fopen("restart","r");

    for(x=0;x<XLENGTH;x++)
    {
        for(y=0;y<YLENGTH;y++)
        {
            for(i=0;i<Q;i++)
            {
                r[x][y][i]=0.0;
                b[x][y][i]=0.0;
            }
        }
    }

    for(x=0;x<XLENGTH;x++)
    {
        for(y=0;y<YLENGTH;y++)
        {
            for(i=0;i<Q;i++)
            {
                fscanf(f4,"%lf\t%lf\t\n",&r[x][y][i],&b[x][y][i]);
                r[x][y][i]=b[x][y][i];
                r[x][y][i]=b[x][y][i];
            }
        }
    }

    for(x=0;x<XLENGTH;x++)
    {
        for(y=0;y<YLENGTH;y++)
        {
        if(x>=xstart && x<=xend)
        { 
            rnode[x][0] = rnode[x][YLENGTH-1]= 1.0;


           //bnode[x][0] = bnode[x][YLENGTH-1]= 1.0;//edit
        }
            
        }
    }

    fflush(f4);
    fclose(f4);
}*/

void	tables	(void)
{ int x = 0;
	 int y = 0;
     int i = 0;


  	for (x = 0; x < XLENGTH; x++)
    {
    for (y = 0; y < YLENGTH; y++)
     {
				for (i = 0; i < Q; i++)
				{   next_x [x][i] = x + cx[i];
		 		if (next_x[x][i] > XLENGTH - 1)
				next_x [x][i] = 0;
	 			if (next_x[x][i] < 0)
				next_x [x][i] = XLENGTH-1;

                next_y [y][i] = y + cy[i];
		  	  	if (next_y[y][i] > YLENGTH - 1)
				next_y [y][i] = 0;
	 		  	if (next_y[y][i] < 0)
				next_y [y][i] = YLENGTH - 1;

				}
     }
   }
    return;
}



void	initialise (void)
{   int        x, y,i, new_x, new_y;
 int        x0,y0,x1;
   double     Center1,Center2,dist_centre,dist_center,ytop,ybottom,slope;
slope = 0.0;

   for (x = 0; x < XLENGTH; x++)
     {
       for (y = 0; y < YLENGTH; y++)
         {
	ytop =yend+slope*(x-xstart) ;
	ybottom=ystart-slope*(x-xstart);
        rnode[x][y] = 0.0;
        bnode[x][y]=0.0;
        if(x>xstart && x<xend)
       {
	if(y<ybottom || y>ytop)
	 { 
           // rnode[x][y]= 1.0;


          bnode[x][y]= 1.0;//edit
        }
      }
        }
      }
      for (x = 0; x < XLENGTH; x++)
      {
          for (y = 0; y < YLENGTH; y++)
          {
	  for (i = 0; i < Q; i++)
           {
            new_x = next_x[x][i];
            new_y = next_y[y][i];

            if( bnode [new_x][new_y]== 1.0 )//edit
            {

           // if(x>xstart && x<xend) edit
            {rnode[x][y]= 1.0;}
            }
           }//edit*/


       for (i = 0; i < Q; i++)
             {

		r[x][y][i]=0;
		b[x][y][i]=0;
              if(bnode[x][y]!=1)//edit check 
            
            {
                  if(x>xfill )
                   {
                       r[x][y][i]=RHO_R*t[i];
                   }
                   else
                   {
                      b[x][y][i]=RHO_B*t[i];
                   }

            }

            }
          }
      }
    return;
}


void	calc_obs ( void )
{   int       x, y, i;
    double  ux_sum, uy_sum;
/////////////////////////////////////////////////omp
#pragma omp parallel
  {////1

  #pragma omp	for private(y, i,ux_sum,uy_sum) //reduction(+:redmass, totmass)////////2
     for (x = 0; x < XLENGTH; x++)
       {
       for (y = 0; y < YLENGTH; y++)
          {
           rho[x][y][0] = rho[x][y][1] = rhot[x][y] = ux[x][y][0] = ux[x][y][1] = uy[x][y][0] = uy[x][y][1] = 0.0;


 	   for (i = 0; i < Q; i++)
             {
                 if (bnode[x][y] !=1)
                {


                rho  [x][y][0]   += r[x][y][i];
                ux   [x][y][0]   += r[x][y][i] * cx[i];
                uy   [x][y][0]   += r[x][y][i] * cy[i];

                rho  [x][y][1]	+= b[x][y][i];
                ux   [x][y][1]   += b[x][y][i] * cx[i];
                uy   [x][y][1]   += b[x][y][i] * cy[i];
                }
     	     }

      ux_sum = ux [x][y][0]*OMEGA_R + ux [x][y][1]*OMEGA_B;/* sum of rho_sigma*ea/tou */
      uy_sum = uy [x][y][0]*OMEGA_R + uy [x][y][1]*OMEGA_B;
      rhot[x][y] =  rho  [x][y][0] +  rho  [x][y][1];/*sum rho_sigma sigma is component */

      if ( rho  [x][y][0] + rho  [x][y][1]!= 0)
          {
         rhoN[x][y] = (rho[x][y][0] - rho[x][y][1])/(rho[x][y][0] + rho[x][y][1]);/* track the interface not to be used in this application*/
           if ( bnode[x][y] == 1.0)
           rhoN[x][y] = 2;
           

            uprime_x[x][y] = (ux_sum)/ ( rho [x][y][0]*OMEGA_R + rho [x][y][1]*OMEGA_B);
            uprime_y[x][y] = (uy_sum)/ ( rho [x][y][0]*OMEGA_R + rho [x][y][1]*OMEGA_B);
//fluid velocity
		U[x][y] = (ux[x][y][0]+ux[x][y][1]+0.5*((bx[x][y][0]+ rx[x][y][0]+ gx[x][y][0])+(bx[x][y][1]+ rx[x][y][1]+ gx[x][y][1])))/(rho[x][y][0]+rho[x][y][1]);
		V[x][y] =( uy[x][y][0]+uy[x][y][1]+0.5*((by[x][y][0]+ ry[x][y][0]+ gy[x][y][0])+(by[x][y][1]+ ry[x][y][1]+ gy[x][y][1])))/(rho[x][y][0]+rho[x][y][1]);
	 
            }
          else
          {
          uprime_x[x][y] = 0 ; uprime_y[x][y] = 0;U[x][y]=0;V[x][y]=0;
          }

          if ( rho [x][y][0] != 0 )
          {
          ux [x][y][0] /= rho [x][y][0];
          }
          else
          {
          ux [x][y][0]= 0;
          }
          if ( rho [x][y][1] != 0 )
          {
          ux [x][y][1] /= rho [x][y][1];
          }
          else
          {
          ux [x][y][1]= 0.0;
          }

          if ( rho [x][y][0] != 0 )
          {
          uy [x][y][0] /= rho [x][y][0];
          }
          else
          {
          uy [x][y][0]= 0;
          }
          if ( rho [x][y][1] != 0 )
          {
          uy [x][y][1] /= rho [x][y][1];
          }
          else
          {
          uy [x][y][1]= 0;
          }



       }
       }
  }
    return;
}

void	project (  )

{   int	        subs,x, y,i, new_x, new_y;

    double	    udotc, u2, f0, rho_N;
    double      omega, by_temp, bx_temp;

#pragma omp parallel
  {////1

  #pragma omp  for private(y)
  for (x = 0; x < XLENGTH; x++)
   {/*force interaction*/
     for (y = 0; y < YLENGTH; y++)

        {
            psi[x][y][0] = rho [x][y][0];
            psi[x][y][1] = rho [x][y][1];
        }

   }

  #pragma omp   for private(y,i,new_x,new_y)
    for (x = 0; x < XLENGTH; x++)
     {/*force interaction*/
       for (y = 0; y < YLENGTH; y++)
         {
           bx [x][y][0] =  by [x][y][0]=  bx [x][y][1] = by [x][y][1] = 0.0;

         for (i = 0; i < Q; i++)
       	    {
            if(bnode[x][y]!=1)
            {
                new_x			= next_x[x][i];
                new_y			= next_y[y][i];
                if(cx[i]*cy[i] !=0)
                {
                    bx [x][y][0]   += t[i]* psi[new_x][new_y][0] * cx[i];/*summation of rho of streamed node * weight * ci*/
                    by [x][y][0]   += t[i]* psi[new_x][new_y][0] * cy[i];

                    bx [x][y][1]   += t[i]* psi[new_x][new_y][1] * cx[i];
                    by [x][y][1]   += t[i]* psi[new_x][new_y][1] * cy[i];
                }
                else
                {
                    bx [x][y][0]   += t[i]* psi[new_x][new_y][0] * cx[i];
                    by [x][y][0]   += t[i]* psi[new_x][new_y][0] * cy[i];

                    bx [x][y][1]   += t[i]* psi[new_x][new_y][1] * cx[i];
                    by [x][y][1]   += t[i]* psi[new_x][new_y][1] * cy[i];
                }

       	    }
       	    }

         }
      }


  #pragma omp   for private(y,bx_temp,by_temp)
  for (x = 0; x < XLENGTH; x++)
     {/*force interaction*/
       for (y = 0; y < YLENGTH; y++)


          {
              if(bnode[x][y]!=1)
              {


                bx_temp =  bx [x][y][1];
                bx [x][y][1]= - G * psi[x][y][1]*bx [x][y][0];
                bx [x][y][0]= - G * psi[x][y][0]*bx_temp;

                by_temp =  by [x][y][1];
                by [x][y][1]= - G * psi[x][y][1]*by [x][y][0];
                by [x][y][0]= - G * psi[x][y][0]*by_temp;
              }

          }

     }


  #pragma omp   for private(y, i,new_x,new_y)
     for (x = 0; x < XLENGTH; x++)
     {/*wetting force*/
       for (y = 0; y < YLENGTH; y++)


          {
           rx [x][y][0] =  ry [x][y][0] = rx [x][y][1] =  ry [x][y][1]= 0.0;

           for (i = 0; i < Q; i++)
       	    {
             new_x			= next_x[x][i];
             new_y			= next_y[y][i];
            if( bnode [new_x][new_y] == 1.0 && bnode[x][y] != 1.0 && x>xstart && x<xend )//edit this if required to x,y of rnode
             //  if( rnode [x][y] == 1.0 && bnode[x][y] != 1.0 && x>xstart && x<xend )//edit
             {
              rx [x][y][0]   += t[i]*  cx[i];//wall weight is 1 for multiple wall ,weight can be changed
              ry [x][y][0]   += t[i]*  cy[i];

              rx [x][y][1]   += t[i]*  cx[i];
              ry [x][y][1]   += t[i]*  cy[i];
             }
            else
            {
             rx [x][y][0]   += t[i]*  cx[i]* 0.0;
             ry [x][y][0]   += t[i]*  cy[i]* 0.0;

             rx [x][y][1]   += t[i]*  cx[i]* 0.0;
             ry [x][y][1]   += t[i]*  cy[i]* 0.0;
              }

            }

          rx [x][y][0] *= - Gads0 * psi[x][y][0];// positive is repulsion
          ry [x][y][0] *= - Gads0 * psi[x][y][0];
          rx [x][y][1] *= - Gads1 * psi[x][y][1];
          ry [x][y][1] *= - Gads1 * psi[x][y][1];

        }

    }

  #pragma omp   for private(y)
  for (x = 0; x < XLENGTH; x++)
     {/*gravity*/
      for (y = 0; y < YLENGTH; y++)
        {


           gx [x][y][0] =  gy [x][y][0] = gx [x][y][1] =  gy [x][y][1]= 0.0;

           gx [x][y][1]= 0 * psi [x][y][1];
           gx [x][y][0]= 0 * psi [x][y][0];

           gy [x][y][1]= 0 * psi [x][y][1];
           gy [x][y][0]= 0 * psi [x][y][0];


        }
     }
/*collision of r and b*/
  #pragma omp   for private(y, i, u2,omega,udotc,f0)
   for (x = 0; x < XLENGTH; x++)
    {
     for (y = 0; y < YLENGTH; y++)
      {

         if ( rho[x][y][0]!= 0 )
          {
            ux_eq0 [x][y] = uprime_x[x][y] +((bx[x][y][0]+ rx[x][y][0]+ gx[x][y][0])/(rho[x][y][0]* OMEGA_R) );
            uy_eq0 [x][y] = uprime_y[x][y] +((by[x][y][0]+ ry[x][y][0]+ gy[x][y][0])/(rho[x][y][0]* OMEGA_R ));

            u2    = ux_eq0[x][y] * ux_eq0[x][y] + uy_eq0[x][y] * uy_eq0[x][y] ;
					omega = OMEGA_R;
/*
	   for (i = 0; i < 1; i++)
             {	udotc  = ux_eq0[x][y] * cx[i] + uy_eq0[x][y] * cy[i];
              if ( bnode[x][y]!= 1.0 )
               {
                f0	   =  t[i] * rho[x][y][0] *(1+3.0 * udotc + 4.5 * udotc * udotc - 1.5 * u2)/*(3.0 - (2.0/RHO_R)+ 3.0 * udotc + 4.5 * udotc * udotc - 1.5 * u2);
                r[x][y][i] += omega * ( f0 - r[x][y][i]); /* collision step

               }
             }*/
           for (i = 0; i < Q; i++)
            {	udotc  = ux_eq0[x][y] * cx[i] + uy_eq0[x][y] * cy[i];
           if ( bnode[x][y]!= 1.0)// && rnode[x][y]!=1.0 )
                {


                f0	   =  t[i] * rho[x][y][0] *(1+3.0 * udotc + 4.5 * udotc * udotc - 1.5 * u2) /*((1.0/RHO_R) + 3.0 * udotc + 4.5 * udotc * udotc - 1.5 * u2)*/;
                r[x][y][i] += omega * ( f0 - r[x][y][i]);
                }

            }
          }
          else
           {
           for (i = 0; i < Q; i++)
             {
              r [x][y][i] += 0;/*no collision */
             }
           }

        if ( rho [x][y][1] != 0 )
           {
             ux_eq1 [x][y] = uprime_x[x][y] +( (bx[x][y][1]+ rx[x][y][1]+ gx[x][y][1])/(rho[x][y][1]* OMEGA_B ));
             uy_eq1 [x][y] = uprime_y[x][y] +( (by[x][y][1]+ ry[x][y][1]+ gy[x][y][1])/(rho[x][y][1] * OMEGA_B ));

             u2    = ux_eq1[x][y] * ux_eq1[x][y] + uy_eq1[x][y] * uy_eq1[x][y];
					omega = OMEGA_B;
/*
	  for (i = 0; i < 1; i++)
            {	udotc  = ux_eq1[x][y] * cx[i] + uy_eq1[x][y] * cy[i];
             if ( bnode[x][y]!= 1.0 )
               {
                f0	   =  t[i] * rho[x][y][1] * (3.0 -(2.0/RHO_B)+ 3.0 * udotc + 4.5 * udotc * udotc - 1.5 * u2);
                b[x][y][i] += omega * ( f0 - b[x][y][i]);

               }
            }*/
         for (i = 0; i < Q; i++)
           {	udotc  = ux_eq1[x][y] * cx[i] + uy_eq1[x][y] * cy[i];
           if ( bnode[x][y]!= 1.0)// && rnode[x][y]!=1.0 )
            {
              f0	   =  t[i] * rho[x][y][1] * ((1.0/RHO_B) + 3.0 * udotc + 4.5 * udotc * udotc - 1.5 * u2);
              b[x][y][i] += omega * ( f0 - b[x][y][i]);

            }
           }
          }
         else
          {
           for (i = 0; i < Q; i++)
            {
             b [x][y][i] += 0;
            }
          }
       }
    }

 }
     return;
}



void	propagate (void)
{	int	x, y, i, new_x, new_y,bounceback;
        double aux, rho0, ru,udotc,u2;

#pragma omp parallel
  {///////////////////////////////1

  #pragma omp for private(y, i, new_x, new_y)
     for (x = 0; x < XLENGTH; x++)
      {
	for (y = 0; y < YLENGTH; y++)
         {


	      for (i = 0; i < Q; i++)
		{	if(bnode[x][y]!=1)
		   {

		    new_x			= next_x[x][i];
			new_y			= next_y[y][i];

			back[x][y][i]	= r  [x][y][i];//post collision
//		back[new_x][new_y][i]	= r  [x][y][i];//post streak

		   }
		}
        }
    }


/*stream*/
#pragma omp for private(y, i,new_x,new_y)
     for (x = 0; x < XLENGTH; x++)
      {
	for (y = 0; y < YLENGTH; y++)
          {

	      for (i = 0; i < Q; i++)
                {
                 if (bnode[x][y] != 1.0)
                    {
			 new_x			= next_x[x][i];
			new_y			= next_y[y][i];


                    r[new_x][new_y][i] = back[x][y][i];
                    }

                }

            }
       }
       /*inlet
      #pragma omp for private(y)
        for (x = 0; x < 1; x++)
     {
        for (y = 0; y < YLENGTH; y++)
     {
        if (bnode[x][y] != 1.0)
        {


        if(x==0)
                        {
                        /* constant pressure
                        ux[x][y][0]=1-(r[x][y][0] + r[x][y][2]+r[x][y][4]+2*(r[x][y][3]+r[x][y][6]+r[x][y][7]))/(RHO_R);
                        r[x][y][1]=r[x][y][3]+2*(RHO_R)*ux[x][y][0]/3;
                        r[x][y][5]=r[x][y][7]-0.5*(r[x][y][2] - r[x][y][4])+(RHO_R)*ux[x][y][0]/6;
                        r[x][y][8]=r[x][y][6]+0.5*(r[x][y][2] - r[x][y][4])+(RHO_R)*ux[x][y][0]/6;*/

                        /*velocity inlet
                        rho0=(r[x][y][0] + r[x][y][2]+r[x][y][4]+2*(r[x][y][3]+r[x][y][6]+r[x][y][7]))/(1-ux0);
                        r[x][y][1]=r[x][y][3]+2*rho0*ux0/3;
                        r[x][y][5]=r[x][y][7]-0.5*(r[x][y][2] - r[x][y][4])+rho0*ux0/6;
                        r[x][y][8]=r[x][y][6]+0.5*(r[x][y][2] - r[x][y][4])+rho0*ux0/6;

                        }


       }
     }
     }*/
     /* outlet condition
     #pragma omp for private(y)
        for (x = XLENGTH-2; x < XLENGTH; x++)
     {
        for (y = 0; y < YLENGTH; y++)
     {
        if (bnode[x][y] != 1.0 )
        {


        if(x==XLENGTH-1)
                        {
                            /*extrapolation of red on other end 2nd order mohammed book
                            r[x][y][3]=2*r[x-1][y][3] - r[x-2][y][3];
                            r[x][y][6]=2*r[x-1][y][6] - r[x-2][y][6];
                            r[x][y][7]=2*r[x-1][y][7] - r[x-2][y][7];
                            r[x][y][3]=r[x-1][y][3] ;
                            r[x][y][6]=r[x-1][y][6] ;
                            r[x][y][7]=r[x-1][y][7] ;

                            /*pressure bc
                        ux[x][y][0]=-1+(r[x][y][0] + r[x][y][2]+r[x][y][4]+2*(r[x][y][1]+r[x][y][5]+r[x][y][8]))/(RHO_R-0.5);
                        r[x][y][3]=r[x][y][1]-2*(RHO_R-0.5)*ux[x][y][0]/3;
                        r[x][y][7]=r[x][y][5]+0.5*(r[x][y][2] - r[x][y][4])-(RHO_R-0.5)*ux[x][y][0]/6;
                        r[x][y][6]=r[x][y][8]-0.5*(r[x][y][2] - r[x][y][4])-(RHO_R-0.5)*ux[x][y][0]/6;*/

                        /*velocity outlet of uxo
                        rho0=-1+(r[x][y][0] + r[x][y][2]+r[x][y][4]+2*(r[x][y][1]+r[x][y][5]+r[x][y][8]))/(1+ux0);
                        r[x][y][3]=r[x][y][1]-2*(RHO_R-0.5)*(ux0)/3;
                        r[x][y][7]=r[x][y][5]+0.5*(r[x][y][2] - r[x][y][4])-(RHO_R-0.5)*(ux0)/6;
                        r[x][y][6]=r[x][y][8]-0.5*(r[x][y][2] - r[x][y][4])-(RHO_R-0.5)*(ux0)/6;
                        }

         }
     }
     }*/


//bounce back
#pragma omp for private(y,new_x,new_y,bounceback,i)
    for (x = 0; x < XLENGTH; x++)
     {
    for (y = 0; y < YLENGTH; y++)
     {
      if ( rnode [x][y] ==1.0 && bnode[x][y]!=1)
       {
      for(i=0;i<Q;i++)
	{	
			 new_x			= next_x[x][i];
			new_y			= next_y[y][i];
			bounceback              =bounce[i];
 	
		if(bnode[new_x][new_y]==1)
		{
			r[x][y][bounceback]=back[x][y][i];
		}
	}
       }
    }
   }


#pragma omp for private(y, i, new_x, new_y)
     for (x = 0; x < XLENGTH; x++)
      {
	for (y = 0; y < YLENGTH; y++)
         {

	    for (i = 0; i < Q; i++)
	     {	if(bnode[x][y]!=1)
            {


	         new_x			= next_x[x][i];
            new_y			= next_y[y][i];

            back[x][y][i]	= b	[x][y][i];
            }
	     }

        }
      }
/*stream*/
#pragma omp for private(y, i,new_x,new_y)
    for (x = 0; x < XLENGTH; x++)
     {
      for (y = 0; y < YLENGTH; y++)
       {

	   for (i = 0; i < Q; i++)
             {
              if (bnode[x][y] != 1.0)
                {
                   new_x			= next_x[x][i];
            new_y			= next_y[y][i];

         b[new_x][new_y][i] = back[x][y][i];

                }
             }

       }
    }
    /*inlet and outlet
#pragma omp for private(y, i)
for (x = 0; x < XLENGTH; x++)
  {
   for (y = 0; y < YLENGTH; y++)
   {
                /*  if(x==XLENGTH-1)
           {
              /* b[x][y][3]=2*b[x+1][y][3] - b[x+2][y][3];
               b[x][y][7]=2*b[x+1][y][7] - b[x+2][y][7];
               b[x][y][6]=2*b[x+1][y][6] - b[x+2][y][6];*/
            /*constant pressure
               ux[x][y][1]=-1+(b[x][y][0] + b[x][y][2]+b[x][y][4]+2*(b[x][y][1]+b[x][y][5]+b[x][y][8]))/RHO_R;
               b[x][y][3]=b[x][y][1]-2*RHO_R*ux[x][y][1]/3;
               b[x][y][7]=b[x][y][5]+0.5*(b[x][y][2] - b[x][y][4])-RHO_R*ux[x][y][1]/6;/*here also RHO_R to maintain same pressure at both ends i.e rho*cs^2 has to bee same
               b[x][y][6]=b[x][y][8]-0.5*(b[x][y][2] - b[x][y][4])-RHO_R*ux[x][y][1]/6;

           }
           if(x==0)extrapolation of red on other end 2nd order mohammed book
           {
               b[x][y][1]=2*b[x+1][y][1] - b[x+2][y][1];
               b[x][y][5]=2*b[x+1][y][5] - b[x+2][y][5];
               b[x][y][8]=2*b[x+1][y][8] - b[x+2][y][8];


           }

   }
  }*/
//bounce back
#pragma omp for private(y,i,new_x,new_y,bounceback)
    for (x = 0; x < XLENGTH; x++)
     {
    for (y = 0; y < YLENGTH; y++)
     {
      if ( rnode [x][y] ==1.0 && bnode[x][y]!=1)
       {
      for(i=0;i<Q;i++)
	{	
			 new_x			= next_x[x][i];
			new_y			= next_y[y][i];
			bounceback              =bounce[i];
 	
		if(bnode[new_x][new_y]==1)
		{
			b[x][y][bounceback]=back[x][y][i];
		}
	}
       }
    }
   }




  





}
}




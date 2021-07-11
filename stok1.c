 /*********************************************************************************************************************************
 Developed by Swathi Padmanabhan, Jaya Prakash, 2021
 Frontiers in Imaging Spectroscopy and Theranostics Lab, IAP, IISc.

 POLARIZED MONTE CARLO FOR MULTIPLE LAYERS

 Built upon Light Monte Carlo programs for Polaarized Light by
 Jessica C. Ramella-Roman, Steve L. Jacques and Scott A. Prahl 2005
 *  
 *  Stok1.c -Meridian Planes MC
 *	Main program for Monte Carlo simulation of photon
 *	travelling into scattering multiple layered media keeping track of 
 *  its status of polarization while updating Stokes Slab geometry.
 *  
 *  by Jessica C. Ramella-Roman
 *
 *  A report in algorithm and implementation of the program:
 *   
 * Paper: Systematic Evaluation of polarized Monte Carlo for multilayered media,2021
 * 
 *	This program is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU General Public License
 *	as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "array.h"
#include "complex.h"
#include "mie.h"
#include "nrutil.h"
#include <time.h>
#include <math.h>


#define ALIVE      	1
#define DEAD       	0
#define	NN      	100
#define THRESHOLD   0.01                /* used in roulette */
#define CHANCE      0.1                 /* used in roulette */
#define COSZERO (1.0-1.0E-12)	
  /* cosine of about 1e-6 rad. */
#define COS90D  1.0E-6		
  /* cosine of about 1.57 - 1e-6 rad. */

#define RandomNum (double) RandomGen(1, 0, NULL)
#define SIGN(x) ((x)>=0 ? 1:-1)
#define InitRandomGen (double) RandomGen(0, 1, NULL)
#define PARTIALREFLECTION    1       

/* Declare Subroutines */
void 	rotSphi(double* S, double phi, double* S2);
double	RandomGen(char Type, long Seed, long *Status);
void 	multS(double* S, double theta, double* S2);
void 	rotateXXYY(double* XX, double* YY,double* ZZ, double phi, double theta, double* XX2, double* YY2,double* ZZ2);
void 	updateU(double* U, double phi, double theta, double* U2);
double  sincos(double *x);

/*Subroutines for multilayers*/
int hitboundary(double *U,double s,double *slabsize,double *mu_sca,double *mu_abs, double s_left,double z,int current_layer);
void criticalangle(int num_layers,double *n_medium,double *crit_ang0,double *crit_ang1);
double Fresnel(int current_layer,double n1,double n2,double ca1,double *cos_ptr);
void crossUporNot(int current_layer,double *U,double z,double *n_medium,double *crit_ang0, double *crit_ang1,short photon_status,int num_layers,double rnd,double *slabsize);
void crossDownorNot(int current_layer,double* U,double z,double* n_medium,double* crit_ang0, double* crit_ang1,short photon_status,int num_layers,double rnd,double *slabsize);
void UpdateSrefl(double* S,double* S2,double *U,int current_layer,double* n_medium);
void UpdateSTrans(double* S,double* S2,double *U,int current_layer,double* n_medium);
void Crossornot(double *U,int current_layer,double z,double *n_medium,double *crit_ang0, double *crit_ang1,double *S,double *S2,int num_layers,short photon_status,double rnd,double *slabsize);
void spin(int current_layer,double I0,double I,double phi,double theta,double*U2,double*U,double *S,double **s11_t,double **s12_t,double **s33_t,double **s43_t,double*S2);

/*************** MAIN ****************************************/	
int main() {
	double pi = 3.1415926535897932384;
	/*For defining the multiple layers*/
	int num_layers; /* number of layers*/
	int current_layer; /*for iterating */
	double *slabsize;  /*For defining slab thickness*/
	double *layer_thick;/*layer thicknesses*/
	double *mu_abs,*mu_sca,*p_size,*n_particle,*n_medium;/* Arrays for various layer parameters*/
	double *crit_ang0,*crit_ang1;/*for calculating cosine of critical angle*/
	double s_new,s_left;
	double *cos_ptr,ca1;/* cos_ptr will act as temp variable in Rfresnel's calculation,ca1 will store the dirction cosine in Z dir*/
    double *tmp; /*temporary array to pass to mie function*/
    double n1,n2;/*refractive indices of incidence and transmission to pass to functions*/


	/* Mie theory stuff */
	double radius,lambda, A;
	long  nangles,i;
	struct complex m;
	
	struct complex s1_t[3][1000];
	struct complex s2_t[3][1000];

	double *mu=NULL;
	double x,qext,qsca,qback,g, rho, vol,dy, dx,hw;
	double *dy_t, *dx_t, *hw_t;
	FILE *target;
	
	double jjj;
	
	/* E field stuff */
	double	phi, theta,I,I0;
	int 	ithedeg;
	double	IT, QT, UT, VT;
	double	IR_1, QR_1, UR_1, VR_1;
	double	**IR, **QR, **UR, **VR;
	double  **IR_trans,**QR_trans,**UR_trans,**VR_trans;
   
	/* Propagation parameters */
	double	y,z;		/* photon position.  x already declared. Also, incrementals & max range. */
	double	s;          /* step sizes. s = -log(RND)/mus [cm] */
	long	i_photon;   /* current photon */
	long	Nphotons;   /* number of photons in simulation */
	short   photon_status;  /*  = ALIVE=1 or DEAD=0 */

	/* other variables */
	double	mua;        /* absorption coefficient [cm^-1] */
	double	mus;        /* scattering coefficient [cm^-1] */
	double	musp;       /* reduced scattering coefficient [cm^-1] */
	double  *albedo=NULL;     /* albedo of tissue */

	/* dummy variables */
	double  rnd;        /* assigned random value 0-1 */
	long     MM;
	double  W,absorb;  /* photon weight */
	/*double  slabsize;*/ 
	int 	j,ix,iy;
	double  cos22,sin22,costheta,sini,cosi;

	/* For calculating the weight matrix */
	double  **Weight,c,temp2,temp1;		/* A weight matrix for recording wt values and temp variables*/
	double x_new,y_new,*xgrid,*ygrid,xnew,ynew; /* for finding new grid positions*/
	double x1,y1;
	int ixnew,iynew;						/*Indices for new grid position*/

	/**** allocate matrices and arrays *******/
	double	*U, *U2;
	double	*S;     	/* */
	double	*S0;     	/* */
	double	*S2;     	/* */
	double	*s11=NULL;
	double	*s12=NULL;
	double	*s33=NULL;
	double	*s43=NULL;
	double	*IQUV;     	/* [I, Q, U, V] Stokes Vector */

	double **s11_t=NULL,**s12_t=NULL,**s33_t=NULL,**s43_t=NULL;
	double	 start_time,finish_time,temp;

	/*For recording diffues reflectance and transmittance*/
	
	double grid_ref[500][200],dr,dz,ird,iad,grid_R[500],grid_T[500],grid_trans[500][200];
	int i_r,i_z,i_a,nr,nz;	/*Reflectance_stokes_diffuse*/
	int na=1;
	double da;
	double grid_sum,grid_sumt,scale1,scale2; /*for scaling the values of reflectance and transmittance*/  
   

start_time = clock();

MM = NN - 1;

U      = new_darray(3);
U2     = new_darray(3);
S      = new_darray(4);
S0     = new_darray(4);
S2     = new_darray(4);/* dummy S*/

IQUV   = new_darray(4);

IR     =   dmatrix(0, MM, 0, MM); /* [0:MM] */
QR     =   dmatrix(0, MM, 0, MM); /* [0:MM] */
UR     =   dmatrix(0, MM, 0, MM); /* [0:MM] */
VR     =   dmatrix(0, MM, 0, MM); /* [0:MM] */
IR_trans = dmatrix(0, MM, 0, MM); /* [0:MM] */
QR_trans = dmatrix(0, MM, 0, MM); /* [0:MM] */
UR_trans = dmatrix(0, MM, 0, MM); /* [0:MM] */
VR_trans = dmatrix(0, MM, 0, MM); /* [0:MM] */

Weight = dmatrix(0, MM, 0, MM);
xgrid = new_darray(MM);
ygrid = new_darray(MM);

/**** end  allocate matrices and arrays *******/
/*Choosing the multilayers parameters*/
num_layers    = 2;
current_layer = 0;

layer_thick   = new_darray(num_layers);
slabsize      = new_darray(num_layers+1);
mu_abs        = new_darray(num_layers);
mu_sca        = new_darray(num_layers);
p_size        = new_darray(num_layers);
n_particle    = new_darray(num_layers);
n_medium      = new_darray(num_layers);
albedo        =new_darray(num_layers);

mu_abs[0]     =0.1;
mu_abs[1]     =0.7;
mu_abs[2]     =1;

p_size[0]        = 8.68/2;
p_size[1]        = 8.72/2;
p_size[2]        = 10.4/2;
n_particle[0]    = 1.55;
n_particle[1]    = 1.55;
n_particle[2]    = 1.55;
n_medium[0]      = 1.4;
n_medium[1]      = 1.3;
n_medium[2]      = 1.45;

tmp			  = new_darray(num_layers);
hw_t		  = new_darray(num_layers);
dx_t		  = new_darray(num_layers);
dy_t		  = new_darray(num_layers);
cos_ptr       = 0;
crit_ang0     =new_darray(num_layers);
crit_ang1	  =new_darray(num_layers);
s_new         =0;
s_left        =0;
n1            =0;
n2            =0;
/* CHOOSE MIE SCATTERING parameters */

lambda 		= 0.850; /* microns */
rho 		= 1.152e-4;/*Dilution 1*/
Nphotons	= 1e6;
mua 		= 0.015; /*�a  */

/* ------------------------*/
nangles 	= 1000;


/* Setup MIE SCATTERING parameters */

mu  = new_darray(nangles);

s11 = new_darray(nangles);
s12 = new_darray(nangles);
s33 = new_darray(nangles);
s43 = new_darray(nangles);
s11_t = dmatrix(0,num_layers,0,nangles);
s12_t = dmatrix(0,num_layers,0,nangles);
s33_t = dmatrix(0,num_layers,0,nangles);
s43_t = dmatrix(0,num_layers,0,nangles);
/*s1_t  = matrix(0,num_layers,0,nangles);
s2_t  = matrix(0,num_layers,0,nangles);*/

/*For defining the scattering matrix for multiple layers*/

for(current_layer=0;current_layer<num_layers;++current_layer){
	for(i=0;i<=nangles;++i){
		s11_t[current_layer][i]=0;
		s12_t[current_layer][i]=0;
		s33_t[current_layer][i]=0;
		s43_t[current_layer][i]=0;
		/*s1_t[current_layer][i]=0+i;
		s2_t[current_layer][i]=999+i;*/
		}
}
/*printf("scattering matrix=%f",s11_t[2][100]);*/
printf("All declared! \n");

for(i=0;i<=nangles;i++)
	{mu[i] = cos(pi*i/nangles);}
/*Calculate mus and mie parameters for all layers*/
for(current_layer=0;current_layer<=num_layers;++current_layer)
{	/*printf("Mie Starting here \n");*/
	struct complex*s1=NULL;
	struct complex*s2=NULL;
	s1  = new_carray(nangles);
	s2  = new_carray(nangles);
    printf("Current layer- %d \n",current_layer);
	m.re = n_particle[current_layer]/n_medium[current_layer];
	/*printf("particle index=%f,medium=%f",n_particle[current_layer],n_medium[current_layer]);*/
	m.im = 0.0;
	tmp[current_layer]    = 2*pi*p_size[current_layer]/(lambda/n_medium[current_layer]);
	
	x = tmp[current_layer]; 
	vol  = 4.0/3*pi*p_size[current_layer]*p_size[current_layer]*p_size[current_layer];
	A    = pi*p_size[current_layer]*p_size[current_layer];
    /*printf("x= %f \n",x);*/

	double *mu_tmp=NULL;
	mu_tmp=new_darray(nangles);
	mu_tmp=mu;
	Mie(x,m,mu_tmp,nangles,s1,s2,&qext,&qsca,&qback,&g); /* <---- Call Mie program ----- */
	mu=mu_tmp;
	
		for(int i=0;i<nangles;i++){
		/*printf("s1_t=%f ,\t, current layer = %d ",s1_t[current_layer][i],current_layer);*/
		s1_t[current_layer][i].re = s1[i].re;
		s1_t[current_layer][i].im = s1[i].im;
		
		/*if(i<5) {printf("s1_t.re=%f \t,s1_t.im=%f \n",s1_t[current_layer][i].re,s1_t[current_layer][i].im);}*/
		
	    s2_t[current_layer][i].re=s2[i].re;
		s2_t[current_layer][i].im=s2[i].im;
		/*if(i<5) {printf("s2_t.re=%f \t,s2_t.im=%f \n",s2_t[current_layer][i].re,s2_t[current_layer][i].im);}*/
		}
	free_carray(s1);
	free_carray(s2);
	
	/*printf("Mie done \n");*/
	
	mu_sca[current_layer] 	= qsca*A*rho*1e4; /* Mus is in cm^-1 */
	printf("musca=%f \n",mu_sca[current_layer]);
	musp 	= mu_sca[current_layer]*(1-g);/* [cm^-1] */
	albedo[current_layer] 	= mu_sca[current_layer]/(mu_sca[current_layer] + mu_abs[current_layer]);
	/*printf("musp=%d \n",musp);*/
	/*printf("albedo=%d \n",albedo[current_layer]);*/
	

	free_darray(mu);


	hw			= 7/mu_sca[0]; /* [cm] , maximum range in x and y for output. */
	dx 			= 2*hw/NN;
	dy 			= 2*hw/NN;
	/*printf("hw=%f \t dx=dy=%f\n",hw,dx);*/
	hw_t[current_layer]=hw;
	dx_t[current_layer]=dx;
	dy_t[current_layer]=dy;
	/*printf("s1_t.re=%f \t,s1_t.im=%f \n",s1_t[0][0].re,s1_t[0][0].im);
	printf("s2_t.re=%f \t,s2_t.im=%f \n",s2_t[0][0].re,s2_t[0][0].im);*/
}
/* End of Calculation of Mie function*/
/*Compute scattering matrix parameters for all the layers using S1 and S2*/
for( int j=0;j<=2;++j){
	printf("\n current layer=%d\n",j);

for( i=0;i<nangles;i++){
	
	s11_t[j][i] = 0.5*cabbs(s2_t[j][i])*cabbs(s2_t[j][i]) + 0.5*cabbs(s1_t[j][i])*cabbs(s1_t[j][i]);
	s12_t[j][i] = 0.5*cabbs(s2_t[j][i])*cabbs(s2_t[j][i]) - 0.5*cabbs(s1_t[j][i])*cabbs(s1_t[j][i]);
	s33_t[j][i] = (cmul(conj(s1_t[j][i]),s2_t[j][i])).re; 
	s43_t[j][i] = (cmul(conj(s1_t[j][i]),s2_t[j][i])).im; 
	/*if(i<5)printf("%5.5f\t %5.5f\t %5.5f\t %5.5f\n",s11_t[j][i],s12_t[j][i],s33_t[j][i],s43_t[j][i]); */
	}
}

/*Setting the slab size based on mu_sca*/
layer_thick[0]	= 0.03;
layer_thick[1]	= 0.04;
layer_thick[2]  = 0.03;

slabsize[0]    = 0;
slabsize[1]    =layer_thick[0];
slabsize[2]    =layer_thick[0] + layer_thick[1];
slabsize[3]    =layer_thick[0] + layer_thick[1]+layer_thick[2];

/*Calculate critical angles for all parameters*/
/*current_layer=1;*/
criticalangle(num_layers,n_medium,crit_ang0,crit_ang1);

printf("crit0=%f,crit1=%f",crit_ang0[1],crit_ang1[1]);

/*Scattering parameters s11 s12 s33 s43*/


 c=-0.98;
			for (iy=0;iy<NN;iy++){
			ygrid[iy]=c;
			xgrid[iy]=c;
			c=c+0.02;}
current_layer=0;

/******** MONTE CARLO *******/
 
	InitRandomGen;

/* LAUNCHNOW*/
	IT=0;/*W*/				
	QT=0;						
	UT=0;
	VT=0;
		
	IR_1=0;/*W*/				
	QR_1=0;						
	UR_1=0;
	VR_1=0;
 		
	temp=0;
	printf("Here its breaking \n"); 	
	for (iy=0; iy<NN; iy++){
				
		for (ix=0; ix<NN; ix++) {
		
			IR[iy][ix] = 0.0;
			
			QR[iy][ix] = 0.0;
			
			UR[iy][ix] = 0.0;
			
			VR[iy][ix] = 0.0;

			Weight[iy][ix] = 0.0;

		}}
  for (jjj = 1; jjj <= 4; jjj++) {

		if (jjj == 1){
		
			S0[0] = 1;	
			S0[1] = 0;
			S0[2] = 0;
			S0[3] = 0;
			printf("launch H\n");}

		if (jjj == 2){
			
			S0[0] = 1;	
			S0[1] =	0;
			S0[2] = 0;
			S0[3] = 0;
			printf("launch V\n");}

		if (jjj == 3){
			
			S0[0] = 1;	
			S0[1] =	0;
			S0[2] = 0;
			S0[3] = 0;
			printf("launch P\n");}


		if (jjj == 4){
			
			S0[0] = 1;	
			S0[1] =	0;
			S0[2] = 0;
			S0[3] = 0;	
			printf("launch R\n");}
	
       		
/* LAUNCH photon */
for (i_photon = 1; i_photon <= Nphotons; i_photon++) {

current_layer=0;
/*printf("Photon launching");*/
/*pencil beam	*/
       
		x = 0.0;
		y = 0.0;
		z = 0.0;
		x_new =0.0;
		y_new =0.0;

/* photon direction cosines */

		U[0] = 0.0;
		U[1] = 0.0;
		U[2] = 1.0;
		

		for (i=0; i<4; i++) S[i] = S0[i]; /* set incident Stokes vector to S0 */
		for (i=0; i<4; i++) S2[i] = 0.0; /* set incident Stokes vector to S0 */
	
		photon_status = ALIVE;
		
		W	= 1; /* photon weight */
	
/********* ALIVE cycle *****************/
		while (photon_status == ALIVE) {
			 if(s_left == 0.0) { /* make a new step. */
             do rnd = RandomNum; 
             while( rnd <= 0.0 ); /* avoid zero. */
	       s = -log(rnd)/(mu_abs[current_layer]+mu_sca[current_layer]);
			/*printf("s without fn=%f",s);*/
              }
           else {	 /* take the leftover. */
	       s = s_left/(mu_abs[current_layer]+mu_sca[current_layer]);
        	s_left = 0.0;
			 /*printf("s after here = %f,s_left =%f",s,s_left);*/
  			}
		  x += U[0]*s;
			y += U[1]*s;
			z += U[2]*s;
			x_new=x;
			y_new=y;
            hw_t[0]=hw;
			dx_t[0]=dx;
			dy_t[0]=dy;

			/**** ABSORB */
		
			absorb = W*(1-albedo[current_layer]);
			W-= absorb;
/*Checks inside medium or not */
			if ( z<=0) {
				phi=atan2(U[1],U[0]);
			
			    rotSphi(S, phi, S2);
			IR_1+=S2[0];
			QR_1+=S2[1];
			UR_1+=S2[2];
			VR_1+=S2[3];
		
			if (x >= -hw)
				ix = (int)(fabs(x + hw)/dx);
					
			if (y >= -hw)
				iy = (int)(fabs(y + hw)/dy);
					
			if (ix > MM) ix = MM; 
					
			if (iy > MM) iy = MM; 
				
			IR[iy][ix] += S2[0];
					
			QR[iy][ix] += S2[1];
					
			UR[iy][ix] += S2[2];
					
			VR[iy][ix] += S2[3]; 

			photon_status = DEAD;
			nr=500;nz=200;
			dr=0.01;
			da=0.5*pi/na;
			ird = sqrt(x_new*x_new+y_new*y_new)/dr;
            if(ird>nr-1) i_r= nr-1;
            else i_r = (int)ird;

			iad = acos(-U[2])/da;
			/*printf("iad = %f \n",iad);*/
            if(iad>na-1) i_a=na-1;
            else i_a = (int)iad;
			/*Defining grid for 'r' as in MCML*/

			/*for (i_r=0;i_r<nr-1;i_r++){
				for (i_a=0;i_a<200;i_a++){*/
                grid_ref[i_r][i_a]+=W;
				grid_sum += grid_ref[i_r][i_a];
				grid_R[i_r] = grid_sum;

	
         	}
		else if ( z>=slabsize[3]) {
			phi=-atan2(U[1],U[0]);
		
		rotSphi(S, phi, S2);

		IT+=S2[0]*W;
		QT+=S2[1]*W;
		UT+=S2[2]*W;
		VT+=S2[3]*W;
		if (x >= -hw)
				ix = (int)(fabs(x + hw)/dx);
					
			if (y >= -hw)
				iy = (int)(fabs(y + hw)/dy);
					
			if (ix > MM) ix = MM; 
					
			if (iy > MM) iy = MM; 
				
			IR[iy][ix] += S2[0];
					
			QR[iy][ix] += S2[1];
					
			UR[iy][ix] += S2[2];
					
			VR[iy][ix] += S2[3]; 

			photon_status = DEAD;
		
		nr=500;nz=200;
			dr=0.01;
			da=0.5*pi/na;
			ird = sqrt(x_new*x_new+y_new*y_new)/dr;
            if(ird>nr-1) i_r= nr-1;
            else i_r = (int)ird;

			iad = acos(-U[2])/da;
			/*printf("iad = %f \n",iad);*/
            if(iad>na-1) i_a=na-1;
            else i_a = (int)iad;
			/*Defining grid for 'r' as in MCML*/

			/*for (i_r=0;i_r<nr-1;i_r++){
				for (i_a=0;i_a<200;i_a++){*/
                grid_trans[i_r][i_a]+=W;
				grid_sumt += grid_trans[i_r][i_a];
				grid_T[i_r] = grid_sumt;
						
		}/*z>slab size*/
		/* Inside medium ; Boundary evaluations and Stokes updation*/
		/*Spin photons, scattering  and move steps*/ 
		if(hitboundary(U,s,slabsize,mu_sca,mu_abs,s_left,z,current_layer)){
		
		/*printf("U here= %f",U[2]);*/
		Crossornot(U,current_layer,z,n_medium,crit_ang0,crit_ang1,S,S2,num_layers,photon_status,rnd,slabsize);
		spin(current_layer,I0,I,phi,theta,U2,U,S,s11_t,s12_t,s33_t,s43_t,S2);}
        else 
		{spin(current_layer,I0,I,phi,theta,U2,U,S,s11_t,s12_t,s33_t,s43_t,S2);}
		

		x_new = x_new +U[0]*s;			
			y_new = y_new+U[1]*s;
			x1 = x_new/0.02;
			xnew=((int)x1*0.02)+0.02;
			
			/*xnew=0.02;*/
			y1= y_new/0.02;
			ynew = ((int)y1*0.02)+0.02;
			/*ynew = 0.020;*/
			
			for(ix=0;ix<NN;ix++){
				temp1=xnew;
				if(xgrid[ix]<=temp1)
				ixnew=ix;	}
			
			for (iy=0;iy<NN;iy++){
				temp2=ynew;
			if ( ygrid[iy]<=temp2) 
				iynew=iy;}
				/*printf("%f \n",temp2);*/
					
			Weight[iynew][ixnew]+=W;
			
/*ROULETTE*/
			rnd=0; while(rnd==0) rnd=RandomNum;
	
			if (W<THRESHOLD){
				if (rnd<=CHANCE)
					W/=CHANCE;
				else photon_status=DEAD;
			}
	
	
	} /* end of single photon launching */

}/* slab size*/
printf("R= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",IR_1/(Nphotons),QR_1/(Nphotons),UR_1/(Nphotons),VR_1/(Nphotons));	
printf("T= %5.5f\t %5.5f\t %5.5f\t %5.5f\n ",IT/(Nphotons),QT/(Nphotons),UT/(Nphotons),VT/(Nphotons));	
		
IT=0;
QT=0;
UT=0;
VT=0;
	
IR_1=0;
QR_1=0;
UR_1=0;
VR_1=0;

if (jjj==1) {

	target = fopen("outHI.dat","w");
	
	for (i=0; i<NN; i++) {
		/*printf("IR=%f",IR[i][0]);*/
		fprintf(target,"%5.5f", IR[i][0]);
		
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", IR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outHQ.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", QR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", QR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	target = fopen("outHU.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", QR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", QR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outHV.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", VR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", VR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("weightH.dat","w");
	
	for (i=0; i<NN; i++) {

		fprintf(target,"%5.5f", Weight[i][0]);
		
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", Weight[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	target = fopen("Reflectance_H.dat","w");
	for (i_r=0; i_r<nr-1; i_r++) {

		fprintf(target,"%5.5f", grid_R[i_r]);
		/*for(i_a=1;i_a<na;i_a++)
		fprintf(target,"\t%5.5f", grid_ref[i_r][i_a]);*/
		fprintf(target,"\n");}
	fclose(target);

	 target = fopen("Transmittance_H.dat","w");
	for (i_r=0; i_r<nr-1; i_r++) {

		fprintf(target,"%5.5f", grid_T[i_r]);
		/*for(i_a=1;i_a<na;i_a++)
		fprintf(target,"\t%5.5f", grid_ref[i_r][i_a]);*/
		fprintf(target,"\n");}
	fclose(target);

	for (i_r=0; i_r<nr-1; i_r++)
	/*for(i_a=1;i_a<na;i_a++) */{
		    grid_R[i_r]=0;
			grid_T[i_r]=0;
			}

 
	for (iy=0; iy<NN; iy++)
		for (ix=0; ix<NN; ix++) {
			IR[iy][ix] = 0.0;
			QR[iy][ix] = 0.0;
			UR[iy][ix] = 0.0;
			VR[iy][ix] = 0.0;
			Weight[iy][ix]=0;
		}
	} 
/*111111111111111111111111111111111111111111111111111111111111111*/

	if (jjj==2) {

/* save data to file */
	target = fopen("outVI.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", IR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", IR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outVQ.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", QR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", QR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outVU.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", UR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", UR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outVV.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", VR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", VR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	target = fopen("weightV.dat","w");
	
	for (i=0; i<NN; i++) {

		fprintf(target,"%5.5f", Weight[i][0]);
		
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", Weight[i][j]);
		fprintf(target,"\n");
		}
		fclose(target);

		target = fopen("Reflectance_V.dat","w");
	for (i_r=0; i_r<nr-1; i_r++) {

		fprintf(target,"%5.5f", grid_R[i_r]);
		/*for(i_a=1;i_a<na;i_a++)*/
		/*fprintf(target,"\t%5.5f", grid_ref[i_r][i_a]);*/
		fprintf(target,"\n");}
	fclose(target);

	 target = fopen("Transmittance_V.dat","w");
	for (i_r=0; i_r<nr-1; i_r++) {

		fprintf(target,"%5.5f ", grid_T[i_r]);
		/*for(i_a=1;i_a<na;i_a++)*/
		/*fprintf(target,"\t%5.5f", grid_ref[i_r][i_a]);*/
		fprintf(target,"\n");}
	fclose(target);

	for (i_r=0; i_r<nr-1; i_r++)
	/*for(i_a=1;i_a<na;i_a++)*/ {
		    grid_R[i_r]=0;
			grid_T[i_r]=0;
			}


		

	for (iy=0; iy<NN; iy++)
		for (ix=0; ix<NN; ix++) {
			IR[iy][ix] = 0.0;
			QR[iy][ix] = 0.0;
			UR[iy][ix] = 0.0;
			VR[iy][ix] = 0.0;
			Weight[iy][ix]=0;
			}
	} 
/* 222222222222222222222222222222222222222222222222222222222222222*/

	if (jjj==3) {

/* save data to file */

	target = fopen("outPI.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", IR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", IR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outPQ.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", QR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", QR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outPU.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", UR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", UR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outPV.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", VR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", VR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	target = fopen("weightP.dat","w");
	
	for (i=0; i<NN; i++) {

		fprintf(target,"%5.5f", Weight[i][0]);
		
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", Weight[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	target = fopen("Reflectance_P.dat","w");
	for (i_r=0; i_r<nr-1; i_r++) {

		fprintf(target,"%5.5f", grid_R[i_r]);
		/*for(i_a=1;i_a<na;i_a++)*/
		/*fprintf(target,"\t%5.5f", grid_ref[i_r][i_a]);*/
		fprintf(target,"\n");}
	fclose(target);

	 target = fopen("Transmittance_P.dat","w");
	for (i_r=0; i_r<nr-1; i_r++) {

		fprintf(target,"%5.5f", grid_T[i_r]);
		/*for(i_a=1;i_a<na;i_a++)*/
		/*fprintf(target,"\t%5.5f", grid_ref[i_r][i_a]);*/
		fprintf(target,"\n");}
	fclose(target);

	for (i_r=0; i_r<nr-1; i_r++)
	for(i_a=1;i_a<na;i_a++) {
		    grid_ref[i_r][i_a]=0;
			grid_trans[i_r][i_a]=0;
			}
	
	for (iy=0; iy<NN; iy++)
		for (ix=0; ix<NN; ix++) {
			IR[iy][ix] = 0.0;
			QR[iy][ix] = 0.0;
			UR[iy][ix] = 0.0;
			VR[iy][ix] = 0.0;
			Weight[iy][ix]=0.0;
		}
	} 
/* 33333333333333333333333333333333333333333333333333333333333333333333*/

	if (jjj==4) {
/* save data to file */
	target = fopen("outRI.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", IR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", IR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	
	target = fopen("outRQ.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", QR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", QR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outRU.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", UR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", UR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);

	target = fopen("outRV.dat","w");
	for (i=0; i<NN; i++) {
		fprintf(target,"%5.5f", VR[i][0]);
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", VR[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);	
	target = fopen("weightR.dat","w");
	
	for (i=0; i<NN; i++) {

		fprintf(target,"%5.5f", Weight[i][0]);
		
		for (j=1; j<NN; j++)
			fprintf(target,"\t%5.5f", Weight[i][j]);
		fprintf(target,"\n");
		}
	fclose(target);
	target = fopen("Reflectance_R.dat","w");
	for (i_r=0; i_r<nr-1; i_r++) {

		fprintf(target,"%5.5f", grid_R[i_r]);
		/*for(i_a=1;i_a<na;i_a++)
		fprintf(target,"\t%5.5f", grid_ref[i_r][i_a]);*/
		fprintf(target,"\n");}
	fclose(target);
	
	 target = fopen("Transmittance_R.dat","w");
	for (i_r=0; i_r<nr-1; i_r++) {

		fprintf(target,"%5.5f", grid_T[i_r]);
		/*for(i_a=1;i_a<na;i_a++)*/
		/*fprintf(target,"\t%5.5f", grid_ref[i_r][i_a]);*/
		fprintf(target,"\n");}
	fclose(target);
		

	}
	/* end of 4 photon launchings */ 
	}

	finish_time= clock();
	printf("Elapsed Time              = %10.2f seconds\n", (double)(finish_time-start_time)/CLOCKS_PER_SEC);

	fflush(NULL);
return 0;
} /* main routine*/ 

/*************** end MAIN ************************************/	

/*************************************************************/	
/* SUBROUTINES */
/**************************************************************************
 *	RandomGen
 *      A random number generator that generates uniformly
 *      distributed random numbers between 0 and 1 inclusive.
 *      The algorithm is based on:
 *      W.H. Press, S.A. Teukolsky, W.T. Vetterling, and B.P.
 *      Flannery, "Numerical Recipes in C," Cambridge University
 *      Press, 2nd edition, (1992).
 *      and
 *      D.E. Knuth, "Seminumerical Algorithms," 2nd edition, vol. 2
 *      of "The Art of Computer Programming", Addison-Wesley, (1981).
 *
 *      When Type is 0, sets Seed as the seed. Make sure 0<Seed<32000.
 *      When Type is 1, returns a random number.
 *      When Type is 2, gets the status of the generator.
 *      When Type is 3, restores the status of the generator.
 *
 *      The status of the generator is represented by Status[0..56].
 *
 *      Make sure you initialize the seed before you get random
 *      numbers.
 ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

double
RandomGen(char Type, long Seed, long *Status){
  static long i1, i2, ma[56];   /* ma[0] is not used. */
  long        mj, mk;
  short       i, ii;

  if (Type == 0) {              /* set seed. */
    mj = MSEED - (Seed < 0 ? -Seed : Seed);
    mj %= MBIG;
    ma[55] = mj;
    mk = 1;
    for (i = 1; i <= 54; i++) {
      ii = (21 * i) % 55;
      ma[ii] = mk;
      mk = mj - mk;
      if (mk < MZ)
        mk += MBIG;
      mj = ma[ii];
    }
    for (ii = 1; ii <= 4; ii++)
      for (i = 1; i <= 55; i++) {
        ma[i] -= ma[1 + (i + 30) % 55];
        if (ma[i] < MZ)
          ma[i] += MBIG;
      }
    i1 = 0;
    i2 = 31;
  } else if (Type == 1) {       /* get a number. */
    if (++i1 == 56)
      i1 = 1;
    if (++i2 == 56)
      i2 = 1;
    mj = ma[i1] - ma[i2];
    if (mj < MZ)
      mj += MBIG;
    ma[i1] = mj;
    return (mj * FAC);
  } else if (Type == 2) {       /* get status. */
    for (i = 0; i < 55; i++)
      Status[i] = ma[i + 1];
    Status[55] = i1;
    Status[56] = i2;
  } else if (Type == 3) {       /* restore status. */
    for (i = 0; i < 55; i++)
      ma[i + 1] = Status[i];
    i1 = Status[55];
    i2 = Status[56];
  } else
    puts("Wrong parameter to RandomGen().");
  return (0);
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC


/************************************************************************************
 *	rotSphi(S,phi,S)
 *		Rotate S by phi [radians] and return as S
 *      multiply S for the rotational matrix of Chandrasekar or Boheren and Hoffman
 *		Uses invtan()
 ****/
void 	rotSphi(double* S, double phi, double* S2) {
	double	cos2phi, sin2phi;

	cos2phi = cos(2*phi);
	sin2phi = sin(2*phi); 
	
	S2[0] = S[0]; 
	S2[1] = S[1]*cos2phi+S[2]*sin2phi; 
	S2[2] = -S[1]*sin2phi+S[2]*cos2phi;  
	S2[3] = S[3]; 
	/*printf("S2[0]=%f",S2[0]);
	printf("S2[1]=%f",S2[1]);
	printf("S2[2]=%f",S2[2]);
	printf("S2[3]=%f",S2[3]);*/
}


/**************************************************************************
 *	updateU(U,U2)
 ****/
void 	updateU(double* U, double phi, double theta, double* U2) {
	double	ux, uy, uz, uxx, uyy, uzz, temp, sintheta, costheta, sinphi, cosphi;
	double 	pi = 3.14159265358979;
	
	ux = U[0];
	uy = U[1];
	uz = U[2];
	
	costheta = cos(theta);
	sintheta = sqrt(1.0 - costheta*costheta); 
	cosphi   = cos(phi);
	if (phi < pi)
		sinphi = sqrt(1.0 - cosphi*cosphi);   
	else
		sinphi = -sqrt(1.0 - cosphi*cosphi);
	
  /* New directional cosines. */
  if (1 - fabs(uz) <= 1.0E-12) {      /* close to perpendicular. */
    uxx = sintheta * cosphi;
    uyy = sintheta * sinphi;
    uzz = costheta * SIGN(uz);   /*  SIGN(x) is faster than division. */
    } 
  else {					/* usually use this option */
    temp = sqrt(1.0 - uz * uz);
    uxx = sintheta * (ux * uz * cosphi - uy * sinphi) / temp + ux * costheta;
    uyy = sintheta * (uy * uz * cosphi + ux * sinphi) / temp + uy * costheta;
    uzz = -sintheta * cosphi * temp + uz * costheta;
    }
  /* Update directional cosines */
  U2[0] = uxx;
  U2[1] = uyy;
  U2[2] = uzz;
}
/* Modified from 10/09/2020                                                                 */
/******************************************************************************************/
/* To check for the boundary conditions and see if its getting 
*   transmitted or reflected to next layers												*/
/*******************************************************************************************/
/* Function for checking the boundary hit*/
int hitboundary(double *U,double s,double *slabsize,double *mu_sca,double *mu_abs, double s_left,double z,int current_layer)
{
   
   double dl_b;   /*distance(length) to the boundary*/
   int hit;		/*return 1 if boundary is hit or return  0*/
  
   double z0=slabsize[current_layer];
   double z1=slabsize[current_layer+1];
  /* printf("z0= %f, \t z1 = %f",z0,z1);*/

/*Distance to the boundary calculation*/
/*printf("z0=%f,z1=%f",z0,z1);*/
if (U[2]>0.0)
	{dl_b=(z1-z)/U[2];
	/*printf("uz<0 \n");*/ }
else if(U[2]<0.0)
    {dl_b=(z0-z)/U[2];
	/*printf("uz>0 \n");*/ }
/*Horizontal crossing*/

if (U[2] !=0.0 && s>dl_b){
	
	double mu_t = mu_abs[current_layer]+mu_sca[current_layer];
	/*printf("mut= %f",mu_t);
	printf("s=%f, dl_b=%f",s,dl_b);*/
	s_left=(s-dl_b)*mu_t;
	/*printf("s_left= %f",s_left);*/
	s = dl_b;
	hit = 1;
}
else 
hit=0;
/*printf("hit=%d \n",hit);*/
return(hit);
}
/*calculating and storing the sleft and the dlb for computing new direction cosines
  before checking the boundary*/

/*s_new = s_left/(mu_abs[current_layer]+mu_sca[current_layer]);*/
/*printf("s=%f",s);*/

 
/*Computing critical angles for all layers and parameters*/
void criticalangle(int num_layers,double* n_medium,double *crit_ang0,double *crit_ang1){
	
	double n1,n2;
	for(int j=0; j<=num_layers; j++)  {
    n1 = n_medium[j];
    n2 = n_medium[j-1];
        crit_ang0[j] = n1>n2 ? 
		sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;
    /*printf("crit angle=%f",crit_ang0[current_layer]);*/
    n2 =n_medium[j+1];
        crit_ang1[j] = n1>n2 ? 
		sqrt(1.0 - n2*n2/(n1*n1)) : 0.0;
	}
	
}

/* Function for calculating Frensel's reflectances*/
double Fresnel(int current_layer,double n1,double n2,double ca1,double *cos_ptr){
double r;
/*printf("ca1=%f",ca1);
printf("n1=%f,n2=%f",n1,n2);*/
if(n1==n2){
	/*printf("ca1=%f",ca1);*/

      *cos_ptr = ca1;
    r = 0.0;
/*printf("r=%f",r);*/
	
}
else if(ca1>COSZERO) {	/** normal incident. **/
    *cos_ptr = ca1;
    r = (n2-n1)/(n2+n1);
    r *= r;
  }
else if(ca1<COS90D)  {	/** very slant. **/
    *cos_ptr = 0;
    r = 1.0;
}
  else  {			  		/** general. **/
    double sa1, sa2;	 /* sine of the incident and transmission angles. */
	double ca2;
    
    sa1 = sqrt(1-ca1*ca1);
    sa2 = n1*sa1/n2;
    if(sa2>=1.0) {	
	  /* double check for total internal reflection. */
      *cos_ptr = 0.0;
      r = 1.0;
    }
    else  {
      double cap, cam;	/* cosines of the sum ap or */
						/* difference am of the two */
						/* angles. ap = a1+a2 */
						/* am = a1 - a2. */
      double sap, sam;	/* sines. */
      
      ca2 = sqrt(1-sa2*sa2);
      
      cap = ca1*ca2 - sa1*sa2; /* c+ = cc - ss. */
	  cam = ca1*ca2 + sa1*sa2; /* c- = cc + ss. */
      sap = sa1*ca2 + ca1*sa2; /* s+ = sc + cs. */
      sam = sa1*ca2 - ca1*sa2; /* s- = sc - cs. */
	  /*printf("cap=%f ,cam=%f , sap=%f,sam =%f",cap,cam,sap,sam);*/
      r = 0.5*sam*sam*(cam*cam+cap*cap)/(sap*sap*cam*cam); 
	  
		/* rearranged for speed. */
    }
  }
  /*printf("r= %f \n",r);*/
  return(r);
}

/*Based on this above function, determine updation of S*/  
/*Updating S after reflection*/
void UpdateSrefl(double* S,double* S2,double *U,int current_layer,double* n_medium){
 double a_i; /*angle of incidence*/
 double a_t; /*angle of transmision*/
 double a_plus,a_minus;/*the sum and difference of a_i and a_t*/  
 /*double **m_array,temp;*/
 double m_array[4][4],temp;
 short tempx,tempy;
 /* m_array =dmatrix(0, 3, 0, 3);*/
  for (tempx=0;tempx<=3;tempx++){
  	for(tempy=0;tempy<=3;tempy++){
		  m_array[tempx][tempy]=0;}}

 
	double ni=n_medium[current_layer];
	double nt;
	if(current_layer==0)
	nt=1;
	else
	nt=n_medium[current_layer-1];
	/*printf("U here=%f",U[2]);*/
 	/*a_i = acos(U[2]);*/
	a_i = acos(-U[2]);
	 /*finding the a_t using Snell's law*/
	a_t = asin((ni/nt)*sin(a_i));
	/*a_t = atan(nt/ni);*/
	/*printf("a_t=%f",a_t);
	printf("a_i=%f",a_i);*/
	a_plus=a_i+a_t;
	a_minus=a_i-a_t;
	m_array[0][0]=pow(cos(a_minus),2)+pow(cos(a_plus),2);
	m_array[0][1]=pow(cos(a_minus),2)-pow(cos(a_plus),2);
	m_array[1][0]=pow(cos(a_minus),2)-pow(cos(a_plus),2);
	m_array[1][1]=pow(cos(a_minus),2)+pow(cos(a_plus),2);
	m_array[2][2]=-(2*cos(a_plus)*cos(a_minus));
	m_array[3][3]=-(2*cos(a_plus)*cos(a_minus));
    temp = 0.5*pow((tan(a_minus)/sin(a_plus)),2);
  /*printf("m_array=%f,temp=%f",m_array[2][2],temp);*/
	S2[0]=(m_array[0][0]*S[0]+m_array[1][0]*S[1])*temp;
	S2[1]=(m_array[1][0]*S[0]+m_array[0][0]*S[1])*temp;
	S2[2]=(m_array[2][2]*S[2])*temp;
	S2[3]=(m_array[3][3]*S[3])*temp;
   /* printf("S2[0]=%f",S2[0]);
	printf("S2[1]=%f",S2[1]);
	printf("S2[2]=%f",S2[2]);
	printf("S2[3]=%f",S2[3]);
	printf("S[0]=%f",S[0]);
	printf("S[1]=%f",S[1]);
	printf("S[2]=%f",S[2]);
	printf("S[3]=%f",S[3]);*/
           
           
		/*return(S2);*/ 
   
}

/*Now to update the Stokes vectors based on the Transmittance from Fresnel's */
void UpdateSTrans(double* S,double* S2,double *U,int current_layer,double* n_medium){
 double a_i; /*angle of incidence*/
 double a_t; /*angle of transmision*/
 double a_plus,a_minus;/*the sum and difference of a_i and a_t*/  
 /*double **m_array,temp;*/
 double m_array[4][4],temp;
 short tempx,tempy;
  /*m_array =dmatrix(0, 3, 0, 3);*/
  for (tempx=0;tempx<=3;tempx++){
  	for(tempy=0;tempy<=3;tempy++){
		  m_array[tempx][tempy]=0;}}


	double ni=n_medium[current_layer];
	double nt;
	if(current_layer==2)
	nt=1;
	else
	nt=n_medium[current_layer+1];
	/*printf("U here=%f",U[2]);*/
 	a_i = acos(U[2]);
	 /*finding the a_t using Snell's law*/
	a_t = asin((ni*sin(a_i))/nt);
	/*a_t = atan(nt/ni);*/
	/*printf("a_t=%f",a_t);
	printf("a_i=%f",a_i);*/
	a_plus=a_i+a_t;
	a_minus=a_i-a_t;
	m_array[0][0]=pow(cos(a_minus),2)+1;
	m_array[0][1]=pow(cos(a_minus),2)-1;
	m_array[1][0]=pow(cos(a_minus),2)-1;
	m_array[1][1]=pow(cos(a_minus),2)+1;
	m_array[2][2]=(2*cos(a_minus));
	m_array[3][3]=(2*cos(a_minus));
    temp = 0.5*(sin(2*a_i)*sin(2*a_t))/pow((sin(a_plus)*cos(a_minus)),2);
 /*printf("m_array=%f,temp=%f",m_array[2][2],temp);*/
	S2[0]=(m_array[0][0]*S[0]+m_array[1][0]*S[1])*temp;
	S2[1]=(m_array[1][0]*S[0]+m_array[0][0]*S[1])*temp;
	S2[2]=(m_array[2][2]*S[2])*temp;
	S2[3]=(m_array[3][3]*S[3])*temp;
/*printf("S2[0]=%f",S2[0]);
	printf("S2[1]=%f",S2[1]);
	printf("S2[2]=%f",S2[2]);
	printf("S2[3]=%f",S2[3]);
	printf("S[0]=%f",S[0]);
	printf("S[1]=%f",S[1]);
	printf("S[2]=%f",S[2]);
	printf("S[3]=%f",S[3]);*/
   

	/*return(S2);*/ 
   
}
/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/
/*Check boundary crossing up or down*/
void crossUporNot(int current_layer,double* U,double z,double *n_medium,double* crit_ang0, double* crit_ang1,short photon_status,int num_layers,double rnd,double *slabsize){
	double *u1 = U; /* z directional cosine. */
	double uz=u1[2];
   double uz1;	/*cosines of transmission alpha. always */
				/* positive. */
    double r=0.0;	/* reflectance */
    double ni = n_medium[current_layer];
	
    double nt=n_medium[current_layer-1];
		/*Get the value of r*/
	if( (-1* uz) <= crit_ang0[current_layer]) 
    {r=1.0;		
	      /* total internal reflection. */
	/*printf("r_tir=%f \n",r);*/}
    else 
	{	if(current_layer==0){
		 /*printf("photon dead here 1");*/
		 photon_status=DEAD;}
		else r = Fresnel(current_layer,ni,nt,(-1*uz),&uz1);
	/*printf("r_fres= %fn \n",r);*/}

	if(rnd > r ) {		/* transmitted to layer-1. */
	
    if(current_layer==0)  {
    
	    uz=-uz1;
     
	  /*printf("photon dead here 2");*/
	   photon_status=DEAD;
    }
    else {
      current_layer--;

      u1[0] *= ni/nt;
      u1[1] *= ni/nt;
      u1[2] = -uz1;
    }
  }
  else 						/* reflected. */
    uz = -uz;
	/*printf("uz=%f",uz);*/
  }
  void crossDownorNot(int current_layer,double* U,double z,double* n_medium,double* crit_ang0, double* crit_ang1,short photon_status,int num_layers,double rnd,double *slabsize){
	double *u1 = U; /* z directional cosine. */
	double uz=u1[2];
    double uz1;	/*cosines of transmission alpha. always */
				/* positive. */
    double r=0.0;	/* reflectance */
    
    double ni = n_medium[current_layer];
    double nt = n_medium[current_layer+1];
	/*printf("r=%f \n",r);*/
	/*if(current_layer==1) printf("current layer=%d",current_layer);*/
	/*Get the value of r*/

	if(  uz <=crit_ang1[current_layer]) 
   { r=1.0;		      /* total internal reflection. */
   /*printf("r_tir=%f \n",r);*/}
    else 	
		if(current_layer==num_layers)

	{ /*printf("current layer=%d",current_layer);*/
		
		photon_status=DEAD;}
		else
		{   /*printf("Breaking");*/
			r = Fresnel(current_layer,ni,nt,uz,&uz1);
			
		    /*printf("r= %f",r);*/
		}

	if(rnd > r) {		/* transmitted to layer-1. */
    if(current_layer==num_layers)  {
      uz = uz1;
	 
	
	   photon_status=DEAD;
    }
    else {
      current_layer++;
	  
      u1[0] *= ni/nt;
      u1[1] *= ni/nt;
      u1[2] = uz1;
    }
  }
  else 						/* reflected. */
    uz = -uz;
	
  }


/*To check cross or not on the whole*/
void Crossornot(double* U,int current_layer,double z,double* n_medium,double* crit_ang0, double* crit_ang1,double* S,double* S2,int num_layers,short photon_status,double rnd, double *slabsize){
 /*printf("U[2]=%f",U[2]);
 printf("current layer= %d",current_layer);*/
 if (U[2]<0.0)
   
   {crossUporNot(current_layer,U,z,n_medium,crit_ang0,crit_ang1,photon_status,num_layers,rnd,slabsize);
    
  UpdateSrefl(S,S2,U,current_layer,n_medium);
   }
 else  
	{crossDownorNot(current_layer,U,z,n_medium,crit_ang0,crit_ang1,photon_status,num_layers,rnd,slabsize);

	UpdateSTrans(S,S2,U,current_layer,n_medium);
	}
	}

	void spin(int current_layer,double I0,double I,double phi,double theta,double*U2,double*U,double *S,double **s11_t,double **s12_t,double **s33_t,double **s43_t,double*S2){
    /* SPIN */
   /*printf("Is now in else loop \n \n");*/
/* REJECTION METHOD to choose azimuthal angle phi and deflection angle theta */
			int ithedeg,i;	
			double pi = 3.1415926535897932384;
			long nangles=1000;
			double costheta,temp,cosi,sini,sin22,cos22;
/*For recording */
		double x1,y1,temp1,temp2,c;
		int ix,iy;
		double x_new,y_new,*xgrid,*ygrid,xnew,ynew;
		int ixnew,iynew;
			do{ theta 	= acos(2*RandomNum-1);       	
			
			    phi = RandomNum*2.0*pi;              
		   							
				I0=s11_t[current_layer][0]*S[0]+s12_t[current_layer][0]*(S[1]*cos(2*phi)+S[2]*sin(2*phi));			
	                        
	 			ithedeg = floor(theta*nangles/pi);  
                                                 
				I=s11_t[current_layer][ithedeg ]*S[0]+s12_t[current_layer][ithedeg]*(S[1]*cos(2*phi)+S[2]*sin(2*phi));
				
			}while(RandomNum*I0>=I);
			
  /*------------------------------------------------------------------------------	
   Scattering : rotate to meridian plane	then scatter															
------------------------------------------------------------------------------*/
	
			updateU(U, phi, theta, U2);  /* update photon trajectory vector */
						
			costheta=cos(theta);

			rotSphi(S, phi, S2);

			S[0]= s11_t[current_layer][ithedeg]*S2[0]+s12_t[current_layer][ithedeg]*S2[1];
				
			S[1]= s12_t[current_layer][ithedeg]*S2[0]+s11_t[current_layer][ithedeg]*S2[1];
	
			S[2]= s33_t[current_layer][ithedeg]*S2[2]+s43_t[current_layer][ithedeg]*S2[3];
				
			S[3]= -s43_t[current_layer][ithedeg]*S2[2]+s33_t[current_layer][ithedeg]*S2[3];

			/*printf("s11_t[current_layer][ithedeg]=%f,\n",s11_t[current_layer][ithedeg]);*/
			temp=(sqrt(1-costheta*costheta)*sqrt(1-U2[2]*U2[2]));
			
			if ( temp==0){
				cosi=0;}
			else{
			
				if ((phi>pi) & (phi<2*pi))
					cosi=(U2[2]*costheta-U[2])/temp;	
				else
					cosi=-(U2[2]*costheta-U[2])/temp;	
				if (cosi<-1) cosi=-1;
				if (cosi>1) cosi=1;
				}
		
			sini = sqrt(1-cosi*cosi);
	
			cos22=2*cosi*cosi-1;
			
			sin22=2*sini*cosi;
			
			S2[0]=S[0];
			
			S2[1]=(S[1]*cos22-S[2]*sin22);
					
			S2[2]=(S[1]*sin22+S[2]*cos22);
					
			S2[3]=S[3];


			S[1]= S2[1]/S2[0];	
			S[2]= S2[2]/S2[0];
			S[3]= S2[3]/S2[0];
			S[0]= 1.0;
           /* printf("S[2]=%f,S2[0]=%f",S[2],S2[0]);*/
   
			
			for (i=0; i<3; i++) U[i] = U2[i]; /* update U */
			
				}

	

/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - solver.c

Content:
  -This file contains the implementation of the Riemann solvers
  
*/


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>

#include "definitions.h"
#include "structures.h"
#include "solvers.h"
#include "closures.h"



void compute_euler_HLLE(t_wall *wall,double *lambda_max){

	int m;
	double WR[6], WL[6]; /**<Auxiliar array of variables rotated for the 1D problem**/
#if ST==2||ST==3 //fluctuation version for source terms	
	double WRe[6], WLe[6];
	double pLe, pRe;
#endif
	double WRprime[6], WLprime[6];
	double uL, uR, vL, vR, wL, wR, pL, pR, HL, HR, cL, cR, gammaL, gammaR;
	double raizrhoR, raizrhoL, sumRaizRho;
	double u_hat, v_hat, w_hat, H_hat, c_hat, gamma_hat;
	double S1, S2, diffS, maxS;
	double FR[5], FL[5];
	double F_star[5];
#if MULTICOMPONENT	
	double phiL, phiR;
#endif

	/**Rotation of array of variables**/

	WR[0]=wall->UR[0];
	WL[0]=wall->UL[0];

    // This is a simplification of the rotation matrix, only valid for cartesian mesh
	WR[1]=wall->UR[1]*wall->nx+wall->UR[2]*wall->ny+wall->UR[3]*wall->nz;
	WL[1]=wall->UL[1]*wall->nx+wall->UL[2]*wall->ny+wall->UL[3]*wall->nz;

	WR[2]=-wall->UR[1]*wall->ny+wall->UR[2]*wall->nx+wall->UR[2]*wall->nz;
	WL[2]=-wall->UL[1]*wall->ny+wall->UL[2]*wall->nx+wall->UL[2]*wall->nz;

      WR[3]=wall->UR[3]*wall->nx+wall->UR[3]*wall->ny-wall->UR[1]*wall->nz;
	WL[3]=wall->UL[3]*wall->nx+wall->UL[3]*wall->ny-wall->UL[1]*wall->nz;

	WR[4]=wall->UR[4];
	WL[4]=wall->UL[4];


	#if ST==2||ST==3
	WRe[0]=wall->URe[0];
	WLe[0]=wall->ULe[0];

	WRe[1]=wall->URe[1]*wall->nx+wall->URe[2]*wall->ny+wall->URe[3]*wall->nz;
	WLe[1]=wall->ULe[1]*wall->nx+wall->ULe[2]*wall->ny+wall->ULe[3]*wall->nz;

	WRe[2]=-wall->URe[1]*wall->ny+wall->URe[2]*wall->nx+wall->URe[2]*wall->nz;
	WLe[2]=-wall->ULe[1]*wall->ny+wall->ULe[2]*wall->nx+wall->ULe[2]*wall->nz;

      WRe[3]=wall->URe[3]*wall->nx+wall->URe[3]*wall->ny-wall->URe[1]*wall->nz;
	WLe[3]=wall->ULe[3]*wall->nx+wall->ULe[3]*wall->ny-wall->ULe[1]*wall->nz;

	WRe[4]=wall->URe[4];
	WLe[4]=wall->ULe[4];
	#endif

	for(m=0;m<5;m++){
		WRprime[m]=WR[m];
		WLprime[m]=WL[m];
	}
	#if ST==2||ST==3
		WRprime[0]=WR[0]-WRe[0];
		WLprime[0]=WL[0]-WLe[0];
		WRprime[2]=WR[2]-WRe[2];
		WLprime[2]=WL[2]-WLe[2];
		WRprime[3]=WR[3]-WRe[3];
		WLprime[3]=WL[3]-WLe[3];
		WRprime[4]=WR[4]-WRe[4];
		WLprime[4]=WL[4]-WLe[4];
      #endif

#if MULTICOMPONENT
	WR[5]=wall->UR[5];
	WL[5]=wall->UL[5];
	phiL=WL[5]/WL[0];
	phiR=WR[5]/WR[0];
	#if MULTI_TYPE==1
		gammaL=phiL;
		gammaR=phiR;
	#else
		gammaL=1.0+1.0/phiL;
		gammaR=1.0+1.0/phiR;
	#endif
#else
	gammaL=_gamma_;
	gammaR=_gamma_;
#endif

	/**Additional variables for the solver**/
	uL=WL[1]/WL[0];
	uR=WR[1]/WR[0];

	vL=WL[2]/WL[0];
	vR=WR[2]/WR[0];

      wL=WL[3]/WL[0];
	wR=WR[3]/WR[0];

	pL=pressure_from_energy(gammaL, WL[4], uL, vL, wL, WL[0], wall->z);
	pR=pressure_from_energy(gammaR, WR[4], uR, vR, wR, WR[0], wall->z);

#if ST==3
	HL=(WL[4]-WL[0]*_g_*wall->z+pL)/WL[0];
	HR=(WR[4]-WR[0]*_g_*wall->z+pR)/WR[0];
#else
	HL=(WL[4]+pL)/WL[0];
	HR=(WR[4]+pR)/WR[0];
#endif

	cL=sqrt(gammaL*pL/WL[0]);

	cR=sqrt(gammaR*pR/WR[0]);

	raizrhoL=sqrt(WL[0]);
	raizrhoR=sqrt(WR[0]);
	sumRaizRho=raizrhoR+raizrhoL;

	/**Hat variables (Roe averages)**/
	u_hat=(uR*raizrhoR+uL*raizrhoL)/sumRaizRho;
	v_hat=(vR*raizrhoR+vL*raizrhoL)/sumRaizRho;
      w_hat=(wR*raizrhoR+wL*raizrhoL)/sumRaizRho;
	H_hat=(HR*raizrhoR+HL*raizrhoL)/sumRaizRho;
#if MULTICOMPONENT
	#if MULTI_TYPE==1
		gamma_hat=1.0+ 1.0/( (phiR*raizrhoR+phiL*raizrhoL)/sumRaizRho );
	#else
		gamma_hat=(gammaR*raizrhoR+gammaL*raizrhoL)/sumRaizRho;
	#endif
#else
	gamma_hat=_gamma_;
#endif

	c_hat=sqrt((gamma_hat-1)*(H_hat-0.5*(u_hat*u_hat+v_hat*v_hat+w_hat*w_hat)));

	/**Physical flux calculation (the F of the eqs.)**/

	FR[0]=WR[1];
	FL[0]=WL[1];

#if ST==2||ST==3
	pRe  =wall->pRe;
	pLe  =wall->pLe;
	FR[1]=WR[1]*uR+(pR-pRe);
	FL[1]=WL[1]*uL+(pL-pLe);
#else
	FR[1]=WR[1]*uR+pR;
	FL[1]=WL[1]*uL+pL;
#endif

	FR[2]=WR[1]*vR;
	FL[2]=WL[1]*vL;

	FR[3]=WR[1]*wR;
	FL[3]=WL[1]*wL;

	FR[4]=uR*(WR[4]+pR);
	FL[4]=uL*(WL[4]+pL);


	/**Wave speed estimation**/

	S1=MIN(uL-cL,u_hat-c_hat);
	S2=MAX(uR+cR, u_hat+c_hat);

	maxS=MAX(ABS(S1),ABS(S2));
	diffS=S2-S1;

	/**HLLE flux calculation**/
	for(m=0;m<5;m++){
		if(S1>=0){
			F_star[m]=FL[m];
		}else if(S2<=0){
			F_star[m]=FR[m];
		}else{
                  F_star[m]=(S2*FL[m]-S1*FR[m]+S1*S2*(WRprime[m]-WLprime[m]))/(diffS);
            }
	}

	/**Inverse rotation of the flux**/
	wall->fR_star[0]=F_star[0]; //Mass is not vectorial
	wall->fR_star[1]=F_star[1]*wall->nx - F_star[2]*wall->ny - F_star[3]*wall->nz;
	wall->fR_star[2]=F_star[1]*wall->ny + F_star[2]*wall->nx + F_star[2]*wall->nz;
	wall->fR_star[3]=F_star[3]*wall->nx + F_star[3]*wall->ny + F_star[1]*wall->nz;
	wall->fR_star[4]=F_star[4]; //Energy is not vectorial

	for(m=0;m<5;m++){
		wall->fL_star[m]=wall->fR_star[m];
	}

	//printf("%14.14e %14.14e %14.14e\n",FR[0],FL[0],F_star[0]);


	*lambda_max=MAX(*lambda_max,maxS);

}

void compute_euler_HLLC(t_wall *wall,double *lambda_max){

	int m;
	double WR[5], WL[5]; /**<Auxiliar array of variables rotated for the 1D problem**/
	double uL, uR, vL, vR, wR, wL, pL, pR, HL, HR, cL, cR;
	double raizrhoR, raizrhoL, sumRaizRho;
	double u_hat, v_hat, w_hat, H_hat, c_hat;
	double S1, S2, maxS, S_star;
      double uK, vK, wK, rhoK, SK, pK, EK, aux;
	double FR[5], FL[5];
	double F_star[5],W_star[5];

	/**Rotation of array of variables**/
	
	WR[0]=wall->UR[0];
	WL[0]=wall->UL[0];

	// This is a simplification of the rotation matrix, only valid for cartesian mesh
	WR[1]=wall->UR[1]*wall->nx+wall->UR[2]*wall->ny+wall->UR[3]*wall->nz;
	WL[1]=wall->UL[1]*wall->nx+wall->UL[2]*wall->ny+wall->UL[3]*wall->nz;

	WR[2]=-wall->UR[1]*wall->ny+wall->UR[2]*wall->nx+wall->UR[2]*wall->nz;
	WL[2]=-wall->UL[1]*wall->ny+wall->UL[2]*wall->nx+wall->UL[2]*wall->nz;
      
	WR[3]=wall->UR[3]*wall->nx+wall->UR[3]*wall->ny-wall->UR[1]*wall->nz;
	WL[3]=wall->UL[3]*wall->nx+wall->UL[3]*wall->ny-wall->UL[1]*wall->nz;
      


	WR[4]=wall->UR[4];
	WL[4]=wall->UL[4];
	
	/**Additional variables for the solver**/
	uL=WL[1]/WL[0];
	uR=WR[1]/WR[0];

	vL=WL[2]/WL[0];
	vR=WR[2]/WR[0];
      
      wL=WL[3]/WL[0];
	wR=WR[3]/WR[0];

	pL=(_gamma_-1.0)*(WL[4]-0.5*WL[0]*(uL*uL+vL*vL+wL*wL));
	pR=(_gamma_-1.0)*(WR[4]-0.5*WR[0]*(uR*uR+vR*vR+wR*wR));
	
	HL=(WL[4]+pL)/WL[0];
	HR=(WR[4]+pR)/WR[0];
	
	cL=sqrt(_gamma_*pL/WL[0]);

	cR=sqrt(_gamma_*pR/WR[0]);
	
	raizrhoL=sqrt(WL[0]);
	raizrhoR=sqrt(WR[0]);
	sumRaizRho=raizrhoR+raizrhoL;

	/**Hat variables (Roe averages)**/
	u_hat=(uR*raizrhoR+uL*raizrhoL)/sumRaizRho;
	v_hat=(vR*raizrhoR+vL*raizrhoL)/sumRaizRho;
      w_hat=(wR*raizrhoR+wL*raizrhoL)/sumRaizRho;
	H_hat=(HR*raizrhoR+HL*raizrhoL)/sumRaizRho;

	c_hat=sqrt((_gamma_-1)*(H_hat-0.5*(u_hat*u_hat+v_hat*v_hat+w_hat*w_hat)));

	/**Physical flux calculation (the F of the eqs.)**/

	FR[0]=WR[1];
	FL[0]=WL[1];

	FR[1]=WR[1]*uR+pR;
	FL[1]=WL[1]*uL+pL;

	FR[2]=WR[1]*vR;
	FL[2]=WL[1]*vL;
      
      FR[3]=WR[1]*wR;
	FL[3]=WL[1]*wL;

	FR[4]=uR*(WR[4]+pR);
	FL[4]=uL*(WL[4]+pL);


	/**Wave speed estimation**/
	
	S1=MIN(uL-cL,u_hat-c_hat);
	S2=MAX(uR+cR, u_hat+c_hat);
      
      //S1=MIN(0.0,S1);
      //S2=MAX(0.0,S2);

	maxS=MAX(ABS(S1),ABS(S2));
      
      S_star= ( pR-pL + WL[1]*(S1-uL) - WR[1]*(S2-uR) ) / ( WL[0]*(S1-uL) - WR[0]*(S2-uR) );

	/**HLLC flux calculation**/
	
      if(S1>=0){
            for(m=0;m<5;m++){
                  F_star[m]=FL[m];
            }
      }else if(S2<=0){
            for(m=0;m<5;m++){
                  F_star[m]=FR[m];
            }
      }else{
            if(S_star<=0){
                  uK=uR ;
                  vK=vR;
                  wK=wR;
                  rhoK=WR[0];
                  SK=S2;
                  pK=pR;
                  EK=WR[4];
            }else{
                  uK=uL; 
                  vK=vL;
                  wK=wL;
                  rhoK=WL[0];
                  SK=S1;
                  pK=pL;
                  EK=WL[4];            
            }
            aux=rhoK*(SK-uK)/(SK-S_star);
            W_star[0]=aux;
            W_star[1]=aux*S_star;
            W_star[2]=aux*vK;
            W_star[3]=aux*wK;
            W_star[4]=aux*( EK/rhoK + (S_star-uK)*(S_star + pK/(rhoK*(SK-uK))) );
            for(m=0;m<5;m++){
                  if(S_star<=0){
                        F_star[m]=FR[m]+S2*(W_star[m]-WR[m]);
                  }else{
                        F_star[m]=FL[m]+S1*(W_star[m]-WL[m]);
                  }
            }
            
      }
	
	    
      /**Inverse rotation of the flux**/
	wall->fR_star[0]=F_star[0]; //Mass is not vectorial
	wall->fR_star[1]=F_star[1]*wall->nx - F_star[2]*wall->ny - F_star[3]*wall->nz;
	wall->fR_star[2]=F_star[1]*wall->ny + F_star[2]*wall->nx + F_star[2]*wall->nz;
	wall->fR_star[3]=F_star[3]*wall->nx + F_star[3]*wall->ny + F_star[1]*wall->nz;
	wall->fR_star[4]=F_star[4]; //Energy is not vectorial
	
      for(m=0;m<5;m++){
            wall->fL_star[m]=wall->fR_star[m];
      }
      
      
	*lambda_max=MAX(*lambda_max,maxS);
		
}

void compute_euler_HLLS(t_wall *wall,double *lambda_max, t_sim *sim){

	int m;
	double WR[6], WL[6], S[6], B[6]; /**<Auxiliar array of variables rotated for the 1D problem**/
	double uL, uR, vL, vR, wL, wR, pL, pR, HL, HR, gammaL, gammaR;
      double pLe, pRe, rhoLe, rhoRe;
	double raizrhoR, raizrhoL, sumRaizRho;
	double u_hat, v_hat, w_hat, H_hat, c_hat, gamma_hat,psi,chi;
	double S1, S2, diffS, maxS;
	double FR[5], FL[5];
	double F_star[5];
#if MULTICOMPONENT	
	double phiL, phiR;
#endif
	/**Rotation of array of variables**/

	WR[0]=wall->UR[0];
	WL[0]=wall->UL[0];

      // This is a simplification of the rotation matrix, only valid for cartesian mesh
	WR[1]=wall->UR[1]*wall->nx+wall->UR[2]*wall->ny+wall->UR[3]*wall->nz;
	WL[1]=wall->UL[1]*wall->nx+wall->UL[2]*wall->ny+wall->UL[3]*wall->nz;

	WR[2]=-wall->UR[1]*wall->ny+wall->UR[2]*wall->nx+wall->UR[2]*wall->nz;
	WL[2]=-wall->UL[1]*wall->ny+wall->UL[2]*wall->nx+wall->UL[2]*wall->nz;

	WR[3]=wall->UR[3]*wall->nx+wall->UR[3]*wall->ny-wall->UR[1]*wall->nz;
	WL[3]=wall->UL[3]*wall->nx+wall->UL[3]*wall->ny-wall->UL[1]*wall->nz;

	WR[4]=wall->UR[4];
	WL[4]=wall->UL[4];

#if MULTICOMPONENT
	WR[5]=wall->UR[5];
	WL[5]=wall->UL[5];
	phiL=WL[5]/WL[0];
	phiR=WR[5]/WR[0];
	#if MULTI_TYPE==1
		gammaL=phiL;
		gammaR=phiR;
	#else
		gammaL=1.0+1.0/phiL;
		gammaR=1.0+1.0/phiR;
	#endif
#else
	gammaL=_gamma_;
	gammaR=_gamma_;
#endif

	/**Additional variables for the solver**/
	uL=WL[1]/WL[0];
	uR=WR[1]/WR[0];

	vL=WL[2]/WL[0];
	vR=WR[2]/WR[0];

      wL=WL[3]/WL[0];
	wR=WR[3]/WR[0];

	pL=(gammaL-1.0)*(WL[4]-0.5*WL[0]*(uL*uL+vL*vL+wL*wL));
	pR=(gammaR-1.0)*(WR[4]-0.5*WR[0]*(uR*uR+vR*vR+wR*wR));

	HL=(WL[4]+pL)/WL[0];
	HR=(WR[4]+pR)/WR[0];

	raizrhoL=sqrt(WL[0]);
	raizrhoR=sqrt(WR[0]);
	sumRaizRho=raizrhoR+raizrhoL;

	/**Hat variables (Roe averages)**/
	u_hat=(uR*raizrhoR+uL*raizrhoL)/sumRaizRho;
	v_hat=(vR*raizrhoR+vL*raizrhoL)/sumRaizRho;
	w_hat=(wR*raizrhoR+wL*raizrhoL)/sumRaizRho;
	H_hat=(HR*raizrhoR+HL*raizrhoL)/sumRaizRho;
#if MULTICOMPONENT
	#if MULTI_TYPE==1
		gamma_hat=1.0+ 1.0/( (phiR*raizrhoR+phiL*raizrhoL)/sumRaizRho );
	#else
		gamma_hat=(gammaR*raizrhoR+gammaL*raizrhoL)/sumRaizRho;
	#endif
#else
	gamma_hat=_gamma_;
#endif

	c_hat=sqrt((gamma_hat-1)*(H_hat-0.5*(u_hat*u_hat+v_hat*v_hat+w_hat*w_hat)));

	/**Physical flux calculation (the F of the eqs.)**/

	FR[0]=WR[1];
	FL[0]=WL[1];

	FR[1]=WR[1]*uR+pR;
	FL[1]=WL[1]*uL+pL;

	FR[2]=WR[1]*vR;
	FL[2]=WL[1]*vL;

      FR[3]=WR[1]*wR;
	FL[3]=WL[1]*wL;

	FR[4]=uR*(WR[4]+pR);
	FL[4]=uL*(WL[4]+pL);


	/**Wave speed estimation**/

	//S1=MIN(uL-cL,u_hat-c_hat); 
	//S2=MAX(uR+cR, u_hat+c_hat);
	S1=u_hat-c_hat;
	S2=u_hat+c_hat;

	maxS=MAX(ABS(S1),ABS(S2));
	diffS=S2-S1;


	/**Source term**/

	pRe   =wall->pRe; //(gamma_hat-1.0)*wall->URe[4];
	pLe   =wall->pLe; //(gamma_hat-1.0)*wall->ULe[4];
	rhoRe = wall->URe[0];
	rhoLe = wall->ULe[0];


	S[0]=0.0;
      if(ABS(wall->nz)>TOL14){ //this is neccessary when considering atm pressures of 1.e5, to keep precision
            S[1]=(WR[0]+WL[0])*(pRe-pLe)/(rhoRe+rhoLe);
      }else{
            S[1]=0.0;
      }
	S[2]=0.0;
	S[3]=0.0;
	S[4]=S[1]*u_hat;

	psi=(rhoRe-rhoLe)*c_hat*c_hat/(pRe-pLe+TOL14);
	chi=0.5*(psi-1.0)*(v_hat*v_hat+w_hat*w_hat);

	B[0]=-psi*S[1]/(S1*S2);
	B[1]=0.0;
	B[2]=-psi*v_hat/(S1*S2)*S[1];
	B[3]=-psi*w_hat/(S1*S2)*S[1];
	B[4]=-(H_hat-u_hat*u_hat+chi)/(S1*S2)*S[1];

	/**HLLS flux calculation**/

	for(m=0;m<5;m++){
		if(S1>=0){
			F_star[m]=FL[m];
		}else if(S2<=0){
			F_star[m]=FR[m]-S[m];
		}else{
                  F_star[m]=(S2*FL[m]-S1*FR[m]+S1*S2*(WR[m]-WL[m])+S1*(S[m]-S2*B[m]))/(diffS);
            }
	}
	/**Inverse rotation of the flux**/
	wall->fL_star[0]=F_star[0]; //Mass is not vectorial
	wall->fL_star[1]=F_star[1]*wall->nx - F_star[2]*wall->ny - F_star[3]*wall->nz;
	wall->fL_star[2]=F_star[1]*wall->ny + F_star[2]*wall->nx + F_star[2]*wall->nz;
      wall->fL_star[3]=F_star[3]*wall->nx + F_star[3]*wall->ny + F_star[1]*wall->nz;
	wall->fL_star[4]=F_star[4]; //Energy is not vectorial

	//fR_star
	for(m=0;m<5;m++){
		if(S1>=0){
			F_star[m]=FL[m]+S[m];
		}else if(S2<=0){
			F_star[m]=FR[m];
		}else{
			F_star[m]=(S2*FL[m]-S1*FR[m]+S1*S2*(WR[m]-WL[m])+S2*(S[m]-S1*B[m]))/(diffS);
		}
	}
	/**Inverse rotation of the flux**/
	wall->fR_star[0]=F_star[0]; //Mass is not vectorial
	wall->fR_star[1]=F_star[1]*wall->nx - F_star[2]*wall->ny - F_star[3]*wall->nz;
	wall->fR_star[2]=F_star[1]*wall->ny + F_star[2]*wall->nx + F_star[2]*wall->nz;
	wall->fR_star[3]=F_star[3]*wall->nx + F_star[3]*wall->ny + F_star[1]*wall->nz;
	wall->fR_star[4]=F_star[4]; //Energy is not vectorial

	/*if(ABS(S[1])>TOL8 ){
	printf(" nz: %14.14e\n",wall->nz);
	printf(" S1,S4: %14.14e %14.14e\n",S[1],S[4]);
      printf(" rho: %14.14e %14.14e %14.14e %14.14e\n",WR[0],WL[0],rhoRe,rhoLe);
	printf(" FL,FR: %14.14e %14.14e\n",FL[1],FR[1]);
	printf(" PLe,PRe: %14.14e %14.14e\n",pLe,pRe);
	printf(" PL,PR: %14.14e %14.14e\n",pL,pR);

	printf(" UR,UL: %14.14e %14.14e\n",wall->UR[1]/wall->UR[0],wall->UL[1]/wall->UL[0]);

	printf(" EL,ER: %14.14e %14.14e\n",WR[4],WL[4]);
	printf(" ELe,ERe: %14.14e %14.14e\n",wall->URe[4],wall->ULe[4]);
	printf(" rho/rhoE: %14.14e\n",(WR[0]+WL[0])/(rhoRe+rhoLe));
	printf(" deltaF, deltaP, delta %14.14e %14.14e  %14.14e \n",FR[1]-FL[1],pRe-pLe,FR[1]-FL[1]-(pRe-pLe) );

	printf(" FL: %14.14e %14.14e\n",wall->fL_star[3],FL[1]);
	printf(" FR: %14.14e %14.14e\n",wall->fR_star[3],FR[1]);
	printf(" dF: %14.14e %14.14e\n",(S2*FL[1]-S1*FR[1]+S2*S[1])-(S2-S1)*FR[1],   FR[1]-FL[1]-S[1]);
	printf(" dF (rho): %14.14e %14.14e\n",wall->fL_star[0]-FL[0],wall->fR_star[0]-FR[0]);
	printf(" Df-S: %14.14e\n",wall->fR_star[3]-wall->fL_star[3]-S[1]);
	printf(" Df-S: %14.14e\n",wall->fR_star[0]-wall->fL_star[0]);
	printf(" Df-S (rho): %14.14e\n",FR[0]-FL[0]-S[0]);
	printf(" Df-S (rho*u): %14.14e\n",FR[1]-FL[1]-S[1]);
	printf(" Df-S (rho*v): %14.14e\n",FR[2]-FL[2]-S[2]);
	printf(" Df-S (rho*w): %14.14e\n",FR[3]-FL[3]-S[3]);
	printf(" Df-S (E): %14.14e\n",FR[4]-FL[4]-S[4]);
	printf(" Du-B (rho): %14.14e %14.14e %14.14e %14.14e\n",WR[0],WL[0],B[0],WR[0]-WL[0]-B[0]);
	printf(" Du-B (rho*u): %14.14e \n",WR[1]-WL[1]-B[1]);
	printf(" Du-B (rho*v): %14.14e \n",WR[2]-WL[2]-B[2]);
	printf(" Du-B (rho*w): %14.14e \n",WR[3]-WL[3]-B[3]);
	printf(" Du-B (E): %14.14e \n",WR[4]-WL[4]-B[4]);
	getchar();
      }*/


	*lambda_max=MAX(*lambda_max,maxS);

}

void compute_transmissive_euler(t_wall *wall, int wp){

	int m;
	double WR[6], WL[6]; /**<Auxiliar array of variables rotated for the 1D problem**/
	double uL, uR, vL, vR, wL, wR, pL, pR, gammaL, gammaR;
	double FR[5], FL[5];
	double F_star[5];
#if MULTICOMPONENT	
	double phiL, phiR;
#endif


	/**Rotation of array of variables**/

	WR[0]=wall->UR[0];
	WL[0]=wall->UL[0];

      // This is a simplification of the rotation matrix, only valid for cartesian mesh
	WR[1]=wall->UR[1]*wall->nx+wall->UR[2]*wall->ny+wall->UR[3]*wall->nz;
	WL[1]=wall->UL[1]*wall->nx+wall->UL[2]*wall->ny+wall->UL[3]*wall->nz;

	WR[2]=-wall->UR[1]*wall->ny+wall->UR[2]*wall->nx+wall->UR[2]*wall->nz;
	WL[2]=-wall->UL[1]*wall->ny+wall->UL[2]*wall->nx+wall->UL[2]*wall->nz;

      WR[3]=wall->UR[3]*wall->nx+wall->UR[3]*wall->ny-wall->UR[1]*wall->nz;
	WL[3]=wall->UL[3]*wall->nx+wall->UL[3]*wall->ny-wall->UL[1]*wall->nz;


	WR[4]=wall->UR[4];
	WL[4]=wall->UL[4];

#if MULTICOMPONENT
	WR[5]=wall->UR[5];
	WL[5]=wall->UL[5];
	phiL=WL[5]/WL[0];
	phiR=WR[5]/WR[0];
	#if MULTI_TYPE==1
		gammaL=phiL;
		gammaR=phiR;
	#else
	gammaL=1.0+1.0/phiL;
	gammaR=1.0+1.0/phiR;
	#endif
#else
	gammaL=_gamma_;
	gammaR=_gamma_;
#endif


	/**Additional variables for the solver**/
	uL=WL[1]/WL[0];
	uR=WR[1]/WR[0];

	vL=WL[2]/WL[0];
	vR=WR[2]/WR[0];

      wL=WL[3]/WL[0];
	wR=WR[3]/WR[0];

	pL=(gammaL-1.0)*(WL[4]-0.5*WL[0]*(uL*uL+vL*vL+wL*wL));
	pR=(gammaR-1.0)*(WR[4]-0.5*WR[0]*(uR*uR+vR*vR+wR*wR));



	/**Physical flux calculation (the F of the eqs.)**/

      FR[0]=WR[1];
	FL[0]=WL[1];

//#if ST==99
//      pRe  =wall->pRe;
//      pLe  =wall->pLe;
//	FR[1]=WR[1]*uR+(pR-pRe);
//	FL[1]=WL[1]*uL+(pL-pLe);
//#else
	FR[1]=WR[1]*uR+pR;
	FL[1]=WL[1]*uL+pL;
//#endif

	FR[2]=WR[1]*vR;
	FL[2]=WL[1]*vL;

      FR[3]=WR[1]*wR;
	FL[3]=WL[1]*wL;

	FR[4]=uR*(WR[4]+pR);
	FL[4]=uL*(WL[4]+pL);


	/**the star flux is equal to the physical flux**/
	for(m=0;m<5;m++){
		if (wp==1 || wp==4 || wp==5){ //bottom and left interfaces, the innercell flux is FR
			F_star[m]=FR[m];
		}else{
			F_star[m]=FL[m];
		}
	}

      /**Inverse rotation of the flux, only valid for cartesian mesh**/
	wall->fR_star[0]=F_star[0]; //Mass is not vectorial
	wall->fR_star[1]=F_star[1]*wall->nx - F_star[2]*wall->ny - F_star[3]*wall->nz;
	wall->fR_star[2]=F_star[1]*wall->ny + F_star[2]*wall->nx + F_star[2]*wall->nz;
      wall->fR_star[3]=F_star[3]*wall->nx + F_star[3]*wall->ny + F_star[1]*wall->nz;
	wall->fR_star[4]=F_star[4]; //Energy is not vectorial

      for(m=0;m<5;m++){
            wall->fL_star[m]=wall->fR_star[m];
      }



}

void compute_solid_euler_hlle(t_wall *wall, double *lambda_max, int wp){

	int m;
	double WR[5], WL[5]; /**<Auxiliar array of variables rotated for the 1D problem**/
	double uL, uR, vL, vR, wL, wR, pL, pR, HL, HR, cL, cR;
#if ST==2||ST==3 //fluctuation version for source terms	
	double pLe, pRe, rhoprimeL, rhoprimeR, EprimeL, EprimeR;
#endif
	double raizrhoR, raizrhoL, sumRaizRho;
	double u_hat, v_hat, w_hat, H_hat, c_hat;
	double S1, S2, diffS, maxS;
	double FR[5], FL[5];
	double F_star[5];

	/**Rotation of array of variables**/

	WR[0]=wall->UR[0];
	WL[0]=wall->UL[0];

      // This is a simplification of the rotation matrix, only valid for cartesian mesh
	WR[1]=wall->UR[1]*wall->nx+wall->UR[2]*wall->ny+wall->UR[3]*wall->nz;
	WL[1]=wall->UL[1]*wall->nx+wall->UL[2]*wall->ny+wall->UL[3]*wall->nz;

	WR[2]=-wall->UR[1]*wall->ny+wall->UR[2]*wall->nx+wall->UR[2]*wall->nz;
	WL[2]=-wall->UL[1]*wall->ny+wall->UL[2]*wall->nx+wall->UL[2]*wall->nz;

      WR[3]=wall->UR[3]*wall->nx+wall->UR[3]*wall->ny-wall->UR[1]*wall->nz;
	WL[3]=wall->UL[3]*wall->nx+wall->UL[3]*wall->ny-wall->UL[1]*wall->nz;



	WR[4]=wall->UR[4];
	WL[4]=wall->UL[4];


      for(m=0;m<5;m++){
		if (wp==1 || wp==4 || wp==5){ //bottom and left interfaces, the innercell state is WR, otherwise is WL
                  if (m==1){
                        WL[m]=-WR[m];
                  }else{
                        WL[m]= WR[m];
                  }
		}else{
                  if (m==1){
                        WR[m]=-WL[m];
                  }else{
                        WR[m]= WL[m];
                  }
		}
	}

	#if ST==2||ST==3
	if (wp==1 || wp==4 || wp==5){ //bottom and left interfaces, the innercell state is WR, otherwise is WL
		pRe = wall->pRe;
		pLe = pRe;
		rhoprimeR=WR[0]-wall->URe[0];
		rhoprimeL=rhoprimeR;
		EprimeR=WR[4]-wall->URe[4];
		EprimeL=EprimeR;
	}else{
		pLe = wall->pLe;
		pRe = pLe;
		EprimeL=WL[4]-wall->ULe[4];
		EprimeR=EprimeL;
		rhoprimeL=WL[0]-wall->ULe[0];
		rhoprimeR=rhoprimeL;
	}



      #endif


	/**Additional variables for the solver**/
	uL=WL[1]/WL[0];
	uR=WR[1]/WR[0];

	vL=WL[2]/WL[0];
	vR=WR[2]/WR[0];

      wL=WL[3]/WL[0];
	wR=WR[3]/WR[0];

	pL=pressure_from_energy(_gamma_, WL[4], uL, vL, wL, WL[0], wall->z);
	pR=pressure_from_energy(_gamma_, WR[4], uR, vR, wR, WR[0], wall->z);

#if ST==3
	HL=(WL[4]-WL[0]*_g_*wall->z+pL)/WL[0];
	HR=(WR[4]-WR[0]*_g_*wall->z+pR)/WR[0];
#else
	HL=(WL[4]+pL)/WL[0];
	HR=(WR[4]+pR)/WR[0];
#endif


	cL=sqrt(_gamma_*pL/WL[0]);
	cR=sqrt(_gamma_*pR/WR[0]);


	raizrhoL=sqrt(WL[0]);
	raizrhoR=sqrt(WR[0]);
	sumRaizRho=raizrhoR+raizrhoL;

	/**Hat variables (Roe averages)**/
	u_hat=(uR*raizrhoR+uL*raizrhoL)/sumRaizRho;
	v_hat=(vR*raizrhoR+vL*raizrhoL)/sumRaizRho;
      w_hat=(wR*raizrhoR+wL*raizrhoL)/sumRaizRho;
	H_hat=(HR*raizrhoR+HL*raizrhoL)/sumRaizRho;

	c_hat=sqrt((_gamma_-1)*(H_hat-0.5*(u_hat*u_hat+v_hat*v_hat+w_hat*w_hat)));

	/**Physical flux calculation (the F of the eqs.)**/

	FR[0]=WR[1];
	FL[0]=WL[1];

#if ST==2||ST==3
	FR[1]=WR[1]*uR+(pR-pRe);
	FL[1]=WL[1]*uL+(pL-pLe);
#else
	FR[1]=WR[1]*uR+pR;
	FL[1]=WL[1]*uL+pL;
#endif

	FR[2]=WR[1]*vR;
	FL[2]=WL[1]*vL;

      FR[3]=WR[1]*wR;
	FL[3]=WL[1]*wL;

	FR[4]=uR*(WR[4]+pR);
	FL[4]=uL*(WL[4]+pL);

#if ST==2||ST==3
	WL[0]=rhoprimeL;
	WR[0]=rhoprimeR;
	WL[4]=EprimeL;
	WR[4]=EprimeR;
#endif




	/**Wave speed estimation**/

	S1=MIN(uL-cL,u_hat-c_hat);
	S2=MAX(uR+cR, u_hat+c_hat);


	maxS=MAX(ABS(S1),ABS(S2));
	diffS=S2-S1;

	/**HLLE flux calculation**/
	for(m=0;m<5;m++){
		if(S1>=0){
			F_star[m]=FL[m];
		}else if(S2<=0){
			F_star[m]=FR[m];
		}else{
                  F_star[m]=(S2*FL[m]-S1*FR[m]+S1*S2*(WR[m]-WL[m]))/(diffS);
            }
	}


	/**Inverse rotation of the flux**/
	wall->fR_star[0]=F_star[0]; //Mass is not vectorial
	wall->fR_star[1]=F_star[1]*wall->nx - F_star[2]*wall->ny - F_star[3]*wall->nz;
	wall->fR_star[2]=F_star[1]*wall->ny + F_star[2]*wall->nx + F_star[2]*wall->nz;
      wall->fR_star[3]=F_star[3]*wall->nx + F_star[3]*wall->ny + F_star[1]*wall->nz;
	wall->fR_star[4]=F_star[4]; //Energy is not vectorial
      //wall->fR_star[5]=0.0; //no solute flux

      for(m=0;m<5;m++){
            wall->fL_star[m]=wall->fR_star[m];
      }


	*lambda_max=MAX(*lambda_max,maxS);

}

void compute_euler_Roe(t_wall *wall,double *lambda_max){

	printf("%s We are working on it. Sorry!\n",WAR);
	abort();
}


void compute_burgers_flux(t_wall *wall,double *lambda_max){

	double Savg;
	double fL,fR;
	double dU;

	fL=wall->UL[0]*wall->UL[0]/2.0;
	fR=wall->UR[0]*wall->UR[0]/2.0;

	Savg=(wall->UL[0]+wall->UR[0])*0.5;
	dU=wall->UR[0]-wall->UL[0];

	wall->fR_star[0]=0.5*(fL+fR - ABS(Savg)*dU);
	wall->fL_star[0]=wall->fR_star[0];

	*lambda_max=MAX(*lambda_max,Savg);

}

void compute_linear_flux(t_wall *wall,double *lambda_max){

	double Savg;
	double fL,fR;
	double dU;

	Savg=wall->vel;

	fL=wall->UL[0]*Savg;
	fR=wall->UR[0]*Savg;

	dU=wall->UR[0]-wall->UL[0];

	wall->fR_star[0]=0.5*(fL+fR - ABS(Savg)*dU);
	wall->fL_star[0]=wall->fR_star[0];

	*lambda_max=MAX(*lambda_max,Savg);

}


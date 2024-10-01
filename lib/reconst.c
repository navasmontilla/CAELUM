/*

Authors:
 - Adrián Navas Montilla
 - Isabel Echeverribar

Copyright (C) 2019-2024 The authors.

File:
  - reconst.c

Content:
  -This file contains all the functions implementing the reconstruction methods for high-order schemes.
  
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
#include "reconst.h"


double weno3R(double *phi){

	double g0, g1;		//gamma optimal weight/	
	double w0, w1;		//WENO weight
#if TYPE_REC < 2
	double b0, b1;		//beta
	double a0, a1;		//alpha
#if TYPE_REC == 1
	double c0, c1;
#endif
#endif
	double UR;

	g0=2.0/3.0;
	g1=1.0/3.0;

#if TYPE_REC == 0

	b0=(phi[1]-phi[0])*(phi[1]-phi[0]);
	b1=(phi[2]-phi[1])*(phi[2]-phi[1]);

	a0=g0/((b0+epsilon)*(b0+epsilon));
	a1=g1/((b1+epsilon)*(b1+epsilon));

	w0=a0/(a0+a1);
	w1=a1/(a0+a1);

#elif TYPE_REC == 1  //TENO

	b0=(phi[1]-phi[0])*(phi[1]-phi[0]);
	b1=(phi[2]-phi[1])*(phi[2]-phi[1]);

	a0=1.0/pow((b0+epsilon2),_Q_);
	a1=1.0/pow((b1+epsilon2),_Q_);

	c0 = a0/(a0 + a1);
      c1 = a1/(a0 + a1);

      c0 = c0 < _CT_ ? 0. : 1.;
      c1 = c1 < _CT_ ? 0. : 1.;

	a0 = g0*c0;
      a1 = g1*c1;

      w0 = a0/(a0 + a1);
      w1 = a1/(a0 + a1);

#else

    w0=g0;
	w1=g1;

#endif

	UR=w0*(0.5*phi[1]+0.5*phi[0])+w1*(-0.5*phi[2]+1.5*phi[1]);


	return UR;
}

double weno3L(double *phi){

	double g0, g1;		//gamma optimal weight/	
	double w0, w1;		//WENO weight
#if TYPE_REC < 2
	double b0, b1;		//beta
	double a0, a1;		//alpha
#if TYPE_REC == 1
	double c0, c1;
#endif
#endif
	double UL;

	g0=1.0/3.0;
	g1=2.0/3.0;

#if TYPE_REC == 0

	b0=(phi[1]-phi[0])*(phi[1]-phi[0]);
	b1=(phi[2]-phi[1])*(phi[2]-phi[1]);

	a0=g0/((b0+epsilon)*(b0+epsilon));
	a1=g1/((b1+epsilon)*(b1+epsilon));

	w0=a0/(a0+a1);
	w1=a1/(a0+a1);

#elif TYPE_REC == 1  //TENO

	b0=(phi[1]-phi[0])*(phi[1]-phi[0]);
	b1=(phi[2]-phi[1])*(phi[2]-phi[1]);

	a0=1.0/pow((b0+epsilon2),_Q_);
	a1=1.0/pow((b1+epsilon2),_Q_);

	c0 = a0/(a0 + a1);
      c1 = a1/(a0 + a1);

      c0 = c0 < _CT_ ? 0. : 1.;
      c1 = c1 < _CT_ ? 0. : 1.;

	a0 = g0*c0;
      a1 = g1*c1;

      w0 = a0/(a0 + a1);
      w1 = a1/(a0 + a1);

#else //UWC

    w0=g0;
	w1=g1;

#endif

	UL=w0*(-0.5*phi[0]+1.5*phi[1])+w1*(0.5*phi[1]+0.5*phi[2]);


	return UL;
}

double weno5R(double *phi){

	double g0, g1, g2;		//gamma optimal weight/
	double w0, w1, w2;		//WENO weight
#if TYPE_REC < 2
	double b0, b1, b2;		//beta
	double a0, a1, a2;		//alpha
#if TYPE_REC == 1
	double c0, c1, c2;
#endif
#endif
	double UR;

	g0=3.0/10.0;
	g1=3.0/5.0;
      g2=1.0/10.0;

#if TYPE_REC == 0    //WENO

	b0=13.0/12.0*(phi[0]-2*phi[1]+phi[2])*(phi[0]-2*phi[1]+phi[2])+0.25*(phi[0]-4*phi[1]+3*phi[2])*(phi[0]-4*phi[1]+3*phi[2]);
	b1=13.0/12.0*(phi[1]-2*phi[2]+phi[3])*(phi[1]-2*phi[2]+phi[3])+0.25*(phi[1]-phi[3])*(phi[1]-phi[3]);
      b2=13.0/12.0*(phi[2]-2*phi[3]+phi[4])*(phi[2]-2*phi[3]+phi[4])+0.25*(3*phi[2]-4*phi[3]+phi[4])*(3*phi[2]-4*phi[3]+phi[4]);


      a0=g0/((b0+epsilon)*(b0+epsilon));
	a1=g1/((b1+epsilon)*(b1+epsilon));
      a2=g2/((b2+epsilon)*(b2+epsilon));


	w0=a0/(a0 + a1 + a2);
	w1=a1/(a0 + a1 + a2);
      w2=a2/(a0 + a1 + a2);

#elif TYPE_REC == 1  //TENO

	b0=13.0/12.0*(phi[0]-2*phi[1]+phi[2])*(phi[0]-2*phi[1]+phi[2])+0.25*(phi[0]-4*phi[1]+3*phi[2])*(phi[0]-4*phi[1]+3*phi[2]);
	b1=13.0/12.0*(phi[1]-2*phi[2]+phi[3])*(phi[1]-2*phi[2]+phi[3])+0.25*(phi[1]-phi[3])*(phi[1]-phi[3]);
      b2=13.0/12.0*(phi[2]-2*phi[3]+phi[4])*(phi[2]-2*phi[3]+phi[4])+0.25*(3*phi[2]-4*phi[3]+phi[4])*(3*phi[2]-4*phi[3]+phi[4]);


	a0=1.0/pow((b0+epsilon2),_Q_);
	a1=1.0/pow((b1+epsilon2),_Q_);
      a2=1.0/pow((b2+epsilon2),_Q_);

	c0 = a0/(a0 + a1 + a2);
      c1 = a1/(a0 + a1 + a2);
      c2 = a2/(a0 + a1 + a2);

      c0 = c0 < _CT_ ? 0. : 1.;
      c1 = c1 < _CT_ ? 0. : 1.;
      c2 = c2 < _CT_ ? 0. : 1.;

	a0 = g0*c0;
      a1 = g1*c1;
      a2 = g2*c2;

      w0 = a0/(a0 + a1 + a2);
      w1 = a1/(a0 + a1 + a2);
      w2 = a2/(a0 + a1 + a2);

#else //UWC

	w0=g0;
	w1=g1;
	w2=g2;


#endif

	UR=w0*(1.0/3.0*phi[2]+5.0/6.0*phi[1]-1.0/6.0*phi[0]) + w1*(-1.0/6.0*phi[3]+5.0/6.0*phi[2]+1.0/3.0*phi[1]) + w2*(1.0/3.0*phi[4]-7.0/6.0*phi[3]+11.0/6.0*phi[2]) ;


	return UR;
}

double weno5L(double *phi){

	double g0, g1, g2;		//gamma optimal weight/
	double w0, w1, w2;		//WENO weight
#if TYPE_REC < 2
	double b0, b1, b2;		//beta
	double a0, a1, a2;		//alpha
#if TYPE_REC == 1
	double c0, c1, c2;
#endif
#endif
	double UL;

	g0=1.0/10.0;
	g1=3.0/5.0;
      g2=3.0/10.0;

#if TYPE_REC == 0

	b0=13.0/12.0*(phi[0]-2*phi[1]+phi[2])*(phi[0]-2*phi[1]+phi[2])+0.25*(phi[0]-4*phi[1]+3*phi[2])*(phi[0]-4*phi[1]+3*phi[2]);
	b1=13.0/12.0*(phi[1]-2*phi[2]+phi[3])*(phi[1]-2*phi[2]+phi[3])+0.25*(phi[1]-phi[3])*(phi[1]-phi[3]);
      b2=13.0/12.0*(phi[2]-2*phi[3]+phi[4])*(phi[2]-2*phi[3]+phi[4])+0.25*(3*phi[2]-4*phi[3]+phi[4])*(3*phi[2]-4*phi[3]+phi[4]);


      a0=g0/((b0+epsilon)*(b0+epsilon));
	a1=g1/((b1+epsilon)*(b1+epsilon));
      a2=g2/((b2+epsilon)*(b2+epsilon));


	w0=a0/(a0 + a1 + a2);
	w1=a1/(a0 + a1 + a2);
      w2=a2/(a0 + a1 + a2);

#elif TYPE_REC == 1  //TENO

	b0=13.0/12.0*(phi[0]-2*phi[1]+phi[2])*(phi[0]-2*phi[1]+phi[2])+0.25*(phi[0]-4*phi[1]+3*phi[2])*(phi[0]-4*phi[1]+3*phi[2]);
	b1=13.0/12.0*(phi[1]-2*phi[2]+phi[3])*(phi[1]-2*phi[2]+phi[3])+0.25*(phi[1]-phi[3])*(phi[1]-phi[3]);
      b2=13.0/12.0*(phi[2]-2*phi[3]+phi[4])*(phi[2]-2*phi[3]+phi[4])+0.25*(3*phi[2]-4*phi[3]+phi[4])*(3*phi[2]-4*phi[3]+phi[4]);


	a0=1.0/pow((b0+epsilon2),_Q_);
	a1=1.0/pow((b1+epsilon2),_Q_);
      a2=1.0/pow((b2+epsilon2),_Q_);

	c0 = a0/(a0 + a1 + a2);
      c1 = a1/(a0 + a1 + a2);
      c2 = a2/(a0 + a1 + a2);

      c0 = c0 < _CT_ ? 0. : 1.;
      c1 = c1 < _CT_ ? 0. : 1.;
      c2 = c2 < _CT_ ? 0. : 1.;

	a0 = g0*c0;
      a1 = g1*c1;
      a2 = g2*c2;

      w0 = a0/(a0 + a1 + a2);
      w1 = a1/(a0 + a1 + a2);
      w2 = a2/(a0 + a1 + a2);

#else

      w0=g0;
	w1=g1;
      w2=g2;

#endif


	UL=w2*(1.0/3.0*phi[2]+5.0/6.0*phi[3]-1.0/6.0*phi[4]) + w1*(-1.0/6.0*phi[1]+5.0/6.0*phi[2]+1.0/3.0*phi[3]) + w0*(1.0/3.0*phi[0]-7.0/6.0*phi[1]+11.0/6.0*phi[2]) ;



	return UL;
}

double weno7R(double *phi){

	double g0, g1, g2, g3;		//gamma optimal weight/
	double w0, w1, w2, w3;		//WENO weight
#if TYPE_REC < 2
	double b0, b1, b2, b3;		//beta
	double a0, a1, a2, a3;		//alpha
#if TYPE_REC == 1
	double c0, c1, c2, c3;
#endif
#endif
	double UR;

	g0=4.0/35.0;
	g1=18.0/35.0;
      g2=12.0/35.0;
      g3=1.0/35.0;


#if TYPE_REC == 0

      b0 = phi[0]*(547.0*phi[0] - 3882.0*phi[1] + 4642.0*phi[2] - 1854.0*phi[3]) + phi[1]*(7043.0*phi[1] - 17246.0*phi[2] + 7042.0*phi[3]) + phi[2]*(11003.0*phi[2] - 9402.0*phi[3]) + phi[3]*2107.0*phi[3];
	b1 = phi[1]*(267.0*phi[1] - 1642.0*phi[2] + 1602.0*phi[3] - 494.0*phi[4]) + phi[2]*(2843.0*phi[2] - 5966.0*phi[3] + 1922.0*phi[4]) + phi[3]*(3443.0*phi[3] - 2522.0*phi[4]) + phi[4]*547.0*phi[4];
	b2 = phi[2]*(547.0*phi[2] - 2522.0*phi[3] + 1922.0*phi[4] - 494.0*phi[5]) + phi[3]*(3443.0*phi[3] - 5966.0*phi[4] + 1602*phi[5]) + phi[4]*(2843.0*phi[4] - 1642*phi[5]) + phi[5]*267.0*phi[5];
	b3 = phi[3]*(2107.0*phi[3] - 9402.0*phi[4] + 7042.0*phi[5] - 1854.0*phi[6]) + phi[4]*(11003.0*phi[4] - 17246.0*phi[5] + 4642.0*phi[6]) + phi[5]*(7043.0*phi[5] - 3882.0*phi[6]) + phi[6]*547.0*phi[6];

      a0=g0/((b0+epsilon)*(b0+epsilon));
	a1=g1/((b1+epsilon)*(b1+epsilon));
      a2=g2/((b2+epsilon)*(b2+epsilon));
      a3=g3/((b3+epsilon)*(b3+epsilon));


	w0=a0/(a0 + a1 + a2 + a3);
	w1=a1/(a0 + a1 + a2 + a3);
      w2=a2/(a0 + a1 + a2 + a3);
      w3=a3/(a0 + a1 + a2 + a3);

#elif TYPE_REC == 1  //TENO

      b0 = phi[0]*(547.0*phi[0] - 3882.0*phi[1] + 4642.0*phi[2] - 1854.0*phi[3]) + phi[1]*(7043.0*phi[1] - 17246.0*phi[2] + 7042.0*phi[3]) + phi[2]*(11003.0*phi[2] - 9402.0*phi[3]) + phi[3]*2107.0*phi[3];
	b1 = phi[1]*(267.0*phi[1] - 1642.0*phi[2] + 1602.0*phi[3] - 494.0*phi[4]) + phi[2]*(2843.0*phi[2] - 5966.0*phi[3] + 1922.0*phi[4]) + phi[3]*(3443.0*phi[3] - 2522.0*phi[4]) + phi[4]*547.0*phi[4];
	b2 = phi[2]*(547.0*phi[2] - 2522.0*phi[3] + 1922.0*phi[4] - 494.0*phi[5]) + phi[3]*(3443.0*phi[3] - 5966.0*phi[4] + 1602*phi[5]) + phi[4]*(2843.0*phi[4] - 1642*phi[5]) + phi[5]*267.0*phi[5];
	b3 = phi[3]*(2107.0*phi[3] - 9402.0*phi[4] + 7042.0*phi[5] - 1854.0*phi[6]) + phi[4]*(11003.0*phi[4] - 17246.0*phi[5] + 4642.0*phi[6]) + phi[5]*(7043.0*phi[5] - 3882.0*phi[6]) + phi[6]*547.0*phi[6];

	a0=1.0/pow((b0+epsilon2),_Q_);
	a1=1.0/pow((b1+epsilon2),_Q_);
      a2=1.0/pow((b2+epsilon2),_Q_);
	a3=1.0/pow((b3+epsilon2),_Q_);

	c0 = a0/(a0 + a1 + a2 + a3);
      c1 = a1/(a0 + a1 + a2 + a3);
      c2 = a2/(a0 + a1 + a2 + a3);
	c3 = a3/(a0 + a1 + a2 + a3);

      c0 = c0 < _CT_ ? 0. : 1.;
      c1 = c1 < _CT_ ? 0. : 1.;
      c2 = c2 < _CT_ ? 0. : 1.;
	c3 = c3 < _CT_ ? 0. : 1.;

	a0 = g0*c0;
      a1 = g1*c1;
      a2 = g2*c2;
	a3 = g3*c3;

      w0 = a0/(a0 + a1 + a2 + a3);
      w1 = a1/(a0 + a1 + a2 + a3);
      w2 = a2/(a0 + a1 + a2 + a3);
	w3 = a3/(a0 + a1 + a2 + a3);

#else //UWC

      w0=g0;
      w1=g1;
      w2=g2;
      w3=g3;

#endif

	UR= w0*(1.0/4.0*phi[3]  + 13.0/12.0*phi[2] - 5.0/12.0*phi[1] + 1.0/12.0*phi[0]) + w1*(-1.0/12.0*phi[4] + 7.0/12.0*phi[3] + 7.0/12.0*phi[2] - 1.0/12.0*phi[1]) + w2*(1.0/12.0*phi[5] - 5.0/12.0*phi[4] + 13.0/12.0*phi[3] + 1.0/4.0*phi[2]) + w3*(-1.0/4.0*phi[6] + 13.0/12.0*phi[5] - 23.0/12.0*phi[4] + 25.0/12.0*phi[3]);


	return UR;
}

double weno7L(double *phi){

	double g0, g1, g2, g3;		//gamma optimal weight/
	double w0, w1, w2, w3;		//WENO weight
#if TYPE_REC < 2
	double b0, b1, b2, b3;		//beta
	double a0, a1, a2, a3;		//alpha
#if TYPE_REC == 1
	double c0, c1, c2, c3;
#endif
#endif
	double UL;

	g0=1.0/35.0;
	g1=12.0/35.0;
      g2=18.0/35.0;
      g3=4.0/35.0;


#if TYPE_REC == 0

      b0 = phi[0]*(547.0*phi[0] - 3882.0*phi[1] + 4642.0*phi[2] - 1854.0*phi[3]) + phi[1]*(7043.0*phi[1] - 17246.0*phi[2] + 7042.0*phi[3]) + phi[2]*(11003.0*phi[2] - 9402.0*phi[3]) + phi[3]*2107.0*phi[3];
	b1 = phi[1]*(267.0*phi[1] - 1642.0*phi[2] + 1602.0*phi[3] - 494.0*phi[4]) + phi[2]*(2843.0*phi[2] - 5966.0*phi[3] + 1922.0*phi[4]) + phi[3]*(3443.0*phi[3] - 2522.0*phi[4]) + phi[4]*547.0*phi[4];
	b2 = phi[2]*(547.0*phi[2] - 2522.0*phi[3] + 1922.0*phi[4] - 494.0*phi[5]) + phi[3]*(3443.0*phi[3] - 5966.0*phi[4] + 1602*phi[5]) + phi[4]*(2843.0*phi[4] - 1642*phi[5]) + phi[5]*267.0*phi[5];
	b3 = phi[3]*(2107.0*phi[3] - 9402.0*phi[4] + 7042.0*phi[5] - 1854.0*phi[6]) + phi[4]*(11003.0*phi[4] - 17246.0*phi[5] + 4642.0*phi[6]) + phi[5]*(7043.0*phi[5] - 3882.0*phi[6]) + phi[6]*547.0*phi[6];

      a0=g0/((b0+epsilon)*(b0+epsilon));
	a1=g1/((b1+epsilon)*(b1+epsilon));
      a2=g2/((b2+epsilon)*(b2+epsilon));
      a3=g3/((b3+epsilon)*(b3+epsilon));

	w0=a0/(a0 + a1 + a2 + a3);
	w1=a1/(a0 + a1 + a2 + a3);
      w2=a2/(a0 + a1 + a2 + a3);
      w3=a3/(a0 + a1 + a2 + a3);

#elif TYPE_REC == 1  //TENO

      b0 = phi[0]*(547.0*phi[0] - 3882.0*phi[1] + 4642.0*phi[2] - 1854.0*phi[3]) + phi[1]*(7043.0*phi[1] - 17246.0*phi[2] + 7042.0*phi[3]) + phi[2]*(11003.0*phi[2] - 9402.0*phi[3]) + phi[3]*2107.0*phi[3];
	b1 = phi[1]*(267.0*phi[1] - 1642.0*phi[2] + 1602.0*phi[3] - 494.0*phi[4]) + phi[2]*(2843.0*phi[2] - 5966.0*phi[3] + 1922.0*phi[4]) + phi[3]*(3443.0*phi[3] - 2522.0*phi[4]) + phi[4]*547.0*phi[4];
	b2 = phi[2]*(547.0*phi[2] - 2522.0*phi[3] + 1922.0*phi[4] - 494.0*phi[5]) + phi[3]*(3443.0*phi[3] - 5966.0*phi[4] + 1602*phi[5]) + phi[4]*(2843.0*phi[4] - 1642*phi[5]) + phi[5]*267.0*phi[5];
	b3 = phi[3]*(2107.0*phi[3] - 9402.0*phi[4] + 7042.0*phi[5] - 1854.0*phi[6]) + phi[4]*(11003.0*phi[4] - 17246.0*phi[5] + 4642.0*phi[6]) + phi[5]*(7043.0*phi[5] - 3882.0*phi[6]) + phi[6]*547.0*phi[6];

	a0=1.0/pow((b0+epsilon2),_Q_);
	a1=1.0/pow((b1+epsilon2),_Q_);
      a2=1.0/pow((b2+epsilon2),_Q_);
	a3=1.0/pow((b3+epsilon2),_Q_);

	c0 = a0/(a0 + a1 + a2 + a3);
      c1 = a1/(a0 + a1 + a2 + a3);
      c2 = a2/(a0 + a1 + a2 + a3);
	c3 = a3/(a0 + a1 + a2 + a3);

      c0 = c0 < _CT_ ? 0. : 1.;
      c1 = c1 < _CT_ ? 0. : 1.;
      c2 = c2 < _CT_ ? 0. : 1.;
	c3 = c3 < _CT_ ? 0. : 1.;

	a0 = g0*c0;
      a1 = g1*c1;
      a2 = g2*c2;
	a3 = g3*c3;

      w0 = a0/(a0 + a1 + a2 + a3);
      w1 = a1/(a0 + a1 + a2 + a3);
      w2 = a2/(a0 + a1 + a2 + a3);
	w3 = a3/(a0 + a1 + a2 + a3);

#else //UWC

      w0=g0;
      w1=g1;
      w2=g2;
      w3=g3;

#endif

      UL = w0*(-1.0/4.0*phi[0] + 13.0/12.0*phi[1] - 23.0/12.0*phi[2] + 25.0/12.0*phi[3]) + w1*(1.0/12.0*phi[1] - 5.0/12.0*phi[2] + 13.0/12.0*phi[3] + 1.0/4.0*phi[4]) + w2*(-1.0/12.0*phi[2] + 7.0/12.0*phi[3]  + 7.0/12.0*phi[4] - 1.0/12.0*phi[5]) + w3*(1.0/4.0*phi[3] + 13.0/12.0*phi[4] - 5.0/12.0*phi[5] + 1.0/12.0*phi[6]);


	return UL;
}


/***************************************************************************
 * Function	: SegmentStress
 * Description	: Stress field of a finite straight dislocation segment
 * 3/11/04 Wei Cai
 * Notes:  regularization constant "a" is used, stress field never singular
 *         refer to matlab code "calsegstr3da.m"
***************************************************************************/
#include <stdio.h>
#include <math.h>
#include "SegmentStress.h"

#define FFACTOR_LMIN        1e-8  /* minimum allowed segment length */
#define FFACTOR_LMIN2       1e-16 /* minimum allowed segment length squared */

/*-------------------------------------------------------------------------
 *
 *      Function:     xvector
 *      Description:  Calculates the cross product of two vectors.
 *
 *      Arguments:
 *          ax, ay, az  Components of first vector
 *          bx, by, bz  Components of second vector
 *          cx, cy, cz  Pointers to locations in which to return to
 *                      the caller the components of the cross product
 *                      of the two vectors.
 *
 *------------------------------------------------------------------------*/
void xvector(real8 ax, real8 ay, real8 az,
             real8 bx, real8 by, real8 bz,
             real8 *cx, real8 *cy, real8 *cz)
{
        *cx =  ay*bz - az*by;
        *cy =  az*bx - ax*bz;
        *cz =  ax*by - ay*bx;

        return;
}

void xunitvector(real8 ax, real8 ay, real8 az,
                        real8 bx, real8 by, real8 bz,
                        real8 *cx, real8 *cy, real8 *cz)
{
     real8 L, dx, dy, dz;
     xvector(ax,ay,az, bx,by,bz, cx,cy,cz);
     L=(*cx)*(*cx)+(*cy)*(*cy)+(*cz)*(*cz);
     if(L<FFACTOR_LMIN2)
     {
         dx=bx+1; dy=by; dz=bz;
         xvector(ax,ay,az, dx,dy,dz, cx,cy,cz);
         L=(*cx)*(*cx)+(*cy)*(*cy)+(*cz)*(*cz);
         if(L<FFACTOR_LMIN2)
         {
             dx=bx; dy=by+1; dz=bz;
             xvector(ax,ay,az, dx,dy,dz, cx,cy,cz);
             L=(*cx)*(*cx)+(*cy)*(*cy)+(*cz)*(*cz);
         }
     }
     L=sqrt(L);
     *cx/=L; *cy/=L; *cz/=L;
}

void Init3x3(double A[3][3])
{
  int i, j;
  for (i = 0; i < 3; i++)
    for (j = 0; j < 3; j++)
      {
        A[i][j] = 0.0;
      }
}

/*-------------------------------------------------------------------------
 *
 *      Function:     Matrix33Mult33
 *      Description:  Multiplies two 3X3 matrices
 *
 *      Arguments:
 *          a    3 X 3 array containing components of the first matrix
 *          b    3 X 3 array containing components of the second matrix
 *          c    3 X 3 array in which to return to the caller the
 *               resulting matrix after multiplying matrix <a> by
 *               matrix <b>.
 *
 *------------------------------------------------------------------------*/
void Matrix33Mult33(real8 a[3][3], real8 b[3][3], real8 c[3][3])
{
        int i, j, k;

        for(i = 0; i < 3; i++) {
            for(j = 0; j < 3; j++) {
                c[i][j] = 0.0;
                for(k = 0; k < 3; k++) {
                    c[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        return;
}

void SegmentStressHor(real8 MU, real8 NU,
                      real8 bx, real8 by, real8 bz,
                      real8 z1, real8 z2,
                      real8 x, real8 z,
                      real8 a,
                      real8 Sigma[3][3] )
{/* segment lies horizontally (0,0,z1) to (0,0,z2)
    field point (x,0,z)  (i.e. y=0)
    Burgers vector (bx,by,bz)
    shear modulus MU, Poisson's ratio NU,
    regularization radius a (dislocation core width),
    return stress: Sigma[3][3]
  */
    real8 s[6];
    real8 L, Lp, L2, Lp2, Ra, Ra2, Ra3, Rap, Rap2, Rap3;
    real8 x2, a2, sunit, invral,invrapl, rhoa2;
    int form;
/*
    printf("%20.10e %20.10e %20.10e %20.10e %20.10e %20.10e %20.10e\n",
            bx,by,bz,z1,z2,x,z);
*/  
    L=z1-z; Lp=z2-z;
    L2=L*L; Lp2=Lp*Lp; x2=x*x; a2=a*a;
    
    Ra2=x2+L2+a2;   Ra=sqrt(Ra2);   Ra3=Ra*Ra2;
    Rap2=x2+Lp2+a2; Rap=sqrt(Rap2); Rap3=Rap*Rap2;

    
    sunit=MU*(1.0/4.0/M_PI/(1-NU));

    if((L>0)&&(Lp>0))
        form=1;
    else if((L<0)&&(Lp<0))
        form=2;
    else
        form=3;

    /* for debug purposes (form 1,2,3 should yield same result) */
    /* form=3; */
    
    switch(form) {
    case(1):
        invral=1/Ra/(Ra+L);
        invrapl=1/Rap/(Rap+Lp);        
        s[0]= by*x*( invrapl*(1-(x2+a2)/Rap2-(x2+a2)*invrapl)
                   - invral *(1-(x2+a2)/Ra2 -(x2+a2)*invral ) );
        s[1]=-bx*x*( invrapl - invral );
        s[2]= by*( (-NU/Rap+(x2+(1-NU)*a2/2.0)/Rap3)
                  -(-NU/Ra +(x2+(1-NU)*a2/2.0)/Ra3 ) );
        s[3]=-by*x*( invrapl*(1+a2/Rap2+a2*invrapl)
                    -invral *(1+a2/Ra2 +a2*invral ) );
        s[4]= (  bx*(NU/Rap-(1-NU)*(a2/2.0)/Rap3)
                -bz*x*(1-NU)*invrapl*(1+(a2/2.0)/Rap2+(a2/2.0)*invrapl) )
             -(  bx*(NU/Ra -(1-NU)*(a2/2.0)/Ra3 )
                -bz*x*(1-NU)*invral *(1+(a2/2.0)/Ra2 +(a2/2.0)*invral ) );
        s[5]= by*x*( (-2*NU*invrapl*(1+(a2/2.0)/Rap2+(a2/2.0)*invrapl)-Lp/Rap3)
                    -(-2*NU*invral *(1+(a2/2.0)/Ra2 +(a2/2.0)*invral )-L /Ra3 ) );
        break;
    case(2):
        invral=1/Ra/(Ra-L);
        invrapl=1/Rap/(Rap-Lp);        
        s[0]=-by*x*( invrapl*(1-(x2+a2)/Rap2-(x2+a2)*invrapl)
                   - invral *(1-(x2+a2)/Ra2 -(x2+a2)*invral ) );
        s[1]= bx*x*( invrapl - invral );
        s[2]= by*( (-NU/Rap+(x2+(1-NU)*a2/2.0)/Rap3)
                  -(-NU/Ra +(x2+(1-NU)*a2/2.0)/Ra3 ) );
        s[3]= by*x*( invrapl*(1+a2/Rap2+a2*invrapl)
                    -invral *(1+a2/Ra2 +a2*invral ) );
        s[4]= (  bx*(NU/Rap-(1-NU)*(a2/2.0)/Rap3)
                +bz*x*(1-NU)*invrapl*(1+(a2/2.0)/Rap2+(a2/2.0)*invrapl) )
             -(  bx*(NU/Ra -(1-NU)*(a2/2.0)/Ra3 )
                +bz*x*(1-NU)*invral *(1+(a2/2.0)/Ra2 +(a2/2.0)*invral ) );
        s[5]= by*x*( ( 2*NU*invrapl*(1+(a2/2.0)/Rap2+(a2/2.0)*invrapl)-Lp/Rap3)
                    -( 2*NU*invral *(1+(a2/2.0)/Ra2 +(a2/2.0)*invral )-L /Ra3 ) );
        break;
    case(3):
        rhoa2=x2+a2;
        s[0]= by*x/rhoa2*( Lp/Rap*(1+rhoa2/Rap2)
                          -L /Ra *(1+rhoa2/Ra2 ) );
        s[1]= bx*x/rhoa2*( Lp/Rap - L/Ra );
        s[2]= by*( (-NU/Rap+x2/Rap3+(1-NU)*(a2/2.0)/Rap3)
                  -(-NU/Ra +x2/Ra3 +(1-NU)*(a2/2.0)/Ra3 ) );
        s[3]= by*x/rhoa2*( Lp/Rap*(1+(a2*2.0)/rhoa2+a2/Rap2)
                          -L /Ra *(1+(a2*2.0)/rhoa2+a2/Ra2 ) );
        s[4]= (  bx*(NU/Rap-(1-NU)*(a2/2.0)/Rap3)
                +bz*x/rhoa2*(1-NU)*Lp/Rap*(1+a2/rhoa2+(a2/2.0)/Rap2) )
             -(  bx*(NU/Ra -(1-NU)*(a2/2.0)/Ra3 )
                +bz*x/rhoa2*(1-NU)*L /Ra *(1+a2/rhoa2+(a2/2.0)/Ra2 ) );
        s[5]= by*x*( ( 2*NU*Lp/rhoa2/Rap*(1+a2/rhoa2+(a2/2.0)/Rap2)-Lp/Rap3)
                    -( 2*NU*L /rhoa2/Ra *(1+a2/rhoa2+(a2/2.0)/Ra2 )-L /Ra3 ) );
        break;
    default:
        fprintf(stderr, "SegmentStressHor: unknown form!");
    }
    Sigma[0][0]=sunit*s[0];
    Sigma[0][1]=sunit*s[1];
    Sigma[0][2]=sunit*s[2];
    Sigma[1][1]=sunit*s[3];
    Sigma[1][2]=sunit*s[4];
    Sigma[2][2]=sunit*s[5];

    Sigma[1][0]=Sigma[0][1];
    Sigma[2][0]=Sigma[0][2];
    Sigma[2][1]=Sigma[1][2];    
}


void SegmentStress(real8 MU, real8 NU,
                   real8 bX, real8 bY, real8 bZ,
                   real8 xA, real8 yA, real8 zA,
                   real8 xB, real8 yB, real8 zB,
                   real8 x0, real8 y0, real8 z0,
                   real8 a,
                   real8 Sigma[3][3] )
{/* segment goes from (xA,yA,zA)->(xB,yB,zB)
    Burgers vector (bX,bY,bZ)
    field point at (x0,y0,z0)
    shear modulus MU, Poisson's ratio NU,
    regularization radius a (dislocation core width),
    return stress: Sigma[3][3]
 */
    real8 tmp[3][3], M[3][3], Mt[3][3];
    real8 xAB, yAB, zAB, xA0, yA0, zA0, rAB;
    real8 e1x, e1y, e1z, e2x, e2y, e2z, e3x, e3y, e3z;
    real8 z1, z2, xr0, zr0, brX, brY, brZ;

    xAB=xB-xA; yAB=yB-yA; zAB=zB-zA;
    xA0=x0-xA; yA0=y0-yA; zA0=z0-zA;

    rAB=sqrt(xAB*xAB+yAB*yAB+zAB*zAB);

    if(rAB<FFACTOR_LMIN) {
		//printf("(%e %e %e) - (%e %e %e)\n",xA,yA,zA,xB,yB,zB);
		Init3x3(Sigma);
		return;
        //Fatal("SegmentStress: rAB < LMIN");
	}

    e1x=xAB/rAB; e1y=yAB/rAB; e1z=zAB/rAB;

    xunitvector(e1x,e1y,e1z,xA0,yA0,zA0,&e3x,&e3y,&e3z);
    xvector(e3x,e3y,e3z,e1x,e1y,e1z,&e2x,&e2y,&e2z);

    /* e1,e2,e3 should all be normalized */

    M[0][0]=e2x; M[0][1]=e3x; M[0][2]=e1x;
    M[1][0]=e2y; M[1][1]=e3y; M[1][2]=e1y;
    M[2][0]=e2z; M[2][1]=e3z; M[2][2]=e1z;

    Mt[0][0]=e2x; Mt[0][1]=e2y; Mt[0][2]=e2z;
    Mt[1][0]=e3x; Mt[1][1]=e3y; Mt[1][2]=e3z;
    Mt[2][0]=e1x; Mt[2][1]=e1y; Mt[2][2]=e1z;

    /*
    printf("%20.10e %20.10e %20.10e\n"
           "%20.10e %20.10e %20.10e\n"
           "%20.10e %20.10e %20.10e\n",
           M[0][0], M[0][1], M[0][2],
           M[1][0], M[1][1], M[1][2],
           M[2][0], M[2][1], M[2][2] );
    */
    
    /* rr21=Mt*r21*/
    z1=0;
    z2=Mt[2][0]*xAB+Mt[2][1]*yAB+Mt[2][2]*zAB;
    /* rr=Mt*rt*/
    xr0=Mt[0][0]*xA0+Mt[0][1]*yA0+Mt[0][2]*zA0;
    zr0=Mt[2][0]*xA0+Mt[2][1]*yA0+Mt[2][2]*zA0;
    /* br=Mt*b*/
    brX=Mt[0][0]*bX+Mt[0][1]*bY+Mt[0][2]*bZ;
    brY=Mt[1][0]*bX+Mt[1][1]*bY+Mt[1][2]*bZ;
    brZ=Mt[2][0]*bX+Mt[2][1]*bY+Mt[2][2]*bZ;
    
    SegmentStressHor(MU,NU,brX,brY,brZ,z1,z2,xr0,zr0,a,Sigma);

    Matrix33Mult33(M,Sigma,tmp);
    Matrix33Mult33(tmp,Mt,Sigma);
    
    return;
}

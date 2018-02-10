/**
 * Definition of generic functions to be used to compute
 * f (ydot=f(y)), g (out=g(y)) and Jac (Jacobian)
 *
 * **/

#include <stdio.h>
#include <math.h>
//#include <mex.h> //MATLAB
#include	<sundials/sundials_types.h> /* definition of type realtype */
#define LH_APP_SQRT_EPS   0.001
#define LH_APP_DERSQRT_MAX   15.81
#define LH_APP_SQRT_C3   -12132
#define LH_APP_SQRT_C2   -171.6
#define LH_APP_SQRT_C1   12.93
#define LH_APP_SQRT_C0   0.0


/** fitting of derivable square root approximation around 0 
 *  if x > eps then sqrt(x)
 *  else if (0< x <=eps) polynom interpolation 
 *  test positive: return -1.0 as error code
 * */

realtype lh_approx_squareroot(realtype x) {
    if (x<0.0) return(-1.0);
    if (x>=LH_APP_SQRT_EPS) return sqrt(x);

    else return (LH_APP_SQRT_C0 + LH_APP_SQRT_C1*x + LH_APP_SQRT_C2*pow(x,2) + LH_APP_SQRT_C3*pow(x,3));
}

realtype lh_approx_der_squareroot(realtype x) {
    if (x<0.0) return(-1.0);
    if (x>=LH_APP_SQRT_EPS) return 1/(2*sqrt(x));

    else {
	    return (LH_APP_SQRT_C1 + 2*LH_APP_SQRT_C2*pow(x,1) + 3*LH_APP_SQRT_C3*pow(x,2));
            printf("DER SQUAREROOT INFO!\n");
    }
}

realtype lh_jacobian_pressure_first_diag(realtype K, realtype T, realtype x1, realtype x2) {
  realtype val, den, num, result;

  val=T*(pow(x1,2)-pow(x2,2));
  den=(2*lh_approx_squareroot(val));
  num= -2*K*x1;
      // first diag is the derivative % x1 of dx1/dt (first term of the jacobian) default is negative
  result=num/den;
  return result;
}

realtype lh_jacobian_pressure_second_diag(realtype K, realtype T, realtype x1, realtype x2) {
  realtype val, den, num, result;

  val=T*(pow(x1,2)-pow(x2,2));
  den=(2*lh_approx_squareroot(val));
  num= 2*K*x2;
  // second diag is the derivative % x2 of dx1/dt (second term of the first line of the jacobian)
  result=num/den;
  return result;
}

realtype lh_jacobian_pressure_main_diag(realtype K1, realtype K2, realtype T1, realtype T2, realtype x1, realtype x2, realtype x3) {
  realtype diag_1, diag_2, result;

  diag_2 = - lh_jacobian_pressure_second_diag(K1, T1, x1, x2);
  diag_1 =  lh_jacobian_pressure_first_diag(K2, T2, x2, x3);
  // main diag is the derivative % xi of the pressure der function=dxi/dt  
//  printf("diag1:%f, diag2:%f\n", diag_1, diag_2);
  result = diag_1 + diag_2;
  return result;
}

realtype lh_compute_valve_K(realtype* K_Table, realtype diam, realtype posdeg, int nx, int ny) {
//should be with #include <mex.h>
//realtype lh_compute_valve_K(mxArray* K_Table, realtype diam, realtype posdeg)
   realtype K1,K2,K;
  // perform a 2D interpolation with (x,y)=(diam, posdeg)
   int i,j;
   realtype x1, x2, y1, y2;
   realtype z1, z2, z3, z4;
   i=0; j=0;
   x1=K_Table[(i+1)*(ny+1)]; 
   x2=K_Table[(i+2)*(ny+1)]; 
   y1=K_Table[j+1];
   y2=K_Table[j+2];  // only to escape warning ...
   // The K_Table contains the x_data in the first column and
   // the y_data in the first row; (first row, first col) elmt=1
   while ((i<nx-1) && (diam >= x1)) {  // search diameter between-points
      x1=K_Table[(i+1)*(ny+1)];
      x2=K_Table[(i+2)*(ny+1)];
      if ( (diam >= x1) && (diam < x2) ) {
         break;
      }
      i++;
   }
   while ((j<ny-1) && (posdeg >= y1)) {
      y1=K_Table[j+1];
      y2=K_Table[j+2];
      if ( (posdeg >= y1) && (posdeg < y2) ) {
         break;
      }
      j++;
   }
 //  printf("i:%d, x1:%f, x2:%f, j:%d, y1:%f, y2:%f\n", i,x1,x2,j,y1,y2);
   if (i==nx)                        // diam outside table, right
      return K_Table[(nx)*(ny+1)+j];
   else if ((i==0) && (diam < x1)) { // diam outside table, left
      return K_Table[(i+1)*(ny+1)+j+1];
   }
   else { 
	   // look-up for the four-between points
      z1=K_Table[(i+1)*(ny+1)+j+1];
      z2=K_Table[(i+2)*(ny+1)+j+1];
      z3=(j<ny)?K_Table[(i+1)*(ny+1)+j+2] : K_Table[(i+1)*(ny+1)+ny];
      z4=(j<ny)?K_Table[(i+2)*(ny+1)+j+2] : K_Table[(i+2)*(ny+1)+ny];

  //   printf("di:¨%f, pos:%f, z1:%f, z2:%f, z3:%f, z4:%f\n", diam,posdeg,z1,z2,z3,z4);

      if (j==ny) {  // then z1==z3 && z2==z4
         K=z1+(diam-x1)/(x2-x1)*(z2-z1);
      }
      else  {
         K1=z1+(posdeg-y1)/(y2-y1)*(z3-z1); 
         K2=z2+(posdeg-y1)/(y2-y1)*(z4-z2);
	 // the final linear inteprolation:
         K=K1+(diam-x1)/(x2-x1)*(K2-K1);	 
      }
      
      return K;
   }
}


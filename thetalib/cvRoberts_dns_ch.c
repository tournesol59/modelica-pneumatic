/*
 * -----------------------------------------------------------------
 * $Revision: 4834 $
 * $Date: 2016-08-01 16:59:05 -0700 (Mon, 01 Aug 2016) $
 * -----------------------------------------------------------------
 * Programmer(s): Scott D. Cohen, Alan C. Hindmarsh and
 *                Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Example problem:
 * 
 * The following is a simple example problem, with the coding
 * needed for its solution by CVODE. The problem is from
 * pneumatic system, following system of equations  
 * y1=press,vol1, P1=fixed pressure source ahead vol1
 * y2=press,vol2, P2=fixed pressure ground downstream vol2
 *
 *    dy1/dt = RT/VOL1*(K1*sqrt((P1^2-y1^2)/T) - K2*sqrt(y1^2-y2^2)/T))
 *    dy2/dt = RT/VOL2*(K2*sqrt(y1^2-y2^2)/T) - K3*sqrt((y2^2-P2^2)/T))
 *
 * on the interval from t = 0.0 to t = 4, with initial
 * conditions: y1 = y2 = 1.0E5. The problem is stiff.
 * While integrating the system, we also use the rootfinding
 * feature to find the points at which y1 = y2 .
 * This program solves the problem with the BDF method,
 * Newton iteration with the CVDENSE dense linear solver, and a
 * user-supplied Jacobian routine.
 * It uses a scalar relative tolerance and a vector absolute
 * tolerance. Output is printed in decades from t = .4 to t = 4
 * Run statistics (optional outputs) are printed at the end.
 * -----------------------------------------------------------------
 */

#include <stdio.h>
#include <unistd.h>
#include <malloc.h>
/* Header files with a description of contents used */

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <cvode/cvode_impl.h>
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <sunmatrix/sunmatrix_dense.h> /* access to dense SUNMatrix            */
#include <sunlinsol/sunlinsol_dense.h> /* access to dense SUNLinearSolver      */
#include <cvode/cvode_direct.h>        /* access to CVDls interface            */
#include <sundials/sundials_dense.h> /* definitions DlsMat DENSE_ELEM */
#include <sundials/sundials_types.h> /* definition of type realtype */

/* Problem constants*/

#define DIAM1 25.0
#define DIAM2 35.0
#define DIAM3 25.0
#define P1 150000.0
#define P2 100000.0
#define VOL1 0.05
#define VOL2 0.025

/* User_data: in case of closed-loop control */
typedef struct {
   realtype u1, u2;
}  *UserData;

/* User-defined vector and matrix accessor macros: Ith, IJth */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   DENSE_ELEM macro in dense.h. DENSE_ELEM numbers rows and columns of a
   dense matrix starting from 0. */

/* These macros are defined in order to write code which exactly matches
   the mathematical problem description given above.

   Ith(v,i) references the ith component of the vector v, where i is in
   the range [1..NEQ] and NEQ is defined below. The Ith macro is defined
   using the N_VIth macro in nvector.h. N_VIth numbers the components of
   a vector starting from 0.

   IJth(A,i,j) references the (i,j)th element of the dense matrix A, where
   i and j are in the range [1..NEQ]. The IJth macro is defined using the
   SM_ELEMENT_D macro in dense.h. SM_ELEMENT_D numbers rows and columns of 
   a dense matrix starting from 0. */

#define Ith(v,i)    NV_Ith_S(v,i-1)         /* Ith numbers components 1..NEQ */
#define IJth(A,i,j) SM_ELEMENT_D(A,i-1,j-1) /* IJth numbers rows,cols 1..NEQ */


/* Problem Constants */

#define NEQ   2                /* number of equations  */
#define Y1    RCONST(1.02E5)      /* initial y components */
#define Y2    RCONST(1.01E5)

#define RTOL  RCONST(1.0e-4)   /* scalar relative tolerance            */
#define ATOL1 RCONST(1.0e-8)   /* vector absolute tolerance components */
#define ATOL2 RCONST(1.0e-8)

#define TIME0    RCONST(0.0)      /* initial time           */
#define TIME1    RCONST(0.00001)      /* first output time      */
#define TPLUS    RCONST(1.0E-5)     /* output time additionner     */
#define TMULT    RCONST(1.005555)
#define NOUT  100000               /* number of output times, here necessary to compute here 1.0 sec */

/* Functions Called by the Solver */

static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data);

static int g(realtype t, N_Vector y, realtype *gout, void *user_data);

static int Jac(realtype t,
               N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

static UserData AllocUserData(void);
static void InitUserData(UserData data);
static void FreeUserData(UserData data);

/* Functions called by f,g and Jac */
realtype lh_approx_squareroot(realtype x);

realtype lh_approx_der_squareroot(realtype x);

realtype lh_jacobian_pressure_first_diag(realtype K, realtype T, realtype x1, realtype x2);

realtype lh_jacobian_pressure_second_diag(realtype K, realtype T, realtype x1, realtype x2);

realtype lh_jacobian_pressure_main_diag(realtype K1, realtype K2, realtype T1, realtype T2, realtype x1, realtype x2, realtype x3);

realtype lh_compute_valve_K(realtype* K_Table, realtype diam, realtype posdeg, int nx, int ny);

/* Private functions to output results */
static void PrintJacobian(SUNMatrix J);
static void PrintOutput(realtype t, realtype y1, realtype y2);
static void PrintRootInfo(int root_f1, int root_f2, int root_f3);

/* Private function to print final statistics */

static void PrintFinalStats(void *cvode_mem);

/* Private function to check function return values */

static int check_flag(void *flagvalue, const char *funcname, int opt);
/*
 *-------------------------------
 * Main Program
 *-------------------------------
 */

int main()
{
  realtype reltol, t, tout;
  N_Vector y, abstol;
  SUNMatrix A;
  SUNLinearSolver LS;
  void *cvode_mem;
  int flag, flagr, iout;
  int rootsfound[3];
  UserData data;
  pid_t child_pid;  
  int file_d[2];
  char *message[10];
  int nbytes;

  y = abstol = NULL;
  A = NULL;
  LS = NULL;
  cvode_mem = NULL;

  data = NULL;

  data = AllocUserData();
  InitUserData(data);
  /* Create serial vector of length NEQ for I.C. and abstol */
  y = N_VNew_Serial(NEQ);
  if (check_flag((void *)y, "N_VNew_Serial", 0)) return(1);
  abstol = N_VNew_Serial(NEQ); 
  if (check_flag((void *)abstol, "N_VNew_Serial", 0)) return(1);

  /* Initialize y */
  Ith(y,1) = Y1;
  Ith(y,2) = Y2;

  /* Set the scalar relative tolerance */
  reltol = RTOL;
  /* Set the vector absolute tolerance */
  Ith(abstol,1) = ATOL1;
  Ith(abstol,2) = ATOL2;

  /* Call CVodeCreate to create the solver memory and specify the 
   * Backward Differentiation Formula and the use of a Newton iteration */
  cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON);
  if (check_flag((void *)cvode_mem, "CVodeCreate", 0)) return(1);
  
  /* Set the pointer to user-defined data */
  flag = CVodeSetUserData(cvode_mem, data);
  if(check_flag(&flag, "CVodeSetUserData", 1)) return(1);

  /* Call CVodeInit to initialize the integrator memory and specify the
   * user's right hand side function in y'=f(t,y), the inital time TIME0, and
   * the initial dependent variable vector y. */
  flag = CVodeInit(cvode_mem, f, TIME0, y);
  if (check_flag(&flag, "CVodeInit", 1)) return(1);

  /* Call CVodeSVtolerances to specify the scalar relative tolerance
   * and vector absolute tolerances */
  flag = CVodeSVtolerances(cvode_mem, reltol, abstol);
  if (check_flag(&flag, "CVodeSVtolerances", 1)) return(1);

  /* Call CVodeRootInit to specify the root function g with 2 components */
  flag = CVodeRootInit(cvode_mem, 2, g);
  if (check_flag(&flag, "CVodeRootInit", 1)) return(1);

  /* Create dense SUNMatrix for use in linear solves */
  A = SUNDenseMatrix(NEQ, NEQ);
  if(check_flag((void *)A, "SUNDenseMatrix", 0)) return(1);

  /* Create dense SUNLinearSolver object for use by CVode */
  LS = SUNDenseLinearSolver(y, A);
  if(check_flag((void *)LS, "SUNDenseLinearSolver", 0)) return(1);

  /* Call CVDlsSetLinearSolver to attach the matrix and linear solver to CVode */
  flag = CVDlsSetLinearSolver(cvode_mem, LS, A);
  if(check_flag(&flag, "CVDlsSetLinearSolver", 1)) return(1);

  /* Set the user-supplied Jacobian routine Jac */
  flag = CVDlsSetJacFn(cvode_mem, Jac);
  if(check_flag(&flag, "CVDlsSetJacFn", 1)) return(1);

  /* In loop, call CVode, print results, and test for error.
     Break out of loop when NOUT preset output times have been reached.  */
  printf(" \n2-pressures pneumatic system\n\n");

  iout = 0;  tout = TIME1;

 pipe(file_d);
 if ( child_pid = fork() == -1) {
  perror("fork error");
  exit(-1);
 }
  while(1) {
   if child_pid == 1) { 
// close file descriptor
     close(file_d[0]);

    //computation:
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);
    if (flag == CV_TOO_MUCH_WORK) {
       printf("INFO: CV too much work: mxstep:%d\n", MXSTEP_DEFAULT);
    }
    PrintOutput(t, Ith(y,1), Ith(y,2));

    if (flag == CV_ROOT_RETURN) {
      flagr = CVodeGetRootInfo(cvode_mem, rootsfound);
      if (check_flag(&flagr, "CVodeGetRootInfo", 1)) return(1);
      PrintRootInfo(rootsfound[0],rootsfound[1],rootsfound[2]);
    }

    if (check_flag(&flag, "CVode", 1)) break;
    if (flag == CV_SUCCESS) {
      iout++;
      if (tout < 1.0E-3) {
         tout *= TMULT;
      }
      else {
	 tout += TPLUS;
      }
    }
    if (iout == NOUT) break;

     send(file_d[1], message, (strlen(message)+1));

   }//end child_pid
   else {
    // parent pid
     close(file_d[1]);

     nbytes = read(file_d[0], message, sizeof(message));

     printf("message received: %s", message);
   }
  } 

  /* Print some final statistics */
  printf("iout:%d\n", iout);
  PrintFinalStats(cvode_mem);

  /* Free y and abstol vectors */
  N_VDestroy_Serial(y);
  N_VDestroy_Serial(abstol);


  /* Free integrator memory */
  CVodeFree(&cvode_mem);

  /* Free the linear solver memory */
  SUNLinSolFree(LS);

  /* Free the matrix memory */
  SUNMatDestroy(A);

  FreeUserData(data);
  return(0);
}


/*
 *-------------------------------
 * Functions called by the solver
 *-------------------------------
 */
static UserData AllocUserData(void)
{
  UserData data;
  data = (UserData) malloc(sizeof *data);
  return (data);
}


/* Load problem constants in data */

static void InitUserData(UserData data)
{
  data->u1 =30.0;
  data->u2 = 0.0;
}

static void FreeUserData(UserData data)
{
  free(data);
}

/*
 * f routine. Compute function f(t,y). 
 */
static int f(realtype t, N_Vector y, N_Vector ydot, void *user_data)
{
  realtype y1, y2, yd1, yd2, val1, val2, val3, K1, K2, K3;
  realtype K_Table[20]={1, 10, 30, 60, 90,  25, 1, 2, 5, 10,  50, 2.5, 5, 12, 20,  100, 5, 12, 25, 45};

  realtype posdeg1, posdeg2, posdeg3;
  realtype Temp=293;
  UserData data;

  data=(UserData)user_data;
  posdeg1=data->u1;
  posdeg2=25; // interpolation not optimized for pos==90 deg=max
  posdeg3=45;

  y1 = Ith(y,1); y2 = Ith(y,2);

  val1=P1*P1-y1*y1;
  val2=y1*y1-y2*y2;
  val3=y2*y2-P2*P2;

  K1=lh_compute_valve_K(K_Table, DIAM1, posdeg1, 3, 4);
  K2=lh_compute_valve_K(K_Table, DIAM2, posdeg2, 3, 4); 
  K3=lh_compute_valve_K(K_Table, DIAM3, posdeg3, 3, 4);

 /*
 * compute derivatives of pressure:
 *    dy1/dt = RT/VOL1*(K1*sqrt((P1^2-y1^2)/T) - K2*sqrt(y1^2-y2^2)/T))
 *    dy2/dt = RT/VOL2*(K2*sqrt(y1^2-y2^2)/T) - K3*sqrt((y2^2-P2^2)/T))
 */
  yd1 = RCONST(Temp*288)/VOL1*(K1*lh_approx_squareroot(val1/Temp)
		  - K2*lh_approx_squareroot(val2/Temp) );
  Ith(ydot,1) = yd1;
  yd2 = RCONST(Temp*288)/VOL2*(K2*lh_approx_squareroot(val2/Temp)
		  - K3*lh_approx_squareroot(val3/Temp) );
  Ith(ydot,2) = yd2;

  return(0);
}

/*
 * g routine. Compute functions g_i(t,y) for i = 0,1. 
 */

static int g(realtype t, N_Vector y, realtype *gout, void *user_data)
{
  realtype y1, y2;
  realtype val1, val2, val3;

  y1 = Ith(y,1); y2 = Ith(y,2);

  val1=P1*P1-y1*y1;
  val2=y1*y1-y2*y2;
  val3=y2*y2-P2*P2;

  gout[0] = lh_approx_squareroot(val1);
  gout[1] = lh_approx_squareroot(val2);
  gout[2] = lh_approx_squareroot(val3);

  return(0);
}

/*
 * Jacobian routine. Compute J(t,y) = df/dy. *
 */

static int Jac( realtype t,
               N_Vector y, N_Vector fy, SUNMatrix J, void *user_data,
               N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)
{
  realtype y1, y2;
  realtype T=300;
  realtype K1, K2, K3;
  realtype posdeg1, posdeg2, posdeg3;
  realtype K_Table[20]={1, 10, 30, 60, 90,  25, 1, 2, 5, 10,  50, 2.5, 5, 12, 20,  100, 5, 12, 25, 45};
  UserData data;

  data=(UserData)user_data;
  posdeg1 = data->u1;
  posdeg2=25;
  posdeg3=45;

  y1 = Ith(y,1); y2 = Ith(y,2);

  K1=lh_compute_valve_K(K_Table, DIAM1, posdeg1, 3, 4);
  K2=lh_compute_valve_K(K_Table, DIAM2, posdeg2, 3, 4);
  K3=lh_compute_valve_K(K_Table, DIAM3, posdeg3, 3, 4);
//  printf("K1:%f, K2:%f, K3:%f\n", K1, K2, K3);
  IJth(J,1,1) = lh_jacobian_pressure_second_diag(K1, T, P1, y1) +
	  lh_jacobian_pressure_first_diag(K2, T, y1, y2);
  IJth(J,1,2) = lh_jacobian_pressure_second_diag(K2, T, y1, y2);
  
  IJth(J,2,1) = -lh_jacobian_pressure_first_diag(K2, T, y1, y2); //default is negative =>*(-1)
  IJth(J,2,2) = lh_jacobian_pressure_main_diag(K2, K3, T, T, y1, y2, P2); 

//  PrintJacobian(J);
  return(0);
}

/*
 *-------------------------------
 * Private helper functions
 *-------------------------------
 */

static void PrintJacobian(SUNMatrix J) {
   int i,j;
   realtype *MatData;

   MatData=SUNDenseMatrix_Data(J); 
   for( i=0; i<NEQ; i++) {
      for (j=0; j<NEQ; j++) {
          printf(" %f", MatData[i*2+j]);
      }
   
      printf("\n");
   }
   return;
}

static void PrintOutput(realtype t, realtype y1, realtype y2)
{
#if defined(SUNDIALS_EXTENDED_PRECISION)
  printf("At t = %0.4Le      y =%14.6Le  %14.6Le\n", t, y1, y2);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
  printf("At t = %0.4e      y =%14.6e  %14.6e\n", t, y1, y2);
#else
  printf("At t = %0.4e      y =%14.6e  %14.6e\n", t, y1, y2);
#endif

  return;
}

static void PrintRootInfo(int root_f1, int root_f2, int root_f3)
{
  printf("    rootsfound[] = %3d %3d %3d\n", root_f1, root_f2, root_f3);

  return;
}

/* 
 * Get and print some final statistics
 */

static void PrintFinalStats(void *cvode_mem)
{
  long int nst, nfe, nsetups, nje, nfeLS, nni, ncfn, netf, nge;
  int flag;

  flag = CVodeGetNumSteps(cvode_mem, &nst);
  check_flag(&flag, "CVodeGetNumSteps", 1);
  flag = CVodeGetNumRhsEvals(cvode_mem, &nfe);
  check_flag(&flag, "CVodeGetNumRhsEvals", 1);
  flag = CVodeGetNumLinSolvSetups(cvode_mem, &nsetups);
  check_flag(&flag, "CVodeGetNumLinSolvSetups", 1);
  flag = CVodeGetNumErrTestFails(cvode_mem, &netf);
  check_flag(&flag, "CVodeGetNumErrTestFails", 1);
  flag = CVodeGetNumNonlinSolvIters(cvode_mem, &nni);
  check_flag(&flag, "CVodeGetNumNonlinSolvIters", 1);
  flag = CVodeGetNumNonlinSolvConvFails(cvode_mem, &ncfn);
  check_flag(&flag, "CVodeGetNumNonlinSolvConvFails", 1);

  flag = CVDlsGetNumJacEvals(cvode_mem, &nje);
  check_flag(&flag, "CVDlsGetNumJacEvals", 1);
  flag = CVDlsGetNumRhsEvals(cvode_mem, &nfeLS);
  check_flag(&flag, "CVDlsGetNumRhsEvals", 1);

  flag = CVodeGetNumGEvals(cvode_mem, &nge);
  check_flag(&flag, "CVodeGetNumGEvals", 1);

  printf("\nFinal Statistics:\n");
  printf("nst = %-6ld nfe  = %-6ld nsetups = %-6ld nfeLS = %-6ld nje = %ld\n",
	 nst, nfe, nsetups, nfeLS, nje);
  printf("nni = %-6ld ncfn = %-6ld netf = %-6ld nge = %ld\n \n",
	 nni, ncfn, netf, nge);
}

/*
 * Check function return value...
 *   opt == 0 means SUNDIALS function allocates memory so check if
 *            returned NULL pointer
 *   opt == 1 means SUNDIALS function returns a flag so check if
 *            flag >= 0
 *   opt == 2 means function allocates memory so check if returned
 *            NULL pointer 
 */

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
  int *errflag;

  /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
  if (opt == 0 && flagvalue == NULL) {
    fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  /* Check if flag < 0 */
  else if (opt == 1) {
    errflag = (int *) flagvalue;
    if (*errflag < 0) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
	      funcname, *errflag);
      return(1); }}

  /* Check if function returned NULL pointer - no memory allocated */
  else if (opt == 2 && flagvalue == NULL) {
    fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
	    funcname);
    return(1); }

  return(0);
}

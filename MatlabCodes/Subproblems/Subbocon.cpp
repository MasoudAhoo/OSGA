/* ===========================================================================
   =============================== SUBBOCON ================================== 
   ===========================================================================             
    _______________________________________________________________________ 
   | Subbocon is a function computing the global maximizer of the rational | 
   | bound-constrained subproblem                                          |
   |            max  -(gamma + <h,x>) / Q(x)   s.t. x \in [xl, xu]         |
   | of OSGA.                                                              |
   |                                                                       |
   | INPUT:                                                                |  
   |       h     : a parameter of the subproblem;                          | 
   |       gamma : a parameter of the subproblem;                          |
   |       x0    : initial point;                                          |
   |       q0    : center of prox-function;                                |
   |       xl    : lower bound on x;                                       |
   |       xu    : upper bound on x;                                       |
   |                                                                       |
   | OUTPUT:                                                               |
   |       u     : the global maximizer of the subproblem;                 |
   |       Egh   : the maximum.                                            |
   |_______________________________________________________________________|

    _______________________________________________________________________ 
   | Reference  : An optimal subgradient algorithm for large-scale bound-  |
   |              constrained convex optimization                          |
   | Authors    : Masoud Ahookhosh, Arnold Neumaier                        |
   | Written by : Masoud Ahookhosh                                         |
   | Date       : March 2014                                               |
   |_______________________________________________________________________|

============================================================================ */
#include <iostream>
#include <cmath> 
#include <algorithm>
#include <climits>
//#include <cstdlib>
//#include <vector>
#include "mex.h"
#include "matrix.h"
#include <omp.h>

using namespace std;

/* -------------------------- Input arguments ----------------------------- */
#define	h_IN	   prhs[0]
#define	gamma_IN   prhs[1]
#define	x0_IN	   prhs[2]
#define	q0_IN	   prhs[3]
#define	xl_IN	   prhs[4]
#define	xu_IN	   prhs[5]
/* ------------------------------------------------------------------------ */

/* -------------------------- Output arguments ---------------------------- */
#define	u_OUT	   plhs[0]
#define	Egh_OUT	   plhs[1]
/* ------------------------------------------------------------------------ */

/* ======================================================================== */
/* ================= Start of the computational routine =================== */
/* ======================================================================== */
static void subboconcpp(
	               double u[],
	               double *Egh,
                       double h[],
	               double gamma,
                       double x0[],
	               double q0,
                       double xl[],
	               double xu[],
                       int    dim
	               )
{
   
   int    n,m,infinity = INT_MAX,i,k,r,s,t;;
   double *lamb  = new double[dim];
   double *lamb1 = new double[dim];
   double *lam   = new double[dim];
     
   /* ------------------ Computing lambl(i) and lambu(i) ------------------ */
   /* n is the length of x0 */   
   n = dim;  
  
   #pragma omp parallel for shared(n,h,x0,xl,xu)//num_threads(8)
   for (int i=0; i<n; i++) 
   {
      double lamb_flag;
      if (h[i] > 0) 
      {
             lamb_flag = (x0[i] - xl[i])/h[i];
      } 
      else if (h[i] < 0) 
      {
             lamb_flag = (x0[i] - xu[i])/h[i];
      }
      else
      {
             lamb_flag = + infinity;
      }
      lamb[i]  = lamb_flag;
      lamb1[i] = lamb_flag;
   }
   /* --------------------------------------------------------------------- */

   /* ----------------- Sorting lamb and constructing lam ----------------- */
   
   /* lamb = sort(lamb,1); */
   sort(lamb1, lamb1+n);
   
   int  q = 0;
   lam[0] = lamb1[0];
   
   bool *dif   = new bool[n-1];
   /*#pragma omp parallel for shared(n,lamb1)  */
   for (int j=0; j<n-1; j++) 
   {
      dif[j] = (lamb1[j+1] != lamb1[j]);
   }

   /*#pragma omp parallel for shared(n,lamb1) firstprivate(q) lastprivate(m)*/
   for (int j=0; j<n; j++) 
   {
      if (dif[j]) 
      {
         q      = q + 1;
         lam[q] = lamb1[j];
      }
   }
 
   if (lam[q] < infinity)
   {
       m      = q+1;
       lam[m] = + infinity;
   }
   else
   {
       m = q;
   }
      
  
   /* --------------------------------------------------------------------- */

   double *E  = new double[m]; 
   double *LH = new double[m];
 
   #pragma omp parallel for shared(m,n,lam,lamb,h,x0,xl,xu) 
   for (k=0; k<m; k++) 
   {
      double hpk, hqk, dk1, sk1, n2pk_x0;
      double ak, bk, ck, dk, sk; 
      double lam_hat, phik, ek1, ek2; 
      double *pk = new double[dim];
      double *qk = new double[dim];   

      hpk     = 0;     
      hqk     = 0;
      n2pk_x0 = 0;
      dk1     = 0;
      sk1     = 0; 

      for (int i=0; i<n; i++) 
      {
         double pk_x0;
         if (lamb[i] <= lam[k]) 
         {
              if (h[i]<0)
              {
                 pk[i] = xu[i];
              }
              else if (h[i]>0)
              {
                 pk[i] = xl[i];
              }
              qk[i] = 0;
         }
         else
         {
              pk[i] = x0[i];
              qk[i] = -h[i];
         }

         hpk     += h[i]*pk[i];
         hqk     += h[i]*qk[i];
         pk_x0    = pk[i]-x0[i];
         n2pk_x0 += pk_x0*pk_x0;
         dk1     += pk_x0*qk[i];
         sk1     += qk[i]*qk[i];
      }
         
      /* ------- Compute ak, bk, ck, dk and sk to construct e(lam) ------- */
      ak = - (gamma + hpk);
      bk = - hqk;
      ck = q0 + 0.5 * n2pk_x0;
      dk = dk1;
      sk = 0.5 * sk1;
      /* ----------------------------------------------------------------- */ 

      /* - Finding the global maximizer of e(lam) for [lam(k), lam(k+1)] - */
      /* --------------------- using Proposition 4 ----------------------- */
      if (bk != 0) 
      {
          double w = ak*ak - bk*(ak*dk - bk*ck)/sk;
          lam_hat  = (-ak + sqrt(w))/bk;
          phik     = bk / (2*sk*lam_hat + dk); 
      } 
      else 
      {
          if (ak > 0) 
          { 
              lam_hat = -dk / (2 * sk);
              phik    = 4*ak*sk / (4*ck*sk - (dk*dk));
          } 
          else if (ak < 0)
          {
              lam_hat = + infinity;
              phik    = 0;
          }
          else
          {
              lam_hat = 0;
              phik    = 0;
          }
      }
      /* ----------------------------------------------------------------- */

      /* --- Checking the feasibility of the solution of Proposition 4 --- */
      if (lam[k] <= lam_hat && lam_hat <= lam[k + 1]) 
      {
          LH[k] = lam_hat;           
          E[k] = phik;
      } 
      else 
      {
          ek1 = (ak+bk*lam[k])/(ck+(dk+sk*lam[k])*lam[k]);
          ek2 = (ak+bk*lam[k+1])/(ck+(dk+sk*lam[k+1])*lam[k+1]);      
          if (ek1 >= ek2) 
          {
              LH[k] = lam[k];
              E[k]  = ek1;               
          } 
          else 
          { 
              LH[k] = lam[k+1];
              E[k] = ek2;
          }
      }
      /* ------------------------------------------------------------------ */
      delete[] pk;
      delete[] qk;

   }/* ------------------------- End of for (k) --------------------------- */

   /*sort(E, E+m);
   double e=E[m-1];*/
   double e  = E[0];
   double lh = LH[0];
   /*int flag = 0;
   #pragma omp parallel for*/
   for (int i=1; i<m; i++) 
   {
       if (E[i] > e) 
       {
          e  = E[i]; 
          lh = LH[i];
       }
   }

   double lambda = lh;
   /*double *u1 = new double[n];
   /*#pragma omp parallel for shared(m,h,x0,xl,xu) firstprivate(lambda)*/
   for (int i=0; i<m; i++) 
   {
       u[i] = x0[i]-lambda*h[i];
       if (u[i] < xl[i]) 
       {
          u[i] = xl[i];
       }
       else if (u[i] > xu[i])
       {
          u[i] = xu[i];
       }
       /*else
       {
          u[i] = u1[1];
       }*/
   }

   /* --------- The global maximizer by searching all m intervals --------- */
   *Egh = e;
   delete[] dif;
   delete[] E;
   delete[] lamb1;
   delete[] lamb;
   delete[] lam; 

   return;
   /* --------------------------------------------------------------------- */

}
/* ======================================================================== */
/* ================== End of the computational routine ==================== */
/* ======================================================================== */


/* ======================================================================== */
/* ===================== Start of the gateway routine ===================== */ 
/* ======================================================================== */
void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
{
   double *Egh, *u;
   double *h, gamma, *x0, q0, *xl, *xu;
   int    dim;
   
   /* ---------------- Check for proper number of arguments --------------- */
   if (nrhs != 6) 
   { 
	      mexErrMsgTxt( "Only 6 inputs allowed for Subbocon.\n" );
   } 
   else if (nlhs > 2) 
   {
	      mexErrMsgTxt( "Too many output arguments.\n" ); 
   }
 
   /* ----------------- Check for proper input arguments ------------------ */
   if ( !mxIsClass(prhs[0],"double") )
   {
        mexErrMsgTxt("1st argument of Subbocon must be double.\n");
   }
   if ( !mxIsClass(prhs[1],"double") )
   {
        mexErrMsgTxt("2st argument of Subbocon must be double.\n");
   }
   if ( !mxIsClass(prhs[2],"double") )
   {
        mexErrMsgTxt("3st argument of Subbocon must be double.\n");
   }
   if ( !mxIsClass(prhs[3],"double") )
   {
        mexErrMsgTxt("4st argument of Subbocon must be double.\n");
   }
   if ( !mxIsClass(prhs[4],"double") )
   {
        mexErrMsgTxt("5st argument of Subbocon must be double.\n");
   }
   if ( !mxIsClass(prhs[5],"double") )
   {
        mexErrMsgTxt("6st argument of Subbocon must be double.\n");
   }

   /* -------------- Create a matrix for the return argument -------------- */ 
   dim     = mxGetM(h_IN);
   u_OUT   = mxCreateDoubleMatrix( dim, 1, mxREAL);
   Egh_OUT = mxCreateDoubleMatrix( 1, 1, mxREAL);

   /* ------------- Assign pointers to the various parameters ------------- */
   h     = mxGetPr (h_IN);
   gamma = mxGetScalar (gamma_IN);
   x0    = mxGetPr (x0_IN);
   q0    = mxGetScalar (q0_IN);
   xl    = mxGetPr (xl_IN);  
   xu    = mxGetPr (xu_IN);

   u     = mxGetPr(u_OUT);
   Egh   = mxGetPr(Egh_OUT);

   /* ------------ Do the actual computations in a subroutine ------------- */
   subboconcpp(u, Egh, h, gamma, x0, q0, xl, xu, dim); 
   return;

}
/* ======================================================================== */
/* ===================== End of the gateway routine ======================= */ 
/* ======================================================================== */



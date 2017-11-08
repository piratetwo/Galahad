#include <limits.h>
#include <float.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#define INT long int
#define INT_INF LONG_MAX
#define INF DBL_MAX

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

#ifndef NULL
#define NULL 0
#endif

typedef struct cg_parameter_struct /* user controlled parameters */
{
    int PrintFinal ;
} cg_parameter ;

void cg_default
(
    cg_parameter   *Parm
)
{
    /* T => print final function value
       F => no printout of final function value */
    Parm->PrintFinal = FALSE ;
}

typedef struct cg_stats_struct /* statistics returned to user */
{
    double               f ; /*function value at solution */
    double           gnorm ; /* max abs component of gradient */
    INT               iter ; /* number of iterations */
    INT            IterSub ; /* number of subspace iterations */
    INT             NumSub ; /* total number subspaces */
    INT              nfunc ; /* number of function evaluations */
    INT              ngrad ; /* number of gradient evaluations */
} cg_stats ;

int cg_descent /*  return status of solution process:
                       0 (convergence tolerance satisfied)
                       1 (change in func <= feps*|f|)
                       2 (total number of iterations exceeded maxit)
                       3 (slope always negative in line search)
                       4 (number of line search iterations exceeds nline)
                       5 (search direction not a descent direction)
                       6 (excessive updating of eps)
                       7 (Wolfe conditions never satisfied)
                       8 (debugger is on and the function value increases)
                       9 (no cost or gradient improvement in
                          2n + Parm->nslow iterations)
                      10 (out of memory)
                      11 (function nan or +-INF and could not be repaired)
                      12 (invalid choice for memory parameter) */
(
    double            *x, /* input: starting guess, output: the solution */
    INT                n, /* problem dimension */
    cg_stats       *Stat, /* structure with statistics (can be NULL) */
    cg_parameter  *UParm, /* user parameters, NULL = use default parameters */
    double      grad_tol, /* StopRule = 1: |g|_infty <= max (grad_tol,
                                           StopFac*initial |g|_infty) [default]
                             StopRule = 0: |g|_infty <= grad_tol(1+|f|) */
    double      (*value) (double *, INT),  /* f = value (x, n) */
    void         (*grad) (double *, double *, INT), /* grad (g, x, n) */
    double    (*valgrad) (double *, double *, INT), /* f = valgrad (g, x, n),
                          NULL = compute value & gradient using value & grad */
    double         *Work  /* NULL => let code allocate memory
                             not NULL => use array Work for required memory
                             The amount of memory needed depends on the value
                             of the parameter memory in the Parm structure.
                             memory > 0 => need (mem+6)*n + (3*mem+9)*mem + 5
                                           where mem = MIN(memory, n)
                             memory = 0 => need 4*n */
)
{
  int status ;
  double f, gnorm ;
  double *g ;
  g = (double *) malloc (n*sizeof (double)) ;

  printf (" Calling dummy cg_descent\n");
  f = value (x, n);
  grad (g, x, n);
  f = valgrad (g, x, n);
  gnorm=fabs(g[0]);

  Stat->f = f ;
  Stat->gnorm = gnorm ;
  Stat->iter = 0 ;
  Stat->IterSub = 0 ;
  Stat->NumSub = 0 ;
  Stat->nfunc = 2 ;
  Stat->ngrad = 2 ;
  printf (" Returning from dummy cg_descent\n");

  status = 0 ;
  return( status ) ;
}

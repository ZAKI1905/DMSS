/*
  LikeLi_Eq class
*/

#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_multiroots.h>

#include <Math/QuantFuncMathMore.h>
#include <Math/Util.h>

#include <Zaki/Math/GSLMultiFWrapper.hpp>
#include <Zaki/Math/GSLMultiFdfWrapper.hpp>

#include "DMSS/LikeLi_Eq.hpp"

// using namespace DMSS ;
//==============================================================

//--------------------------------------------------------------
/// Default Constructor
LikeLi_Eq::LikeLi_Eq() {}

//--------------------------------------------------------------
/// F(x) = 0 Equation
int LikeLi_Eq::Eq(const gsl_vector* x, gsl_vector* in_f)
{
  
  const double sig_str = gsl_vector_get(x,0);
  double val    = 0 ;
  double sig_bg_rate  = 0 ;

  for(size_t i=0 ; i<obs_set.size() ; ++i)
    for (size_t j = 0; j < obs_set[i].size(); j++)
  {

    sig_bg_rate = sig_str*sig_shape_hist[i][ j+1 ] + gsl_vector_get(x,i+1) ;

    /// safe evaluation of log(x) with a protections against negative or zero argument 
    val += obs_set[i][j].val*ROOT::Math::Util::EvalLog(sig_bg_rate) - sig_bg_rate ;
  }




  const double x0 = gsl_vector_get(x,0);
  const double x1 = gsl_vector_get(x,1);

  // Simple system of equations
  // f_1(x,y) = a*x - 4
  // f_2(x,y) = b*x + y*x
  // const double y0 = pars.a * x0  - 4 ;
  // const double y1 = pars.b * x0 + x1 * x0 ;

  // Rosenbrock system of equations
  // f_1(x,y) = a * (1 - x)
  // f_2(x,y) = b * (y - x^2)
  // const double y0 = pars.a * (1 - x0);
  // const double y1 = pars.b * (x1 - x0 * x0);

  // gsl_vector_set (in_f, 0, y0);
  // gsl_vector_set (in_f, 1, y1);

  return GSL_SUCCESS ;
}

//--------------------------------------------------------------
/// F(x) = 0 Equation Solver
int LikeLi_Eq::SolveEq()
{
  const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrid ;
  gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (T, 2) ;

  printf ("Using '%s' solver...\n", gsl_multiroot_fsolver_name(s));


  Zaki::Math::GSLMultiFWrapper<LikeLi_Eq, 
    int (LikeLi_Eq::*)(const gsl_vector*, gsl_vector*)> 
  func(this, &LikeLi_Eq::Eq, 2) ;

  gsl_multiroot_function F = static_cast<gsl_multiroot_function> (func) ; 

  // gsl_set_error_handler_off() ;
  gsl_vector *r_guess = gsl_vector_alloc (2);

  gsl_vector_set (r_guess, 0, root_guess[0]);
  gsl_vector_set (r_guess, 1, root_guess[1]);

  gsl_multiroot_fsolver_set (s, &F, r_guess);

  int status;
  int iter = 0, max_iter = 1000;
  // gsl_vector* r ;
  do
  {
    iter++;
    status = gsl_multiroot_fsolver_iterate (s);

    PrintState(iter, s);
    
    if (status)   /* check if solver is stuck */
        break;

    status = gsl_multiroot_test_residual (s->f, 1e-7) ;
  }
  while (status == GSL_CONTINUE && iter < max_iter);

  gsl_multiroot_fsolver_free (s);
  gsl_vector_free(r_guess) ;
  
  return status ;
}
//--------------------------------------------------------------
// Prints the status
void LikeLi_Eq::PrintState(const size_t& iter, gsl_multiroot_fsolver * s)
{
  printf ("iter = %3zu x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
}
//--------------------------------------------------------------
// Prints the status
void LikeLi_Eq::PrintState(const size_t& iter, gsl_multiroot_fdfsolver * s)
{
  printf ("iter = %3zu x = % .3f % .3f "
          "f(x) = % .3e % .3e\n",
          iter,
          gsl_vector_get (s->x, 0),
          gsl_vector_get (s->x, 1),
          gsl_vector_get (s->f, 0),
          gsl_vector_get (s->f, 1));
}
//--------------------------------------------------------------
/// Jacobian
int LikeLi_Eq::Jacob(const gsl_vector* x, gsl_matrix* J)
{
  const double x0 = gsl_vector_get(x,0);
  const double x1 = gsl_vector_get(x,1);

  // gsl_matrix_set (J, 0, 0, pars.a );
  gsl_matrix_set (J, 0, 1, 0      );
  // gsl_matrix_set (J, 1, 0, pars.b + x1 );
  gsl_matrix_set (J, 1, 1, x0     );

  return 0 ;
}
//--------------------------------------------------------------
/// { F(x), F'(x) }
int LikeLi_Eq::EqFdf(const gsl_vector* x, gsl_vector* f, gsl_matrix* J) 
{

  const double x0 = gsl_vector_get(x,0);
  const double x1 = gsl_vector_get(x,1);

  // const double y0 = pars.a * x0  - 4 ;
  // const double y1 = pars.b * x0 + x1 * x0 ;

  // gsl_vector_set (f, 0, y0);
  // gsl_vector_set (f, 1, y1);

  // gsl_matrix_set (J, 0, 0, pars.a );
  gsl_matrix_set (J, 0, 1, 0      );
  // gsl_matrix_set (J, 1, 0, pars.b + x1 );
  gsl_matrix_set (J, 1, 1, x0     );

  return 0 ;
}

//--------------------------------------------------------------
int LikeLi_Eq::SolveEqJacob()
{
  const gsl_multiroot_fdfsolver_type * T = gsl_multiroot_fdfsolver_gnewton;
  gsl_multiroot_fdfsolver * s = gsl_multiroot_fdfsolver_alloc (T, 2);

  printf ("Using '%s' solver...\n", gsl_multiroot_fdfsolver_name(s));

  Zaki::Math::GSLMultiFdfWrapper<LikeLi_Eq, 
    int (LikeLi_Eq::*)(const gsl_vector*, gsl_vector*), 
    int (LikeLi_Eq::*)(const gsl_vector*, gsl_matrix*),
    int (LikeLi_Eq::*)(const gsl_vector*, gsl_vector*, gsl_matrix*)
    > 

  fdf(this, &LikeLi_Eq::Eq, &LikeLi_Eq::Jacob, 
      &LikeLi_Eq::EqFdf, 2) ;

  gsl_multiroot_function_fdf Fdf = static_cast<gsl_multiroot_function_fdf> (fdf) ; 

  gsl_vector *r_guess = gsl_vector_alloc (2);

  gsl_vector_set (r_guess, 0, root_guess[0]);
  gsl_vector_set (r_guess, 1, root_guess[1]);

  gsl_multiroot_fdfsolver_set(s, &Fdf, r_guess) ;

  // Iterations begin
  int status ;
  int iter = 0, max_iter = 1000;
  do
  {
    iter++ ;
    status = gsl_multiroot_fdfsolver_iterate(s) ;

    // ................................................
    // Error handling
    if (status == GSL_EBADFUNC)
    {
      std::cout <<  "Error: The iteration encountered a singular" 
                    " point where the function or its derivative"
                    " evaluated to Inf or NaN.\n" ;
    } else if (status == GSL_EZERODIV)
    {
      std::cout <<  "Error: The derivative of the function vanished"
                    " at the iteration point, preventing the algorithm"
                    " from continuing without a division by zero.\n" ;
    }
    // ................................................

    PrintState(iter, s);

    if (status)
      break;

    status = gsl_multiroot_test_residual (s->f, 1e-7);

  } 
  while ( status == GSL_CONTINUE && iter < max_iter );
  


  gsl_multiroot_fdfsolver_free(s) ;
  gsl_vector_free(r_guess) ;

  return status ;
}
// --------------------------------------------------------------
// void LikeLi_Eq::SetParams(const Params& in_pars)
// {
//   pars = in_pars ;
// }

//==============================================================

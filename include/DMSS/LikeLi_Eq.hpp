#ifndef DMSS_LikeLi_Eq_H
#define DMSS_LikeLi_Eq_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

#include "DMSS/Bin.hpp"
#include "DMSS/Prog.hpp"

//==============================================================
// namespace DMSS
// {
//==============================================================
class LikeLi_Eq : public Prog
{
  
  double root_guess[2] = {-10, -5} ;
  //--------------------------------------------------------------
  public:

    /// Default Constructor
    LikeLi_Eq() ;

    /// Destructor
    ~LikeLi_Eq() {} ;

    int SolveEq() ;
    int SolveEqJacob() ;

    int Eq(const gsl_vector* x, gsl_vector* f) ;
    int Jacob(const gsl_vector* x, gsl_matrix* J) ;
    int EqFdf(const gsl_vector* x, gsl_vector* f, gsl_matrix* J) ;

    void AddObsCounts(const Bin&) ;
    void AddSigShape(const TH1F&) ;

    void PrintState(const size_t& iter, gsl_multiroot_fsolver * s) ;
    void PrintState(const size_t& iter, gsl_multiroot_fdfsolver * s) ;

  private:
    double theErrorDef = 1 ;
    std::vector<std::vector<Zaki::Math::Quantity> > obs_set ;
    std::vector<TH1F> sig_shape_hist ;

};

//==============================================================
// } // DMSS namespace
//==============================================================
#endif /*DMSS_LikeLi_Eq_H*/

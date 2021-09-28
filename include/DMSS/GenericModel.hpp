#ifndef DMSS_GenericModel_H
#define DMSS_GenericModel_H

#include "DMSS/Model.hpp"

//==============================================================
class GenericModel : public Model
{

  //--------------------------------------------------------------
  public:

    GenericModel()  ;
    ~GenericModel() ;

    //............................................
    // Plotters
    //............................................
    // Plot 

    // Overriding the base class (Model) method
    // Threshold case
    double ContFuncThresh_LG(double, double) override ;
    double ContFuncThresh_MG(double, double) override ;
    double ContFuncThresh_LM(double, double) override ;

    // Boosted case
    double ContFuncBoost_LG(double, double) override ;
    double ContFuncBoost_MG(double, double) override ;
    double ContFuncBoost_LM(double, double) override ;

  //--------------------------------------------------------------
  private:
    
    virtual GenericModel* IClone() const override ;
};

//==============================================================
#endif /*DMSS_DarkPhoton_H*/

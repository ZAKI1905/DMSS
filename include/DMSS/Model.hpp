#ifndef DMSS_Model_H
#define DMSS_Model_H

#include <map>

#include "DMSS/Prog.hpp"
#include "DMSS/Data.hpp"

//==============================================================
enum class PlotMode
{
  Ldec_Gann = 0, Mdm_Gann, Ldec_Mdm
};

//==============================================================
struct PlotOptions
{ 
  // width and height
  double w,h ;

  // X-axis and Y-axis labels
  std::string xLabel, yLabel;

  // connect the points or not
  bool joined;
};

//==============================================================
struct PlotGrid
{
  Zaki::Math::Grid2D grid ;
  PlotOptions option ;
  bool fixed = false ;
};

//==============================================================
class Model : public Prog
{
  //--------------------------------------------------------------
  public:
    //..............................................
    ~Model() ; 
    // ***********  
    // Constructors
    // ***********  
    /// Default
    Model() ;  
    
    /// Via name
    Model(const std::string&) ;
       
    // Via name & dimension
    // Model(const std::string& in_name, const size_t in_dim) 
    //   : name(in_name) {  SetDim(in_dim) ;} 

    std::unique_ptr<Model> Clone() const
    {
      return std::unique_ptr<Model>(IClone()) ;
    }
    //..............................................

    // ***********  
    // Setters
    // ***********  
    virtual void SetDMMass(const double) ;
    virtual void SetDMMass(const std::vector<double>&) ;

    void SetECut(const double) ;
    void SetTGrid(const PlotMode&, const Zaki::Math::Grid2D&, const PlotOptions&, const bool fixed) ;
    void SetBGrid(const PlotMode&, const Zaki::Math::Grid2D&,  const PlotOptions&, const bool fixed) ;
    void SetDecayProd(const std::vector<DMSolarSignal::Particle>&);

    /**
     Initializing the model for example for fixing the relic abundance, etc.
     */
    virtual void Init() ;

    void SetActiveBin(const Bin&) ;
    //..............................................

    // ***********  
    // Getters
    // ***********  
    double GetDMMass() const ;
    std::vector<double> GetDMMassSet() const ;
    double GetECut() const ;
    std::vector<DMSolarSignal::Particle> GetDecayProd() const ;
    
    PlotGrid* GetTGrid(const PlotMode&) ;
    PlotGrid* GetBGrid(const PlotMode&) ;
    
    virtual double GetSpectrum(Bin&) ;

    /**
     If the grid values (in probing the parameter space)
     are fixed with respect to e.g. the dark matter mass
     if yes, it will save a lot of time in generating
     threshold plots. If no, the grid has to be updated
     for each contour.
    */
    // bool HasFixedGrid(const PlotMode&) ;

    // Override the Prog class print method
    // void Print() const override ;

    /**
     A function of 2 variables that will be used
     for making contour plots using CONFIND
     which is called from 'Analysis' class
     and via the child class only
     */
    
    ///  Threshold case, in (L_dec, G_ann) plane
    virtual double ContFuncThresh_LG(double, double) = 0 ;

    ///  Threshold case, in (M_dm, G_ann) plane
    virtual double ContFuncThresh_MG(double, double) = 0 ;

    ///  Threshold case, in (L_dec, M_dm) plane
    virtual double ContFuncThresh_LM(double, double) = 0 ;

    ///  Boosted case, in (L_dec, G_ann) plane
    virtual double ContFuncBoost_LG(double, double) = 0 ;

    ///  Boosted case, in (M_dm, G_ann) plane
    virtual double ContFuncBoost_MG(double, double) = 0 ;

    ///  Boosted case, in (L_dec, M_dm) plane
    virtual double ContFuncBoost_LM(double, double) = 0 ;

    /// Gets the bin that is under consideration (active)
    Bin GetActiveBin() const ;
  //--------------------------------------------------------------
  protected:
    
    virtual Model* IClone() const = 0 ;
    
    /// Dark matter mass (GeV)
    double mDM = -1 ;
    
    /// Energy cut (GeV)
    double energy_cut = 50 ;

    /// Boosted grid
    std::map<PlotMode, PlotGrid> boos_grids ;

    /// Threshold grid
    std::map<PlotMode, PlotGrid> thresh_grids ;

    /// A set of decay products
    std::vector<DMSolarSignal::Particle> decay_products ;

    /// Dark matter mass set
    std::vector<double> DM_mass_set ;

    Bin active_bin ;

};
//==============================================================

//==============================================================
#endif /*DMSS_Model_H*/

#ifndef DMSS_SatBundle_H
#define DMSS_SatBundle_H

#include "DMSS/Satellite.hpp"

//==============================================================
class SatBundle : public Prog
{
  //--------------------------------------------------------------
  private:
  std::shared_ptr<Satellite> m_SatPtr;
  std::vector<int> m_BinPeriods  ;

  //--------------------------------------------------------------
  public:
  /// Default Constructor
  SatBundle() : Prog("SatBundle", true), m_SatPtr(nullptr) { }

  /// Constructor via satellite pointer
  SatBundle(std::shared_ptr<Satellite> in_SatPtr) 
  : Prog("SatBundle", true), m_SatPtr(in_SatPtr) {}

  /// Constructor via satellite pointer & the bin periods
  SatBundle(std::shared_ptr<Satellite> in_SatPtr, 
            const std::vector<int>& in_BinPeriods) 
  : Prog("SatBundle", true), m_SatPtr(in_SatPtr), 
    m_BinPeriods(in_BinPeriods)
  {
    if ( m_SatPtr->GetData().GetBins().size() != m_BinPeriods.size() )
    {
      Z_LOG_ERROR(" The bin periods set size doesn't match the dataset size.") ;
    }
  }

  /// Copy Constructor
  SatBundle(const SatBundle& other) 
  : Prog("SatBundle", true), m_BinPeriods(other.m_BinPeriods)
  {
    if(other.m_SatPtr)
      m_SatPtr = other.m_SatPtr->Clone() ;
  }

  /// Assignment Operator
  SatBundle& operator=(const SatBundle&) ;

  /// Arrow '->' Operator
  std::shared_ptr<Satellite> operator->() 
  {
    return m_SatPtr;
  } 

  void SetBinPeriod(const std::vector<int>&) ;

  /// Gets the satellite pointer
  std::shared_ptr<Satellite> GetSatellite() ;

  /// Gets the bin periods
  std::vector<int> GetBinPeriod() const ;

  /// Gets the timed i-th energy bin with the j-th time binning 
  Bin GetTimeBin(const size_t& e_i, const size_t& t_j) ;
};

//==============================================================
#endif /*DMSS_SatBundle_H*/

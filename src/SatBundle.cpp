/*
  SatBundle class

*/

#include "DMSS/SatBundle.hpp"

//==============================================================
//--------------------------------------------------------------
SatBundle& SatBundle::operator=(const SatBundle& other)
{
  if(this == &other) return *this ;
  
  m_BinPeriods = other.m_BinPeriods ;
  if(other.m_SatPtr)
    m_SatPtr = other.m_SatPtr->Clone() ;
  
  return *this ;
}
//--------------------------------------------------------------
void SatBundle::SetBinPeriod(const std::vector<int>& per)
{
  if(!m_SatPtr)
  {
    Z_LOG_ERROR("Add a satellite first!") ;
    return ;
  }

  // if ( m_SatPtr->GetData().GetBins().size() != per.size() )
  // {
  //   Z_LOG_ERROR(" The bin periods set size doesn't match the dataset size.") ;
  // }

  m_BinPeriods = per ;
}

//--------------------------------------------------------------
std::vector<int> SatBundle::GetBinPeriod() const
{
  return m_BinPeriods ;
}

//--------------------------------------------------------------
/// Returns the satellite pointer
std::shared_ptr<Satellite> SatBundle::GetSatellite()
{
  return m_SatPtr ;
}

//--------------------------------------------------------------
Bin SatBundle::GetTimeBin(const size_t& b_idx, const size_t& t_idx)
{
  Bin b = m_SatPtr->GetData().GetBins()[b_idx];

  m_SatPtr->TimeBin(m_BinPeriods[t_idx]);

  size_t chops = m_SatPtr->GetSigShape(b.GetECenter().val).GetSize() -2  ;
  b.divide(chops);

  return b ;
}
//--------------------------------------------------------------

//==============================================================

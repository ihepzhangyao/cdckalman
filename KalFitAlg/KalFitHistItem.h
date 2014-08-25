// 
// Description: class for ntuple items
//

#ifndef KALFITHISTITEM_H
#define KALFITHISTITEM_H

#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/NTuple.h"
#include "AIDA/IHistogram1D.h"
#include "AIDA/IHistogram2D.h"

NTuple::Tuple* m_ntTest;  // test drift dist in Kalman yzhang 2014-05-16  yzhang 2014-05-16  yzhang 2014-05-16  yzhang 2014-05-16 
NTuple::Item<double> m_ddRec,m_ddKal,m_ddLayer,m_fiRec,m_fiKal,m_ddFiTerm,m_ddFltLenTerm,m_ddFltlen,m_ddWire,m_ddMass,m_ddTruth,m_ddFltlenKal,m_dchi2R,m_dchi2L,m_dchi2min,m_chi2Sum,m_chi2Add,m_dp,m_lastdp,m_pCurrent,m_pTruthCurrent;
NTuple::Item<long> m_iTurn,m_iFit,m_iKeep;
#endif 

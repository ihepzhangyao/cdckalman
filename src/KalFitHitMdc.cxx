//-----------------------------------------------------------------------
// File from KalFit module
//
// Filename : KalFitHitMdc.cc
//------------------------------------------------------------------------
// Description : Hit in Mdc 
//------------------------------------------------------------------------
// Modif :
//------------------------------------------------------------------------
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
	    typedef HepGeom::Point3D<double> HepPoint3D;
	#endif
#include "KalFitAlg/KalFitHitMdc.h"
#include "KalFitAlg/KalFitWire.h"
#include "MdcRecEvent/RecMdcHit.h"

// constructor
KalFitHitMdc::KalFitHitMdc(int id, int LR, double tdc, double dist[2],
		       double erdist[2], KalFitWire* wire, RecMdcHit* rechitptr):
  LR_(LR), tdc_(tdc), wire_(wire), rechitptr_(rechitptr), chi2_(0), chi2_back_(0), id_(id) {
  dist_[0]=dist[0];
  dist_[1]=dist[1];
  erdist_[0]=erdist[0];
  erdist_[1]=erdist[1];
}

KalFitHitMdc::KalFitHitMdc(int id, int LR, double tdc, double dist[2],
                       double erdist[2], KalFitWire* wire):
LR_(LR), tdc_(tdc), wire_(wire), chi2_(0), chi2_back_(0), id_(id) {

  dist_[0]=dist[0];
  dist_[1]=dist[1];
  erdist_[0]=erdist[0];
  erdist_[1]=erdist[1];
}


// destructor
KalFitHitMdc::~KalFitHitMdc(void){}

int KalFitHitMdc::isolated(void){
  
  return 1;
}

#include "KalFitAlg/KalFitHitMdc.h"
#include "KalFitAlg/KalFitHelixSeg.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector/ThreeVector.h"
KalFitHelixSeg::KalFitHelixSeg(KalFitHitMdc* hit,HepPoint3D pivot,HepVector a,
 HepSymMatrix Ea) :
Helix(pivot,a,Ea),
hitmdc_(hit),
a_pre_fwd_(5, 0),
a_pre_bwd_(5, 0),
a_filt_fwd_(5, 0),
a_filt_bwd_(5, 0),
Ea_pre_fwd_(5, 0),
Ea_filt_fwd_(5, 0),
Ea_pre_bwd_(5, 0),
Ea_filt_bwd_(5, 0)
{
 lr_ =0; 
 residual_exclude_ = 0.0;
 residual_include_ = 0.0;
 doca_exclude_ = 0.0;
 doca_include_ = 0.0;
 tof_ = 0.0;
 dt_ = 0.0;
 dd_ = 0.0;
}

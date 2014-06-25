//-----------------------------------------------------------------------
// File from KalFit module
//
// Filename : KalFitWire.cc
//------------------------------------------------------------------------
// Description : Description of a Mdc Wire
//------------------------------------------------------------------------
// Modif :
//------------------------------------------------------------------------
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "KalFitAlg/KalFitWire.h"
#include "KalFitAlg/KalFitLayer_Mdc.h"
#include "KalFitAlg/KalFitSuper_Mdc.h"

const double KalFitWire::A[NREGION] = { 8.5265E-7, 
				      1.1368E-6, 
				      1.2402E-6};
const double KalFitWire::F[NLAYER] = { 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
				    47.2,46.6,47.0,46,45,46.7,46,45,44.6,
				    45.9,45.9,45.3,44.5,44.0,46.0,45.6,
				    44.2,43.6,44.0,43.5,43.1,42.9,42.4,
				    42.7,42.5,42.2,41,42,42,42,41.4,42.5,
				    43.0,43.2,43.4,43.4};
const double L49_2 = 220.4*220.4;

// constructor
KalFitWire::KalFitWire(const int localID,
		   const KalFitLayer_Mdc & layer,
		   const HepPoint3D & fwd,
		   const HepPoint3D & bck,
		   KalFitWire * const vt,
		   unsigned int geoID,
		   unsigned int stereo):
  localId_(localID), fwd_(fwd), bck_(bck), layer_(layer), 
  state_(WireHitInvalid), geoID_(geoID), stereo_(stereo), 
  ddl_(0), ddr_(0), distance_(0) {

  int layer_ID = layer_.layerId();
  xyPosition_ = 0.5 * (fwd_+bck_);
  //yzhang delete for CDC 2014-04-10 
  //if (!stereo_){
  //  dx_ = 0;
  //  dy_ = 0;
  //  x_ = xyPosition_.x();
  //  y_ = xyPosition_.y();
  //} else {
    dx_ = (double)fwd_.x() - (double)bck_.x();
    dy_ = (double)fwd_.y() - (double)bck_.y();
    x_ = bck_.x();
    y_ = bck_.y();
  //}
  // Wire sag coeff
  Hep3Vector wire_;
  wire_ = (CLHEP::Hep3Vector)fwd - (CLHEP::Hep3Vector)bck;
  lzx_ = sqrt(wire_.z()*wire_.z()+wire_.x()*wire_.x());

  // problem with layer_ and its id inside sagcoef ?!!
  //  sagcoef();
  if (layer_ID < 3)
    A_ = A[0];
  else {
    if (layer_ID< 14)
      A_ = A[1];
    else if (layer_ID < 50){
      double f_current(F[layer_ID]);
      double l(wire_.mag());
      A_ = A[2]*F[49]*F[49]*L49_2/(f_current*f_current*l*l);
    } else
      std::cout << "*** PROBLEM WIRE !!!!! " << std::endl;
  }

  // Neighbor :
  neighbor_[0] = (KalFitWire *) innerLeft(vt);
  neighbor_[1] = (KalFitWire *) innerRight(vt);

  neighbor_[2] = (KalFitWire *) left();
  neighbor_[3] = (KalFitWire *) right();

  neighbor_[4] = (KalFitWire *) outerLeft(vt);
  neighbor_[5] = (KalFitWire *) outerRight(vt);

}

KalFitWire::KalFitWire():
  localId_(0), fwd_(0,0,0), bck_(0,0,0), layer_(*(KalFitLayer_Mdc *)NULL), 
  state_(WireHitInvalid), x_(0), y_(0), dx_(0), dy_(0)
{
  // Neighbor :
  neighbor_[0] = NULL;
  neighbor_[1] = NULL;
  neighbor_[2] = NULL;
  neighbor_[3] = NULL;
  neighbor_[4] = NULL;
  neighbor_[5] = NULL;
}

// destructor
KalFitWire::~KalFitWire(void){}

// wire sag coefficient determination :
void KalFitWire::sagcoef(void){

  Hep3Vector wire_;
  wire_ =(CLHEP::Hep3Vector)fwd_ - (CLHEP::Hep3Vector)bck_;
  int layer_ID = layer_.layerId();

  if (layer_ID < 3){
    A_ = A[0];
  } else {
    if (layer_ID< 14){
      A_ = A[1];
    } else if (layer_ID < 50){
      double f_current(F[layer_ID]);
      double l(wire_.mag());
      A_ = A[2]*F[49]*F[49]*L49_2/(f_current*f_current*l*l);
    } else {
      std::cout << "*** PROBLEM WIRE !!!!! " << std::endl;
    }
  }
}

const KalFitWire * KalFitWire::left(void) const{
  const KalFitWire * tmp = this;
  if (!localId_) tmp += layer_.superLayer().nWire();
  return --tmp;
}

const KalFitWire * KalFitWire::right(void) const{
  const KalFitWire * tmp = this;
  if (localId_==layer_.superLayer().localMaxId()){
    tmp -= layer_.superLayer().nWire();
  }
  return ++tmp;
}

const KalFitWire * KalFitWire::innerLeft(KalFitWire * const vtWire) const{

  if (!layer_.localLayerId()) return vtWire;
  const KalFitWire * tmp = this;
  if (layer_.offset()&1){
    tmp -= layer_.superLayer().nWire();
  }else{
    if (localId_) tmp -= layer_.superLayer().nWire();
    --tmp;
  }
  return tmp; 
}

const KalFitWire * KalFitWire::innerRight(KalFitWire * const vtWire) const{
   if (!layer_.localLayerId()) return vtWire;
  const KalFitWire * tmp = this;
  if (layer_.offset()&1){
    if (localId_==layer_.superLayer().localMaxId()){
      tmp -= layer_.superLayer().nWire();
    }
    tmp -= (layer_.superLayer().nWire()-1);
  }else{
    tmp -= layer_.superLayer().nWire();
  }
  return tmp; 
}

const KalFitWire * KalFitWire::outerLeft(KalFitWire * const vtWire) const{
  if (layer_.localLayerId()==layer_.superLayer().layerMaxId()) return vtWire;
  const KalFitWire * tmp = this;
  if (layer_.offset()&1){
    tmp += layer_.superLayer().nWire();
  }else{
    if (!localId_) tmp += layer_.superLayer().nWire();
    tmp += (layer_.superLayer().nWire()-1);
  }
  return tmp;  
}

const KalFitWire * KalFitWire::outerRight(KalFitWire * const vtWire) const{
  if (layer_.localLayerId()==layer_.superLayer().layerMaxId()) return vtWire;
  const KalFitWire * tmp = this;
  if (layer_.offset()&1){
    if (localId_^layer_.superLayer().localMaxId()){
      tmp += layer_.superLayer().nWire();
    }
    tmp++;
  }else{
    tmp += layer_.superLayer().nWire();
  }
  return tmp;
}

// Check left and right side to see if there is hit valid
void KalFitWire::chk_left_and_right(void) {
  if (((**(neighbor_+2)).state_&WireHit)&&
      ((**(neighbor_+3)).state_&WireHit)){
    state_ |= WireHitInvalid;
    (**(neighbor_+2)).state_ |= WireHitInvalid;
    (**(neighbor_+3)).state_ |= WireHitInvalid;
  }
}

 

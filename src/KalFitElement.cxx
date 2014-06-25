//-----------------------------------------------------------------------
// File from KalFit module
//
// Filename : KalFitElement.cc
//------------------------------------------------------------------------
// Description :
// Element is a base class which represents an element of detector 
//------------------------------------------------------------------------
// Modif :
//------------------------------------------------------------------------
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
	    typedef HepGeom::Point3D<double> HepPoint3D;
	#endif
//#include "TrackUtil/Helix.h"
#include "KalFitAlg/helix/Helix.h"
	    
#include "KalFitAlg/KalFitMaterial.h"
#include "KalFitAlg/KalFitTrack.h"
#include "KalFitAlg/KalFitElement.h"

int KalFitElement::loss_;
int KalFitElement::muls_;

void KalFitElement::updateTrack(KalFitTrack& track, int index)
{

  HepPoint3D x;
  double path = intersect(track, x);

  //cout<<__FILE__<<" "<<__LINE__<<" in updateTrack: path= "<<path<<" intersect x "<<x<<endl;
  //cout<<__FILE__<<" "<<__LINE__<<" in updateTrack path= "<<path<<" intersect x "<<x<<" change pivot to pivot_numf x"<<endl;
  
  if(path > 0){
    // move pivot

    //std::cout<<"KalFitElement: track helix1= "<<track.a()<<std::endl;
    //std::cout<<" KalFitTrack::numf_: "<<KalFitTrack::numf_<<std::endl;
    //std::cout<<" KalFitTrack::muls_: "<<muls_<<std::endl;
    //std::cout<<" KalFitTrack::loss_: "<<loss_<<std::endl;
    
    track.pivot_numf(x);
  
    //std::cout<<"KalFitElement: track helix2= "<<track.a()<<std::endl;
  
    
    // multiple scattering and energy loss
    int index_element(index);
    if (index_element==0) index_element=1;
    if(muls_) track.ms(path, *material_, index_element);
    if(loss_) track.eloss(path, *material_, index_element);

    //cout<<"KalFitElement: track helix3= "<<track.a()<<endl;
  }
  //cout<<"KalfitElement: track helix2= "<<track.a()<<endl;
}

//yzhang add for CDC 2014-05-02 
void KalFitElement::updateTrackBack(KalFitTrack& track, int index)
{

  HepPoint3D x;
  double path = intersectBack(track, x);

  //cout<<__FILE__<<" "<<__LINE__<<" in updateTrackBack path= "<<path<<" intersect x "<<x<<" change pivot to x"<<endl;
  
  if(path > 0){
    // move pivot

    //std::cout<<"KalFitElement: track helix1= "<<track.a()<<std::endl;
    //std::cout<<" KalFitTrack::numf_: "<<KalFitTrack::numf_<<std::endl;
    //std::cout<<" KalFitTrack::muls_: "<<muls_<<std::endl;
    //std::cout<<" KalFitTrack::loss_: "<<loss_<<std::endl;
    
    //yzhang change TEMP 2014-05-02 
    //track.pivot_numf(x);
    track.pivot(x,true);//yzhang phi not +2pi
  
    //std::cout<<"KalFitElement: track helix2= "<<track.a()<<std::endl;
  
    
    // multiple scattering and energy loss
    int index_element(index);
    if (index_element==0) index_element=1;
    if(muls_) track.ms(path, *material_, index_element);
    if(loss_) track.eloss(path, *material_, index_element);

    //cout<<"KalFitElement: track helix3= "<<track.a()<<endl;
  }
  //cout<<"KalfitElement: track helix2= "<<track.a()<<endl;
}


void KalFitElement::updateTrack_rphi(KalFitTrack& track, int index)
{
  HepPoint3D x;
  double path = intersect(track, x);
  if(path > 0){
    // move pivot
    track.pivot_numf(x);

    // multiple scattering and energy loss
    if(muls_) track.ms(path, *material_, index);
    if(loss_) track.eloss(path, *material_, index);
  }
}

void KalFitElement::updateTrack_alreadyfound(KalFitTrack& track,
					   int index)
{
  HepPoint3D x;
  double path = intersect(track, x);
  if(path > 0){
    // move pivot
    track.pivot_numf(x);

    // multiple scattering and energy loss
    if(muls_) track.ms(path, *material_, index);
    if(loss_) track.eloss(path, *material_, index);
  }
}


void KalFitElement::asso_rphi(Lpav& circ, KalFitTrack& track)
{}

void KalFitElement::asso_rphi(KalFitTrack& track)
{}

void KalFitElement::muls(int i)
{
  muls_ = i;
}

void KalFitElement::loss(int i)
{
  loss_ = i;
}

int KalFitElement::muls(void)
{
  return muls_;
}

int KalFitElement::loss(void)
{
  return loss_;
}

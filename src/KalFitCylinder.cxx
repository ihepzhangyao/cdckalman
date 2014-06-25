//-----------------------------------------------------------------------
// File from KalFit module
//
// Filename : KalFitCylinder.cc
//------------------------------------------------------------------------
// Description :
// Cylinder is an Element whose shape is a cylinder.
//------------------------------------------------------------------------
// Modif :
//------------------------------------------------------------------------
#include <float.h>
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
//#include "TrackUtil/Helix.h"
#include "KalFitAlg/helix/Helix.h"

#include "KalFitAlg/KalFitTrack.h"
#include "KalFitAlg/KalFitMaterial.h"
#include "KalFitAlg/KalFitCylinder.h"

//yzhang add 2014-05-02 
double KalFitCylinder::intersectBack(const KalFitTrack& track,
			       HepPoint3D& x) const
{
  double dPhi[4];
  dPhi[0] = track.intersect_cylinder_back(ri_);
  if(dPhi[0] == 0) return -1;
  dPhi[1] = track.intersect_cylinder_back(ro_);
  if(dPhi[1] == 0) return -1;
  dPhi[2] = track.intersect_xy_plane(zf_);
  dPhi[3] = track.intersect_xy_plane(zb_);
  //std::cout<<__FILE__<<"  "<<__LINE__<<" intersectBack dphi "<<dPhi[0]<<" "<<dPhi[1]<<" "<<dPhi[2]<<" "<<dPhi[3]<<std::endl;

  int n[2];
  int j = 0;
  for(int i = 0; i < 4 && j < 2; i++){
    HepPoint3D xx = track.x(dPhi[i]);
    if(isInside(xx)) n[j++] = i;
  }
  if(j < 2) return -1;

  x = track.x((dPhi[n[0]] + dPhi[n[1]]) * .5);

  double tanl = track.tanl();
  //std::cout<<__FILE__<<"   "<<__LINE__<< " KalFitCylinder: track radius"<<track.radius()<<" dphi0 " <<dPhi[n[0]]<<" dphi1 "<<dPhi[n[1]]<<" tanl "<<tanl<<endl;
  return fabs(track.radius() * (dPhi[n[0]] - dPhi[n[1]])
             * sqrt(1 + tanl * tanl));
}

double KalFitCylinder::intersect(const KalFitTrack& track,
			       HepPoint3D& x) const
{
  double dPhi[4];
  dPhi[0] = track.intersect_cylinder(ro_);
  if(dPhi[0] == 0) return -1;
  dPhi[1] = track.intersect_cylinder(ri_);
  if(dPhi[1] == 0) return -1;
  dPhi[2] = track.intersect_xy_plane(zf_);
  dPhi[3] = track.intersect_xy_plane(zb_);

 //std::cout<<__FILE__<<"  "<<__LINE__<<" intersect dphi "<<dPhi[0]<<" "<<dPhi[1]<<" "<<dPhi[2]<<" "<<dPhi[3]<<std::endl;
  int n[2];
  int j = 0;
  for(int i = 0; i < 4 && j < 2; i++){
    HepPoint3D xx = track.x(dPhi[i]);
    if(isInside(xx)) n[j++] = i;
  }
  if(j < 2) return -1;

  x = track.x((dPhi[n[0]] + dPhi[n[1]]) * .5);

  double tanl = track.tanl();
  //cout<<"KalFitCylinder: track radius"<<track.radius()<<" dphi0 "
  //    <<dPhi[n[0]]<<" dphi1 "<<dPhi[n[1]]<<" tanl "<<tanl<<endl;
  return fabs(track.radius() * (dPhi[n[0]] - dPhi[n[1]])
             * sqrt(1 + tanl * tanl));
}


double KalFitCylinder::intersect(const KalFitTrack& track,
                                   HepPoint3D& x, const HepPoint3D& point) const
{

  const double ro = sqrt(point.x()*point.x()+point.y()*point.y());

  //std::cout<<" ro: "<<ro<<std::endl;
  
  double dPhi[4];
  dPhi[0] = track.intersect_cylinder(ro);
  if(dPhi[0] == 0) return -1;
  dPhi[1] = track.intersect_cylinder(ri_);
  if(dPhi[1] == 0) return -1;
  dPhi[2] = track.intersect_xy_plane(zf_);
  dPhi[3] = track.intersect_xy_plane(zb_);

  //for(int ii=0; ii<4; ii++)
    //std::cout<<"dPhi["<<ii<<"]"<<dPhi[ii]<<std::endl;
  
  int n[2];
  int j = 0;
  for(int i = 0; i < 4 && j < 2; i++){
    HepPoint3D xx = track.x(dPhi[i]);
    if(isInside(xx)) n[j++] = i;
  }

  if(j < 2) return -1;

  x = track.x((dPhi[n[0]] + dPhi[n[1]]) * .5);
  
  double tanl = track.tanl();

  return fabs(track.radius() * (dPhi[n[0]] - dPhi[n[1]])
      * sqrt(1 + tanl * tanl));
 }




bool KalFitCylinder::isInside(const HepPoint3D& x) const
{
  double r = x.perp();
  double z = x.z();
  //std::cout<<"r: "<<r<<" z: "<<z<<" ri: "<<ri_<<" ro: "<<ro_<<" zb_: "<<zb_<<"zf: "<<zf_<<std::endl;
  
  return (r >= ri_ - FLT_EPSILON &&
          r <= ro_ + FLT_EPSILON &&
          z >= zb_ - FLT_EPSILON &&
          z <= zf_ + FLT_EPSILON);
}


bool KalFitCylinder::isInside2(const HepPoint3D& x) const
{
  double r = x.perp();
  double z = x.z();   
  //std::cout<<"r: "<<r<<" z: "<<z<<" ri: "<<ri_<<" ro: "<<ro_<<" zb_: "<<zb_<<"zf: "<<zf_<<std::endl;

  return (r <= ro_ + FLT_EPSILON &&
          z >= zb_ - FLT_EPSILON &&
          z <= zf_ + FLT_EPSILON);
}


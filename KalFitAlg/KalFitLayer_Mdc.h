//
//  Directly inspisrated from FTLayer of Kakuno-san
//
//
#ifndef _DEFINE_LAYER_Mdc_H_
#define _DEFINE_LAYER_Mdc_H_

#include <math.h>
#include <iostream>
#include "KalFitAlg/KalFitSuper_Mdc.h"
//class KalFitSuper_Mdc;

class KalFitLayer_Mdc{
public:
  /// constructor
  KalFitLayer_Mdc(const KalFitSuper_Mdc & super,
		const double radius, const double stereoAngle, 
		const double zf, const double zb, const double offset,
		const int layerID, const int localLayerID);
  
  /// destructor
  ~KalFitLayer_Mdc(){};

public: // Selectors
  /// returns local-layer ID
  const int localLayerId(void) const;

  /// returns layer ID
  const int layerId(void) const;

  /// returns tangent of slant angle
  const double tanSlant(void) const;

  /// returns r form origin
  const double r(void) const;

  /// returns z of forward end-plate
  const double zf(void) const;

  /// returns z of backward end-plate
  const double zb(void) const;

  /// returns offset of numbering(local ID)
  const int offset(void) const;

  /// returns super-layer
  const KalFitSuper_Mdc & superLayer(void) const;

  /// returns limit of "d" for stereo layer
  const double limit(void) const;

  /// returns z for "d" in r-phi plane
  double z(const double d) const;

  /// returns cell size
  double csize(void) const;

private: // private data members
  const double _radius;
  const double _tanSlant;
  const double _zf;				// z of forward endplate
  const double _zb;				// z of backward endplate
  const int _layerId;
  const int _localLayerId;
  const int _offset;
  const KalFitSuper_Mdc & _superLayer;
};

//----------------------------------
#ifdef KalFitLayer_Mdc_NO_INLINE
#define inline 
#else
#undef inline
#define KalFitLayer_Mdc_INLINE_DEFINE_HERE
#endif

#ifdef KalFitLayer_Mdc_INLINE_DEFINE_HERE

inline
KalFitLayer_Mdc::KalFitLayer_Mdc(const KalFitSuper_Mdc & super,
			     const double radius, const double stereoAngle, 
			     const double zf, const double zb, const double offset,
			     const int layerID, const int localLayerID)  : _radius(radius),
  _tanSlant(1./tan(stereoAngle)),
  _zf(zf),
  _zb(zb),
  _layerId(layerID),
  _localLayerId(localLayerID),
  _offset((int)(2.0*offset)),
  _superLayer(super)
{
}

inline
const int KalFitLayer_Mdc::layerId(void) const { return _layerId;}
inline
const int KalFitLayer_Mdc::localLayerId(void) const { return _localLayerId;}
inline
const double KalFitLayer_Mdc::tanSlant(void) const { return _tanSlant;}

inline
const double KalFitLayer_Mdc::r(void) const { return _radius;}

inline
const double KalFitLayer_Mdc::zf(void) const { return _zf;}

inline
const double KalFitLayer_Mdc::zb(void) const { return _zb;}

inline
const double KalFitLayer_Mdc::limit(void) const {  return (double)(_zf-_zb)/_tanSlant;}
  
inline
double 
KalFitLayer_Mdc::z(const double d) const
{
  return (double)_zb+d*_tanSlant;
}

inline
const int
KalFitLayer_Mdc::offset(void) const
{
  return _offset;
}

inline
const KalFitSuper_Mdc & 
KalFitLayer_Mdc::superLayer(void) const
{
  return _superLayer;
}

inline
double
KalFitLayer_Mdc::csize(void) const
{
  return 2*M_PI*_radius/_superLayer.nWire();
}

#endif

#undef inline

#endif 


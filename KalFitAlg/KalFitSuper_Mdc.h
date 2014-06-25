#ifndef KalFitSuper_Mdc_FLAG_
#define KalFitSuper_Mdc_FLAG_

#include "KalFitAlg/KalFitList.h"


class KalFitWire;
class KalFitSuper_Mdc{
public:
  /// Constructors and destructor
  KalFitSuper_Mdc(const int firstWireID,const int NWire,
			     const int firstLayerID,const int NLayer,
			     const int superLayerID);
  ~KalFitSuper_Mdc();

public:
  /// clear object
  void clear(void);


public: // Selectors

  /// returns number of wires
  const int nWire(void) const;

  /// returns number of layers
  const int nLayer(void) const;

  /// returns layer max ID
  const int layerMaxId(void) const;

  /// returns local max ID
  const int localMaxId(void) const;

  /// returns super-layer ID
  const int superLayerId(void) const;

  /// append a wireHit to the list of hits of the superlayer
  void appendHit(KalFitWire *);

private: //static data members
  static const unsigned int _neighborsMask[6];

private: // private data members
  const int _superLayerId;
  const int _firstWireId;
  const int _Nwire;
  const int _firstLayerId;
  const int _Nlayer;
  KalFitList<KalFitWire *>& _wireHits;
  KalFitList<KalFitWire *>& _singleHits;
};

#ifdef KalFitSuper_Mdc_NO_INLINE
#define inline
#else
#undef inline
#define KalFitSuper_Mdc_INLINE_DEFINE_HERE
#endif

#ifdef KalFitSuper_Mdc_INLINE_DEFINE_HERE

inline void KalFitSuper_Mdc::appendHit(KalFitWire *h){
  _wireHits.append(h);
}

#endif
#undef inline

#endif /* KalFitSuper_Mdc_FLAG_ */

//
#include "KalFitAlg/KalFitList.h"
#include "KalFitAlg/KalFitLayer_Mdc.h"
#include "KalFitAlg/KalFitSuper_Mdc.h"
#include "KalFitAlg/KalFitWire.h"

//...Globals...
const unsigned int
KalFitSuper_Mdc::_neighborsMask[6] = {452,420,340,172,92,60};
//                          [0] = 452 = FTWireHitAppended + FTWireNeighbor345
//                          [1] = 420 = FTWireHitAppended + FTWireNeighbor245
//     definitions of       [2] = 340 = FTWireHitAppended + FTWireNeighbor135
//    _neighborsMask[]      [3] = 172 = FTWireHitAppended + FTWireNeighbor024
//                          [4] = 92  = FTWireHitAppended + FTWireNeighbor013
//                          [5] = 60  = FTWireHitAppended + FTWireNeighbor012

KalFitSuper_Mdc::KalFitSuper_Mdc(const int firstWireID,const int NWire,
			   const int firstLayerID,const int NLayer,
			   const int superLayerID)
  : _superLayerId(superLayerID),
    _firstWireId(firstWireID),
    _Nwire(NWire),
    _firstLayerId(firstLayerID),
    _Nlayer(NLayer),
    _wireHits(*(new KalFitList<KalFitWire *>(500))),
    _singleHits(*(new KalFitList<KalFitWire *>(100)))
{

}

KalFitSuper_Mdc::~KalFitSuper_Mdc()
{
  clear();
  delete &_wireHits;
  delete &_singleHits;
}

void
KalFitSuper_Mdc::clear(void)
{
  if(_wireHits.length()){
    register KalFitWire ** hptr = _wireHits.firstPtr();
    KalFitWire ** const last = _wireHits.lastPtr();
    do {(**hptr).state(WireHitInvalid);}while((long)(hptr++)^(long)last);
    _wireHits.clear();
  }
}

const int
KalFitSuper_Mdc::nWire(void) const{
  return _Nwire;
}


const int
KalFitSuper_Mdc::nLayer(void) const{
  return _Nlayer;
}


const int
KalFitSuper_Mdc::localMaxId(void) const{
  return (_Nwire - 1);
}


const int
KalFitSuper_Mdc::layerMaxId(void) const{
  return (_Nlayer - 1);
}


const int
KalFitSuper_Mdc::superLayerId(void) const{
  return _superLayerId;
}


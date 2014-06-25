#include "KalFitAlg/KalFitAlg.h"
#include "KalFitAlg/KalFitTrack.h"
#include "KalFitAlg/KalFitWire.h"
//#include "EvTimeEvent/RecEsTime.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/Bootstrap.h"
#include "MagneticField/MagneticFieldSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/Algorithm.h"
#include "MdcRawEvent/MdcDigi.h"
#include "RawEvent/RawDataUtil.h"
#include "EventModel/Event.h"
#include "Identifier/MdcID.h"
#include "Identifier/Identifier.h"
#include "CLHEP/Matrix/SymMatrix.h"

double KalFitTrack::EventT0_ = 0.;
HepSymMatrix KalFitTrack::initMatrix_(5,0);
const MdcCalibFunSvc* KalFitTrack::CalibFunSvc_ = 0;
const IMagneticFieldSvc* KalFitTrack::MFSvc_ = 0;
IMdcGeomSvc* KalFitTrack::iGeomSvc_ = 0;
MdcDigiCol* KalFitTrack::mdcDigiCol_ = 0;
int KalFitTrack::tprop_ = 1;


void KalFitTrack::setInitMatrix(HepSymMatrix m)
{ 
	initMatrix_ = m;
}

HepSymMatrix KalFitTrack::getInitMatrix() const
{
	return initMatrix_ ;
}

//------------------set  event  start time-----------
void  KalFitTrack::setEventStartTime(double eventstart ) 
{
	EventT0_ = eventstart;
}

double  KalFitTrack::getEventStartTime(void ) const 
{
	//------------------get  event  start time-----------
	return EventT0_ ;
}

double KalFitTrack::getDriftTime(KalFitHitMdc& hitmdc , double toftime) const
{
	int layerid = hitmdc.wire().layer().layerId();
	int wireid = hitmdc.wire().geoID();
	double zhit = hitmdc.rechitptr()->getZhit();//cm
        double Q = hitmdc.rechitptr()->getAdc();
	double timewalk = CalibFunSvc_->getTimeWalk(layerid, Q);
	double timeoffset = CalibFunSvc_->getT0(wireid);
	double eventt0 = getEventStartTime(); 
	double tp = CalibFunSvc_->getTprop(layerid,zhit*10.);//zhit cm to mm
	double rawtime = RawDataUtil::MdcTime(hitmdc.rechitptr()->getTdc());

	double drifttime = rawtime - eventt0 - toftime - timewalk - timeoffset - tp;

	if(debug_ == 4 ) {
	  //std::cout<<"eventt0 "<<eventt0<<" toftime "<<toftime<<" timewalk "<<timewalk<<" timeoffset "<<timeoffset<<" tp "<<tp<<" drifttime "<<drifttime<<" rechit driftT "<<hitmdc.rechitptr()->getDriftT()<<std::endl;
	}
	if(drifttime_choice_ == 0){
	  return drifttime;
	}else if(drifttime_choice_ == 1){
	  // use the driftT caluculated by track-finding
	  return hitmdc.rechitptr()->getDriftT(); 
	}
}


// attention , the unit of the driftdist is mm
double KalFitTrack::getDriftDist(KalFitHitMdc& hitmdc, double drifttime, int lr) const
{
	int layerid = hitmdc.wire().layer().layerId();
	int  cellid = MdcID::wire(hitmdc.rechitptr()->getMdcId());
	//if(debug_ == 4 ){
	//	std::cout<<"the cellid is .."<<cellid<<std::endl;
	//} 
	double entrangle = hitmdc.rechitptr()->getEntra();

	//std::cout<<" entrangle: "<<entrangle<<std::endl;
	if(driftdist_choice_ == 0){
	  return CalibFunSvc_->driftTimeToDist(drifttime, layerid, cellid,  lr, entrangle);
	}else if(driftdist_choice_ == 1){
	  return fabs(hitmdc.rechitptr()->getDriftDistRight()); 
	}
}

// attention , the unit of the sigma is mm
double  KalFitTrack::getSigma(KalFitHitMdc& hitmdc, double tanlam, int lr, double dist) const
{ 
	int layerid = hitmdc.wire().layer().layerId();
	double entrangle = hitmdc.rechitptr()->getEntra();
	double z = hitmdc.rechitptr()->getZhit();
	double Q = hitmdc.rechitptr()->getAdc();
	double sigma = CalibFunSvc_->getSigma(layerid, lr, dist, entrangle, tanlam, z*10. , Q );
	if(debug_ == 4 ){
	  //cout<<"layerid "<<layerid<<" lr "<<lr<<" dist "<<dist<<" entra "<<entrangle<<" tanlam "<<tanlam<<" z "<<z<<" Q "<<Q<<" sigma from calib "<<sigma<<endl;//wangll
	}
	return sigma;
}

void KalFitTrack::setMdcCalibFunSvc(const MdcCalibFunSvc*  calibsvc)
{
	CalibFunSvc_ = calibsvc;
}

void KalFitTrack::setMagneticFieldSvc(IMagneticFieldSvc* mf)
{
	MFSvc_ = mf;
	if(MFSvc_==0) cout<<"KalFitTrack2:: Could not load MagneticFieldSvc!"<<endl;
} 

void KalFitTrack::setIMdcGeomSvc(IMdcGeomSvc*  igeomsvc)
{
	iGeomSvc_ = igeomsvc;
}

void KalFitTrack::setMdcDigiCol(MdcDigiCol*  digicol)
{
	mdcDigiCol_ = digicol;
}


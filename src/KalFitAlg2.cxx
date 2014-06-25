#include "KalFitAlg/KalFitAlg.h"
#include "KalFitAlg/KalFitTrack.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "MdcGeomSvc/IMdcGeomSvc.h"
#include <fstream>



 IMdcGeomSvc* KalFitAlg::imdcGeomSvc_ = 0;

/*void KalFitAlg::setInitMatrix(HepSymMatrix m)
 {
 initMatrix_ = m;

 }*/

 /*
 void KalFitAlg::setMaterial_Mdc(void)
{
  if(debug_ == 4) cout << " Mdc cylinder material file :\n";
  //  string  fn("$HOME/geomdc_material.dat");
  if(debug_ == 4) cout << " "  << matfile_ << std::endl;
  ifstream fin(matfile_.c_str() ); 
  if(!fin){
    perror(matfile_.c_str());
    cout<<"KalFitAlg::setMaterial.......file reading error"<<endl;
    exit(1);
  }

  double z, a, i, rho, x0;

  for (int j = 0; j < 6 &&
	 fin >> z >> a >> i >> rho >> x0; j++){
     KalFitMaterial  mat(z, a, i, rho, x0);
     _BesKalmanFitMaterials.push_back(mat);
  };
  fin.close();
}
*/


/*
void KalFitAlg::setCylinder_Mdc(void)
{
  double radius, thick, length , z0;

  /// innerwall of inner drift chamber
  radius = _BesKalmanFitTubs[0].GetInnerRadius()/(cm);
  thick  = _BesKalmanFitTubs[0].GetOuterRadius()/(cm) - _BesKalmanFitTubs[0].GetInnerRadius()/(cm);
  length = 2.0*_BesKalmanFitTubs[0].GetZHalfLength()/(cm); 
  z0     = 0.0;
  std::cout<<"innerwall: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl;
  KalFitCylinder innerwallCylinder(&(_BesKalmanFitMaterials[1]), radius, thick, length , z0);
  _BesKalmanFitWalls.push_back(innerwallCylinder);
  
  
  /// outer air, be attention the calculation of the radius and thick of the air cylinder is special 
  radius = _BesKalmanFitTubs[2].GetOuterRadius()/(cm);
  thick  = _BesKalmanFitTubs[0].GetInnerRadius()/(cm) - _BesKalmanFitTubs[2].GetOuterRadius()/(cm);
  length = 2.0*_BesKalmanFitTubs[0].GetZHalfLength()/(cm);
  z0     = 0.0;
  std::cout<<"outer air: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl;
  KalFitCylinder outerAirCylinder(&(_BesKalmanFitMaterials[4]), radius, thick, length , z0);
  _BesKalmanFitWalls.push_back(outerAirCylinder);
	   
 
  /// outer beryllium layer   
  radius = _BesKalmanFitTubs[2].GetInnerRadius()/(cm);
  thick  = _BesKalmanFitTubs[2].GetOuterRadius()/(cm) - _BesKalmanFitTubs[2].GetInnerRadius()/(cm);
  length = 2.0*_BesKalmanFitTubs[2].GetZHalfLength()/(cm);
  z0     = 0.0;
  std::cout<<"outer Be: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl;
  KalFitCylinder outerBeCylinder(&(_BesKalmanFitMaterials[2]), radius, thick, length , z0);
  _BesKalmanFitWalls.push_back(outerBeCylinder);
  

  ///  helium gas  
  radius = _BesKalmanFitTubs[3].GetInnerRadius()/(cm);
  thick  = _BesKalmanFitTubs[3].GetOuterRadius()/(cm) - _BesKalmanFitTubs[3].GetInnerRadius()/(cm);
  length = 2.0*_BesKalmanFitTubs[3].GetZHalfLength()/(cm);
  z0     = 0.0;
  std::cout<<"He gas: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl;
  KalFitCylinder heCylinder(&(_BesKalmanFitMaterials[3]), radius, thick, length , z0);
  _BesKalmanFitWalls.push_back(heCylinder);

    
  /// inner beryllium layer   
  radius = _BesKalmanFitTubs[4].GetInnerRadius()/(cm);
  thick  = _BesKalmanFitTubs[4].GetOuterRadius()/(cm) - _BesKalmanFitTubs[4].GetInnerRadius()/(cm);
  length = 2.0*_BesKalmanFitTubs[4].GetZHalfLength()/(cm);
  z0     = 0.0;
  std::cout<<"inner Be: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl;
  KalFitCylinder innerBeCylinder(&(_BesKalmanFitMaterials[2]), radius, thick, length , z0);
  _BesKalmanFitWalls.push_back(innerBeCylinder);

  /// inner air
  radius = 0.;
  thick  = _BesKalmanFitTubs[4].GetInnerRadius()/(cm);
  length = 2.0*_BesKalmanFitTubs[4].GetZHalfLength()/(cm);
  z0     = 0.0;
  std::cout<<"inner air: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl;
  KalFitCylinder innerAirCylinder(&(_BesKalmanFitMaterials[4]), radius, thick, length , z0);
  _BesKalmanFitWalls.push_back(innerAirCylinder);
  
}
*/


void KalFitAlg::set_Mdc(void)
{
   MdcGeomSvc * mdcGeomSvc = dynamic_cast<MdcGeomSvc* >(imdcGeomSvc_);
  if(! mdcGeomSvc) {
       std::cout<<"ERROR OCCUR when dynamic_cast in KalFitAlg2.cxx ..!!"<<std::endl;
   }

  if(debug_ == 4)  {
     cout << "KalFitAlg:: MDC initialisation " << std::endl;
   }
  _wire = NULL;
  _layer = NULL; 
  _superLayer = NULL;

  const int Nwire = mdcGeomSvc->getWireSize();
  const int Nlyr =  mdcGeomSvc->getLayerSize();
  const int Nsup =  mdcGeomSvc->getSuperLayerSize();

  if (!Nwire || !Nlyr || !Nsup){
    cout << ".....MdcGeom Objects are missing !! " << std::endl;
    exit(-1); 
  }
      
  if (!_wire) 
    _wire = (KalFitWire *) malloc((Nwire+1) * sizeof(KalFitWire));
  if (!_layer) 
    _layer = (KalFitLayer_Mdc *) malloc(Nlyr * sizeof(KalFitLayer_Mdc));
  if (!_superLayer) 
    _superLayer = (KalFitSuper_Mdc *) malloc(Nsup * sizeof(KalFitSuper_Mdc));

  if (!_wire || !_layer || !_superLayer){  
    std::cerr << "KalFitAlg::Cannot allocate geometries" << std::endl;
    std::cerr << "JOB will stop" << std::endl;
    exit(-1);
  }
   
  int superLayerID = 0;
  int layerID = 0;
  int localLayerID = 0;
  int localWireID = 0; 
  int localID = 0;
  int wireID;

  MdcGeoLayer * layer_back = NULL;
  MdcGeoSuper * superLayer_back = NULL;
  int k = 0;
  int Nlayer[12];
  int Nlocal[12];
  int NlocalWireID[43];

  for (wireID = 0;wireID <= Nwire; wireID++) {
    MdcGeoLayer * layer = (wireID==Nwire) ? NULL : mdcGeomSvc->Wire(wireID)->Lyr();
    if (layer != layer_back) {
      layer_back = layer;
      MdcGeoSuper * superLayer = (wireID==Nwire) ? NULL : mdcGeomSvc->Layer(layerID)->Sup();
      if (superLayer != superLayer_back) {
        superLayer_back = superLayer;
        Nlayer[k] = localLayerID;
        Nlocal[k] = localID;
        localLayerID = 0;
        k++;
      }
      NlocalWireID[layerID] = localWireID;
      localID = 0;
      localWireID = 0;
      layerID++;
      localLayerID++;
    }
    localID++;
    localWireID++;
  }

  superLayerID = -1;
  layerID = -1;
  localLayerID = 0;
  localID = 0;
  layer_back = NULL;
  superLayer_back = NULL;
  for (wireID = 0;wireID < Nwire; wireID++) {
    MdcGeoLayer * layer = (wireID==Nwire) ? 
                          NULL : mdcGeomSvc->Wire(wireID)->Lyr();
    if (layer != layer_back){
      layer_back = layer;
      MdcGeoSuper * superLayer = (wireID==Nwire) ? 
                                 NULL : mdcGeomSvc->Layer(layerID+1)->Sup();
      if (superLayer != superLayer_back){
        superLayer_back = superLayer;
        // initialize super-layer
        superLayerID++;
        new(_superLayer+superLayerID) KalFitSuper_Mdc(wireID,
                                                   Nlocal[superLayerID+1],
                                                   layerID+1,
                                                   Nlayer[superLayerID+1],
                                                   superLayerID);
        localLayerID=0;
      }
      // initialize layer
      layerID++;
      double slantWithSymbol = (mdcGeomSvc->Layer(layerID)->Slant())
                              *(mdcGeomSvc->Layer(layerID)->Sup()->Type());
      new(_layer+layerID) KalFitLayer_Mdc(_superLayer[superLayerID],
                                  0.1*layer->Radius(), (layer->Slant())*(layer->Sup()->Type()),
                                  0.1*(layer->Length()/2), 
                                  0.1*(-layer->Length()/2), layer->Offset(),
				  layerID, localLayerID++ );
      localID = 0;
    }
    // initialize wire

    const MdcGeoWire * wire = (wireID==Nwire) ? NULL : mdcGeomSvc->Wire(wireID); 
    HepPoint3D fwd(0.1*wire->Backward());
    HepPoint3D bck(0.1*wire->Forward());


    //yzhang changed for cdc 2014-04-04 all stereo layer
    //if (superLayerID == 2 || superLayerID == 3  ||
    //    superLayerID == 4 || superLayerID == 9 ||
    //    superLayerID == 10) {     // axial wire
      //new(_wire+wireID) KalFitWire(localID++,_layer[layerID],
        			   //fwd,bck, _wire+Nwire,wire->Id(),0);
        		      
    //} else {                      // stereo wire
      new(_wire+wireID) KalFitWire(localID++,_layer[layerID],
        			   fwd,bck, _wire+Nwire,wire->Id(),1);
    //}
  }
 
  // make virtual wire object for the pointer of boundary's neighbor
  new(_wire+Nwire) KalFitWire();

  if (debug_ == 4) cout << "MDC geometry reading done" << std::endl;  

  return;
}

void KalFitAlg::setCalibSvc_init()
 {
  IMdcCalibFunSvc* imdcCalibSvc; 
  MsgStream log(msgSvc(), name());
  StatusCode sc = service ("MdcCalibFunSvc", imdcCalibSvc);
  m_mdcCalibFunSvc_ = dynamic_cast<MdcCalibFunSvc*>(imdcCalibSvc);
  if ( sc.isFailure() ){
    log << MSG::FATAL << "Could not load MdcCalibFunSvc!" << endreq;
  }
 KalFitTrack::setMdcCalibFunSvc( m_mdcCalibFunSvc_);
 }

void KalFitAlg::setGeomSvc_init()
 {
  IMdcGeomSvc*  imdcGeomSvc;
  MsgStream log(msgSvc(), name());
  StatusCode sc = service ("MdcGeomSvc", imdcGeomSvc);
//  m_mdcGeomSvc_ = dynamic_cast<MdcGeomSvc*>(imdcGeomSvc);
  if ( sc.isFailure() ){
    log << MSG::FATAL << "Could not load MdcGeomSvc!" << endreq;
  }
 imdcGeomSvc_ = imdcGeomSvc;

 KalFitTrack::setIMdcGeomSvc( imdcGeomSvc);
 }

 /*void KalFitAlg::getEventStarTime()
 {
 double t0=0.;
  MsgStream log(msgSvc(), name());
 // t_t0 = -1;
 // t_t0Stat = -1;
  SmartDataPtr<EvTimeCol> evtimeCol(eventSvc(),"/Event/Recon/EvTimeCol");
  if (evtimeCol) {
    EvTimeCol::iterator iter_evt = evtimeCol->begin();
    t0 =  (*iter_evt)->getTest()*1.e-9;
 // t_t0 = (*iter_evt)->getTest();
 // t_t0Stat = (*iter_evt)->getStat();
  }else{
    log << MSG::WARNING << "Could not find EvTimeCol" << endreq;
  }

 KalFitTrack::setT0(t0);

 }
*/

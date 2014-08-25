//-----------------------------------------------------------------------
// Description : Main file of the module KalFit in charge of :
// 1/ Refit of the Mdc tracks using Kalman filter
// 2/ Backward filter (smoothing)
// 3/ and also several mass hypothesis, multiple scattering, energy loss,
//    non unif mag field treatment, wire sag effect...
//------------------------------------------------------------------------
// Modif :
//------------------------------------------------------------------------
#include <cstdio>
#include <fstream>
#include <string.h>
#include <map>
#include <vector>
#include <algorithm>
#include "CLHEP/Geometry/Vector3D.h"
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Matrix.h"
#include "GaudiKernel/MsgStream.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/IMessageSvc.h"
#include "GaudiKernel/IDataManagerSvc.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/IJobOptionsSvc.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/StatusCode.h"
#include "GaudiKernel/ContainedObject.h"
#include "GaudiKernel/SmartRef.h"
#include "GaudiKernel/SmartRefVector.h"
#include "GaudiKernel/ObjectVector.h"

//#include "TrackUtil/Helix.h"
#include "KalFitAlg/helix/Helix.h"

#include "KalFitAlg/lpav/Lpav.h"
// Bfield:
#include "KalFitAlg/coil/Bfield.h"
#include "KalFitAlg/KalFitTrack.h"
#include "KalFitAlg/KalFitHitMdc.h"
#include "KalFitAlg/KalFitHelixSeg.h"
#include "KalFitAlg/KalFitElement.h"
#include "KalFitAlg/KalFitAlg.h"
#include "McTruth/McParticle.h"
#include "McTruth/MdcMcHit.h" 
#include "EventModel/EventHeader.h"
#include "EvTimeEvent/RecEsTime.h"
#include "ReconEvent/ReconEvent.h"
#include "MdcRawEvent/MdcDigi.h"
#include "MdcRecEvent/RecMdcHit.h"
#include "MdcRecEvent/RecMdcTrack.h"    
#include "MdcRecEvent/RecMdcKalHelixSeg.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "MdcGeomSvc/MdcGeomSvc.h"
#include "MagneticField/IMagneticFieldSvc.h"
#include "MagneticField/MagneticFieldSvc.h"

#include "VertexFit/IVertexDbSvc.h"

#include "Identifier/Identifier.h"                                         
#include "Identifier/MdcID.h"       
#include "GaudiKernel/IPartPropSvc.h"
#include "GaudiKernel/INTupleSvc.h"
#include "KalFitAlg/KalFitHistItem.h"
#include "MdcPrintSvc/IMdcPrintSvc.h"

using CLHEP::HepVector; 
using CLHEP::Hep3Vector;
using CLHEP::HepMatrix;
using CLHEP::HepSymMatrix;

using namespace Event;
using namespace KalmanFit;



// Radius of the inner  wall of mdc
const double KalFitAlg::RIW = 6.35;

/// Constructor
KalFitAlg::KalFitAlg(const std::string& name, ISvcLocator* pSvcLocator):
  Algorithm(name, pSvcLocator), m_mdcCalibFunSvc_(0),m_MFSvc_(0),
  _wire(0), _layer(0), _superLayer(0),
  pathl_(1), wsag_(4), back_(1), pT_(0.0), lead_(2), mhyp_(31),
  pe_cut_(2.0), pmu_cut_(2.0), ppi_cut_(2.0), pk_cut_(2.0), pp_cut_(2.0),   
  muls_(1), loss_(1), enhance_(0),drifttime_choice_(0),driftdist_choice_(0),choice_(0),
  fac_h1_(1),fac_h2_(1),fac_h3_(1),fac_h4_(1),fac_h5_(1),
  i_back_(-1), debug_kft_(0), debug_(0), ntuple_(0),eventno(-1),
  Tds_back_no(0),m_nt1(0),m_nt2(0),m_nt3(0),m_nt4(0),m_nt5(0),
  iqual_back_(1),tprop_(1),m_usevtxdb(0),m_dangcut(10),m_dphicut(10),
  nTotalTrks(0)
{
  declareProperty("dchi2cut_all",dchi2cut_all = -1);
  declareProperty("dchi2cut_fltLen",dchi2cut_fltLen);
  declareProperty("gain1",gain1_ = 1.);
  declareProperty("gain2",gain2_ = 1.);
  declareProperty("gain3",gain3_ = 1.);
  declareProperty("gain4",gain4_ = 1.);
  declareProperty("gain5",gain5_ = 1.);
  declareProperty("fitnocut",fitnocut_ = 5);
  declareProperty("inner_steps",inner_steps_ = 3);
  declareProperty("outer_steps",outer_steps_ = 3);
  declareProperty("choice",choice_ = 0);
  declareProperty("numfcor",numfcor_ = 0);
  declareProperty("numf",numf_ = 0);
  declareProperty("steplev",steplev_ = 2);
  declareProperty("usage",usage_=0);
  declareProperty("i_front",i_front_=1);
  declareProperty("lead",lead_=2);
  declareProperty("mhyp",mhyp_=31);//yzhang add
  declareProperty("muls",muls_=1);
  declareProperty("loss",loss_=1);
  declareProperty("matrixg",matrixg_=100.0);
  declareProperty("lr",lr_=1);  
  declareProperty("debug_kft",debug_kft_=0);
  declareProperty("debug",debug_=0);  
  declareProperty("ntuple",ntuple_=0);  
  declareProperty("activeonly",activeonly_=0);  
  declareProperty("matfile",matfile_="geomdc_material.dat"); 
  declareProperty("cylfile",cylfile_="geomdc_cylinder.dat"); 
  declareProperty("dchi2cutf",dchi2cutf_=1000.0);
  declareProperty("dchi2cuts",dchi2cuts_=1000.0);
  declareProperty("pt",pT_=0.0);
  declareProperty("pe_cut",pe_cut_=2.0);
  declareProperty("pmu_cut",pmu_cut_=2.0);
  declareProperty("ppi_cut",ppi_cut_=2.0);
  declareProperty("pk_cut",pk_cut_=2.0);
  declareProperty("pp_cut",pp_cut_=2.0); 
  declareProperty("wsag",wsag_=4);
  declareProperty("back",back_=1);
  declareProperty("i_back",i_back_=1);
  declareProperty("tofflag",tofflag_=1);
  declareProperty("tofflagFltlenFromRec",tofflag_fltLen_=0);
  declareProperty("tof_hyp",tof_hyp_=1);
  declareProperty("resolution",resolution_=1);
  declareProperty("fstrag",fstrag_=0.9);
  declareProperty("drifttime_choice",drifttime_choice_=0);
  declareProperty("driftdist_choice",driftdist_choice_=0);
  declareProperty("tprop",tprop_=1);
  declareProperty("pt_cut",pt_cut_= 0.2);
  declareProperty("theta_cut",theta_cut_= 0.8);
  declareProperty("usevtxdb",m_usevtxdb= 0);
  declareProperty("cosmicflag",m_csmflag= 0);
  declareProperty("dangcut",m_dangcut=10.);
  declareProperty("dphicut",m_dphicut=10.);
  declareProperty("testCDC", testCDC=1);
  declareProperty("iniTrkPar", iniTrkPar=1);
  declareProperty("nPivotCnv", m_nPivotCnv=2);//yzhang 2014-04-11 
  declareProperty("pivotCnvSelector", m_pivotCnvSelector=1);//yzhang 2014-04-11 
  declareProperty("fix2TurnTrackInnerWall", m_fix2TurnTrackInnerWall=0);//yzhang 2014-05-02
  declareProperty("fixInnerWall", m_fixInnerWall=true);//yzhang 2014-05-03 
  declareProperty("t0", m_t0=0);//yzhang 2014-05-30 
  for(int i=0; i<5; i++) nFailedTrks[i]=0;
}

// destructor
KalFitAlg::~KalFitAlg(void)
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO <<" Start cleaning, delete  Mdc geometry objects" << endreq;
  clean();
  log << MSG::INFO << "End cleaning " << endreq;
}

void KalFitAlg::clean(void)
{
  // delete all Mdc objects :
  _BesKalmanFitWalls.clear();
  _BesKalmanFitMaterials.clear();
  if(_wire)free(_wire);
  if(_layer)free(_layer);
  if(_superLayer)free(_superLayer);
}

// initialization
StatusCode KalFitAlg::initialize()
{
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in initize()" 
    << "KalFit> Initialization for current run " << endreq;
  log << MSG::INFO << "Present Parameters: muls: " << muls_ <<"  loss:  "
    << loss_ <<" matrixg "<< matrixg_ <<" lr "<< lr_ 
    << " debug " << debug_ << " ntuple " << ntuple_
    << " activeonly "<< activeonly_ << endreq;

  KalFitTrack::LR(lr_);
  KalFitTrack::resol(resolution_);
  KalFitTrack::Tof_correc_ = tofflag_;
  KalFitTrack::Tof_correc_fltLenFromRec_ = tofflag_fltLen_;
  KalFitTrack::tofall_ = tof_hyp_;
  KalFitTrack::chi2_hitf_ = dchi2cutf_;
  KalFitTrack::chi2_hits_ = dchi2cuts_;
  KalFitTrack::factor_strag_ = fstrag_;
  KalFitTrack::debug_ = debug_kft_;
  KalFitTrack::drifttime_choice_ = drifttime_choice_;
  KalFitTrack::driftdist_choice_ = driftdist_choice_;
  KalFitTrack::steplev_ = steplev_;
  KalFitTrack::numfcor_ = numfcor_;
  KalFitTrack::inner_steps_ = inner_steps_;
  KalFitTrack::outer_steps_ = outer_steps_;
  KalFitTrack::tprop_ = tprop_;

  setDchisqCut();
  // Delete properly the geometry objects if already existing
  clean();
  log << MSG::INFO << ".....building Mdc " << endreq;

  // Set objects and material properties for Mdc detector :
  //setMaterial_Mdc();
  //setCylinder_Mdc();

  setBesFromGdml();
  // initialize the MdcCalibFunSvc 
  setCalibSvc_init(); 
  //
  //  // initialize the MdcGeomSvc
  //  setGeomSvc_init();
  //getEventStarTime();
  //  // Wires, Layers and SuperLayers of Mdc :
  //  set_Mdc();
  // define histograms and Ntuples
  hist_def();


  IMagneticFieldSvc* IMFSvc; 
  StatusCode sc = service ("MagneticFieldSvc",IMFSvc); 
  if(sc!=StatusCode::SUCCESS) { 
    log << MSG::ERROR << "Unable to open Magnetic field service"<<endreq; 
  }
  KalFitTrack::setMagneticFieldSvc(IMFSvc);

  // Nominal magnetic field :
  if (KalFitTrack::numfcor_){
    KalFitTrack::Bznom_ = (IMFSvc->getReferField())*10000; //unit is KGauss
    if(0 == KalFitTrack::Bznom_)   KalFitTrack::Bznom_ = -10;

    if(4 == debug_){
      std::cout<<" initialize, referField from MagneticFieldSvc: "<< (IMFSvc->getReferField())*10000 <<std::endl;
      std::cout<<" magnetic field: "<<KalFitTrack::Bznom_<<std::endl;
    }

  }

  // Print out of the status of the flags :
  if (ntuple_) 
    log << MSG::INFO <<" ntuple out, the option is  "<< ntuple_ <<endreq;
  if (debug_ >0 ) {
    cout << "KalFitAlg> DEBUG open,Here is the important Parameters :\n";
    cout << " Leading particule with mass hyp = " << lead_ << std::endl;
    cout << " mhyp = " << mhyp_ << std::endl;
    cout << "===== Effects taking into account : " << std::endl;
    cout << " - multiple scattering = " << muls_ << std::endl;
    cout << " - energy loss = " << loss_ << std::endl;
    if (KalFitTrack::strag_)
      cout << " - straggling for the energy loss " << std::endl;
    cout << " - nominal Bz value = " << KalFitTrack::Bznom_ << std::endl;

    if (KalFitTrack::numf_ > 19)
      cout << " - non uniform magnetic field treatment " 
	<< KalFitTrack::numf_ << std::endl;
    cout << " - wire sag correction = " << wsag_ << std::endl;
    cout << " - Tof correction = " << KalFitTrack::Tof_correc_ << std::endl;
    cout << "   Tof correction fltLen = " << KalFitTrack::Tof_correc_fltLenFromRec_ << std::endl;
    cout << " - chi2 cut for each hit = " << KalFitTrack::chi2_hitf_ 
      << std::endl << " is applied after " << KalFitTrack::nmdc_hit2_ 
      << " hits included " << std::endl;

    if (back_){
      cout << " Backward filter is on with a pT cut value = " << pT_ << endl;
    }
    if(debug_ == 4) cout << " pathl = " << pathl_ << std::endl;

    if (KalFitTrack::LR_==1)
      cout << " Decision L/R from MdcRecHit " << std::endl;
  }

  KalFitElement::muls(muls_);
  KalFitElement::loss(loss_);
  KalFitTrack::lead(lead_);
  KalFitTrack::back(back_);
  KalFitTrack::numf(numf_);
  // Get the Particle Properties Service
  IPartPropSvc* p_PartPropSvc;
  static const bool CREATEIFNOTTHERE(true);
  StatusCode PartPropStatus = service("PartPropSvc", p_PartPropSvc, CREATEIFNOTTHERE);
  if (!PartPropStatus.isSuccess() || 0 == p_PartPropSvc) {
    log << MSG::WARNING << " Could not initialize Particle Properties Service" << endreq;
    return StatusCode::SUCCESS; 
  }
  m_particleTable = p_PartPropSvc->PDT();

  //Initailize MdcPrintSvc
  IMdcPrintSvc* imdcPrintSvc;
  sc = service ("MdcPrintSvc", imdcPrintSvc);
  m_mdcPrintSvc = dynamic_cast<MdcPrintSvc*> (imdcPrintSvc);
  if ( sc.isFailure() ){
    log << MSG::FATAL << "Could not load MdcPrintSvc!" << endreq;
    return StatusCode::FAILURE;
  }

  return StatusCode::SUCCESS; 
}

StatusCode KalFitAlg::finalize()
{
  MsgStream log(msgSvc(), name());
  log << MSG::DEBUG<<"KalFitAlg:: nTotalTrks = "<<nTotalTrks<<endreq;
  log << MSG::DEBUG<<"        e: "<<nFailedTrks[0]<<" failed, "<<nTotalTrks-nFailedTrks[0]<<" successed"<<endreq;
  log << MSG::DEBUG<<"       mu: "<<nFailedTrks[1]<<" failed, "<<nTotalTrks-nFailedTrks[1]<<" successed"<<endreq;
  log << MSG::DEBUG<<"       pi: "<<nFailedTrks[2]<<" failed, "<<nTotalTrks-nFailedTrks[2]<<" successed"<<endreq;
  log << MSG::DEBUG<<"        K: "<<nFailedTrks[3]<<" failed, "<<nTotalTrks-nFailedTrks[3]<<" successed"<<endreq;
  log << MSG::DEBUG<<"        p: "<<nFailedTrks[4]<<" failed, "<<nTotalTrks-nFailedTrks[4]<<" successed"<<endreq;

  cout<<"KalFitAlg:: nTotalTrks = "<<nTotalTrks<<endl;
  cout<<"        e: "<<nFailedTrks[0]<<" failed, "<<nTotalTrks-nFailedTrks[0]<<" successed"<<endl;
  cout<<"       mu: "<<nFailedTrks[1]<<" failed, "<<nTotalTrks-nFailedTrks[1]<<" successed"<<endl;
  cout<<"       pi: "<<nFailedTrks[2]<<" failed, "<<nTotalTrks-nFailedTrks[2]<<" successed"<<endl;
  cout<<"        K: "<<nFailedTrks[3]<<" failed, "<<nTotalTrks-nFailedTrks[3]<<" successed"<<endl;
  cout<<"        p: "<<nFailedTrks[4]<<" failed, "<<nTotalTrks-nFailedTrks[4]<<" successed"<<endl;

  return StatusCode::SUCCESS;                                                  
}                       

// begin run setting
StatusCode KalFitAlg::beginRun( )
{
  MsgStream log(msgSvc(), name());  
  log << MSG::INFO << "in beginRun()" << endreq;
  log << MSG::INFO << "Present Parameters: muls: " << muls_ <<"  loss:  "
    << " activeonly "<< activeonly_ << endreq;

  // initialize the MdcGeomSvc
  setGeomSvc_init();
  // Wires, Layers and SuperLayers of Mdc :
  set_Mdc();

  IMagneticFieldSvc* IMFSvc;
  StatusCode sc = service ("MagneticFieldSvc",IMFSvc);
  if(sc!=StatusCode::SUCCESS) {
    log << MSG::ERROR << "Unable to open Magnetic field service"<<endreq;
  }

  // Nominal magnetic field :
  if (KalFitTrack::numfcor_){
    KalFitTrack::Bznom_ = (IMFSvc->getReferField())*10000; //unit is KGauss
    if(0 == KalFitTrack::Bznom_)   KalFitTrack::Bznom_ = -10;

    if(4 == debug_){  
      std::cout<<" beginRun, referField from MagneticFieldSvc:"<< (IMFSvc->getReferField())*10000 <<std::endl;
      std::cout<<" magnetic field: "<<KalFitTrack::Bznom_<<std::endl;
    }
  }

  return StatusCode::SUCCESS;  
}




// hist_def function
void KalFitAlg::hist_def ( void )
{
  if(ntuple_&1) {               
    NTuplePtr  ntTest(ntupleSvc(),"FILE104/test");
    StatusCode status;
    if ( ntTest ) m_ntTest = ntTest;                                                    
    else {                                                                     
      m_ntTest= ntupleSvc()->book("FILE104/test",CLID_ColumnWiseTuple,"Test");
      if ( m_ntTest )  { 
	status = m_ntTest->addItem("ddRec",m_ddRec);
	status = m_ntTest->addItem("ddKal",m_ddKal);
	status = m_ntTest->addItem("ddTruth",m_ddTruth);
	status = m_ntTest->addItem("layer",m_ddLayer);
	status = m_ntTest->addItem("wire",m_ddWire);
	status = m_ntTest->addItem("fiRec",m_fiRec);
	status = m_ntTest->addItem("fiKal",m_fiKal);
	status = m_ntTest->addItem("fiTerm",m_ddFiTerm);
	status = m_ntTest->addItem("fltLenTerm",m_ddFltLenTerm);
	status = m_ntTest->addItem("fltLen",m_ddFltlen);
	status = m_ntTest->addItem("fltLenKal",m_ddFltlenKal);
	status = m_ntTest->addItem("dchi2R",m_dchi2R);
	status = m_ntTest->addItem("dchi2L",m_dchi2L);
	status = m_ntTest->addItem("dchi2min",m_dchi2min);
	status = m_ntTest->addItem("iTurn",m_iTurn);
	status = m_ntTest->addItem("iFit",m_iFit);
	status = m_ntTest->addItem("mass",m_ddMass);
	status = m_ntTest->addItem("chi2Sum",m_chi2Sum);
	status = m_ntTest->addItem("chi2Add",m_chi2Add);
	status = m_ntTest->addItem("dp",m_dp);
	status = m_ntTest->addItem("lastdp",m_lastdp);
	status = m_ntTest->addItem("pCurrent",m_pCurrent);
	status = m_ntTest->addItem("pTruthCurrent",m_pTruthCurrent);
      }
    }

    NTuplePtr  nt1(ntupleSvc(),"FILE104/n101");
    if ( nt1 ) m_nt1 = nt1;                                                    
    else {                                                                     
      m_nt1= ntupleSvc()->book("FILE104/n101",CLID_ColumnWiseTuple,"KalFit");
      if ( m_nt1 )  { 

	status = m_nt1->addItem("trackid",m_trackid);
	status = m_nt1->addItem("stat",5,2,m_stat);
	status = m_nt1->addItem("ndf",5,2,m_ndf);
	status = m_nt1->addItem("chisq",5,2,m_chisq);
	status = m_nt1->addItem("length",5,m_length);
	status = m_nt1->addItem("tof",5,m_tof);
	status = m_nt1->addItem("nhits",5,m_nhits);
	status = m_nt1->addItem("zhelix",5,m_zhelix);
	status = m_nt1->addItem("zhelixe",5,m_zhelixe);
	status = m_nt1->addItem("zhelixmu",5,m_zhelixmu);
	status = m_nt1->addItem("zhelixk",5,m_zhelixk);
	status = m_nt1->addItem("zhelixp",5,m_zhelixp);
	status = m_nt1->addItem("zptot",m_zptot);
	status = m_nt1->addItem("zptote",m_zptote);
	status = m_nt1->addItem("zptotmu",m_zptotmu);
	status = m_nt1->addItem("zptotk",m_zptotk);
	status = m_nt1->addItem("zptotp",m_zptotp);

	status = m_nt1->addItem("zpt",m_zpt);
	status = m_nt1->addItem("zpte",m_zpte);
	status = m_nt1->addItem("zptmu",m_zptmu);
	status = m_nt1->addItem("zptk",m_zptk);
	status = m_nt1->addItem("zptp",m_zptp);

	status = m_nt1->addItem("fptot",m_fptot);
	status = m_nt1->addItem("fptote",m_fptote);
	status = m_nt1->addItem("fptotmu",m_fptotmu);
	status = m_nt1->addItem("fptotk",m_fptotk);
	status = m_nt1->addItem("fptotp",m_fptotp);
	status = m_nt1->addItem("fpt",m_fpt);
	status = m_nt1->addItem("fpte",m_fpte);
	status = m_nt1->addItem("fptmu",m_fptmu);
	status = m_nt1->addItem("fptk",m_fptk);
	status = m_nt1->addItem("fptp",m_fptp);

	status = m_nt1->addItem("lptot",m_lptot);
	status = m_nt1->addItem("lptote",m_lptote);
	status = m_nt1->addItem("lptotmu",m_lptotmu);
	status = m_nt1->addItem("lptotk",m_lptotk);
	status = m_nt1->addItem("lptotp",m_lptotp);
	status = m_nt1->addItem("lpt",m_lpt);
	status = m_nt1->addItem("lpte",m_lpte);
	status = m_nt1->addItem("lptmu",m_lptmu);
	status = m_nt1->addItem("lptk",m_lptk);
	status = m_nt1->addItem("lptp",m_lptp);

	status = m_nt1->addItem("zsigp",m_zsigp);
	status = m_nt1->addItem("zsigpe",m_zsigpe);
	status = m_nt1->addItem("zsigpmu",m_zsigpmu);
	status = m_nt1->addItem("zsigpk",m_zsigpk);
	status = m_nt1->addItem("zsigpp",m_zsigpp);
	status = m_nt1->addItem("fhelix",5,m_fhelix);
	status = m_nt1->addItem("fhelixe",5,m_fhelixe);
	status = m_nt1->addItem("fhelixmu",5,m_fhelixmu);
	status = m_nt1->addItem("fhelixk",5,m_fhelixk);
	status = m_nt1->addItem("fhelixp",5,m_fhelixp);
	status = m_nt1->addItem("lhelix",5,m_lhelix);
	status = m_nt1->addItem("lhelixe",5,m_lhelixe);
	status = m_nt1->addItem("lhelixmu",5,m_lhelixmu);
	status = m_nt1->addItem("lhelixk",5,m_lhelixk);
	status = m_nt1->addItem("lhelixp",5,m_lhelixp);
	if(ntuple_&32) {
	  status = m_nt1->addItem("zerror",15,m_zerror);
	  status = m_nt1->addItem("zerrore",15,m_zerrore);
	  status = m_nt1->addItem("zerrormu",15,m_zerrormu);
	  status = m_nt1->addItem("zerrork",15,m_zerrork);
	  status = m_nt1->addItem("zerrorp",15,m_zerrorp);
	  status = m_nt1->addItem("ferror",15,m_ferror);
	  status = m_nt1->addItem("ferrore",15,m_ferrore);
	  status = m_nt1->addItem("ferrormu",15,m_ferrormu);
	  status = m_nt1->addItem("ferrork",15,m_ferrork);
	  status = m_nt1->addItem("ferrorp",15,m_ferrorp);
	  status = m_nt1->addItem("lerror",15,m_lerror);
	  status = m_nt1->addItem("lerrore",15,m_lerrore);
	  status = m_nt1->addItem("lerrormu",15,m_lerrormu);
	  status = m_nt1->addItem("lerrork",15,m_lerrork);
	  status = m_nt1->addItem("lerrorp",15,m_lerrorp);
	}
	if((ntuple_&16)&&(ntuple_&1)) {
	  status = m_nt1->addItem("evtid",m_evtid);
	  status = m_nt1->addItem("mchelix",5,m_mchelix);
	  status = m_nt1->addItem("mcptot",m_mcptot);
	  status = m_nt1->addItem("mcpid",m_mcpid);
	}
	if( status.isFailure() ) cout<<"Ntuple1 add item failed!"<<endl;       
      }                                                                        
    }
  }

  if(ntuple_&4) {  
    NTuplePtr  nt2(ntupleSvc(),"FILE104/n102");
    StatusCode status2;
    if ( nt2 ) m_nt2 = nt2;                                                    
    else {                                                                     
      m_nt2= ntupleSvc()->book("FILE104/n102",CLID_ColumnWiseTuple,"KalFitComp");
      if ( m_nt2 )  {
	status2 = m_nt2->addItem("delx",m_delx);
	status2 = m_nt2->addItem("dely",m_dely);
	status2 = m_nt2->addItem("delz",m_delz);
	status2 = m_nt2->addItem("delthe",m_delthe);
	status2 = m_nt2->addItem("delphi",m_delphi);
	status2 = m_nt2->addItem("delp",m_delp);                         
	status2 = m_nt2->addItem("delpx",m_delpx);
	status2 = m_nt2->addItem("delpy",m_delpy);
	status2 = m_nt2->addItem("delpz",m_delpz);

	if( status2.isFailure() ) cout<<"Ntuple2 add item failed!"<<endl; 
      }
    }
  }

  if(ntuple_&2) {                                            
    NTuplePtr  nt3(ntupleSvc(),"FILE104/n103");
    StatusCode status3;
    if ( nt3 ) m_nt3 = nt3;                                                    
    else {                                                                     
      m_nt3= ntupleSvc()->book("FILE104/n103",CLID_ColumnWiseTuple,"PatRec");
      if ( m_nt3 )  {
	status3 = m_nt3->addItem("trkhelix",5,m_trkhelix);
	status3 = m_nt3->addItem("trkptot",m_trkptot);
	if(ntuple_&32) {
	  status3 = m_nt3->addItem("trkerror",15,m_trkerror);
	  status3 = m_nt3->addItem("trksigp",m_trksigp);
	}
	status3 = m_nt3->addItem("trkndf",m_trkndf);
	status3 = m_nt3->addItem("trkchisq",m_trkchisq);
	if( status3.isFailure() ) cout<<"Ntuple3 add item failed!"<<endl; 
      }
    }
  }  

  if(ntuple_&4) {                                            
    NTuplePtr  nt4(ntupleSvc(),"FILE104/n104");
    StatusCode status4;
    if ( nt4 ) m_nt4 = nt4;                                                    
    else {                                                                     
      m_nt4= ntupleSvc()->book("FILE104/n104",CLID_ColumnWiseTuple,"PatRecComp");
      if ( m_nt4 )  {
	status4 = m_nt4->addItem("trkdelx",m_trkdelx);
	status4 = m_nt4->addItem("trkdely",m_trkdely);
	status4 = m_nt4->addItem("trkdelz",m_trkdelz);
	status4 = m_nt4->addItem("trkdelthe",m_trkdelthe);
	status4 = m_nt4->addItem("trkdelphi",m_trkdelphi);
	status4 = m_nt4->addItem("trkdelp",m_trkdelp);                         
	if( status4.isFailure() ) cout<<"Ntuple4 add item failed!"<<endl; 
      }
    }
  }
  if(ntuple_&8) {                                            
    NTuplePtr  nt5(ntupleSvc(), "FILE104/n105");
    StatusCode status5;
    if ( nt5 ) m_nt5 = nt5;                                                    
    else {                                                                     
      m_nt5= ntupleSvc()->book("FILE104/n105",CLID_ColumnWiseTuple,"KalFitdChisq");
      if ( m_nt5 )  {
	status5 = m_nt5->addItem("dchi2",m_dchi2);
	status5 = m_nt5->addItem("masshyp",m_masshyp);
	status5 = m_nt5->addItem("residual_estim",m_residest);
	status5 = m_nt5->addItem("residual",m_residnew);
	status5 = m_nt5->addItem("layer",m_layer);
	status5 = m_nt5->addItem("kaldr",m_anal_dr);
	status5 = m_nt5->addItem("kalphi0",m_anal_phi0);
	status5 = m_nt5->addItem("kalkappa",m_anal_kappa);
	status5 = m_nt5->addItem("kaldz",m_anal_dz);
	status5 = m_nt5->addItem("kaltanl",m_anal_tanl);
	status5 = m_nt5->addItem("dr_ea",m_anal_ea_dr);
	status5 = m_nt5->addItem("phi0_ea",m_anal_ea_phi0);
	status5 = m_nt5->addItem("kappa_ea",m_anal_ea_kappa);
	status5 = m_nt5->addItem("dz_ea",m_anal_ea_dz);
	status5 = m_nt5->addItem("tanl_ea",m_anal_ea_tanl);
	if( status5.isFailure() ) cout<<"Ntuple5 add item failed!"<<endl; 
      }
    }
    NTuplePtr  nt6(ntupleSvc(),"FILE104/n106");
    StatusCode status6;
    if ( nt6 ) m_nt6 = nt6;             
    else {                            
      m_nt6= ntupleSvc()->book("FILE104/n106",CLID_ColumnWiseTuple,"kal seg");
      if ( m_nt6 )  {
	status6 = m_nt6->addItem("docaInc",m_docaInc);
	status6 = m_nt6->addItem("docaExc",m_docaExc);
	status6 = m_nt6->addItem("tdr",m_tdrift);
	status6 = m_nt6->addItem("layerid", m_layerid);
	status6 = m_nt6->addItem("event", m_eventNo);
	status6 = m_nt6->addItem("residualInc",m_residualInc);
	status6 = m_nt6->addItem("residualExc",m_residualExc);
	status6 = m_nt6->addItem("lr",m_lr);
	status6 = m_nt6->addItem("dd",m_dd);
	status6 = m_nt6->addItem("y",m_yposition);

	if( status6.isFailure() ) cout<<"Ntuple6 add item failed!"<<endl;
      }
    }
  }
}



void KalFitAlg::setDchisqCut()
{
  int layid = 0;
  if(dchi2cut_all>0) {
    for(layid = 0; layid<43; layid++) {
      KalFitTrack::dchi2cutf_anal[layid] = dchi2cut_all;
      KalFitTrack::dchi2cutf_calib[layid] = dchi2cut_all;
    }
    return;
  }
}

StatusCode KalFitAlg::execute() {
  MsgStream log(msgSvc(), name());
  log << MSG::INFO << "in execute()" << endreq;
  //std::cout<<"begin to deal with EVENT  ..."<<(++eventno)<<std::endl;
  for(int i=0; i<5; i++) iqual_front_[i] = 1;
  iqual_back_ = 1;

  /// 
  IMagneticFieldSvc* IMFSvc; 
  StatusCode sc = service ("MagneticFieldSvc",IMFSvc);
  if(sc!=StatusCode::SUCCESS) {
    log << MSG::ERROR << "Unable to open Magnetic field service"<<endreq; 
  }

  // Nominal magnetic field : 
  if (KalFitTrack::numfcor_){
    KalFitTrack::Bznom_ = (IMFSvc->getReferField())*10000; //unit is KGauss
    if(0 == KalFitTrack::Bznom_)   KalFitTrack::Bznom_ = -10;

    if(4 == debug_){
      std::cout<<" execute, referField from MagneticFieldSvc: "<< (IMFSvc->getReferField())*10000 <<std::endl;
      std::cout<<" magnetic field: "<<KalFitTrack::Bznom_<<std::endl;
    }
  }


  MdcGeomSvc* const geosvc = dynamic_cast<MdcGeomSvc*>(imdcGeomSvc_);
  if(!geosvc) {
    std::cout<<"ERROR OCCUR when dynamic_cast in KalFitAlg::execute ...!!"<<std::endl; 
  }

  SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");
  if (!eventHeader) {
    log << MSG::WARNING << "Could not find Event Header" << endreq;
    return StatusCode::FAILURE;
  }
  eventNo = eventHeader->eventNumber();
  runNo = eventHeader->runNumber();
  if(runNo>0) wsag_=4;//yzhang FIXME for data
  else wsag_=0;

  double t0=m_t0;//for chargedgeantino 2014-05-30 TEMP
  //double t0=0.;
  //SmartDataPtr<RecEsTimeCol> estimeCol(eventSvc(),"/Event/Recon/RecEsTimeCol");
  //if (estimeCol && estimeCol->size()) {
  //	RecEsTimeCol::iterator iter_evt = estimeCol->begin();
  //	t0 =  (*iter_evt)->getTest();
  //	// t0Stat = (*iter_evt)->getStat();
  //}else{
  //	log << MSG::WARNING << "Could not find EvTimeCol" << endreq;
  //	return StatusCode::SUCCESS;
  //}


  if(debug_==4) {
    std::cout<<"in KalFitAlg , we get the event start time = "<<t0<<std::endl;
  } 
  KalFitTrack::setEventStartTime(t0);

  SmartDataPtr<MdcDigiCol> mdcDigiCol(eventSvc(),"/Event/Digi/MdcDigiCol");
  if (sc!=StatusCode::SUCCESS) {
    log << MSG::FATAL << "Could not find MdcDigiCol!" << endreq;
    return StatusCode::SUCCESS;
  }
  KalFitTrack::setMdcDigiCol(mdcDigiCol);

  //yzhang debug 2014-07-10 
  for(int iii = 0; iii<20;iii++){ for(int jjj = 0; jjj<300;jjj++){ pTruthHit[iii][jjj]=-999;}}
  //yzhang 2014-07-10 for debug
  SmartDataPtr<Event::MdcMcHitCol> mdcMcHitCol(eventSvc(),"/Event/MC/MdcMcHitCol");
  if (!mdcMcHitCol) {
    log << MSG::WARNING << "Could not find event" << endreq;
  }    
  Event::MdcMcHitCol::iterator iterMcHit = mdcMcHitCol->begin();
  for(; iterMcHit != mdcMcHitCol->end(); ++iterMcHit) {
    if((*iterMcHit)->getPositionFlag()!=-999) continue;
    Identifier mdcid = (*iterMcHit)->identify();
    int layerId = MdcID::layer(mdcid);
    int wireId = MdcID::wire(mdcid);
    double p = sqrt((*iterMcHit)->getPositionX()*(*iterMcHit)->getPositionX() + (*iterMcHit)->getPositionY()*(*iterMcHit)->getPositionY()+(*iterMcHit)->getPositionZ()*(*iterMcHit)->getPositionZ());
    pTruthHit[layerId][wireId] = p;
    //double pt = sqrt((*iterMcHit)->getPositionX()*(*iterMcHit)->getPositionX() + (*iterMcHit)->getPositionY()*(*iterMcHit)->getPositionY());
  }

  int nHits=0; int nSt=0;
  double maxTrackLength = -999;
  Hep3Vector trueMomFinal;
  HepPoint3D positionFinal;
  HitRefVec  hitRefVec;

  // register RecMdcTrack and MdcRecHit collection 
  bool mcstat = false;//true;

  SmartDataPtr<McParticleCol> mcPartCol(eventSvc(),"/Event/MC/McParticleCol");
  if (!mcPartCol) {
    log << MSG::WARNING << "Could not find McParticle" << endreq;
    mcstat = false;
  }

  if(mcstat) {
    McParticleCol::iterator i_mcTrk = mcPartCol->begin();
    int i_mc=0;
    for (;i_mcTrk != mcPartCol->end(); i_mcTrk++) {
      if(!(*i_mcTrk)->primaryParticle()) continue;
      const HepLorentzVector& mom((*i_mcTrk)->initialFourMomentum());
      const HepLorentzVector& pos = (*i_mcTrk)->initialPosition();
      log << MSG::DEBUG << "MCINFO:particleId=" << (*i_mcTrk)->particleProperty()
	<< " theta=" << mom.theta() <<" phi="<< mom.phi()
	<<" px="<< mom.px() <<" py="<< mom.py() <<" pz="<< mom.pz()
	<< endreq;
      double charge = 0.0;
      int pid = (*i_mcTrk)->particleProperty();
      int mother_id = ((*i_mcTrk)->mother()).particleProperty();
      if( pid >0 ) { 
	charge = m_particleTable->particle( pid )->charge();
      } else if ( pid <0 ) {
	charge = m_particleTable->particle( -pid )->charge();
	charge *= -1;
      } else {
	log << MSG::WARNING << "wrong particle id, please check data" <<endreq;
      }
      HepPoint3D pos2(pos.x(),pos.y(),pos.z());
      Hep3Vector mom2(mom.px(),mom.py(),mom.pz());
      if(debug_ == 4) {
	cout<<endl<<"MC: "<<i_mc<<endl;
	cout<<"pos="<<pos2<<endl;
	cout<<"mom="<<mom2<<endl;
	cout<<"charge="<<charge<<endl;
      }

      Helix mchelix(pos2, mom2, charge);
      log << MSG::DEBUG << "charge of the track " << charge << endreq;
      if( debug_ == 4) cout<< "MC helix: "<<mchelix.a()<<endl;
      mchelix.pivot( HepPoint3D(0,0,0) );
      for( int j =0; j<5; j++) {
	myMCHelix[j] = mchelix.a()[j];
      }
      if(debug_ == 4) {
	cout<<"MC pid, mother_id = "<<pid<<", "<<mother_id<<endl;
	cout<<"              p4 = "<<(*i_mcTrk)->initialFourMomentum()<<endl;
	cout<<"        start from "<<(*i_mcTrk)->initialPosition()<<endl;
	cout<<"            end at "<<(*i_mcTrk)->finalPosition()<<endl;
	cout<<"            MC Helix: "<<mchelix.a()<<endl;
      }
      i_mc++;
    }	
  }



  if((ntuple_&16)&&(ntuple_&1)) {
    // McTruth infor,Retrieve MC track truth
    // bool mcstat = true;
    // more research needed ...

    m_evtid = eventHeader->eventNumber();
    bool mcstat = true;

    SmartDataPtr<McParticleCol> mcPartCol(eventSvc(),"/Event/MC/McParticleCol");
    if (!mcPartCol) {
      log << MSG::WARNING << "Could not find McParticle" << endreq;
      mcstat = false;
    }

    if(mcstat) {
      McParticleCol::iterator i_mcTrk = mcPartCol->begin();
      for (;i_mcTrk != mcPartCol->end(); i_mcTrk++) {
	if(!(*i_mcTrk)->primaryParticle()) continue;
	const HepLorentzVector& mom((*i_mcTrk)->initialFourMomentum());
	const HepLorentzVector& pos = (*i_mcTrk)->initialPosition();
	log << MSG::DEBUG << "MCINFO:particleId=" << (*i_mcTrk)->particleProperty()
	  << " theta=" << mom.theta() <<" phi="<< mom.phi()
	  <<" px="<< mom.px() <<" py="<< mom.py() <<" pz="<< mom.pz()
	  << endreq;
	double charge = 0.0;
	int pid = (*i_mcTrk)->particleProperty();
	if( pid >0 ) { 
	  charge = m_particleTable->particle( pid )->charge();
	} else if ( pid <0 ) {
	  charge = m_particleTable->particle( -pid )->charge();
	  charge *= -1;
	} else {
	  log << MSG::WARNING << "wrong particle id, please check data" <<endreq;
	}
	HepPoint3D pos2(pos.x(),pos.y(),pos.z());
	Hep3Vector mom2(mom.px(),mom.py(),mom.pz());

	Helix mchelix(pos2, mom2, charge);
	log << MSG::DEBUG << "charge of the track " << charge << endreq;
	if( debug_ == 4) cout<< "helix: "<<mchelix.a()<<endl;
	mchelix.pivot( HepPoint3D(0,0,0) );
	for( int j =0; j<5; j++) {
	  m_mchelix[j] = mchelix.a()[j];
	}
	m_mcpid = pid;
	m_mcptot = sqrt(1+pow(m_mchelix[4],2))/m_mchelix[2];
	//std::cout<<__FILE__<<"   "<<__LINE__<<" mchelix.a()  "<<mchelix.a()<< " m_mcptot "<<m_mcptot<<std::endl;
      }	
    }
  }

  Identifier mdcid;

  //retrieve  RecMdcTrackCol from TDS
  SmartDataPtr<RecMdcTrackCol> newtrkCol(eventSvc(),"/Event/Recon/RecMdcTrackCol");
  if (!newtrkCol) {
    log << MSG::FATAL << "Could not find RecMdcTrackCol" << endreq;
    return( StatusCode::SUCCESS);
  }
  log << MSG::INFO << "Begin to make MdcRecTrkCol and MdcRecWirhitCol"<<endreq;

  vector<MdcRec_trk>* mtrk_mgr = MdcRecTrkCol::getMdcRecTrkCol();
  mtrk_mgr->clear(); 
  vector<MdcRec_trk_add>* mtrkadd_mgr = MdcRecTrkAddCol::getMdcRecTrkAddCol();
  mtrkadd_mgr->clear();    
  vector<MdcRec_wirhit>* mhit_mgr = MdcRecWirhitCol::getMdcRecWirhitCol();
  mhit_mgr->clear();

  double trkx1,trkx2,trky1,trky2,trkz1,trkz2,trkthe1,trkthe2,trkphi1,trkphi2,trkp1,trkp2,trkr1,trkr2,trkkap1,trkkap2,trktanl1,trktanl2; 
  Hep3Vector csmp3[2];
  double csmphi[2];
  int status_temp=0;
  RecMdcTrackCol::iterator iter_trk = newtrkCol->begin();
  for(int kj = 1; iter_trk != newtrkCol->end(); iter_trk++,kj++) {
    if(kj<3){
      csmp3[kj-1]=(*iter_trk)->p3();
      csmphi[kj-1] = (*iter_trk)->phi();
    }
    if(ntuple_&2) {
      //check trackcol, track level
      for( int j = 0, ij = 0; j<5; j++) {
	m_trkhelix[j] = (*iter_trk)->helix()[j];
	if(ntuple_&32) {
	  for(int k=0; k<=j; k++,ij++) {
	    m_trkerror[ij] = (*iter_trk)->err()[j][k];
	  }
	}
      }
      m_trkptot = sqrt(1+pow(m_trkhelix[4],2))/m_trkhelix[2];
      if(ntuple_&32){
	m_trksigp = sqrt(pow((m_trkptot/m_trkhelix[2]),2)*m_trkerror[5]+
	    pow((m_trkhelix[4]/m_trkptot),2)*pow((1/m_trkhelix[2]),4)*m_trkerror[14]-
	    2*m_trkhelix[4]*m_trkerror[12]*pow((1/m_trkhelix[2]),3));
      }
      m_trkndf =  (*iter_trk)->ndof();
      m_trkchisq =  (*iter_trk)->chi2();

      if (debug_ == 4) cout<<"Ea from RecMdcTrackCol..." <<(*iter_trk)->err()<<endl;

      StatusCode sc3 = m_nt3->write();
      if( sc3.isFailure() ) cout<<"Ntuple3 filling failed!"<<endl;
    }
    //end of track level check and prepare evt check
    //end prepare	

    log << MSG::DEBUG << "retrieved MDC tracks:"
      << "   Nhits " <<(*iter_trk)->getNhits()
      << "   Nster " <<(*iter_trk)->nster() <<endreq;
    // so ,use this to get the hits vector belong to this track ... 
    HitRefVec gothits = (*iter_trk)->getVecHits();

    MdcRec_trk* rectrk = new MdcRec_trk;

    rectrk->id = (*iter_trk)->trackId();
    rectrk->chiSq = (*iter_trk)->chi2();
    rectrk->ndf = (*iter_trk)->ndof();
    rectrk->fiTerm = (*iter_trk)->getFiTerm();
    rectrk->nhits = (*iter_trk)->getNhits();
    rectrk->nster = (*iter_trk)->nster();
    rectrk->nclus = 0;
    rectrk->stat = (*iter_trk)->stat();
    status_temp = (*iter_trk)->stat();
    MdcRec_trk_add* trkadd = new MdcRec_trk_add;
    trkadd->id = (*iter_trk)->trackId(); 
    trkadd->quality = 0;
    trkadd->kind = 1;
    trkadd->decision = 0;
    trkadd->body = rectrk;
    rectrk->add = trkadd;

    for ( int i=0; i<5; i++) {
      rectrk->helix[i] = (*iter_trk)->helix()[i];
      if( i<3 ) rectrk->pivot[i] = (*iter_trk)->getPivot()[i];
      for( int j = 1; j<i+2;j++) {
	rectrk->error[i*(i+1)/2+j-1] = (*iter_trk)->err()(i+1,j);
      }
    }
    //std::sort(gothits.begin(),  gothits.end(),  order_rechits);
    HitRefVec::iterator it_gothit = gothits.begin();
    for( ; it_gothit != gothits.end(); it_gothit++) {

      //yzhang delete 2014-04-21 
      //if( (*it_gothit)->getStat() != 1 ) {
      //	if(activeonly_) {
      //		log<<MSG::WARNING<<"this hit is not used in helix fitting!"<<endreq; 
      //		cout<<"this hit is not used in helix fitting!"<<endl; 
      //		continue;
      //	}
      //}

      log << MSG::DEBUG << "retrieved hits in MDC tracks:"
	<< "   hits DDL " <<(*it_gothit)->getDriftDistLeft()
	<< "   hits DDR " <<(*it_gothit)->getDriftDistRight()
	<< "   error DDL " <<(*it_gothit)->getErrDriftDistLeft()
	<< "   error DDR " <<(*it_gothit)->getErrDriftDistRight()
	<< "   id of hit "<<(*it_gothit)->getId() 
	<< "   track id of hit "<<(*it_gothit)->getTrkId()
	<< "   hits ADC " <<(*it_gothit)->getAdc() << endreq;

      MdcRec_wirhit* whit = new MdcRec_wirhit;
      whit->id = (*it_gothit)->getId();
      whit->ddl = (*it_gothit)->getDriftDistLeft();
      whit->ddr = (*it_gothit)->getDriftDistRight();
      whit->erddl = (*it_gothit)->getErrDriftDistLeft();
      whit->erddr = (*it_gothit)->getErrDriftDistRight();
      whit->pChiSq = (*it_gothit)->getChisqAdd();
      whit->lr = (*it_gothit)->getFlagLR();
      whit->stat = (*it_gothit)->getStat();
      mdcid = (*it_gothit)->getMdcId();
      int layid =  MdcID::layer(mdcid);
      int localwid = MdcID::wire(mdcid);
      int w0id = geosvc->Layer(layid)->Wirst();
      int wid = w0id + localwid;
      //if(debug_==4)cout
      //	<< "(" << layid
      //	<< "," << localwid
      //	<< ") lr from MDC: "<<whit->lr
      //	<< endl;

      const MdcGeoWire * const wirgeo = geosvc->Wire(wid); 

      //std::cout<<"the track id of *it_gothit... "<<(*it_gothit)->getTrackId()<<std::endl; 
      whit->rechitptr = *it_gothit;
      whit->geo = wirgeo;
      whit->dat = 0;
      whit->trk = rectrk;
      whit->tdc = (*it_gothit)->getTdc();
      whit->adc= (*it_gothit)->getAdc();
      rectrk->hitcol.push_back(whit);
      mhit_mgr->push_back(*whit);
    }
    mtrk_mgr->push_back(*rectrk);
    mtrkadd_mgr->push_back(*trkadd);

    delete rectrk;
    delete trkadd;
  }

  // check trkcol: evt level
  if(ntuple_&4) {
    m_trkdelx = trkx1 - trkx2;
    m_trkdely = trky1 - trky2;
    m_trkdelz = trkz1 - trkz2;
    m_trkdelthe = trkthe1 + trkthe2;
    m_trkdelphi = trkphi1- trkphi2;
    m_trkdelp = trkp1 - trkp2;
    StatusCode sc4 = m_nt4->write();
    if( sc4.isFailure() ) cout<<"Ntuple4 filling failed!"<<endl;
  }

  if(debug_ == 4) { 
    std::cout<<"before refit,ntrk,nhits,nadd..."<<mtrk_mgr->size()
      <<"********"<<mhit_mgr->size()<<"****"<<mtrkadd_mgr->size()<<endl;
  }
  // Actual fitter procedure :
  if(usage_ == 0) kalman_fitting_anal();

  log << MSG::DEBUG <<"after kalman_fitting(),but in execute...."<<endreq;
  clearTables();   


  MdcID mdcId;
  SmartDataPtr<RecMdcKalTrackCol> recmdckaltrkCol(eventSvc(),"/Event/Recon/RecMdcKalTrackCol");
  //cout<<"------------------------ new event ---------------------"<<endl;
  //cout<<"recmdckaltrkCol->size()="<<recmdckaltrkCol->size()<<endl;
  //cout<<"--------------------------------------------------------"<<endl;
  RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
  RecMdcKalTrackCol::iterator KalTrk= recmdckaltrkCol->begin();
  int i_trk=0;
  for(;KalTrk !=recmdckaltrkCol->end();KalTrk++){
    //cout<<"*** track "<<i_trk++<<" ***"<<endl;
    for(int hypo=0; hypo<5; hypo++)
    {
      if((*KalTrk)->getStat(0,hypo)==1) nFailedTrks[hypo]++;
    }
    HelixSegRefVec gothelixsegs = (*KalTrk)->getVecHelixSegs();
    HelixSegRefVec::iterator iter_hit = gothelixsegs.begin();
    //if(iter_hit == gothelixsegs.end())cout<<"iter_hit == gothelixsegs.end()"<<endl;
    int nhitofthistrk=0;
    for( ; iter_hit != gothelixsegs.end(); iter_hit++){
      nhitofthistrk++;
      //cout<<"layerId: "<<(*iter_hit)->getLayerId()<<endl;
      //cout<<"Identifier: "<<(*iter_hit)->getMdcId()<<endl;
      //cout<<"layerId: "<<mdcId.layer((*iter_hit)->getMdcId())<<endl;
      //cout<<"getDT: "<<(*iter_hit)->getDT()<<endl;
    }
    iter_hit=gothelixsegs.begin();
    //for(int m=0; m<nhitofthistrk/5;m++){
    //	identifier = (*iter_hit) -> getMdcId();
    //}
  }

  return StatusCode::SUCCESS;
}

// Fill TDS:
void KalFitAlg::fillTds(MdcRec_trk& TrasanTRK, KalFitTrack& track, 
    RecMdcKalTrack* trk , int l_mass) {

  //yzhang delete change pivot 2014-07-03 
  //HepPoint3D IP(0,0,0);
  //track.pivot(IP);

  // Fit quality
  int iqual(1);
  int trasster = TrasanTRK.nster, trakster = track.nster(),
      trasax(TrasanTRK.nhits-trasster), trakax(track.nchits()-trakster);
  if (TrasanTRK.nhits-track.nchits()>fitnocut_ || TrasanTRK.helix[2]*track.a()[2]<0){

    if(TrasanTRK.nhits-track.nchits()>fitnocut_) {
      std::cout<<__FILE__<<"   "<<__LINE__<<" run "<<runNo<<" evt "<<eventNo <<" hypo "<<l_mass<<" nHitMdc="<<TrasanTRK.nhits<<" nHitKal="<<track.nchits()<<" fit failed by nhits diff ="<<TrasanTRK.nhits-track.nchits()<<" > "<<fitnocut_<<endl;
    }
    if(TrasanTRK.helix[2]*track.a()[2]<0){
      std::cout<<__FILE__<<"   "<<__LINE__<<" run "<<runNo<<" evt "<<eventNo<<" hypo "<<l_mass<<" fit failed by charge sign reversed mdc track helix[2] "<<TrasanTRK.helix[2] <<" kal helix[2] "<<track.a()[2]<<endl;
    }
    iqual = 0;
  }


  if (debug_ == 4) {
    std::cout<<" iqual  "<<iqual<<" nhit diff "<<TrasanTRK.nhits-track.nchits()<<">?"<<fitnocut_<<" || TrasanTRK.helix[2]*track.a()[2] "<< TrasanTRK.helix[2]*track.a()[2]<<"<?0"<<std::endl;
    cout<< "trasster trakster trasax trakax  TrasK      trackK  iqual"<<endl
      <<trasster<<"        "<<trakster<<"        "<<trasax<<"        "<<trakax
      <<"   "<<TrasanTRK.helix[2]<<"   "<<track.a()[2]<<"   "<<iqual<<endl;  
    cout<<"FillTds> track.chiSq..."<<track.chiSq()<<" nchits "<<track.nchits()
      <<" nster "<<track.nster()<<" iqual "<<iqual<<" track.Ea "<< track.Ea()<<endl; 

    cout<<"fillTds>.....track.Ea[2][2] "<<track.Ea()[2][2]<<endl;
    cout << " TRASAN stereo = " << trasster
      << " and KalFitTrack = " << trakster << std::endl;
    //cout << " TRASAN axial = " << trasax 
    //<< " and KalFitTrack = " << trakax << std::endl;

    if (!iqual) {
      cout << "...there is a problem during fit !! " << std::endl;
      if (trasster-trakster>5) 
	cout << " because stereo " << trasster-trakster << std::endl;
      //if (trasax-trakax >5)	
      //cout << " because axial " << std::endl;
      if (TrasanTRK.helix[2]*track.a()[2]<0)
	cout << " because kappa sign " << std::endl;
    }
    //std::cout<<__FILE__<<"   "<<__LINE__<<" nhits >5? "<<track.nchits()<<" nst>1? "<<track.nster() << " chi2>0? "<<track.chiSq()<<" Ea>0? iqual>0 "<<iqual<<std::endl;
  }
  if (debug_ >0){
    if(track.nchits()<=5) {cout<<"evt "<<eventNo<<" hypo "<<l_mass<<" fit failed by nhits="<<track.nchits()<<" >5"<<endl;}
    if(track.nster()<=1) {cout<<"evt "<<eventNo<<" hypo "<<l_mass<<" fit failed by nst="<<track.nster()<<" >1"<<endl;}
    if(track.chiSq()<=0) {cout<<"evt "<<eventNo<<" hypo "<<l_mass<<" fit failed by chi2="<<track.chiSq()<<" >0"<<endl;}
    if(!(track.Ea()[0][0] > 0 && track.Ea()[1][1] > 0 && track.Ea()[2][2] > 0 && track.Ea()[3][3] > 0 && track.Ea()[4][4] > 0)){cout<<"evt "<<eventNo<<" hypo "<<l_mass<<" fit failed by error matrix >0"<<endl;}
  }
  // Protection : if any problem, we keep the original information !!!!
  if (track.nchits() > 5 && track.nster() > 1 && //track.nchits()-track.nster() > 2 && 
      track.chiSq() > 0 &&
      track.Ea()[0][0] > 0 && track.Ea()[1][1] > 0 && 
      track.Ea()[2][2] > 0 && track.Ea()[3][3] > 0 && 
      track.Ea()[4][4] > 0 && iqual) {

    cout<<__FILE__<<" "<<__LINE__<<" run "<<runNo<<" evt "<<eventNo <<" hypo "<<l_mass<<" fit success nHitKal="<< track.nchits() <<" nHitMdc "<<TrasanTRK.nhits <<" nHitDiff "<<TrasanTRK.nhits-track.nchits()<<" fp "<<track.momentum().mag()<<" pivot "<<track.pivot()<<endl;//yzhang debug
    if(debug_ == 4) cout<<"fillTds>.....going on "<<endl;
    trk->setStat(0,0,l_mass);
    if(debug_ >0) cout<<"evt "<<eventNo<<" mass[" <<l_mass<<"] filter success  "<<std::endl;
    trk->setMass(track.mass(),l_mass);

    // chisq & ndf 
    trk->setChisq(track.chiSq(),0,l_mass);
    trk->setNdf(track.nchits()-5,0,l_mass);
    trk->setNhits(track.nchits(),l_mass);
    //trkid
    trk->setTrackId(TrasanTRK.id);

    //std::cout<<__FILE__<<"   "<<__LINE__<<"  setFHelix "<<track.a()<<" l_mass "<<l_mass<<std::endl;
    trk->setFHelix(track.a(),l_mass);
    trk->setFError(track.Ea(),l_mass);


  } else {

    if(debug_) cout<<"ALARM: FillTds Not refit with KalFilter!!!"<<endl;
    // NOT refit with Kalman filter :
    trk->setStat(1,0,l_mass);
    if(debug_ >0 ) cout<<"evt "<< eventNo <<" mass[" <<l_mass<<"] filter failed"<<std::endl;
    trk->setMass(KalFitTrack::mass(l_mass),l_mass);
    // chisq & ndf (0 : filter ; 1 : smoother;)
    trk->setChisq(TrasanTRK.chiSq,0,l_mass);
    trk->setNdf(TrasanTRK.nhits-5,0,l_mass);
    //trkid
    trk->setTrackId(TrasanTRK.id);
    // nhits  
    trk->setNhits(TrasanTRK.nhits,l_mass);
    double a_trasan[5], ea_trasan[15];
    for( int i =0 ; i <5; i++){
      a_trasan[i] = TrasanTRK.helix[i];
    }
    for( int j =0 ; j <15; j++){
      ea_trasan[j] = TrasanTRK.error[j];
    }
    trk->setFHelix(a_trasan, l_mass);
    trk->setFError(ea_trasan,l_mass);
  }
}


void KalFitAlg::fillTds_ip(MdcRec_trk& TrasanTRK, KalFitTrack& track, 
    RecMdcKalTrack* trk,   int l_mass)
{
  HepPoint3D IP(0,0,0);
  track.pivot(IP);

  if (debug_ == 4&& l_mass==lead_) {
    cout << "fillTds_IP>......"<<endl;
    cout << " dr = " << track.a()[0] 
      << ", Er_dr = " << sqrt(track.Ea()[0][0]) << std::endl;
    cout << " phi0 = " << track.a()[1] 
      << ", Er_phi0 = " << sqrt(track.Ea()[1][1]) << std::endl;
    cout << " PT = " << 1/track.a()[2] 
      << ", Er_kappa =" << sqrt(track.Ea()[2][2]) << std::endl;
    cout << " dz = " << track.a()[3] 
      << ", Er_dz = " << sqrt(track.Ea()[3][3]) << std::endl;
    cout << " tanl = " << track.a()[4] 
      << ", Er_tanl = " << sqrt(track.Ea()[4][4]) << std::endl;
  }

  if (TrasanTRK.nhits-track.nchits()>fitnocut_ ||
      TrasanTRK.helix[2]*track.a()[2]<0)
    iqual_front_[l_mass] = 0;

  if(4==debug_) {
    cout<<"iqual_front_["<<l_mass<<"] = "<<iqual_front_[l_mass]<<endl;
    cout<<"TrasanTRK.nhits="<<TrasanTRK.nhits<<endl;
    cout<<"TrasanTRK.helix[2]="<<TrasanTRK.helix[2]<<endl;
    cout<<"track.nchits()="<<track.nchits()<<" > 5"<<endl;
    cout<<"track.nster()="<<track.nster()<<" > 1"<<endl;
    cout<<"straight wires ="<<track.nchits()-track.nster()<<" > 2 skipped"<<endl;
  }


  //yzhang skip axial hits number cut
  if (track.nchits() > 5 && track.nster() > 1 &&	//track.nchits()-track.nster() > 2 && 
      track.chiSq() > 0 &&
      track.Ea()[0][0] > 0 && track.Ea()[1][1] > 0 &&
      track.Ea()[2][2] > 0 && track.Ea()[3][3] > 0 &&
      track.Ea()[4][4] > 0 && iqual_front_[l_mass]) {
    if(4==debug_) {
      cout<<"evt "<<eventNo<<" kalman fitting succeed!"<<endl;
      cout<<"kalman track parames:"<<track.a()<<endl;
      cout<<"MDC track parames:"<<endl;
      for(int i = 0; i < 5; i++) { cout<<TrasanTRK.helix[i]<<" "<<endl; }
      cout<<"track err mat "<<track.Ea()<<endl;
    }

    // fill track information
    double dr    = track.a()[0];
    double phi0  = track.a()[1];
    double kappa = track.a()[2];
    double dz    = track.a()[3];
    double tanl  = track.a()[4];
    int nLayer = track.nLayerUsed();
    trk->setNlayer(nLayer, l_mass);

    // vertex of the track
    double vx = dr*cos(phi0); 
    double vy = dr*sin(phi0);
    double vz = dz;

    // see Belle note148 for the formulas
    // initial momentum of the track
    if(0==kappa) kappa = 10e-10;
    double px = -sin(phi0)/fabs(kappa);
    double py = cos(phi0)/fabs(kappa);
    double pz = tanl/fabs(kappa);

    trk->setX(vx, l_mass);
    trk->setY(vy, l_mass);
    trk->setZ(vz, l_mass);
    trk->setPx(px, l_mass);
    trk->setPy(py, l_mass);
    trk->setPz(pz, l_mass);

    const HepPoint3D poca(dr*cos(phi0),dr*sin(phi0),dz); 
    trk->setPoca(poca,l_mass);

    trk->setZHelix(track.a(),l_mass);
    trk->setZError(track.Ea(),l_mass);


    //set charge
    int charge=0;
    if (kappa > 0.0000000001)
      charge = 1;
    else if (kappa < -0.0000000001)
      charge = -1;
    trk->setCharge(charge,l_mass);

    //set theta
    double ptot = sqrt(px*px+py*py+pz*pz);
    trk->setTheta(acos(pz/ptot),l_mass);

    //yzhang add set mom
    setRecMdcKalTrackMom(trk, track.a(),l_mass);
  } else{
    //cout<<"copy Mdc Helix in fillTds_ip()"<<endl;

    // fill track information
    double dr    = TrasanTRK.helix[0];
    double phi0  = TrasanTRK.helix[1];
    double kappa = TrasanTRK.helix[2];
    double dz    = TrasanTRK.helix[3];
    double tanl  = TrasanTRK.helix[4];

    double vx = dr*cos(phi0); 
    double vy = dr*sin(phi0);
    double vz = dz;

    if(0==kappa) kappa = 10e-10;
    double px = -sin(phi0)/fabs(kappa);
    double py = cos(phi0)/fabs(kappa);
    double pz = tanl/fabs(kappa);

    trk->setX(vx, l_mass);
    trk->setY(vy, l_mass);
    trk->setZ(vz, l_mass);

    trk->setPx(px, l_mass);
    trk->setPy(py, l_mass);
    trk->setPz(pz, l_mass);

    const HepPoint3D poca(dr*cos(phi0),dr*sin(phi0),dz); 

    trk->setPoca(poca,l_mass);				 
    //trk->setZHelix(TrasanTRK.helix,l_mass);
    //trk->setZError(TrasanTRK.error,l_mass);
    double a_trasan[5], ea_trasan[15];
    for( int i =0 ; i <5; i++){
      a_trasan[i] = TrasanTRK.helix[i];
    }
    for( int j =0 ; j <15; j++){
      ea_trasan[j] = TrasanTRK.error[j];
    }
    trk->setZHelix(a_trasan,l_mass);
    trk->setZError(ea_trasan,l_mass);

    //set charge
    int charge=0;
    if (kappa > 0.0000000001)
      charge = 1;
    else if (kappa < -0.0000000001)
      charge = -1;
    trk->setCharge(charge,l_mass);

    //set theta
    double ptot = sqrt(px*px+py*py+pz*pz);
    trk->setTheta(acos(pz/ptot),l_mass);

    //cout<<"MdcRec_trk: ID = "<<TrasanTRK.id<<endl;

    SmartDataPtr<RecMdcTrackCol> mdcTrkCol(eventSvc(),"/Event/Recon/RecMdcTrackCol");
    //int nMdcTrk = mdcTrkCol.size();
    //cout<<"number of Mdc Tracks: "<<nMdcTrk<<endl;
    RecMdcTrackCol::iterator iter_mdcTrk = mdcTrkCol->begin();
    bool findMdcTrk = false;
    for(; iter_mdcTrk != mdcTrkCol->end(); iter_mdcTrk++) {
      if(TrasanTRK.id==(*iter_mdcTrk)->trackId()) {
	findMdcTrk = true;
	break;
      }
    }
    int nLayer = (*iter_mdcTrk)->nlayer();
    trk->setNlayer(nLayer, l_mass);

    //yzhang add set mom
    setRecMdcKalTrackMom(trk, track.a(),l_mass);

  }

  if(4==debug_) {
    RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
    std::cout<<"fillTds_ip p:"<<sqrt(trk->px()*trk->px()+trk->py()*trk->py()+trk->pz()*trk->pz())<<" px: "<<trk->px()<<" py: "<<trk->py()<<" pz: "<<trk->pz()<<std::endl;
    std::cout<<"vx: "<<trk->x()<<" vy: "<<trk->y()<<" vz: "<<trk->z()<<std::endl;
  }
}

void KalFitAlg::sameas(RecMdcKalTrack* trk, int l_mass, int imain)
{
  // note: for this function,only imain==lead(2) considered
  //std::cout<<"BEGINNING THE sameas() function ..."<<std::endl;
  trk->setMass(trk->getMass(imain), l_mass);
  trk->setLength(trk->getLength(imain), l_mass);
  trk->setTof(trk->getTof(imain), l_mass);
  trk->setNhits(trk->getNhits(imain), l_mass);

  for(int jj = 0; jj<2; jj++) {
    trk->setStat(trk->getStat(jj,imain), jj, l_mass);
    trk->setChisq(trk->getChisq(jj, l_mass), jj, l_mass);
    trk->setNdf(trk->getChisq(jj, l_mass), jj, l_mass);
  }
  trk->setLHelix(trk->getFHelix(),l_mass);
  trk->setLError(trk->getFError(),l_mass);
}

void KalFitAlg::filter_fwd_anal(KalFitTrack& track, int l_mass, int way, HepSymMatrix& Eakal) {
  //cout<<__LINE__<<" filter_fwd_anal helix: "<<track.a()<<" pivot "<<track.pivot()<<endl;

  // cout<<"**********************"<<endl;//wangll
  // cout<<"filter pid type "<<l_mass<<endl;//wangll
  // retrieve Mdc  geometry information

  IMdcGeomSvc* igeomsvc;
  StatusCode sc = Gaudi::svcLocator()->service("MdcGeomSvc", igeomsvc);
  if(sc==StatusCode::FAILURE) cout << "GeoSvc failing!!!!!!!SC=" << sc << endl; 
  MdcGeomSvc* geomsvc = dynamic_cast<MdcGeomSvc*>(igeomsvc);
  if(!geomsvc){
    std::cout<<"ERROR OCCUR when dynamic_cast in KalFitTrack.cxx !!"<<std::endl;
  }

  Hep3Vector x0inner(track.pivot());
  HepVector pos_old(3,0);
  double r0kal_prec(0);
  int  nhits_read(0);
  int nhit = track.HitsMdc().size();
  if(debug_ == 4) 
  {
    cout<<"filter_fwd.......... nhit="<<nhit<<endl;
    cout<<"x0inner = "<<x0inner<<endl;
    cout<<"B of track: "<<track.bFieldZ()<<endl;
  }
  double lastHitTrackLength = -999;
  double lastHitLayer = -999;
  for( int i=0 ; i < nhit; i++ ) {
    KalFitHitMdc& HitMdc = track.HitMdc(i);
    // veto on some hits :
    if (HitMdc.chi2()<0) continue;    
    const KalFitWire& Wire = HitMdc.wire();
    int layerf = Wire.layer().layerId();

    bool b_secondTurn = false;
    double trackLength = HitMdc.rechitptr()->getFltLen();//cm
    if(i>0 && (lastHitLayer==layerf) && fabs(trackLength - lastHitTrackLength)> 10.){//diff 10cm think it cross inner wall
      b_secondTurn = true;
      if(debug_==4)std::cout<<__FILE__<<"   "<<__LINE__<<" "<<i<<" hits  trackLength "<<trackLength<<" lastHitTrackLength "<<lastHitTrackLength<<" b_secondTurn is true"<<std::endl;
    }
    if(m_fix2TurnTrackInnerWall){
      //bool turnOverPi = false;
      //yzhang add 2014-05-01 
      //fix innerwall ,layer No same and track length cut 10cm hit TEMP FIXME
      if(debug_ == 4) cout<<endl<<"for inner wall: layer "<<layerf<<" lastHitLayer "<<lastHitLayer<<" trackLength "<<trackLength<<" lastHitTrackLength "<<lastHitTrackLength<<" trackLenghtDiff "<<trackLength-lastHitTrackLength<<endl;
      if(b_secondTurn){
	if(debug_ == 4){
	  double phi0 = track.a()[1];
	  double pxy = fabs(1./track.a()[2]);
	  double tanl = track.a()[4];
	  double px = pxy * (-sin(phi0));
	  double py = pxy * cos(phi0);
	  double pz = pxy * tanl;
	  double ptot = sqrt(px*px+py*py+pz*pz);
	  cout<<"before fix innerwall once p"<<ptot<<" helix:"<<track.a()<<", Ea:"<<track.Ea()<<endl;
	  cout<<endl<<" hit ("<<layerf<<","<<Wire.localId()<<") fix inner wall: layer "<<layerf<<" lastHitLayer "<<lastHitLayer<<" trackLength "<<trackLength<<" lastHitTrackLength "<<lastHitTrackLength<<endl;
	}
	innerwall(track, l_mass, way); 
	if(debug_ == 4){
	  double phi0 = track.a()[1];
	  double pxy = fabs(1./track.a()[2]);
	  double tanl = track.a()[4];
	  double px = pxy * (-sin(phi0));
	  double py = pxy * cos(phi0);
	  double pz = pxy * tanl;
	  double ptot = sqrt(px*px+py*py+pz*pz);
	  cout<<"after fix innerwall once p"<<ptot<<" helix:"<<track.a()<<", Ea:"<<track.Ea()<<endl;
	}
	//HepPoint3D IP(0,0,0);//FIXME for track not from (0,0,0) target 
	//track.pivot(IP,b_secondTurn);
	if(debug_==4) std::cout<<__FILE__<<"   "<<__LINE__<<" change pivot by fltlen diff: "<<trackLength-lastHitTrackLength<<std::endl;
	track.pivotByFltLen(trackLength-lastHitTrackLength);

	if(debug_==4) std::cout<<__FILE__<<"   "<<__LINE__<<" before fix innerwall pivot by fltLen "<<track.pivot()<<" helix: "<<track.a()<<std::endl;
	innerwall(track, l_mass, way*-1); 
	if(debug_ == 4) {
	  double phi0 = track.a()[1];
	  double pxy = fabs(1./track.a()[2]);
	  double tanl = track.a()[4];
	  double px = pxy * (-sin(phi0));
	  double py = pxy * cos(phi0);
	  double pz = pxy * tanl;
	  double ptot = sqrt(px*px+py*py+pz*pz);
	  cout<<"after fix innerwall twice helix:"<<track.a()<<", Ea:"<<track.Ea()<<track.pivot()<<" p "<<ptot<<endl;
	}
	//turnOverPi = true;
      }
      lastHitTrackLength = trackLength;
      lastHitLayer = layerf;
    }

    int wireid = Wire.geoID();
    nhits_read++;
    HepPoint3D fwd(Wire.fwd());
    HepPoint3D bck(Wire.bck());
    Hep3Vector wire = (CLHEP::Hep3Vector)fwd -(CLHEP::Hep3Vector)bck;
    Helix work = *(Helix*)&track;
    work.ignoreErrorMatrix();
    work.pivot((fwd + bck) * .5);//yzhang delete temp
    if(debug_ == 4) std::cout<<__FILE__<<"   "<<__LINE__<<" change pivot of work to "<< (fwd + bck) * .5<< " after that pivot: "<<work.pivot()<<std::endl;
    //turnOverPi = false;

    if(debug_ == 4) {
      cout<<endl<<"===filter hit No."<<i<<" ("<<layerf<<","<<Wire.localId()<<")"<<endl;
      //cout<<"layerf,cellid wireid = "<<layerf<<", "<<Wire.localId()<<" "<<wireid<<endl;
      cout<<"helix:"<<track.a()<<", Ea:"<<track.Ea()<<endl;
      std::cout<<" work.pivot: (fwd + bck) * .5: "<<(fwd + bck)*.5<<std::endl;
      std::cout<<" track.x(0): "<<track.x(0) <<" track.dz "<<track.dz()<<" tanl "<<track.tanl()<<" curv "<<track.curv()<<" pivot "<<track.pivot()<<std::endl;
      std::cout<<" work.x(0): "<<work.x(0)<<" work.dz "<<work.dz()<<" tanl "<<work.tanl()<<" curv "<<work.curv()<<" pivot "<<work.pivot()<<std::endl;
      std::cout<<" bck: "<<bck<<std::endl;
    }

    //HepPoint3D workZ(work.x(0).z());
    //if(b_secondTurn) workZ = work.x(M_PI*2.);//yzhang 2014-06-18 
    HepPoint3D x0kal = (work.x(0).z() - bck.z())/ wire.z() * wire + bck;
    //HepPoint3D x0kal = (workZ - bck.z())/ wire.z() * wire + bck;

    if(4 == debug_) std::cout<<" x0kal before sag: "<<x0kal<<std::endl;

    // Modification to take account of the wire sag :
    if (wsag_ == 4){
      Hep3Vector result;
      const MdcGeoWire* geowire = geomsvc->Wire(wireid); 
      double tension = geowire->Tension();
      //std::cout<<" tension: "<<tension<<std::endl;
      double zinit(x0kal.z()), lzx(Wire.lzx());
      // double A(Wire.Acoef());
      double A = 47.35E-6/tension;
      double Zp = (zinit - bck.z())*lzx/wire.z();

      if(4 == debug_){
	std::cout<<" sag in filter_fwd_anal: "<<std::endl;
	std::cout<<" x0kal.x(): "<<std::setprecision(10)<<x0kal.x()<<std::endl;
	std::cout<<"zinit: "<<zinit<<" bck.z(): "<<bck.z()<<std::endl;
	std::cout<<" wire.x()*(zinit-bck.z())/wire.z(): "<<std::setprecision(10)
							   <<(wire.x()*(zinit-bck.z())/wire.z())<<std::endl;  
	std::cout<<"bck.x(): "<<std::setprecision(10)<<bck.x()<<std::endl;
	std::cout<<" wire.x()*(zinit-bck.z())/wire.z() + bck.x(): "<<std::setprecision(10)
								     <<(wire.x()*(zinit-bck.z())/wire.z() + bck.x())<<std::endl;
      }

      result.setX(wire.x()*(zinit-bck.z())/wire.z() + bck.x());
      result.setY((A*(Zp-lzx)+wire.y()/lzx)*Zp+bck.y());
      result.setZ((A*(2*Zp-lzx)*lzx+wire.y())/wire.z());

      wire.setX(wire.x()/wire.z());
      wire.setY(result.z());
      wire.setZ(1);
      x0kal.setX(result.x());
      x0kal.setY(result.y());
    }

    if(4 == debug_) std::cout<<" x0kal after sag: "<<x0kal<<std::endl;

    // If x0kal is after the inner wall and x0kal_prec before :
    double r0kal = x0kal.perp();

    // change PIVOT :
    double pathl(0);

    track.pivot_numf(x0kal, pathl);
    if(debug_ == 4) {
      cout<<"after pivot_numf:  pathl:"<<pathl<<"\nhelix:"<<track.a()<<", Ea:"<<track.Ea()<<endl;
    }

    if (nhits_read == 1) { 
      track.Ea(Eakal);
    } else {
      if(KalFitElement::muls()) track.msgasmdc(pathl, way);
      if(debug_ == 4)	cout<<"after muls"<<":  helix:"<<track.a()<<", Ea:"<<track.Ea()<<endl;
      if(KalFitElement::loss()) track.eloss(pathl, _BesKalmanFitMaterials[0], way);
      if(debug_ == 4)	cout<<"after Eloss"<<":  helix:"<<track.a()<<", Ea:"<<track.Ea()<<endl;
    }

    double dtracknew = 0.;
    double dtrack = 0.;
    double dtdc = 0.;
    // Add info hit wire :
    if(fabs(track.kappa())>0&&fabs(track.kappa())<1000.0&&fabs(track.tanl())<7.02) {
      if(debug_==4) cout<<" track parameter valid "<<endl;
      Hep3Vector meas = track.momentum(0).cross(wire).unit();
      double diff_chi2 = track.chiSq();
      //yzhang delete following 4 lines 2014-06-24 
      //Hep3Vector IP(0,0,0);
      //Helix work_bef = *(Helix*)&track;
      //work_bef.ignoreErrorMatrix();
      //work_bef.pivot(IP);
      int inext(-1);
      if (i+1<nhit)
	for( unsigned k=i+1 ; k < nhit; k++ )
	  if (!(track.HitMdc(k).chi2()<0)) {
	    inext = (signed) k;
	    break;
	  }
      double dchi2 = -1.0;

      if(debug_ == 4) {
	cout<<"inext "<<inext<<",meas "<<meas<<",way "<<way<<",dch2 "<<dchi2<<",dtrack "<<dtrack<<",dtracknew "<<dtracknew<<",dtdc "<<dtdc<<endl;
      }

      if(inext==1) track.fltLenTerm(HitMdc.rechitptr()->getFltLen());//cm//yzhang add 2014-07-01 
      double chi2 = track.update_hits(HitMdc,inext,meas,way,dchi2,dtrack,dtracknew,dtdc,m_csmflag);
      if(debug_ == 4) {
	cout<<"after update_hits hit:("<<Wire.layer().layerId()<<","<<Wire.localId()<<") "<<endl;
	cout<<"dchi2="<<dchi2<<endl;
	cout<<"track.a()="<<track.a()<<endl;
	cout<<"track.Ea()="<<track.Ea()<<endl;
	cout<<"track.pivot()="<<track.pivot()<<endl;
      }


      /// get the doca from another other independent method    

      /*
	 std::cout<<" step0: "<<std::endl;
	 KalFitTrack temp2(track);
	 std::cout<<" step1: "<<std::endl;

	 Helix       temp3(track.pivot(),track.a(),track.Ea());
	 Helix       temp4(track.pivot(),track.a(),track.Ea());

	 std::cout<<" step2: "<<std::endl;
	 double doca25 = temp2.approach(HitMdc, false);
	 std::cout<<" step3: "<<std::endl;

	 temp2.pivot(IP);
	 std::cout<<" a2: "<<temp2.a()<<std::endl;

	 std::cout<<" step4: "<<std::endl;

	 double doca26 = temp3.approach(HitMdc, false);
	 std::cout<<" another doca2.6: "<<doca26<<std::endl;

	 temp3.pivot(IP);
	 std::cout<<" a3: "<<temp3.a()<<std::endl;

	 temp4.bFieldZ(-10);
	 temp4.pivot(IP);
	 std::cout<<" a4: "<<temp4.a()<<std::endl;

	 std::cout<<" step5: "<<std::endl;

	 double doca1 = track.approach(HitMdc, false);
	 double doca2 = temp2.approach(HitMdc, false);
	 double doca3 = temp3.approach(HitMdc, false);
	 double doca4 = temp4.approach(HitMdc, false);

	 std::cout<<" dtrack: "<<dtrack<<std::endl;
	 std::cout<<" another doca1: "<<doca1<<std::endl;
	 std::cout<<" another doca2: "<<doca2<<std::endl;
	 std::cout<<" another doca2.5: "<<doca25<<std::endl;
	 std::cout<<" another doca3: "<<doca3<<std::endl;
	 std::cout<<" another doca4: "<<doca4<<std::endl;
       */  


      if( dchi2 <0 ) {
	std::cout<<" ... ERROR OF dchi2... "<<std::endl;
      }

      if (ntuple_&8) {
	m_dchi2 = dchi2;
	m_masshyp = l_mass;
	m_residest = dtrack-dtdc;
	m_residnew = dtracknew -dtdc;
	m_layer =  Wire.layer().layerId();
	Helix worktemp = *(Helix*)&track;
	m_anal_dr = worktemp.a()[0];
	m_anal_phi0 = worktemp.a()[1]; 
	m_anal_kappa = worktemp.a()[2];
	m_anal_dz = worktemp.a()[3]; 
	m_anal_tanl = worktemp.a()[4]; 
	m_anal_ea_dr = worktemp.Ea()[0][0];
	m_anal_ea_phi0 = worktemp.Ea()[1][1];
	m_anal_ea_kappa = worktemp.Ea()[2][2];
	m_anal_ea_dz = worktemp.Ea()[3][3];
	m_anal_ea_tanl = worktemp.Ea()[4][4];
	StatusCode sc5 = m_nt5->write();
	if( sc5.isFailure() ) cout<<"Ntuple5 filling failed!"<<endl;     
      }
      //yzhang delete follwing 3 ines 2014-06-24 
      //Helix work_aft = *(Helix*)&track;
      //work_aft.ignoreErrorMatrix();
      //work_aft.pivot(IP);
      diff_chi2 = chi2 - diff_chi2; 
      HitMdc.chi2(diff_chi2);
    }// if track parameter valid
    else{
      if(debug_==4) cout<< "("<<Wire.layer().layerId()<<","<<Wire.localId()<<") "<<" track parameter NOT valid "<<endl;
    }
    r0kal_prec = r0kal;
  }// end of hits looping
}//end of filter_fwd_anal()



// Take account of the inner wall (forward filter) :
void KalFitAlg::innerwall(KalFitTrack& track, int l_mass, int way){

  if(!m_fixInnerWall) {
    if(debug_==4) std::cout<<__FILE__<<"   "<<__LINE__<<" skip innerwall fix  "<<std::endl;
    return;
  }
  if (debug_ ==4) cout<<"innerwall..... way "<< way<<" _BesKalmanFitWalls.size() "<<_BesKalmanFitWalls.size()<<endl;
  if(way>0){
    if (debug_ == 4) cout<<" update the track!(filter) from outter to inner"<<endl;
    for(int i = 0; i < _BesKalmanFitWalls.size(); i++) {
      if (debug_ == 4) cout<<"Wall "<<i<<" update the track!(filter) track.helix: "<<track.a()<<endl;
      _BesKalmanFitWalls[i].updateTrack(track, i);
    }
  }else{
    if (debug_ == 4) cout<<" update the track!(filter) from inner to outter"<<endl;
    for(int i = _BesKalmanFitWalls.size()-1; i>=0; i--) {
      if (debug_ == 4) cout<<"Wall "<<i<<" update the track!(filter) track.helix: "<<track.a()<<endl;
      _BesKalmanFitWalls[i].updateTrackBack(track, i);
    }
  }
}


// Use the information of trasan and refit the best tracks
void KalFitAlg::kalman_fitting_anal(void) {

  //cout<<"kalman_fitting_anal deal with a new event"<<endl;//wangll

  MsgStream log(msgSvc(), name());
  double Pt_threshold(0.3);
  Hep3Vector IP(0,0,0);
  vector<MdcRec_trk>* mdcMgr = MdcRecTrkCol::getMdcRecTrkCol();
  vector<MdcRec_trk_add>* mdc_addMgr = MdcRecTrkAddCol::getMdcRecTrkAddCol();
  vector<MdcRec_wirhit>* whMgr = MdcRecWirhitCol::getMdcRecWirhitCol();    

  // Table Manager
  if ( !&whMgr ) return;

  // Get reduced chi**2 of Mdc track :
  int ntrk = mdcMgr->size();
  double* rPt = new double[ntrk];
  int* rOM = new int[ntrk];
  unsigned int* order = new unsigned int[ntrk];
  unsigned int* rCont = new unsigned int[ntrk];
  unsigned int* rGen = new unsigned int[ntrk];

  int index = 0;
  for(vector<MdcRec_trk>::iterator it  = mdcMgr->begin(),
      end = mdcMgr->end(); it != end; it++) {
    // Pt
    rPt[index] = 0;
    if (it->helix[2])
      rPt[index] = 1 / fabs(it->helix[2]);
    if(debug_ == 4) cout<<"rPt...[ "<<index<<" ]...."<< rPt[index] <<endl;
    if(rPt[index] < 0) rPt[index] = DBL_MAX;
    // Outermost layer 
    std::vector<MdcRec_wirhit*> pt = it->hitcol ;
    //if(debug_ == 4) cout<<"pt size:  "<< pt.size()<<endl;
    int outermost(-1);
    for (vector<MdcRec_wirhit*>::iterator ii = pt.end()-1;
	ii !=pt.begin()-1; ii--) {
      int lyr((*ii)->geo->Lyr()->Id());
      if (outermost < lyr) outermost = lyr;
      //if(debug_ == 4) cout<<"outmost:  "<<outermost<<"   lyr:  "<<lyr<<endl;
    }
    rOM[index] = outermost;
    order[index] = index;
    ++index;
  }
  if(debug_==4) cout<<"Sort Mdc tracks by Pt"<<endl;
  // Sort Mdc tracks by Pt
  for (int j, k = ntrk - 1; k >= 0; k = j){
    j = -1;
    for(int i = 1; i <= k; i++)
      if(rPt[i - 1] < rPt[i]){
	j = i - 1;
	std::swap(order[i], order[j]);
	std::swap(rPt[i], rPt[j]);
	std::swap(rOM[i], rOM[j]);
	std::swap(rCont[i], rCont[j]);
	std::swap(rGen[i], rGen[j]);
      }
  }
  delete [] rPt;

  //
  int newcount(0);
  //check whether  Recon  already registered
  DataObject *aReconEvent;
  eventSvc()->findObject("/Event/Recon",aReconEvent);
  if(!aReconEvent) {
    // register ReconEvent Data Object to TDS;
    ReconEvent* recevt = new ReconEvent;
    StatusCode sc = eventSvc()->registerObject("/Event/Recon",recevt );
    if(sc!=StatusCode::SUCCESS) {
      log << MSG::FATAL << "Could not register ReconEvent" <<endreq;
      return;
    }
  }
  RecMdcKalTrackCol* kalcol = new RecMdcKalTrackCol;   

  IDataProviderSvc* evtSvc = NULL;
  Gaudi::svcLocator()->service("EventDataSvc", evtSvc);
  if (evtSvc) {
    log << MSG::INFO << "makeTds:event Svc has been found" << endreq;
  } else {
    log << MSG::FATAL << "makeTds:Could not find eventSvc" << endreq;
    return;
    //return StatusCode::SUCCESS;
  }

  StatusCode kalsc;
  IDataManagerSvc *dataManSvc;
  dataManSvc= dynamic_cast<IDataManagerSvc*>(evtSvc);
  DataObject *aKalTrackCol;
  evtSvc->findObject("/Event/Recon/RecMdcKalTrackCol",aKalTrackCol);
  if(aKalTrackCol != NULL) {
    dataManSvc->clearSubTree("/Event/Recon/RecMdcKalTrackCol");
    evtSvc->unregisterObject("/Event/Recon/RecMdcKalTrackCol");
  }

  kalsc = evtSvc->registerObject("/Event/Recon/RecMdcKalTrackCol", kalcol);
  if( kalsc.isFailure() ) {
    log << MSG::FATAL << "Could not register RecMdcKalTrack" << endreq;
    //return StatusCode::SUCCESS;
    return;
  }
  log << MSG::INFO << "RecMdcKalTrackCol registered successfully!" <<endreq;

  DataObject *aKalHelixSegCol;
  evtSvc->findObject("/Event/Recon/RecMdcKalHelixSegCol", aKalHelixSegCol);
  if(aKalHelixSegCol != NULL){
    dataManSvc->clearSubTree("/Event/Recon/RecMdcKalHelixSegCol");
    evtSvc->unregisterObject("/Event/Recon/RecMdcKalHelixSegCol");
  }
  //RecMdcKalHelixSegCol *helixsegcol = new RecMdcKalHelixSegCol;  
  RecMdcKalHelixSegCol *segcol =new RecMdcKalHelixSegCol;
  kalsc = evtSvc->registerObject("/Event/Recon/RecMdcKalHelixSegCol", segcol);
  //kalsc = evtSvc->registerObject("/Event/Recon/RecMdcKalHelixSegCol", helixsegcol);
  if( kalsc.isFailure()){
    log<< MSG::FATAL << "Could not register RecMdcKalHelixSeg" <<endreq;
    return;
    //return StatusCode::SUCCESS;
  }
  log << MSG::INFO << "RecMdcKalHelixSegCol register successfully!" <<endreq;


  //make RecMdcKalTrackCol
  log << MSG::INFO << "beginning to make RecMdcKalTrackCol" <<endreq;     
  // Loop over tracks given by trasan :
  for(int l = 0; l < ntrk; l++) {

    //cout<<"----------------"<<endl;//wangll
    //cout<<"track "<<l<<" : "<<endl;//wangll

    nTotalTrks++;

    for(int i=0; i<5; i++) iqual_front_[i] = 1; // wangll 2010-11-01

    //    m_timer[3]->start();
    MdcRec_trk& TrasanTRK = *(mdcMgr->begin() + order[l]);     
    MdcRec_trk_add& TrasanTRK_add = *(mdc_addMgr->begin()+order[l]);

    // Reject the ones with quality != 0 
    int trasqual = TrasanTRK_add.quality;
    if(debug_ == 4) cout<<"kalman_fitting>trasqual... "<<trasqual<<endl; 
    if (trasqual) continue;

    newcount++;
    if (debug_ == 4)
      cout << "******* KalFit NUMBER : " << newcount << std::endl;      

    // What kind of KalFit ? 
    int type(0);
    if ((TrasanTRK_add.decision & 32) == 32 ||
	(TrasanTRK_add.decision & 64) == 64)      type = 1;

    // Initialisation : (x, a, ea)
    HepPoint3D x(TrasanTRK.pivot[0], TrasanTRK.pivot[1], TrasanTRK.pivot[2]);
    //cout<<"TrasanTRK.pivot: "<<x<<endl;

    HepVector a(5);
    if(debug_ == 4) cout<<"Mdc track helix: "<<endl;
    for(int i = 0; i < 5; i++) {
      a[i] = TrasanTRK.helix[i];
      if(debug_ == 4) cout<<"a["<<i<<"] = "<<a[i]<<endl;
    }
    if(iniTrkPar==-1) {
      if(debug_ == 4) cout<<"==initialize helix from truth par:"<<endl;
      for(int i = 0; i < 5; i++) {
	std::cout<<__FILE__<<"   "<<__LINE__<<" myMCHelix  "<<i<<" "<<myMCHelix[i]<<std::endl;
	a[i] = myMCHelix[i];
	if(debug_ == 4) cout<<"a["<<i<<"] = "<<a[i]<<endl;
      }
    }


    HepSymMatrix ea(5);
    for(int i = 0, k = 0; i < 5; i++) {
      for(int j = 0; j <= i; j++) {
	ea[i][j] = matrixg_*TrasanTRK.error[k++];
	ea[j][i] = ea[i][j];
      }
    }
    double fiTerm = TrasanTRK.fiTerm;
    int way(1);
    // Prepare the track found :
    KalFitTrack track_lead = KalFitTrack(x, a, ea, lead_, 0, 0);

    track_lead.bFieldZ(KalFitTrack::Bznom_);

    // Mdc Hits 
    int inlyr(999), outlyr(-1);
    int* rStat = new int[43];
    for(int irStat=0;irStat<43;++irStat)rStat[irStat]=0;
    std::vector<MdcRec_wirhit*> pt=TrasanTRK.hitcol;
    int hit_in(0);
    if(debug_ == 4) cout<<"*********Pt size(TrasanTRK.hitcol.size)****"<< pt.size()<<endl;
    //// Number of hits/layer 
    //int Num[43] = {0};
    //for (vector<MdcRec_wirhit*>::iterator ii = pt.end()-1;
    //    ii != pt.begin()-1; ii--) {
    //  Num[(*ii)->geo->Lyr()->Id()]++;
    //}
    int hit_asso(0);
    int outmostLayer=-1;
    for (vector<MdcRec_wirhit*>::iterator ii = pt.end()-1; ii != pt.begin()-1; ii--) { 
      hit_asso++;

      hit_in++;
      MdcRec_wirhit & rechit = **ii;
      double dist[2] = {rechit.ddl, rechit.ddr};
      double erdist[2] = {rechit.erddl, rechit.erddr};
      const MdcGeoWire* geo = rechit.geo;

      int lr_decision(0);
      if (KalFitTrack::LR_ == 1){
	if (rechit.lr==2 || rechit.lr==0) lr_decision=-1;
	//	if (rechit.lr==0) lr_decision=-1;
	else if (rechit.lr==1) lr_decision=1;
      } 

      int ind(geo->Lyr()->Id());

      //  ATTENTION HERE!!!
      track_lead.appendHitsMdc( KalFitHitMdc(rechit.id, lr_decision, rechit.tdc, dist, erdist, _wire+(geo->Id()), rechit.rechitptr));
      int layer = geo->Layer();//yzhang add for debug 2014-07-10 
      int cell = geo->Cell();//yzhang add for debug 2014-07-10 
      track_lead.HitMdc(hit_asso-1).pTruth(pTruthHit[layer][cell]);//yzhang add for debug 2014-07-10 
      // inner/outer layer :
      rStat[ind]++;
      if (inlyr>ind) inlyr = ind;
      if (outlyr<ind) outlyr = ind;
    }
    if (debug_ == 4) 
      cout << "**** NUMBER OF Mdc HITS (TRASAN) = " << hit_asso << std::endl;

    // Empty layers :
    int empty_between(0), empty(0);
    for (int i= inlyr; i <= outlyr; i++)
      if (!rStat[i]) empty_between++;
    empty = empty_between+inlyr+(42-outlyr);
    delete [] rStat;

    // RMK high momentum track under study, probably not neeeded...
    //track_lead.order_wirhit(1);
    //track_lead.order_hits();

    if(debug_ == 4) {
      std::cout<<" hits after sorting: nhits="<<track_lead.HitsMdc().size()<<std::endl;    
      for(std::vector<KalFitHitMdc>::iterator it_hit = track_lead.HitsMdc().begin(); it_hit!=track_lead.HitsMdc().end(); it_hit++){
	std::cout<<"("<<it_hit->id()<<": "<<it_hit->wire().layer().layerId()<<","<<it_hit->wire().localId()<<";"<<it_hit->rechitptr()->getFlagLR()<<") rec ddr "<<it_hit->rechitptr()->getDriftDistRight()<<" rec dt "<<it_hit->rechitptr()->getDriftT()<<endl;
	//std::cout<<" the id of this hits after sorting in PatRec is "<<it_hit->id()<<std::endl;    
	//std::cout<<" the layerid of the hit is "<<it_hit->wire().layer().layerId()<<std::endl; 
	//std::cout<<" the cellid of this wire is "<<it_hit->wire().localId()<<std::endl;
      }
      std::cout<<std::endl;
    }

    track_lead.type(type);
    unsigned int nhit = track_lead.HitsMdc().size();
    if (!nhit && debug_ == 4) {
      cout << " ATTENTION TRACK WITH ONLY HITS " << nhit << std::endl;
      continue;
    }

    // Initialisation :
    double  KalFitst(0), KalFitax(0), KalFitschi2(0);
    // Move to the outer hit : 

    Hep3Vector outer_pivot(track_lead.x(fiTerm));

    track_lead.pivot(outer_pivot);

    track_lead.bFieldZ(KalFitTrack::Bznom_);
    //yzhang delete lead 2014-06-25
    //// attention best_chi2 reinitialize !!! 
    //if (nhit>=3 && !KalFitTrack::LR_) 
    //  start_seed(track_lead, lead_, way, TrasanTRK);
    HepSymMatrix Eakal(5,0);

    /// very low momentum and very big costheta angle, use special initial error matrix
    double costheta = track_lead.a()[4] / sqrt(1.0 + track_lead.a()[4]*track_lead.a()[4]); 
    if( (1.0/fabs(track_lead.a()[2]) < pt_cut_ ) && (fabs(costheta)> theta_cut_) ) {
      choice_ = 6;
    }

    /// chose different initial error matrix
    init_matrix(choice_,TrasanTRK, Eakal);


    if (debug_ == 4){
      cout << "from Mdc Pattern Recognition: " << std::endl;
      HepPoint3D IP(0,0,0);
      Helix work(track_lead.pivot(), 
	  track_lead.a(),
	  track_lead.Ea());
      work.pivot(IP);
      cout << " dr = " << work.a()[0] 
	<< ", Er_dr = " << sqrt(work.Ea()[0][0]) << std::endl;
      cout << " phi0 = " << work.a()[1] 
	<< ", Er_phi0 = " << sqrt(work.Ea()[1][1]) << std::endl;
      cout << " PT = " << 1/work.a()[2] 
	<< ", Er_kappa = " << sqrt(work.Ea()[2][2]) << std::endl;
      cout << " dz = " << work.a()[3] 
	<< ", Er_dz = " << sqrt(work.Ea()[3][3]) << std::endl;
      cout << " tanl = " << work.a()[4] 
	<< ", Er_tanl = " << sqrt(work.Ea()[4][4]) << std::endl;
    }

    //yzhang 2014-06-25 close lead particle fitting
    //filter_fwd_anal(track_lead, lead_, way, Eakal);

    // std::cout<<" step3, track_lead.Ea: "<<track_lead.Ea()<<std::endl;
    //track_lead.update_forMdc();

    //HepPoint3D IP(0,0,0);    
    //if (debug_ == 4) {
    //  cout << " Mdc FIRST KALMAN FIT " << std::endl;
    //  Helix work(track_lead.pivot(), 
    //      track_lead.a(),
    //      track_lead.Ea());
    //  work.pivot(IP);
    //  cout << " dr = " << work.a()[0] 
    //    << ", Er_dr = " << sqrt(work.Ea()[0][0]) << std::endl;
    //  cout << " phi0 = " << work.a()[1] 
    //    << ", Er_phi0 = " << sqrt(work.Ea()[1][1]) << std::endl;
    //  cout << " PT = " << 1/work.a()[2] 
    //    << ", Er_kappa = " << sqrt(work.Ea()[2][2]) << std::endl;
    //  cout << " dz = " << work.a()[3] 
    //    << ", Er_dz = " << sqrt(work.Ea()[3][3]) << std::endl;
    //  cout << " tanl = " << work.a()[4] 
    //    << ", Er_tanl = " << sqrt(work.Ea()[4][4]) << std::endl;
    //}

    // fill TDS
    RecMdcKalTrack* kaltrk = new RecMdcKalTrack;
    // Complete the track (other mass assumption, backward) and 
    complete_track(TrasanTRK, TrasanTRK_add, track_lead, kaltrk,kalcol,segcol,1);
  }


  //IDataProviderSvc* eventSvc = NULL;
  //Gaudi::svcLocator()->service("EventDataSvc", eventSvc);
  //if (eventSvc) {
  //  log << MSG::INFO << "makeTds:event Svc has been found" << endreq;
  //} else {
  //  log << MSG::FATAL << "makeTds:Could not find eventSvc" << endreq;
  //  return ;
  //}

  //StatusCode kalsc;
  //IDataManagerSvc *dataManSvc;
  //dataManSvc= dynamic_cast<IDataManagerSvc*>(eventSvc);
  //DataObject *aKalTrackCol;
  //eventSvc->findObject("/Event/Recon/RecMdcKalTrackCol",aKalTrackCol);
  //if(aKalTrackCol != NULL) {
  //  dataManSvc->clearSubTree("/Event/Recon/RecMdcKalTrackCol");
  //  eventSvc->unregisterObject("/Event/Recon/RecMdcKalTrackCol");
  //}

  //kalsc = eventSvc->registerObject("/Event/Recon/RecMdcKalTrackCol", kalcol);
  //if( kalsc.isFailure() ) {
  //  log << MSG::FATAL << "Could not register RecMdcKalTrack" << endreq;
  //  return ; 
  //}
  //log << MSG::INFO << "RecMdcKalTrackCol registered successfully!" <<endreq;

  //StatusCode segsc;
  ////check whether the RecMdcKalHelixSegCol has been already registered
  //DataObject *aRecKalSegEvent;
  //eventSvc->findObject("/Event/Recon/RecMdcKalHelixSegCol", aRecKalSegEvent);
  //if(aRecKalSegEvent!=NULL) {
  //  //then unregister RecMdcKalHelixSegCol
  //  segsc = eventSvc->unregisterObject("/Event/Recon/RecMdcKalHelixSegCol");
  //  if(segsc != StatusCode::SUCCESS) {
  //    log << MSG::FATAL << "Could not unregister RecMdcKalHelixSegCol collection" << endreq;
  //    return;
  //  }
  //}

  //segsc = eventSvc->registerObject("/Event/Recon/RecMdcKalHelixSegCol", segcol);
  //if( segsc.isFailure() ) {
  //  log << MSG::FATAL << "Could not register RecMdcKalHelixSeg" << endreq;
  //  return;
  //}
  //log << MSG::INFO << "RecMdcKalHelixSegCol registered successfully!" <<endreq;


  double x1(0.),x2(0.),y1(0.),y2(0.),z1(0.),z2(0.),the1(0.),the2(0.),phi1(0.),phi2(0.),p1(0.),p2(0.);
  double r1(0.),r2(0.),kap1(999.),kap2(999.),tanl1(0),tanl2(0.); 
  double px1(0.),px2(0.),py1(0.),py2(0.),pz1(0.),pz2(0.),charge1(0.),charge2(0.);

  //check the result:RecMdcKalTrackCol   
  SmartDataPtr<RecMdcKalTrackCol> kaltrkCol(eventSvc(),"/Event/Recon/RecMdcKalTrackCol");
  if (!kaltrkCol) { 
    log << MSG::FATAL << "Could not find RecMdcKalTrackCol" << endreq;
    return;
  }else{
    log << MSG::INFO <<__LINE__ << "Begin to check RecMdcKalTrackCol"<<endreq; 
    RecMdcKalTrackCol::iterator iter_trk = kaltrkCol->begin();
    for( int jj=1; iter_trk != kaltrkCol->end(); iter_trk++,jj++) {
      if(debug_>0){
	cout << "retrieved MDC Kalmantrack:"
	  << "Track Id: " << (*iter_trk)->getTrackId()
	  << " nHits "<< (*iter_trk)->getNhits(0)
	  << " Mass of the fit: "<< (*iter_trk)->getMass(0)
	  << " Length of the track: "<< (*iter_trk)->getLength(0)
	  << " Tof of the track: "<< (*iter_trk)->getTof(0)
	  << " Chisq of the fit: "<< (*iter_trk)->getChisq(0,0)
	  << "  "<< (*iter_trk)->getChisq(1,0) 
	  << " Ndf of the fit: "<< (*iter_trk)->getNdf(0,0)
	  << "  "<< (*iter_trk)->getNdf(1,0)
	  << " Kappa " << (*iter_trk)->getZHelix()[0]
	  << " fHelix Kappa" << (*iter_trk)->getFHelix()[0]
	  << " zHelix " << (*iter_trk)->getZHelix()
	  << endl;

	HepVector zHelixE = (*iter_trk)->getZHelixE();
	HepVector fHelixE = (*iter_trk)->getFHelixE();
	cout << "zHelixE " << zHelixE << "fHelixE " << fHelixE <<endl;
      }
      //for( int i = 0; i<43; i++) {
      //  log << MSG::DEBUG << "retrieved pathl["<<i<<"]= "
      //    << (*iter_trk)->getPathl(i) << " ";
      //  log<<endreq;
      //}



      if(ntuple_&1) {
	//cout<<"fill ntuple 1"<<endl;
	m_trackid = (*iter_trk)->getTrackId();
	//cout<<"fill ntuple 2"<<endl;

	for( int jjj =0, iii=0; jjj<5; jjj++) {

	  m_length[jjj] = (*iter_trk)->getLength(jjj);
	  m_tof[jjj] = (*iter_trk)->getTof(jjj);
	  m_nhits[jjj] = (*iter_trk)->getNhits(jjj);
	  m_zhelix[jjj] = (*iter_trk)->getZHelix()[jjj];
	  m_zhelixe[jjj] = (*iter_trk)->getZHelixE()[jjj];
	  m_zhelixmu[jjj] = (*iter_trk)->getZHelixMu()[jjj];
	  m_zhelixk[jjj] = (*iter_trk)->getZHelixK()[jjj];
	  m_zhelixp[jjj] = (*iter_trk)->getZHelixP()[jjj];
	  if(debug_>0)std::cout<<__FILE__<<"   "<<__LINE__<<" getFHelixE  "<<jjj<<" "<<(*iter_trk)->getFHelixE()[jjj]<<std::endl;
	  m_fhelix[jjj] = (*iter_trk)->getFHelix()[jjj];
	  m_fhelixe[jjj] = (*iter_trk)->getFHelixE()[jjj];
	  m_fhelixmu[jjj] = (*iter_trk)->getFHelixMu()[jjj];
	  m_fhelixk[jjj] = (*iter_trk)->getFHelixK()[jjj];
	  m_fhelixp[jjj] = (*iter_trk)->getFHelixP()[jjj];
	  m_lhelix[jjj] = (*iter_trk)->getLHelix()[jjj];
	  m_lhelixe[jjj] = (*iter_trk)->getLHelixE()[jjj];
	  m_lhelixmu[jjj] = (*iter_trk)->getLHelixMu()[jjj];
	  m_lhelixk[jjj] = (*iter_trk)->getLHelixK()[jjj];
	  m_lhelixp[jjj] = (*iter_trk)->getLHelixP()[jjj];

	  if(ntuple_&32) {
	    for(int kk=0; kk<=jjj; kk++,iii++) {
	      m_zerror[iii] = (*iter_trk)->getZError()[jjj][kk];
	      m_zerrore[iii] = (*iter_trk)->getZErrorE()[jjj][kk];
	      m_zerrormu[iii] = (*iter_trk)->getZErrorMu()[jjj][kk];
	      m_zerrork[iii] = (*iter_trk)->getZErrorK()[jjj][kk];
	      m_zerrorp[iii] = (*iter_trk)->getZErrorP()[jjj][kk];
	      m_ferror[iii] = (*iter_trk)->getFError()[jjj][kk];
	      m_ferrore[iii] = (*iter_trk)->getFErrorE()[jjj][kk];
	      m_ferrormu[iii] = (*iter_trk)->getFErrorMu()[jjj][kk];
	      m_ferrork[iii] = (*iter_trk)->getFErrorK()[jjj][kk];
	      m_ferrorp[iii] = (*iter_trk)->getFErrorP()[jjj][kk];
	      m_lerror[iii] = (*iter_trk)->getLError()[jjj][kk];
	      m_lerrore[iii] = (*iter_trk)->getLErrorE()[jjj][kk];
	      m_lerrormu[iii] = (*iter_trk)->getLErrorMu()[jjj][kk];
	      m_lerrork[iii] = (*iter_trk)->getLErrorK()[jjj][kk];
	      m_lerrorp[iii] = (*iter_trk)->getLErrorP()[jjj][kk];
	    }
	  }
	}
	//cout<<"fill ntuple 3"<<endl;

	// RootConversion changed in BOSS6.0, so use thefollowing:
	m_chisq[0][0] = (*iter_trk)->getChisq(0,0);
	m_chisq[1][0] = (*iter_trk)->getChisq(0,1);
	m_chisq[2][0] = (*iter_trk)->getChisq(0,2);
	m_chisq[3][0] = (*iter_trk)->getChisq(0,3);
	m_chisq[4][0] = (*iter_trk)->getChisq(0,4);
	m_chisq[0][1] = (*iter_trk)->getChisq(1,0);
	m_chisq[1][1] = (*iter_trk)->getChisq(1,1);
	m_chisq[2][1] = (*iter_trk)->getChisq(1,2);
	m_chisq[3][1] = (*iter_trk)->getChisq(1,3);
	m_chisq[4][1] = (*iter_trk)->getChisq(1,4);

	m_ndf[0][0] = (*iter_trk)->getNdf(0,0);
	m_ndf[1][0] = (*iter_trk)->getNdf(0,1);
	m_ndf[2][0] = (*iter_trk)->getNdf(0,2);
	m_ndf[3][0] = (*iter_trk)->getNdf(0,3);
	m_ndf[4][0] = (*iter_trk)->getNdf(0,4);
	m_ndf[0][1] = (*iter_trk)->getNdf(1,0);
	m_ndf[1][1] = (*iter_trk)->getNdf(1,1);
	m_ndf[2][1] = (*iter_trk)->getNdf(1,2);
	m_ndf[3][1] = (*iter_trk)->getNdf(1,3);
	m_ndf[4][1] = (*iter_trk)->getNdf(1,4);

	m_stat[0][0] = (*iter_trk)->getStat(0,0);
	m_stat[1][0] = (*iter_trk)->getStat(0,1);
	m_stat[2][0] = (*iter_trk)->getStat(0,2);
	m_stat[3][0] = (*iter_trk)->getStat(0,3);
	m_stat[4][0] = (*iter_trk)->getStat(0,4);
	m_stat[0][1] = (*iter_trk)->getStat(1,0);
	m_stat[1][1] = (*iter_trk)->getStat(1,1);
	m_stat[2][1] = (*iter_trk)->getStat(1,2);
	m_stat[3][1] = (*iter_trk)->getStat(1,3);
	m_stat[4][1] = (*iter_trk)->getStat(1,4);
	//cout<<"fill ntuple 4"<<endl;

	m_fptot = sqrt(1+pow(m_fhelix[4],2))/m_fhelix[2];
	m_fptote = sqrt(1+pow(m_fhelixe[4],2))/m_fhelixe[2];
	m_fptotmu = sqrt(1+pow(m_fhelixmu[4],2))/m_fhelixmu[2];
	m_fptotk = sqrt(1+pow(m_fhelixk[4],2))/m_fhelixk[2];
	m_fptotp = sqrt(1+pow(m_fhelixp[4],2))/m_fhelixp[2];

	m_lptot = sqrt(1+pow(m_lhelix[4],2))/m_lhelix[2];
	m_lptote = sqrt(1+pow(m_lhelixe[4],2))/m_lhelixe[2];
	m_lptotmu = sqrt(1+pow(m_lhelixmu[4],2))/m_lhelixmu[2];
	m_lptotk = sqrt(1+pow(m_lhelixk[4],2))/m_lhelixk[2];
	m_lptotp = sqrt(1+pow(m_lhelixp[4],2))/m_lhelixp[2];

	m_lpt = 1/m_lhelix[2];
	m_lpte = 1/m_lhelixe[2];
	m_lptmu = 1/m_lhelixmu[2];
	m_lptk = 1/m_lhelixk[2];
	m_lptp = 1/m_lhelixp[2];

	m_fpt = 1/m_fhelix[2];
	m_fpte = 1/m_fhelixe[2];
	m_fptmu = 1/m_fhelixmu[2];
	m_fptk = 1/m_fhelixk[2];
	m_fptp = 1/m_fhelixp[2];
	//cout<<"fill ntuple 5"<<endl;

	if(debug_ >= 3){   
	  std::cout<<"                                           "<<std::endl;
	  std::cout<<"in file Kalman_fitting_anal ,the  m_fpt is .."<<m_fpt<<std::endl;
	  std::cout<<"in file Kalman_fitting_anal ,the  m_fpte is .."<<m_fpte<<std::endl;
	  std::cout<<"in file Kalman_fitting_anal ,the  m_fptmu is .."<<m_fptmu<<std::endl;
	  std::cout<<"in file Kalman_fitting_anal ,the  m_fptk is .."<<m_fptk<<std::endl;
	  std::cout<<"in file Kalman_fitting_anal ,the  m_fptp is .."<<m_fptp<<std::endl;
	}

	m_zpt = 1/m_zhelix[2];
	m_zpte = 1/m_zhelixe[2];
	m_zptmu = 1/m_zhelixmu[2];
	m_zptk = 1/m_zhelixk[2];
	m_zptp = 1/m_zhelixp[2];

	if(debug_ >= 3) {
	  std::cout<<"in file Kalman_fitting_anal ,the  m_zpt is .."<<m_zpt<<std::endl;
	  std::cout<<"in file Kalman_fitting_anal ,the  m_zpte is .."<<m_zpte<<std::endl; 
	  std::cout<<"in file Kalman_fitting_anal ,the  m_zptmu is .."<<m_zptmu<<std::endl;
	  std::cout<<"in file Kalman_fitting_anal ,the  m_zptk is .."<<m_zptk<<std::endl;   
	  std::cout<<"in file Kalman_fitting_anal ,the  m_zptp is .."<<m_zptp<<std::endl;                                           
	}
	m_zptot = sqrt(1+pow(m_zhelix[4],2))/m_zhelix[2];
	m_zptote = sqrt(1+pow(m_zhelixe[4],2))/m_zhelixe[2];
	m_zptotmu = sqrt(1+pow(m_zhelixmu[4],2))/m_zhelixmu[2];
	m_zptotk = sqrt(1+pow(m_zhelixk[4],2))/m_zhelixk[2];
	m_zptotp = sqrt(1+pow(m_zhelixp[4],2))/m_zhelixp[2];
	//cout<<"fill ntuple 6"<<endl;

	if(debug_ >= 3) {
	  std::cout<<"in file Kalman_fitting_anal ,the  m_zptot is .."<<m_zptot<<std::endl;
	  std::cout<<"in file Kalman_fitting_anal ,the  m_zptote is .."<<m_zptote<<std::endl;   
	  std::cout<<"in file Kalman_fitting_anal ,the  m_zptotmu is .."<<m_zptotmu<<std::endl;
	  std::cout<<"in file Kalman_fitting_anal ,the  m_zptotk is .."<<m_zptotk<<std::endl;   
	  std::cout<<"in file Kalman_fitting_anal ,the  m_zptotp is .."<<m_zptotp<<std::endl;
	}

	if(ntuple_&32) {
	  m_zsigp = sqrt(pow((m_zptot/m_zhelix[2]),2)*m_zerror[5]+
	      pow((m_zhelix[4]/m_zptot),2)*pow((1/m_zhelix[2]),4)*m_zerror[14]-
	      2*m_zhelix[4]*m_zerror[12]*pow((1/m_zhelix[2]),3));
	  m_zsigpe = sqrt(pow((m_zptote/m_zhelixe[2]),2)*m_zerrore[5]+
	      pow((m_zhelixe[4]/m_zptote),2)*pow((1/m_zhelixe[2]),4)*m_zerrore[14]-
	      2*m_zhelixe[4]*m_zerrore[12]*pow((1/m_zhelixe[2]),3));
	  m_zsigpmu = sqrt(pow((m_zptotmu/m_zhelixmu[2]),2)*m_zerrormu[5]+
	      pow((m_zhelixmu[4]/m_zptotmu),2)*pow((1/m_zhelixmu[2]),4)*m_zerrormu[14]-
	      2*m_zhelixmu[4]*m_zerrormu[12]*pow((1/m_zhelixmu[2]),3));
	  m_zsigpk = sqrt(pow((m_zptotk/m_zhelixk[2]),2)*m_zerrork[5]+
	      pow((m_zhelixk[4]/m_zptotk),2)*pow((1/m_zhelixk[2]),4)*m_zerrork[14]-
	      2*m_zhelixk[4]*m_zerrork[12]*pow((1/m_zhelixk[2]),3));
	  m_zsigpp = sqrt(pow((m_zptotp/m_zhelixp[2]),2)*m_zerrorp[5]+
	      pow((m_zhelixp[4]/m_zptotp),2)*pow((1/m_zhelixp[2]),4)*m_zerrorp[14]-
	      2*m_zhelixp[4]*m_zerrorp[12]*pow((1/m_zhelixp[2]),3));
	}

	StatusCode sc1 = m_nt1->write();
	if( sc1.isFailure() ) cout<<"Ntuple1 filling failed!"<<endl;     
      }

      if(ntuple_&4) {
	if(jj == 1) {
	  phi1 = (*iter_trk)->getFFi0();
	  r1 = (*iter_trk)->getFDr();
	  z1 = (*iter_trk)->getFDz();	
	  kap1 = (*iter_trk)->getFCpa();	
	  tanl1 = (*iter_trk)->getFTanl();	
	  charge1 = kap1/fabs(kap1);
	  x1 = r1*cos(phi1);
	  y1 = r1*sin(phi1);
	  p1 = sqrt(1+tanl1*tanl1)/kap1;
	  the1 = M_PI/2-atan(tanl1);
	  px1 = -sin(phi1)/fabs(kap1);
	  py1 = cos(phi1)/fabs(kap1);
	  pz1= tanl1/fabs(kap1);

	} else if(jj == 2) {
	  phi2 = (*iter_trk)->getFFi0();
	  r2 = (*iter_trk)->getFDr();
	  z2 = (*iter_trk)->getFDz();	
	  kap2 = (*iter_trk)->getFCpa();	
	  tanl2 = (*iter_trk)->getFTanl();	
	  charge2 = kap2/fabs(kap2);
	  x2 = r1*cos(phi2);
	  y2 = r1*sin(phi2);
	  p2 = sqrt(1+tanl2*tanl2)/kap1;
	  the2 = M_PI/2-atan(tanl2);
	  px2 = -sin(phi2)/fabs(kap2);
	  py2 = cos(phi2)/fabs(kap2);
	  pz2= tanl2/fabs(kap2);
	}
      }

      if(ntuple_&4) {
	m_delx = x1 - x2;
	m_dely = y1 - y2;
	m_delz = z1 - z2;
	m_delthe = the1 + the2;
	m_delphi = phi1- phi2;
	m_delp = p1 - p2;
	m_delpx = charge1*fabs(px1) + charge2*fabs(px2);
	m_delpy = charge1*fabs(py1) + charge2*fabs(py2);
	m_delpz = charge1*fabs(pz1) + charge2*fabs(pz2);

	StatusCode sc2 = m_nt2->write();
	if( sc2.isFailure() ) cout<<"Ntuple2 filling failed!"<<endl;      
      } 
    }//end of loop over RecMdcKalTrack
  }

  delete [] order;
  delete [] rCont;
  delete [] rGen;
  delete [] rOM;

  if (debug_ == 4)
    cout << "Kalfitting finished " << std::endl;
}//end of kalman_fitting_anal


void KalFitAlg::complete_track(MdcRec_trk& TrasanTRK, MdcRec_trk_add& TrasanTRK_add, KalFitTrack& track_lead, RecMdcKalTrack*  kaltrk, RecMdcKalTrackCol* kalcol,RecMdcKalHelixSegCol *segcol,int flagsmooth) {
  static int nmass = KalFitTrack::nmass();
  int way(1);
  MsgStream log(msgSvc(), name());

  //yzhang delete lead 2014-06-25 
  //KalFitTrack track_first(track_lead);
  //KalFitTrack track_ip(track_lead);
  //innerwall(track_ip, lead_, way); 
  //// Fill Tds
  //fillTds_ip(TrasanTRK, track_ip, kaltrk, lead_);
  //// std::cout<<"track_first nster"<<track_first.nster()<<std::endl;

  //fillTds_lead(TrasanTRK, track_first, kaltrk, lead_);
  //// Refit for different mass assumptions :  
  //double pp = track_lead.momentum().mag();

  if(!(i_front_<0)){

    for(int l_mass = 0, flg = 1; l_mass < nmass; l_mass++, flg <<= 1) {

      if (!(mhyp_ & flg)) continue;
      //yzhang delete lead 2014-06-25 
      //      if (l_mass == lead_) continue;

      if(debug_ == 4) cout<<"complete_track..REFIT ASSUMPION " << l_mass << endl;
      // Initialisation :
      double chiSq = 0;
      int nhits = 0;

      // Initialisation : (x, a, ea)
      HepPoint3D x_trasan(TrasanTRK.pivot[0], TrasanTRK.pivot[1], TrasanTRK.pivot[2]);
      HepVector a_trasan(5);
      for(int i = 0; i < 5; i++)
	a_trasan[i] = TrasanTRK.helix[i];

      HepSymMatrix ea_trasan(5);
      for(int i = 0, k = 0; i < 5; i++) {
	for(int j = 0; j <= i; j++) {
	  ea_trasan[i][j] = matrixg_*TrasanTRK.error[k++];
	  ea_trasan[j][i] = ea_trasan[i][j];
	}
      }

      KalFitTrack track(x_trasan,a_trasan, ea_trasan, l_mass, chiSq, nhits);
      track.HitsMdc(track_lead.HitsMdc());	


      double fiTerm = TrasanTRK.fiTerm;
      if(debug_ == 4)
      {
	cout<<"before pivot: "<<track.pivot()<<endl;
	cout<<"helix: "<<track.a()<<", Ea: "<<track.Ea()<<endl;
      }
      if(m_pivotCnvSelector==0){
	//yzhang change pivot m_nPivotCnv times
	for (int i=1;i<=m_nPivotCnv;i++){
	  HepPoint3D pivotTemp(track.x(1./((double)m_nPivotCnv)*fiTerm));
	  track.pivot(pivotTemp);
	  if(debug_==4)cout<<"after pivot"<<i<<": "<<track.pivot()<<endl;
	}
      }else if(m_pivotCnvSelector==1){
	HepPoint3D pivot1(track.x(0.5*fiTerm));
	HepPoint3D pivot2(track.x(fiTerm));
	//track.pivot(track.x(fiTerm));
	track.pivot(pivot1);
	cout<<"after pivot1: "<<track.pivot()<<endl;
	track.pivot(pivot2);
	cout<<"after pivot2: "<<track.pivot()<<endl;
      }else if(m_pivotCnvSelector==2){
	if(debug_==4) std::cout<<__FILE__<<"   "<<__LINE__<<" yzhang not convert pivot  "<<std::endl;
      }
      if(debug_ == 4) {
	cout<<"after pivot: "<<track.pivot()<<endl;
	cout<<"helix: "<<track.a()<<", Ea: "<<track.Ea()<<endl;
      }
      HepSymMatrix Eakal(5,0);

      double costheta = track.a()[4] / sqrt(1.0 + track.a()[4]*track.a()[4]);
      if( (1.0/fabs(track.a()[2]) < pt_cut_ ) && (fabs(costheta)> theta_cut_) ) {
	choice_ = 6;
      }

      init_matrix(choice_,TrasanTRK, Eakal);
      filter_fwd_anal(track, l_mass, way, Eakal);

      // Fill tds
      fillTds(TrasanTRK, track, kaltrk, l_mass);     

      ///fill tds  with results got at (0,0,0)  
      KalFitTrack track_z(track);
      innerwall(track_z, l_mass, way);     
      fillTds_ip(TrasanTRK, track_z, kaltrk, l_mass);  
    }
  } // end of //end of if (!(i_front<0))


  //// Refit with an enhancement of the error matrix at Mdc level :
  //if (enhance_) {

  //  HepPoint3D x_first(0, 0, 0);
  //  HepVector a_first(kaltrk->getFHelix());
  //  HepSymMatrix ea_first(kaltrk->getFError());
  //  HepVector fac(5);
  //  fac[0]=fac_h1_; fac[1]=fac_h2_; fac[2]=fac_h3_; fac[3]=fac_h4_; fac[4]=fac_h5_;
  //  for(int i = 0; i < 5; i++)
  //    for(int j = 0; j <= i; j++)
  //      ea_first[i][j] = fac[i]*fac[j]*ea_first[i][j];
  //  KalFitTrack track(x_first, a_first, ea_first, 2, 0, 0);
  //}

  //// Backward filter :  
  //KalFitTrack track_back(track_lead);
  //if (debug_ == 4) {
  //  cout << " Backward fitting flag:" << back_<< endl;
  //  cout << "track_back pivot " << track_back.pivot() 
  //    << " track_lead kappa " << track_lead.a()[2]
  //    <<endl;
  //}

  //if (back_ && track_lead.a()[2] != 0 && 
  //    1/fabs(track_lead.a()[2]) > pT_) {
  //  track_back.HitsMdc(track_lead.HitsMdc());

  //  if (KalFitTrack::tofall_) {
  //    double p_kaon(0), p_proton(0);
  //    if (!(kaltrk->getStat(0,3))) {
  //      p_kaon = 1 / fabs(kaltrk->getZHelixK()[2]) * 
  //        sqrt(1 + kaltrk->getZHelixK()[4]*kaltrk->getZHelixK()[4]);
  //      track_back.p_kaon(p_kaon);
  //    } else {
  //      p_kaon = 1 / fabs(track_back.a()[2]) *
  //        sqrt(1 + track_back.a()[4]*track_back.a()[4]);
  //      track_back.p_kaon(p_kaon);	  
  //    }
  //    if (!(kaltrk->getStat(0,4))) {
  //      p_proton = 1 / fabs(kaltrk->getZHelixP()[2]) * 
  //        sqrt(1 + kaltrk->getZHelixP()[4]*kaltrk->getZHelixP()[4]);
  //      track_back.p_proton(p_proton);
  //    } else {
  //      p_proton = 1 / fabs(track_back.a()[2]) *
  //        sqrt(1 + track_back.a()[4]*track_back.a()[4]);
  //      track_back.p_proton(p_proton);
  //    }
  //  }

  //yzhang close smoother 2014-06-30 
  //if (!(i_back_<0)) {
  //  //cout<<" *** in smoothing process ***"<<endl;
  //  for(int l_mass = 0; l_mass < nmass; l_mass++) {
  //    //cout<<" --- in hypothesis "<<l_mass<<" :"<<endl;
  //    KalFitTrack track_seed(track_back);
  //    track_seed.chgmass(l_mass);     
  //    /*cout<<"---------------"<<endl;//wangll
  //      cout<<"smooth track "<<l_mass<<endl;//wangll
  //      cout<<"  pivot :"<<track_seed.pivot()<<endl;//wangll
  //      cout<<"  helix :"<<track_seed.a()<<endl;//wangll
  //     */
  //    smoother_anal(track_seed, -way); 
  //    // if( usage_ == 1) smoother_calib(track_seed, -way);
  //    //cout<<"fillTds_back 1"<<endl;
  //    // fill TDS  for backward filter :
  //    fillTds_back(track_seed, kaltrk, TrasanTRK, l_mass, segcol, 1);
  //    //cout<<"nHits: "<<kaltrk->getVecHelixSegs().size()<<endl;
  //  }
  //} else {
  //  smoother_anal(track_back, -way);
  //  fillTds_back(track_back, kaltrk, TrasanTRK, lead_, segcol, 1);
  //}
  //}
  /*
  // Take care of the pointers (use lead. hyp results by default)
  for(int pid = 0; pid < nmass; pid++) {
  if (pid == lead_) continue;
  if (kaltrk->getStat(1,pid)) 
  sameas(kaltrk, pid, lead_);
  }
   */

  //check: before register into TDS
  if(debug_>0){
    cout<< "registered MDC Kalmantrack:"<<__LINE__
      << " Track Id: " << kaltrk->getTrackId()
      << " Mass of the fit: "<< kaltrk->getMass(0)
      << " Length of the track: "<< kaltrk->getLength(0)
      << " Tof of the track: "<< kaltrk->getTof(0)
      << " Chisq of the fit: "<< kaltrk->getChisq(0,0)
      << " "<< kaltrk->getChisq(1,0)
      << " Ndf of the fit: "<< kaltrk->getNdf(0,0)
      << " "<< kaltrk->getNdf(1,0)<<endl;
    HepVector zHelix = kaltrk->getZHelixE();
    HepVector fHelix = kaltrk->getFHelixE();
    cout << "zHelix " << zHelix << "fHelix " << fHelix<<endl;
  }

  kalcol->push_back(kaltrk);
  if(debug_>0){
    std::cout<<__FILE__<<"   "<<__LINE__<<" printRecMdcKalTrackCol  "<<std::endl;
    m_mdcPrintSvc->printRecMdcKalTrackCol();
  }
  track_lead.HitsMdc().clear();
}


void KalFitAlg::init_matrix(MdcRec_trk& trk, HepSymMatrix& Eakal )
{
  for ( int i=0; i<5; i++) {
    for( int j = 1; j<i+2;j++) {
      Eakal(i+1,j) = matrixg_*(trk.error[i*(i+1)/2+j-1]);
      Eakal(j,i+1) = Eakal(i+1,j);
    }
  }

  if (debug_ == 4) cout<<"initialised Ea.. "<<Eakal<<endl;
}



void KalFitAlg::init_matrix(int k, MdcRec_trk& trk, HepSymMatrix& Eakal )
{
  if(0==k){
    for ( int i=0; i<5; i++) {
      for( int j = 1; j<i+2;j++) {
	Eakal(i+1,j) = matrixg_*(trk.error[i*(i+1)/2+j-1]);
	Eakal(j,i+1) = Eakal(i+1,j);
      }
      Eakal(1,1) = Eakal(1,1)* gain1_;
      Eakal(2,2) = Eakal(2,2)* gain2_;
      Eakal(3,3) = Eakal(3,3)* gain3_;
      Eakal(4,4) = Eakal(4,4)* gain4_;
      Eakal(5,5) = Eakal(5,5)* gain5_;
    }
  }
  //  HepSymMatrix ea_temp(5);
  //    for(int i = 0, k = 0; i < 5; i++) {
  //      for(int j = 0; j <= i; j++) {
  //        ea_temp[i][j] = matrixg_*trk.error[k++];
  //        ea_temp[j][i] = ea_temp[i][j];
  //      }
  //    }

  if(1==k){
    Eakal(1,1) = .2;
    Eakal(2,2) = .001;
    Eakal(3,3) = 1.0;
    Eakal(4,4) = 10.0;
    Eakal(5,5) = 0.002;
  }

  if(2==k){
    Eakal(1,1) = .2;
    Eakal(2,2) = 0.1;
    Eakal(3,3) = 1.0;
    Eakal(4,4) = 25.0;
    Eakal(5,5) = 0.10;
  }


  if(3==k){
    Eakal(1,1) = .2;
    Eakal(2,2) = .001;
    Eakal(3,3) = 1.0;
    Eakal(4,4) = 25.0;
    Eakal(5,5) = 0.10;
  }

  if(4==k){
    Eakal(1,1) = .2;
    Eakal(2,2) = .01;
    Eakal(3,3) = 0.01;
    Eakal(4,4) = 1.;
    Eakal(5,5) = .01;
  }

  if(5==k) {
    Eakal(1,1) = 2.;
    Eakal(2,2) = 0.1;
    Eakal(3,3) = 1.;
    Eakal(4,4) = 20.;
    Eakal(5,5) = 0.1;
  }

  if(6==k) {
    Eakal(1,1) = 0.01;
    Eakal(2,2) = 0.01;
    Eakal(3,3) = 0.01;
    Eakal(4,4) = 100.;
    Eakal(5,5) = 0.5;
  }

  if(k!=0){
    Eakal(3,3) = 0.2;
    Eakal(1,1) = 1;
    Eakal(4,4) = 1;
  }

  if (debug_ == 4) cout<<"initialised Eakal.. "<<Eakal<<endl;
}




void KalFitAlg::start_seed(KalFitTrack& track, int lead_, int way, MdcRec_trk& TrasanTRK)
{
  if (debug_ == 4) 
    cout << "start_seed begin... " << std::endl;
  // keep initial helix parameters
  Hep3Vector x_init(track.pivot());
  HepSymMatrix Ea_init(5,0);
  Ea_init = track.Ea();
  HepVector a_init(5);
  a_init = track.a();

  // LR assumption :
  unsigned int nhit_included(10);
  int LR[8][3] = { 
    {1,1,1},
    {1,1,-1},
    {1,-1,1},
    {1,-1,-1},
    {-1,1,1},
    {-1,1,-1},
    {-1,-1,1},
    {-1,-1,-1}
  };

  unsigned int nhit = track.HitsMdc().size();
  double chi2_min(DBL_MAX);
  int i_min(-1);
  for (int ktrial = 0; ktrial < 8; ktrial++) {

    // Come back to trasan seed :
    track.pivot(x_init);
    track.a(a_init);
    track.Ea(Ea_init);

    track.chiSq(0);
    track.chiSq_back(0);
    track.nchits(0);
    track.nster(0);
    track.ndf_back(0);

    HepSymMatrix Eakal(5,0);

    init_matrix(choice_,TrasanTRK, Eakal);
    // initialize the Mdc hits :
    for( unsigned i=0 ; i < nhit; i++ ){
      KalFitHitMdc& HitMdc = track.HitMdc(i);
      int lr_decision(0);
      if (i<3) lr_decision = LR[ktrial][i];
      HitMdc.LR(lr_decision);
      if (i<nhit_included)
	HitMdc.chi2(0);
      else 
	HitMdc.chi2(-1);
    }
    // Mdc fit the ... first hits :

    if(usage_==0) filter_fwd_anal(track, lead_, way, Eakal);
    way=999;

    // Check the chi2 
    if (debug_ == 4) 
      cout << "---- Result for " << ktrial << " case : chi2 = " << track.chiSq()
	<< ", nhits included = " << track.nchits() << " for nhits available = "
	<< nhit << std::endl;

    if (track.chiSq() < chi2_min && 
	(track.nchits() == nhit_included || track.nchits() == nhit)){
      chi2_min = track.chiSq();
      i_min = ktrial;
    }
  }

  // Initialize to the best solution :
  if (debug_ == 4) 
    cout << "*** i_min = " << i_min << " with a chi2 = " << chi2_min << std::endl;

  for( unsigned i=0 ; i < nhit; i++ ){
    KalFitHitMdc& HitMdc = track.HitMdc(i);
    int lr_decision(0);
    if (i_min >= 0 && i < 3) 
      lr_decision = LR[i_min][i];
    HitMdc.LR(lr_decision);
    HitMdc.chi2(0);
    HitMdc.chi2_back(0);
  }
  track.pivot(x_init);
  track.a(a_init);
  track.Ea(Ea_init);
  track.chiSq(0);
  track.chiSq_back(0);
  track.nchits(0);
  track.nster(0);
  track.ndf_back(0);

  // For debugging purpose :
  if (debug_ == 4) {
    for( unsigned i=0 ; i < 3; i++ ){
      KalFitHitMdc& HitMdc = track.HitMdc(i);
      cout << " LR(" << i << ") = " << HitMdc.LR()
	<< ", stereo = " << HitMdc.wire().stereo() 
	<< ", layer = " << HitMdc.wire().layer().layerId() << std::endl;
    }
  }
}

// this function is added to clear tables after processing each event
// to avoid memory leak,because of the usage of MdcTables etc.
//                  Apr. 2005
void KalFitAlg::clearTables( ) {

  if(debug_ == 4) cout<<"Begining to clear Tables ...."<<endl;
  vector<MdcRec_trk>* mdcMgr = MdcRecTrkCol::getMdcRecTrkCol();                 
  vector<MdcRec_trk_add>* mdc_addMgr = MdcRecTrkAddCol::getMdcRecTrkAddCol(); 
  vector<MdcRec_wirhit>* whMgr = MdcRecWirhitCol::getMdcRecWirhitCol();  
  vector<MdcRec_trk>::iterator tit=mdcMgr->begin();
  for( ; tit != mdcMgr->end(); tit++) {
    vector<MdcRec_wirhit*>::iterator vit= tit->hitcol.begin() ;
    for(; vit != tit->hitcol.end(); vit++) {
      delete (*vit);
    }
  }

  mdcMgr->clear();
  mdc_addMgr->clear();
  whMgr->clear();
}

bool  KalFitAlg::order_rechits(const SmartRef<RecMdcHit>& m1, const SmartRef<RecMdcHit>& m2)   {
  return  MdcID::layer(m1->getMdcId()) > MdcID::layer(m2->getMdcId());
}

void KalFitAlg::setRecMdcKalTrackMom(RecMdcKalTrack* kalTk,const HepVector helix, int l_mass){
  int m = l_mass;
  int charge;
  double pxy(0.),px(0.),py(0.),pz(0.),ptot(0.),dr(0.),phi0(0.),kappa(0.),dz(0.),tanl(0.),x(0.),y(0.),z(0.),vx0(0.),vy0(0.),vz0(0.);
  dr    = helix[0];
  phi0  = helix[1];
  kappa = helix[2];
  dz    = helix[3];
  tanl  = helix[4];

  x = dr * cos(phi0);
  y = dr * sin(phi0);
  z = dz;
  //  std::cout<<"x: "<<x<<"y: "<<y<<"z: "<<z<<std::endl;

  if (kappa > 0.0000000001)
    charge = 1;
  else if (kappa < -0.0000000001)
    charge = -1;
  else
    charge = 0;

  if(kappa!=0) pxy = 1.0/fabs(kappa);
  else pxy = 0;

  px = pxy * (-sin(phi0));
  py = pxy * cos(phi0);
  pz = pxy * tanl;
  ptot = sqrt(px*px+py*py+pz*pz);
  // cout<<"x, y, z, charge, px, py, pz, ptot: "<<x<<" , "<<y<<" , "<<z<<" , "<<charge<<" , "<<px<<" , "<<py<<" , "<<pz<<" , "<<ptot<<endl;
  //kalTk->setCharge(charge,m);
  kalTk->setPxy(pxy,m);
  kalTk->setPx(px,m);
  kalTk->setPy(py,m);
  kalTk->setPz(pz,m);
  kalTk->setP(ptot,m);
  kalTk->setTheta(acos(pz/ptot),m);
  kalTk->setPhi(atan2(py,px),m);
  kalTk->setX(x,m);
  kalTk->setY(y,m);
  kalTk->setZ(z,m);
  kalTk->setR(sqrt(x*x+y*y),m);
  /// the next three function need more consideration, temporarily

}//end of setRecMdcKalTrackMom()

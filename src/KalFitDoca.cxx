


#include "KalFitAlg/helix/Helix.h"
#include "KalFitAlg/KalFitHitMdc.h"
#include "KalFitAlg/KalFitWire.h"
//#include "KalFitAlg/KalFitTrack.h"


using namespace  KalmanFit;
// this approach function is used to calculate poca beteween from the helix and 
// the wire , well it also can deal with the wire sag.  
// this function transpalted from belle trasan packages.

double 
Helix::approach(KalFitHitMdc& hit, bool doSagCorrection) const {
    //...Cal. dPhi to rotate...
    //const TMDCWire & w = * link.wire();
    HepPoint3D positionOnWire, positionOnTrack;
    HepPoint3D pv = pivot();
    //std::cout<<"the track pivot in approach is "<<pv<<std::endl;
    HepVector Va = a();
    //std::cout<<"the track parameters is "<<Va<<std::endl;
    HepSymMatrix Ma = Ea(); 
    //std::cout<<"the error matrix is "<<Ma<<std::endl;

    Helix _helix(pv, Va ,Ma);
    //_helix.pivot(IP);
    const KalFitWire& w = hit.wire();
    Hep3Vector Wire = w.fwd() - w.bck(); 
    //xyPosition(), returns middle position of a wire. z componet is 0.
    //w.xyPosition(wp);
    double wp[3]; 
    wp[0] = 0.5*(w.fwd() + w.bck()).x();   
    wp[1] = 0.5*(w.fwd() + w.bck()).y();
    wp[2] = 0.;
    double wb[3]; 
    //w.backwardPosition(wb);
    wb[0] = w.bck().x();
    wb[1] = w.bck().y();
    wb[2] = w.bck().z();
    double v[3];
    v[0] = Wire.unit().x();
    v[1] = Wire.unit().y();
    v[2] = Wire.unit().z();
    //std::cout<<"Wire.unit() is "<<Wire.unit()<<std::endl; 
    
    //...Sag correction...
    /* if (doSagCorrection) {
	HepVector3D dir = w.direction();
	HepPoint3D xw(wp[0], wp[1], wp[2]);
	HepPoint3D wireBackwardPosition(wb[0], wb[1], wb[2]);
	w.wirePosition(link.positionOnTrack().z(),
		       xw,
		       wireBackwardPosition,
		       dir);
	v[0] = dir.x();
	v[1] = dir.y();
	v[2] = dir.z();
	wp[0] = xw.x();
	wp[1] = xw.y();
	wp[2] = xw.z();
	wb[0] = wireBackwardPosition.x();
	wb[1] = wireBackwardPosition.y();
	wb[2] = wireBackwardPosition.z();
     }
    */
    //...Cal. dPhi to rotate...
    const HepPoint3D & xc = _helix.center();

    //std::cout<<" helix center: "<<xc<<std::endl;
    
    double xt[3]; _helix.x(0., xt);
    double x0 = - xc.x();
    double y0 = - xc.y();
    double x1 = wp[0] + x0;
    double y1 = wp[1] + y0;
    x0 += xt[0];
    y0 += xt[1];
    double dPhi = atan2(x0 * y1 - y0 * x1, x0 * x1 + y0 * y1);
    
    //std::cout<<"dPhi is "<<dPhi<<std::endl;

    //...Setup...
    double kappa = _helix.kappa();
    double phi0 = _helix.phi0();

    //...Axial case...
  /*  if (!w.stereo()) {
	positionOnTrack = _helix.x(dPhi);
	HepPoint3D x(wp[0], wp[1], wp[2]);
	x.setZ(positionOnTrack.z());
	 positionOnWire = x;
     //link.dPhi(dPhi);
     std::cout<<" in axial wire : positionOnTrack is "<<positionOnTrack
              <<" positionOnWire is "<<positionOnWire<<std::endl;
         return (positionOnTrack - positionOnWire).mag();
    }
   */
    double firstdfdphi = 0.;
    static bool first = true;
    if (first) {
//	extern BelleTupleManager * BASF_Histogram;
//	BelleTupleManager * m = BASF_Histogram;
//	h_nTrial = m->histogram("TTrack::approach nTrial", 100, 0., 100.);
    }
//#endif

    //...Stereo case...
    double rho = Helix::ConstantAlpha / kappa;
    double tanLambda = _helix.tanl();
    static HepPoint3D x;
    double t_x[3];
    double t_dXdPhi[3];
    const double convergence = 1.0e-5;
    double l;
    unsigned nTrial = 0;
    while (nTrial < 100) {

//	x = link.positionOnTrack(_helix->x(dPhi));
        positionOnTrack = _helix.x(dPhi);
        x = _helix.x(dPhi);
	t_x[0] = x[0];
	t_x[1] = x[1];
	t_x[2] = x[2];

        l = v[0] * t_x[0] + v[1] * t_x[1] + v[2] * t_x[2]
	    - v[0] * wb[0] - v[1] * wb[1] - v[2] * wb[2];

	double rcosPhi = rho * cos(phi0 + dPhi);
	double rsinPhi = rho * sin(phi0 + dPhi);
	t_dXdPhi[0] =   rsinPhi;
	t_dXdPhi[1] = - rcosPhi;
	t_dXdPhi[2] = - rho * tanLambda;

	//...f = d(Distance) / d phi...
	double t_d2Xd2Phi[2];
	t_d2Xd2Phi[0] = rcosPhi;
	t_d2Xd2Phi[1] = rsinPhi;

	//...iw new...
	double n[3];
	n[0] = t_x[0] - wb[0];
	n[1] = t_x[1] - wb[1];
	n[2] = t_x[2] - wb[2];
	
	double a[3];
	a[0] = n[0] - l * v[0];
	a[1] = n[1] - l * v[1];
	a[2] = n[2] - l * v[2];
	double dfdphi = a[0] * t_dXdPhi[0]
	    + a[1] * t_dXdPhi[1]
	    + a[2] * t_dXdPhi[2];

//#ifdef TRKRECO_DEBUG
	if (nTrial == 0) {
//	    break;
	    firstdfdphi = dfdphi;
	}

	//...Check bad case...
	if (nTrial > 3) {
	    std::cout<<" bad case, calculate approach ntrial = "<<nTrial<< std::endl;
	}
//#endif

	//...Is it converged?...
	if (fabs(dfdphi) < convergence)
	    break;

	double dv = v[0] * t_dXdPhi[0]
	    + v[1] * t_dXdPhi[1]
	    + v[2] * t_dXdPhi[2];
	double t0 = t_dXdPhi[0] * t_dXdPhi[0]
	    + t_dXdPhi[1] * t_dXdPhi[1]
	    + t_dXdPhi[2] * t_dXdPhi[2];
	double d2fd2phi = t0 - dv * dv
	    + a[0] * t_d2Xd2Phi[0]
	    + a[1] * t_d2Xd2Phi[1];
//	    + a[2] * t_d2Xd2Phi[2];

	dPhi -= dfdphi / d2fd2phi;

//  	cout << "nTrial=" << nTrial << endl;
//    	cout << "iw f,df,dphi=" << dfdphi << "," << d2fd2phi << "," << dPhi << endl;

	++nTrial;
    }
    //...Cal. positions...
    positionOnWire[0] = wb[0] + l * v[0];
    positionOnWire[1] = wb[1] + l * v[1];
    positionOnWire[2] = wb[2] + l * v[2];
   
    //std::cout<<" positionOnTrack is "<<positionOnTrack
    //         <<" positionOnWire is "<<positionOnWire<<std::endl;
    return (positionOnTrack - positionOnWire).mag();

    //link.dPhi(dPhi);

    // #ifdef TRKRECO_DEBUG
    // h_nTrial->accumulate((float) nTrial + .5);
    // #endif
    // return nTrial;
}


double
Helix::approach(HepPoint3D pfwd, HepPoint3D pbwd,  bool doSagCorrection) const 
{
// ...Cal. dPhi to rotate...
// const TMDCWire & w = * link.wire();
    HepPoint3D positionOnWire, positionOnTrack;
    HepPoint3D pv = pivot();
    HepVector Va = a();
    HepSymMatrix Ma = Ea(); 

    Helix _helix(pv, Va ,Ma);
    Hep3Vector Wire;
    Wire[0] = (pfwd - pbwd).x();
    Wire[1] = (pfwd - pbwd).y();
    Wire[2] = (pfwd - pbwd).z(); 
// xyPosition(), returns middle position of a wire. z componet is 0.
// w.xyPosition(wp);
    double wp[3]; 
    wp[0] = 0.5*(pfwd + pbwd).x();   
    wp[1] = 0.5*(pfwd + pbwd).y();
    wp[2] = 0.;
    double wb[3]; 
// w.backwardPosition(wb);
    wb[0] = pbwd.x();
    wb[1] = pbwd.y();
    wb[2] = pbwd.z();
    double v[3];
    v[0] = Wire.unit().x();
    v[1] = Wire.unit().y();
    v[2] = Wire.unit().z();
// std::cout<<"Wire.unit() is "<<Wire.unit()<<std::endl; 

// ...Sag correction...
    /* if (doSagCorrection) {
	HepVector3D dir = w.direction();
	HepPoint3D xw(wp[0], wp[1], wp[2]);
	HepPoint3D wireBackwardPosition(wb[0], wb[1], wb[2]);
	w.wirePosition(link.positionOnTrack().z(),
		       xw,
		       wireBackwardPosition,
		       dir);
	v[0] = dir.x();
	v[1] = dir.y();
	v[2] = dir.z();
	wp[0] = xw.x();
	wp[1] = xw.y();
	wp[2] = xw.z();
	wb[0] = wireBackwardPosition.x();
	wb[1] = wireBackwardPosition.y();
	wb[2] = wireBackwardPosition.z();
     }
    */
    // ...Cal. dPhi to rotate...
    const HepPoint3D & xc = _helix.center();
    double xt[3];
     _helix.x(0., xt);
    double x0 = - xc.x();
    double y0 = - xc.y();
    double x1 = wp[0] + x0;
    double y1 = wp[1] + y0;
    x0 += xt[0];
    y0 += xt[1];
    //std::cout<<" x0 is: "<<x0<<" y0 is: "<<y0<<std::endl;
    //std::cout<<" x1 is: "<<x1<<" y1 is: "<<y1<<std::endl;
    //std::cout<<" xt[0] is: "<<xt[0]<<" xt[1] is: "<<xt[1]<<std::endl;

    double dPhi = atan2(x0 * y1 - y0 * x1, x0 * x1 + y0 * y1);
    //std::cout<<" x0 * y1 - y0 * x1 is: "<<(x0 * y1 - y0 * x1)<<std::endl;
    //std::cout<<" x0 * x1 + y0 * y1 is: "<<(x0 * x1 + y0 * y1)<<std::endl;
    //std::cout<<" before loop dPhi is "<<dPhi<<std::endl;
    //...Setup...
    double kappa = _helix.kappa();
    double phi0 = _helix.phi0();

    //...Axial case...
  /*  if (!w.stereo()) {
	positionOnTrack = _helix.x(dPhi);
	HepPoint3D x(wp[0], wp[1], wp[2]);
	x.setZ(positionOnTrack.z());
	 positionOnWire = x;
     //link.dPhi(dPhi);
     std::cout<<" in axial wire : positionOnTrack is "<<positionOnTrack
              <<" positionOnWire is "<<positionOnWire<<std::endl;
         return (positionOnTrack - positionOnWire).mag();
    }
   */
    double firstdfdphi = 0.;
    static bool first = true;
    if (first) {
//	extern BelleTupleManager * BASF_Histogram;
//	BelleTupleManager * m = BASF_Histogram;
//	h_nTrial = m->histogram("TTrack::approach nTrial", 100, 0., 100.);
    }
//#endif

    //...Stereo case...
    double rho = Helix::ConstantAlpha / kappa;
    double tanLambda = _helix.tanl();
    static HepPoint3D x;
    double t_x[3];
    double t_dXdPhi[3];
    const double convergence = 1.0e-5;
    double l;
    unsigned nTrial = 0;
    while (nTrial < 100) {

// x = link.positionOnTrack(_helix->x(dPhi));
        positionOnTrack = _helix.x(dPhi);
        x = _helix.x(dPhi);
	t_x[0] = x[0];
	t_x[1] = x[1];
	t_x[2] = x[2];

        l = v[0] * t_x[0] + v[1] * t_x[1] + v[2] * t_x[2]
	    - v[0] * wb[0] - v[1] * wb[1] - v[2] * wb[2];

	double rcosPhi = rho * cos(phi0 + dPhi);
	double rsinPhi = rho * sin(phi0 + dPhi);
	t_dXdPhi[0] =   rsinPhi;
	t_dXdPhi[1] = - rcosPhi;
	t_dXdPhi[2] = - rho * tanLambda;

	//...f = d(Distance) / d phi...
	double t_d2Xd2Phi[2];
	t_d2Xd2Phi[0] = rcosPhi;
	t_d2Xd2Phi[1] = rsinPhi;

	//...iw new...
	double n[3];
	n[0] = t_x[0] - wb[0];
	n[1] = t_x[1] - wb[1];
	n[2] = t_x[2] - wb[2];
	
	double a[3];
	a[0] = n[0] - l * v[0];
	a[1] = n[1] - l * v[1];
	a[2] = n[2] - l * v[2];
	double dfdphi = a[0] * t_dXdPhi[0]
	    + a[1] * t_dXdPhi[1]
	    + a[2] * t_dXdPhi[2];

	if (nTrial == 0) {
//	    break;
	    firstdfdphi = dfdphi;
	}

	//...Check bad case...
	if (nTrial > 3) {
//	    std::cout<<" BAD CASE!!, calculate approach ntrial = "<<nTrial<< endl;
	}
	//...Is it converged?...
	if (fabs(dfdphi) < convergence)
	    break;

	double dv = v[0] * t_dXdPhi[0]
	    + v[1] * t_dXdPhi[1]
	    + v[2] * t_dXdPhi[2];
	double t0 = t_dXdPhi[0] * t_dXdPhi[0]
	    + t_dXdPhi[1] * t_dXdPhi[1]
	    + t_dXdPhi[2] * t_dXdPhi[2];
	double d2fd2phi = t0 - dv * dv
	    + a[0] * t_d2Xd2Phi[0]
	    + a[1] * t_d2Xd2Phi[1];
        //  + a[2] * t_d2Xd2Phi[2];

	dPhi -= dfdphi / d2fd2phi;
	++nTrial;
    }
    //std::cout<<" dPhi is: "<<dPhi<<std::endl;
    //...Cal. positions...
    positionOnWire[0] = wb[0] + l * v[0];
    positionOnWire[1] = wb[1] + l * v[1];
    positionOnWire[2] = wb[2] + l * v[2];

    //std::cout<<"wb[0] is: "<<wb[0]<<" l is: "<<l<<" v[0] is: "<<v[0]<<std::endl;
    //std::cout<<"wb[1] is: "<<wb[1]<<" v[1] is: "<<v[1]<<std::endl;
    //std::cout<<"wb[2] is: "<<wb[2]<<" v[2] is: "<<v[2]<<std::endl;

    //std::cout<<" positionOnTrack is "<<positionOnTrack
    //         <<" positionOnWire is "<<positionOnWire<<std::endl;

    return (positionOnTrack - positionOnWire).mag();
    // link.dPhi(dPhi);
    // return nTrial;
}

//
// Bfield class
//
// 27-Mar-1999 : KUNIYA Toshio
//   Enabled QCS compornent
//
// 21-Feb-1999 : KUNIYA Toshio
//   Keeping comatibility, Bfield class is modified.
//   No longer fortran common block is used for bfield map.
//   Access functions are prepared for fortran call.
//


#include "KalFitAlg/coil/Bfield.h"


#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Geometry/Point3D.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
 typedef HepGeom::Point3D<double> HepPoint3D;
#endif



using CLHEP::HepVector; 
using CLHEP::Hep3Vector;
using CLHEP::HepMatrix;
using CLHEP::HepSymMatrix;


// prototype of file scoped function to call fortran routine
// ... read B field map with ID = imap from a file
/*wangdy
extern "C" {
  extern void geo_coil_readmap_
    (int *imap, double (*bz)[399], double (*br)[399], double (*bphi)[399]);
  extern void geo_coil_readqcsrmap_
    (double (*bzqr)[163], double (*brqr)[163], double (*bphiqr)[163]);
  extern void geo_coil_readqcslmap_
    (double (*bzrl)[51][52], double (*brql)[51][52], double (*bphiqr)[51][52]);
}
 */ 


Bfield *
Bfield::_field[200] = {0};   // All of 200 elements are initialized.

Bfield *
Bfield::getBfield(int imap) {
  if (! _field[imap]) _field[imap] = new Bfield(imap);
  return _field[imap];
}


Bfield::Bfield(int imap) {
  std::cout << std::endl;
  std::cout << "***********************************************" << std::endl;
  std::cout << "           Bfield class  MAP ID = " << imap << std::endl;
  std::cout << "    #### R < 174 cm, -152 < Z < 246 cm ####    " << std::endl;
  std::cout << "                C++ version 1.00               " << std::endl;
  std::cout << "***********************************************" << std::endl;

  const double uniformBz[10] = {0, 15, 14.5, 14, 13, 12.5, 12, 11, 10, 15.5};

  //...initialization
  for(int i=0; i<175; i++)
    for(int j=0; j<399; j++) {
      _Bz[i][j] = 0;
      _Br[i][j] = 0;
      _Bz[i][j] = 0;
      if(i<101 && j<163) {
	_BzQR[i][j] = 0;
	_BrQR[i][j] = 0;
	_BphiQR[i][j] = 0;
      }
    }
  for(int i=0; i<17; i++)
    for(int j=0; j<51; j++)
      for(int k=0; k<52; k++) {
	_BzQL[i][j][k] = 0;
	_BrQL[i][j][k] = 0;
	_BphiQL[i][j][k] =0;
      }
  
  //...
  _fieldID  = imap;

  //...read B field map 

  if(imap<10){
    //
    // uniform B field map
    //
    m_Bx = 0.;
    m_By = 0.;
    m_Bz = uniformBz[imap];
    m_Bfld.setX((double) m_Bx);
    m_Bfld.setY((double) m_By);
    m_Bfld.setZ((double) m_Bz);
    std::cout << "Bfield class >> creating uniform B field with id = " << imap;
    std::cout << ", (Bx,By,Bz) = "<<m_Bfld<<std::endl;
  } else {
    //
    // non-uniform B field map
    //
    std::cout << "Bfield class >> loading non-uniform B field map" << std::endl;
/*wangdy
     geo_coil_readmap_(&imap, _Bz, _Br, _Bphi);
    if( _fieldID == 21 ) {
      std::cout << "Bfield class >> loading QCS" << std::endl;
      geo_coil_readqcsrmap_(_BzQR,_BrQR, _BphiQR);
      geo_coil_readqcslmap_(_BzQL,_BrQL, _BphiQL);
    }
 */ 
    updateCache(0., 0., 0.);
  }
  std::cout << std::endl;

}

//Bfield::~Bfield(){};

const Hep3Vector &
Bfield::fieldMap(double x, double y, double z) const {

  if(_fieldID > 10){
    if(x != m_x || y != m_y || z != m_z)  updateCache(x, y, z);
  }

  return m_Bfld;
}

const Hep3Vector &
Bfield::fieldMap(const HepPoint3D &xyz) const{

  if(_fieldID > 10){
    double x = xyz.x();
    double y = xyz.y();
    double z = xyz.z();
    if(x != m_x || y != m_y || z != m_z)  updateCache(x, y, z);
  }

  return m_Bfld;
}

void
Bfield::fieldMap(double *position, double *field) {
  if(_fieldID > 10){
    double x = position[0];
    double y = position[1];
    double z = position[2];
    if(x != m_x || y != m_y || z != m_z)  updateCache(x, y, z);
  }
  field[0] = m_Bx;
  field[1] = m_By;
  field[2] = m_Bz;
  return;
}

double
Bfield::bx(double x, double y, double z) const {

  if(_fieldID > 10){
    if(x != m_x || y != m_y || z != m_z)  updateCache(x, y, z);
  }

  return m_Bx;
}

double
Bfield::by(double x, double y, double z) const {

  if(_fieldID > 10){
    if(x != m_x || y != m_y || z != m_z)  updateCache(x, y, z);
  }

  return m_By;
}

double
Bfield::bz(double x, double y, double z) const {

  if(_fieldID > 10){
    if(x != m_x || y != m_y || z != m_z)  updateCache(x, y, z);
  }

  return m_Bz;
}

double
Bfield::bx(const HepPoint3D &xyz) const{

  if(_fieldID > 10){
    double x = xyz.x();
    double y = xyz.y();
    double z = xyz.z();
    if(x != m_x || y != m_y || z != m_z)  updateCache(x, y, z);
  }

  return m_Bx;
}

double
Bfield::by(const HepPoint3D &xyz) const{

  if(_fieldID > 10){
    double x = xyz.x();
    double y = xyz.y();
    double z = xyz.z();
    if(x != m_x || y != m_y || z != m_z)  updateCache(x, y, z);
  }

  return m_By;
}

double
Bfield::bz(const HepPoint3D &xyz) const{

  if(_fieldID > 10){
    double x = xyz.x();
    double y = xyz.y();
    double z = xyz.z();
    if(x != m_x || y != m_y || z != m_z)  updateCache(x, y, z);
  }
  return m_Bz;
}

void
Bfield::updateCache(double x, double y, double z) const{

  // this function is only for non-uniform B field

  if( _fieldID <= 10 ) return;

      double PI = 3.14159;

      //...
      double  r   = (double)sqrt((double)x*(double)x + (double)y*(double)y);
      double  phi = (double)atan2((double)y, (double)x);

      //... [cm] --> [mm]
      double zmm = z * 10.;
      double rmm = r * 10.;
      //... make index
      int tz = (int) (( zmm + 1520.)/10.);
      int tr = (int) (rmm/10.);

      //... 
      double bz = 0., br = 0., bphi = 0.;

      if(zmm >= -1520. && zmm < 2460. && rmm >= 0. && rmm < 1740.){
	if(_Bz[tr][tz] && _Bz[tr][tz+1]){
       	    double pz = (zmm + 1520.)/10.- (double)tz;
            double pr =  rmm/10.- (double)tr;

            //...
            bz = (_Bz[tr][tz]*(1.- pz)+_Bz[tr][tz+1]*pz)*(1.-pr)+
                 (_Bz[tr+1][tz]*(1.-pz)+_Bz[tr+1][tz+1]*pz)*pr;
            //...
            br = (_Br[tr][tz]*(1.-pz)+_Br[tr][tz+1]*pz)*(1.-pr)+
                 (_Br[tr+1][tz]*(1.-pz)+_Br[tr+1][tz+1]*pz)*pr;
            //...
            bphi = (_Bphi[tr][tz]*(1.-pz)+_Bphi[tr][tz+1]*pz)*(1.-pr)+
                   (_Bphi[tr+1][tz]*(1.-pz)+_Bphi[tr+1][tz+1]*pz)*pr;

	    if(_fieldID == 21) {
	      //
	      // QCS Right
	      //
	      if( zmm>=800. && zmm < 2420. && rmm < 1000. ) {
		int tqz = (int)( (zmm-800.)/10. );
		bz += (((_BzQR[tr][tqz]*(1.-pz)+_BzQR[tr][tqz+1]*pz)*(1.-pr)
		      +(_BzQR[tr+1][tqz]*(1.-pz)+_BzQR[tr+1][tqz+1]*pz)*pr)
		      *(double)sin((double)(2.*phi+2./180.*(double)PI)));
		br += (((_BrQR[tr][tqz]*(1.-pz)+_BrQR[tr][tqz+1]*pz)*(1.-pr)
		      +(_BrQR[tr+1][tqz]*(1.-pz)+_BrQR[tr+1][tqz+1]*pz)*pr)
		      *(double)sin((double)(2.*phi+2./180.*(double)PI)));
		bphi += (((_BphiQR[tr][tqz]*(1.-pz)
			+_BphiQR[tr][tqz+1]*pz)*(1.-pr)
			+(_BphiQR[tr+1][tqz]*(1.-pz)
			+_BphiQR[tr+1][tqz+1]*pz)*pr)
			*(double)cos((double)(2.*phi+2./180.*(double)PI)));
	      }
	      //
	      // QCS Left
	      //
	      if(zmm<=-500. && zmm>-1520. && rmm<1000.) {
		int tqz = (int)((-zmm-500.)/20.);
		int tqr = (int)(tr/2.);
		pz = (pz+(double)( tz-(int)(tz/2.)*2. ))/2.;
		pr = ( pr + (double)(tr-tqr*2) )/2.;
		double f = 1.;
		//		if( phi < (PI/2.) && phi >= (-PI/2.) ) {
		//		  phi = PI-phi;
		//		  f =-1.;
		//		} else if( phi < -PI/2. ) {
		//		  phi = 2.*PI+phi;
		//		}
		//		double pphi = ( phi-PI/2. )/(11.25/180.*PI);
		double phi_tmp = phi;
		if( phi_tmp < (PI/2.) && phi_tmp >= (-PI/2.) ) {
		  phi_tmp = PI-phi_tmp;
		  f =-1.;
		} else if( phi_tmp < -PI/2. ) {
		  phi_tmp = 2.*PI+phi_tmp;
		}
		double pphi = ( phi_tmp-PI/2. )/(11.25/180.*PI);
		int tphi = (int)pphi;
		pphi -= (double)tphi;
		if (tphi >= 16) tphi -= 16;
		
		bz += f*
		  (((_BzQL[tphi][tqr][tqz]*(1.-pz)+_BzQL[tphi][tqr][tqz+1]*pz)
		  *(1-pr)+(_BzQL[tphi][tqr+1][tqz]*(1.-pz)
		  +_BzQL[tphi][tqr+1][tqz+1]*pz)*pr)*(1.-pphi)
		  +((_BzQL[tphi+1][tqr][tqz]*(1.-pz)
		  +_BzQL[tphi+1][tqr][tqz+1]*pz)*(1.-pr)
		  +(_BzQL[tphi+1][tqr+1][tqz]*(1.-pz)
		  +_BzQL[tphi+1][tqr+1][tqz+1]*pz)*pr)*pphi);
		br += f*
		  (((_BrQL[tphi][tqr][tqz]*(1.- pz)
		  +_BrQL[tphi][tqr][tqz+1]*pz)*(1.-pr)
		  +(_BrQL[tphi][tqr+1][tqz]*(1.-pz)
		  +_BrQL[tphi][tqr+1][tqz+1]*pz)*pr)*(1-pphi)
		  +((_BrQL[tphi+1][tqr][tqz]*(1.- pz)
		  +_BrQL[tphi+1][tqr][tqz+1]*pz)*(1.-pr)
		  +(_BrQL[tphi+1][tqr+1][tqz]*(1.-pz)
		  +_BrQL[tphi+1][tqr+1][tqz+1]*pz)*pr)*pphi);
		bphi += f*
		  (((_BphiQL[tphi][tqr][tqz]*(1.- pz)
		  +_BphiQL[tphi][tqr][tqz+1]*pz)*(1.-pr)
		  +(_BphiQL[tphi][tqr+1][tqz]*(1.-pz)
		  +_BphiQL[tphi][tqr+1][tqz+1]*pz)*pr)*(1.-pphi)
		  +((_BphiQL[tphi+1][tqr][tqz]*(1.- pz)
		  +_BphiQL[tphi+1][tqr][tqz+1]*pz)*(1.-pr)
		  +(_BphiQL[tphi+1][tqr+1][tqz]*(1.-pz)
		  +_BphiQL[tphi+1][tqr+1][tqz+1]*pz)*pr)*pphi);
	      }
	    }
         }  else {
           bz=0.;
           br=0.;
           bphi=0.;
         }
      } else if(zmm>=-1520. && zmm<=2460. && rmm<=0. && rmm<=1740.){
	if(tz == 246) tz=tz-1;
	if(tr == 174) tr=tr-1;
	double pz= (zmm + 1520.)/10.- (double)tz;
	double pr=  rmm/10.- (double)tr;
	bz   = (_Bz[tr][tz]*(1.- pz)+_Bz[tr][tz+1]*pz)*(1.-pr)+
	       (_Bz[tr+1][tz]*(1.-pz)+_Bz[tr+1][tz+1]*pz)*pr;

	br   = (_Br[tr][tz]*(1.-pz)+_Br[tr][tz+1]*pz)*(1.-pr)+
	       (_Br[tr+1][tz]*(1.-pz)+_Br[tr+1][tz+1]*pz)*pr;

	bphi = (_Bphi[tr][tz]*(1.-pz)+_Bphi[tr][tz+1]*pz)*(1.-pr)+
	       (_Bphi[tr+1][tz]*(1.-pz)+_Bphi[tr+1][tz+1]*pz)*pr;
      } else {
	bz =0.;
	br =0.;
	bphi =0.;
      }
      

      //... Set B field
      double Bmag_xy = (double)sqrt(br*br + bphi*bphi);
      double Bphi_rp = (double)atan2( (double)bphi, (double)br );
      m_Bx = Bmag_xy * (double)cos((double)phi + Bphi_rp)/1000.;
      m_By = Bmag_xy * (double)sin((double)phi + Bphi_rp)/1000.;   
      //m_Bx = br * (double)cos((double)phi)/1000.;
      //m_By = br * (double)sin((double)phi)/1000.;   
      m_Bz = bz/1000.;
      m_x = x;
      m_y = y;
      m_z = z;
      m_Bfld.setX((double) m_Bx);
      m_Bfld.setY((double) m_By);
      m_Bfld.setZ((double) m_Bz);
}


//
//... Fortran callable access functions to Bfield class.
//

// initialize bfield map
extern "C"
void init_bfield_(int *imap) {
  // create Bfiled map ID = imap
  //  Bfield *thisMap = Bfield::getBfield(*imap);
  (void) Bfield::getBfield(*imap);
  // It is OK even though 'thisMap' losts its scope at here.
  // Because address of field map is not deleted from memory
  // due to static linkaged Bfield class.
  return;
}

// retrieving B field corresponding to the POSition
extern "C"
void get_bfield_(int *imap, double *pos, double *field, int *error) {
  Bfield *thisMap = Bfield::getBfield(*imap);
  // std::cout << " > accessing Bfield class from fortran routine." << std::endl;
  if( thisMap != 0 ) {
    thisMap->fieldMap(pos,field);
    *error = 0;
  }
  else *error = 1;
  return;
}

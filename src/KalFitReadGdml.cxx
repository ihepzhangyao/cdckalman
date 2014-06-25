

#include "G4Geo/MdcG4Geo.h"
#include "G4Geo/BesG4Geo.h"
#include "KalFitAlg/KalFitAlg.h"
#include "KalFitAlg/KalFitTrack.h"
#include "G4Material.hh"
#include "G4Tubs.hh"
#include "GDMLProcessor.hh"


void KalFitAlg::setBesFromGdml(void){

	int i(0);
	double Z(0.),A(0.),Ionization(0.),Density(0.),Radlen(0.);

	G4LogicalVolume *logicalMdc = 0;
	MdcG4Geo* aMdcG4Geo = new MdcG4Geo();
	logicalMdc = aMdcG4Geo->GetTopVolume();   

	/// mdcgas
	G4Material* mdcMaterial = logicalMdc->GetMaterial();  

	for(i=0; i<mdcMaterial->GetElementVector()->size(); i++){
		Z += (mdcMaterial->GetElement(i)->GetZ())*
			(mdcMaterial->GetFractionVector()[i]);
		A += (mdcMaterial->GetElement(i)->GetA())*
			(mdcMaterial->GetFractionVector()[i]);
	}
	Ionization = mdcMaterial->GetIonisation()->GetMeanExcitationEnergy();
	Density = mdcMaterial->GetDensity()/(g/cm3);
	Radlen = mdcMaterial->GetRadlen();
	std::cout<<"mdcgas: Z: "<<Z<<" A: "<<(A/(g/mole))<<" Ionization: "<<(Ionization/eV)<<" Density: "<<Density<<" Radlen: "<<Radlen<<std::endl;
	KalFitMaterial FitMdcMaterial(Z,A/g/mole,Ionization/eV,Density,Radlen/10.); 
	_BesKalmanFitMaterials.push_back(FitMdcMaterial);
	KalFitTrack::mdcGasRadlen_ = Radlen/10.;

	/// inner wall shield fiml1 Al by wangll 2012-09-07
	G4LogicalVolume* innerWallFilm1Volume = const_cast<G4LogicalVolume*>(GDMLProcessor::GetInstance()->GetLogicalVolume("LogicalMdcInnerFilm1"));
	G4Material* innerWallFilm1Material = innerWallFilm1Volume->GetMaterial();
	G4Tubs* innerwallFilm1Tub = dynamic_cast<G4Tubs*>(innerWallFilm1Volume->GetSolid());

	Z = 0.;
	A = 0.;
	for(i=0; i<innerWallFilm1Material->GetElementVector()->size(); i++){    
		Z += (innerWallFilm1Material->GetElement(i)->GetZ())*
			(innerWallFilm1Material->GetFractionVector()[i]); 
		A += (innerWallFilm1Material->GetElement(i)->GetA())*	
			(innerWallFilm1Material->GetFractionVector()[i]);
	}

	Ionization = innerWallFilm1Material->GetIonisation()->GetMeanExcitationEnergy();
	Density = innerWallFilm1Material->GetDensity()/(g/cm3);
	Radlen = innerWallFilm1Material->GetRadlen();
	std::cout<<"Mdc innerwall Film1, Al: Z: "<<Z<<" A: "<<(A/(g/mole))<<" Ionization: "<<(Ionization/eV)<<" Density: "<<Density<<" Radlen: "<<Radlen<<std::endl;
	KalFitMaterial FitInnerwallFilm1Material(Z,A/g/mole,Ionization/eV,Density,Radlen/10.);
	_BesKalmanFitMaterials.push_back(FitInnerwallFilm1Material);


	/// inner wall CarbonFiber by wll 2012-09-06
	G4LogicalVolume* innerwallVolume = const_cast<G4LogicalVolume*>(GDMLProcessor::GetInstance()->GetLogicalVolume("logicalMdcSegment2"));
	G4Material* innerwallMaterial = innerwallVolume->GetMaterial();
	G4Tubs* innerwallTub = dynamic_cast<G4Tubs*>(innerwallVolume->GetSolid());

	Z = 0.;
	A = 0.;
	for(i=0; i<innerwallMaterial->GetElementVector()->size(); i++){    
		Z += (innerwallMaterial->GetElement(i)->GetZ())*
			(innerwallMaterial->GetFractionVector()[i]); 
		A += (innerwallMaterial->GetElement(i)->GetA())*	
			(innerwallMaterial->GetFractionVector()[i]);
	}

	Ionization = innerwallMaterial->GetIonisation()->GetMeanExcitationEnergy();
	Density = innerwallMaterial->GetDensity()/(g/cm3);
	Radlen = innerwallMaterial->GetRadlen();
	std::cout<<"Mdc innerwall, Al: Z: "<<Z<<" A: "<<(A/(g/mole))<<" Ionization: "<<(Ionization/eV)<<" Density: "<<Density<<" Radlen: "<<Radlen<<std::endl;
	KalFitMaterial FitInnerwallMaterial(Z,A/g/mole,Ionization/eV,Density,Radlen/10.);
	_BesKalmanFitMaterials.push_back(FitInnerwallMaterial);

	/// inner wall shield film0 Al by wangll 2012-09-07
	G4LogicalVolume* innerWallFilm0Volume = const_cast<G4LogicalVolume*>(GDMLProcessor::GetInstance()->GetLogicalVolume("LogicalMdcInnerFilm0"));
	G4Material* innerWallFilm0Material = innerWallFilm0Volume->GetMaterial();
	G4Tubs* innerwallFilm0Tub = dynamic_cast<G4Tubs*>(innerWallFilm0Volume->GetSolid());

	Z = 0.;
	A = 0.;
	for(i=0; i<innerWallFilm0Material->GetElementVector()->size(); i++){    
		Z += (innerWallFilm0Material->GetElement(i)->GetZ())*
			(innerWallFilm0Material->GetFractionVector()[i]); 
		A += (innerWallFilm0Material->GetElement(i)->GetA())*	
			(innerWallFilm0Material->GetFractionVector()[i]);
	}

	Ionization = innerWallFilm0Material->GetIonisation()->GetMeanExcitationEnergy();
	Density = innerWallFilm0Material->GetDensity()/(g/cm3);
	Radlen = innerWallFilm0Material->GetRadlen();
	std::cout<<"Mdc innerwall Film0, Al: Z: "<<Z<<" A: "<<(A/(g/mole))<<" Ionization: "<<(Ionization/eV)<<" Density: "<<Density<<" Radlen: "<<Radlen<<std::endl;
	KalFitMaterial FitInnerwallFilm0Material(Z,A/g/mole,Ionization/eV,Density,Radlen/10.);
	_BesKalmanFitMaterials.push_back(FitInnerwallFilm0Material);

	///////////////////////////////////////////////////////////////////////////////////////////////////  
	G4LogicalVolume *logicalBes = 0;
	BesG4Geo* aBesG4Geo = new BesG4Geo();
	logicalBes = aBesG4Geo->GetTopVolume();

	/// air
	G4LogicalVolume* logicalAirVolume = const_cast<G4LogicalVolume*>(GDMLProcessor::GetInstance()->GetLogicalVolume("logicalWorld"));
	G4Material* airMaterial = logicalAirVolume->GetMaterial();
	Z = 0.;
	A = 0.;
	for(i=0; i<airMaterial->GetElementVector()->size(); i++){
		Z += (airMaterial->GetElement(i)->GetZ())*
			(airMaterial->GetFractionVector()[i]);
		A += (airMaterial->GetElement(i)->GetA())*
			(airMaterial->GetFractionVector()[i]);
	}

	Ionization = airMaterial->GetIonisation()->GetMeanExcitationEnergy();
	Density = airMaterial->GetDensity()/(g/cm3);
	Radlen = airMaterial->GetRadlen();
	std::cout<<"air: Z: "<<Z<<" A: "<<(A/(g/mole))<<" Ionization: "<<(Ionization/eV)<<" Density: "<<Density<<" Radlen: "<<Radlen<<std::endl;
	KalFitMaterial FitAirMaterial(Z,A/g/mole,Ionization/eV,Density,Radlen/10.);
	_BesKalmanFitMaterials.push_back(FitAirMaterial);

	/// outer beryllium pipe 
	G4LogicalVolume* logicalOuterBeVolume = const_cast<G4LogicalVolume*>(GDMLProcessor::GetInstance()->GetLogicalVolume("logicalouterBe"));
	G4Material* outerBeMaterial = logicalOuterBeVolume->GetMaterial();

	G4Tubs* outerBeTub = dynamic_cast<G4Tubs*>(logicalOuterBeVolume->GetSolid());
	Z = 0.;
	A = 0.;
	for(i=0; i<outerBeMaterial->GetElementVector()->size(); i++){   
		Z += (outerBeMaterial->GetElement(i)->GetZ())*
			(outerBeMaterial->GetFractionVector()[i]);   
		A += (outerBeMaterial->GetElement(i)->GetA())*   
			(outerBeMaterial->GetFractionVector()[i]);
	}
	Ionization =  outerBeMaterial->GetIonisation()->GetMeanExcitationEnergy();
	Density = outerBeMaterial->GetDensity()/(g/cm3);
	Radlen = outerBeMaterial->GetRadlen();
	std::cout<<"outer beryllium: Z: "<<Z<<" A: "<<(A/(g/mole))<<" Ionization: "<<(Ionization/eV)<<" Density: "<<Density<<" Radlen: "<<Radlen<<std::endl;
	KalFitMaterial FitOuterBeMaterial(Z,A/g/mole,Ionization/eV,Density,Radlen/10.);
	_BesKalmanFitMaterials.push_back(FitOuterBeMaterial);


	/// cooling oil 
	G4LogicalVolume* logicalOilLayerVolume = const_cast<G4LogicalVolume*>(GDMLProcessor::GetInstance()->GetLogicalVolume("logicaloilLayer"));
	G4Material* oilLayerMaterial = logicalOilLayerVolume->GetMaterial();
	G4Tubs* oilLayerTub = dynamic_cast<G4Tubs*>(logicalOilLayerVolume->GetSolid());

	Z = 0.;
	A = 0.;
	for(i=0; i<oilLayerMaterial->GetElementVector()->size(); i++){        
		Z += (oilLayerMaterial->GetElement(i)->GetZ())*
			(oilLayerMaterial->GetFractionVector()[i]);             
		A += (oilLayerMaterial->GetElement(i)->GetA())*             
			(oilLayerMaterial->GetFractionVector()[i]);
	}
	Ionization = oilLayerMaterial->GetIonisation()->GetMeanExcitationEnergy();
	Density = oilLayerMaterial->GetDensity()/(g/cm3);
	Radlen = oilLayerMaterial->GetRadlen();
	std::cout<<"cooling oil: Z: "<<Z<<" A: "<<(A/(g/mole))<<" Ionization: "<<(Ionization/eV)<<" Density: "<<Density<<" Radlen: "<<Radlen<<std::endl;
	KalFitMaterial FitOilLayerMaterial(Z,A/g/mole,Ionization/eV,Density,Radlen/10.);
	_BesKalmanFitMaterials.push_back(FitOilLayerMaterial);


	/// inner beryllium pipe 
	G4LogicalVolume* logicalInnerBeVolume = const_cast<G4LogicalVolume*>(GDMLProcessor::GetInstance()->GetLogicalVolume("logicalinnerBe"));
	G4Material* innerBeMaterial = logicalInnerBeVolume->GetMaterial();
	G4Tubs* innerBeTub = dynamic_cast<G4Tubs*>(logicalInnerBeVolume->GetSolid());
	Z = 0.;
	A = 0.;
	for(i=0; i<innerBeMaterial->GetElementVector()->size(); i++){
		Z += (innerBeMaterial->GetElement(i)->GetZ())*
			(innerBeMaterial->GetFractionVector()[i]);
		A += (innerBeMaterial->GetElement(i)->GetA())*
			(innerBeMaterial->GetFractionVector()[i]);
	}

	Ionization = innerBeMaterial->GetIonisation()->GetMeanExcitationEnergy();
	Density = innerBeMaterial->GetDensity()/(g/cm3);
	Radlen = innerBeMaterial->GetRadlen();
	std::cout<<"inner beryllium: Z: "<<Z<<" A: "<<(A/(g/mole))<<" Ionization: "<<(Ionization/eV)<<" Density: "<<Density<<" Radlen: "<<Radlen<<std::endl;
	KalFitMaterial FitInnerBeMaterial(Z,A/g/mole,Ionization/eV,Density,Radlen/10.);
	_BesKalmanFitMaterials.push_back(FitInnerBeMaterial);


	/// gold
	G4LogicalVolume* logicalGoldLayerVolume = const_cast<G4LogicalVolume*>(GDMLProcessor::GetInstance()->GetLogicalVolume("logicalgoldLayer"));
	G4Material* goldLayerMaterial = logicalGoldLayerVolume->GetMaterial();
	G4Tubs* goldLayerTub = dynamic_cast<G4Tubs*>(logicalGoldLayerVolume->GetSolid());

	Z = 0.;
	A = 0.;
	for(i=0; i<goldLayerMaterial->GetElementVector()->size(); i++){
		Z += (goldLayerMaterial->GetElement(i)->GetZ())*
			(goldLayerMaterial->GetFractionVector()[i]);
		A += (goldLayerMaterial->GetElement(i)->GetA())*
			(goldLayerMaterial->GetFractionVector()[i]);
	}
	Ionization = goldLayerMaterial->GetIonisation()->GetMeanExcitationEnergy();
	Density = goldLayerMaterial->GetDensity()/(g/cm3);
	Radlen = goldLayerMaterial->GetRadlen();
	std::cout<<"gold layer: Z: "<<Z<<" A: "<<(A/(g/mole))<<" Ionization: "<<(Ionization/eV)<<" Density: "<<Density<<" Radlen: "<<Radlen<<std::endl;
	KalFitMaterial FitGoldLayerMaterial(Z,A/g/mole,Ionization/eV,Density,Radlen/10.);
	_BesKalmanFitMaterials.push_back(FitGoldLayerMaterial);
	
	if(testCDC>0)
	{
		// CDC gas (He:C4H10=90:10) by wuchen 2013-06-23
		Z = 3.88432;
		A = 7.67152*g/mole;
		Ionization = 48.6658*eV;
		Density = 0.000403321;//*g/cm3
		Radlen = 1414900*mm;
		std::cout<<"CDC gas: Z: "<<Z<<" A: "<<(A/(g/mole))<<" Ionization: "<<(Ionization/eV)<<" Density: "<<Density<<" Radlen: "<<Radlen<<std::endl;
		KalFitMaterial FitCdcMaterial(Z,A/g/mole,Ionization/eV,Density,Radlen/10.); 
		_BesKalmanFitMaterials.clear();
		_BesKalmanFitMaterials.push_back(FitCdcMaterial);// 0
		KalFitTrack::mdcGasRadlen_ = Radlen/10.;
		
		/// CDC inner wall shield (Carbon Fiber) by wuchen 2013-06-23
		Z = 6.5633;
		A = 13.1278*g/mole;
		Ionization = 83.4461*eV;
		Density = 1.57;// *g/cm3
		Radlen = 253.835*mm;
		std::cout<<"CDC innerwall : Z: "<<Z<<" A: "<<(A/(g/mole))<<" Ionization: "<<(Ionization/eV)<<" Density: "<<Density<<" Radlen: "<<Radlen<<std::endl;
		KalFitMaterial FitInnerwallFilm1Material(Z,A/g/mole,Ionization/eV,Density,Radlen/10.);
		_BesKalmanFitMaterials.push_back(FitInnerwallFilm1Material);// 1
	}


	/// now construct the cylinders
	double radius, thick, length , z0;


	/// film1 of the innerwall of inner drift chamber
	radius = innerwallFilm1Tub->GetInnerRadius()/(cm);
	thick  = innerwallFilm1Tub->GetOuterRadius()/(cm) - innerwallFilm1Tub->GetInnerRadius()/(cm);
	length = 2.0*innerwallFilm1Tub->GetZHalfLength()/(cm);
	z0     = 0.0;
	std::cout<<"innerwallFilm1: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl;
	KalFitCylinder innerwallFilm1Cylinder(&_BesKalmanFitMaterials[1], radius, thick, length , z0);
	_BesKalmanFitWalls.push_back(innerwallFilm1Cylinder);


	/// innerwall of inner drift chamber
	radius = innerwallTub->GetInnerRadius()/(cm);
	thick  = innerwallTub->GetOuterRadius()/(cm) - innerwallTub->GetInnerRadius()/(cm);
	length = 2.0*innerwallTub->GetZHalfLength()/(cm);
	z0     = 0.0;
	std::cout<<"innerwall: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl;
	KalFitCylinder innerwallCylinder(&_BesKalmanFitMaterials[2], radius, thick, length , z0);
	_BesKalmanFitWalls.push_back(innerwallCylinder);

	/// film0 of the innerwall of inner drift chamber
	radius = innerwallFilm0Tub->GetInnerRadius()/(cm);
	thick  = innerwallFilm0Tub->GetOuterRadius()/(cm) - innerwallFilm0Tub->GetInnerRadius()/(cm);
	length = 2.0*innerwallFilm0Tub->GetZHalfLength()/(cm);
	z0     = 0.0;
	std::cout<<"innerwallFilm0: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl;
	KalFitCylinder innerwallFilm0Cylinder(&_BesKalmanFitMaterials[3], radius, thick, length , z0);
	_BesKalmanFitWalls.push_back(innerwallFilm0Cylinder);

	/// outer air, be attention the calculation of the radius and thick of the air cylinder is special 
	radius = outerBeTub->GetOuterRadius()/(cm);
	thick  = innerwallTub->GetInnerRadius()/(cm) - outerBeTub->GetOuterRadius()/(cm);
	length = 2.0*innerwallTub->GetZHalfLength()/(cm);
	z0     = 0.0;
	std::cout<<"outer air: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl;
	KalFitCylinder outerAirCylinder(&_BesKalmanFitMaterials[4], radius, thick, length , z0);
	_BesKalmanFitWalls.push_back(outerAirCylinder);

	/// outer Beryllium layer
	radius = outerBeTub->GetInnerRadius()/(cm);
	thick  = outerBeTub->GetOuterRadius()/(cm) - outerBeTub->GetInnerRadius()/(cm);
	length = 2.0*outerBeTub->GetZHalfLength()/(cm);
	z0     = 0.0;
	std::cout<<"outer Be: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl; 
	KalFitCylinder outerBeCylinder(&_BesKalmanFitMaterials[5], radius, thick, length , z0);
	_BesKalmanFitWalls.push_back(outerBeCylinder);

	/// oil layer
	radius = oilLayerTub->GetInnerRadius()/(cm);
	thick  = oilLayerTub->GetOuterRadius()/(cm) - oilLayerTub->GetInnerRadius()/(cm);
	length = 2.0*oilLayerTub->GetZHalfLength()/(cm);
	z0     = 0.0;
	std::cout<<"oil layer: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl; 
	KalFitCylinder oilLayerCylinder(&_BesKalmanFitMaterials[6], radius, thick, length , z0);
	_BesKalmanFitWalls.push_back(oilLayerCylinder);

	/// inner Beryllium layer
	radius = innerBeTub->GetInnerRadius()/(cm);
	thick  = innerBeTub->GetOuterRadius()/(cm) - innerBeTub->GetInnerRadius()/(cm);
	length = 2.0*innerBeTub->GetZHalfLength()/(cm);
	z0     = 0.0;
	std::cout<<"inner Be: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl; 
	KalFitCylinder innerBeCylinder(&_BesKalmanFitMaterials[7], radius, thick, length , z0);
	_BesKalmanFitWalls.push_back(innerBeCylinder);

	/// gold layer
	radius = goldLayerTub->GetInnerRadius()/(cm);
	thick  = goldLayerTub->GetOuterRadius()/(cm) - goldLayerTub->GetInnerRadius()/(cm);
	length = 2.0*goldLayerTub->GetZHalfLength()/(cm);
	z0     = 0.0;
	std::cout<<"gold layer: "<<" radius: "<<radius<<" thick:"<<thick<<" length: "<<length<<std::endl; 
	KalFitCylinder goldLayerCylinder(&_BesKalmanFitMaterials[8], radius, thick, length , z0);
	_BesKalmanFitWalls.push_back(goldLayerCylinder);

	if(testCDC>0) {
	  /// CDC inner wall by wuchen 2013-06-23
	  radius = 50.16; // cm
	  thick  = 0.04; // cm
	  length = 150; // cm
	  z0     = 0.0; // cm
	  std::cout<<"CDC inner wall: "<<" radius: "<<radius<<" thick:"<<thick<<" outerRadius:"<<radius+thick<<" length: "<<length<<std::endl;
	  KalFitCylinder innerWallCDC(&_BesKalmanFitMaterials[1], radius, thick, length , z0);
	  _BesKalmanFitWalls.clear();
	  _BesKalmanFitWalls.push_back(innerWallCDC);
	}
}

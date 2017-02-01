//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id$
//
// G4TheoFSGenerator
//
// 20110307  M. Kelsey -- Add call to new theTransport->SetPrimaryProjectile()
//		to provide access to full initial state (for Bertini)
// 20110805  M. Kelsey -- Follow change to G4V3DNucleus::GetNucleons()

#include "G4DynamicParticle.hh"
#include "G4TheoFSGenerator.hh"
#include "G4ReactionProductVector.hh"
#include "G4ReactionProduct.hh"
#include "G4IonTable.hh"
#include <iomanip>

G4TheoFSGenerator::G4TheoFSGenerator(const G4String& name)
    : G4HadronicInteraction(name)
    , theTransport(0), theHighEnergyGenerator(0)
    , theQuasielastic(0), fP0(0)
 {
 theParticleChange = new G4HadFinalState;
}

G4TheoFSGenerator::~G4TheoFSGenerator()
{
  delete theParticleChange;
}

void G4TheoFSGenerator::ModelDescription(std::ostream& outFile) const
{
  outFile << GetModelName() <<" consists of a " << theHighEnergyGenerator->GetModelName()
		  << " string model and a stage to de-excite the excited nuclear fragment."
		  << "\n<p>"
		  << "The string model simulates the interaction of\n"
          << "an incident hadron with a nucleus, forming \n"
          << "excited strings, decays these strings into hadrons,\n"
          << "and leaves an excited nucleus. \n"
          << "<p>The string model:\n";
  theHighEnergyGenerator->ModelDescription(outFile);
  outFile <<"\n<p>";
  theTransport->PropagateModelDescription(outFile);
}


G4HadFinalState * G4TheoFSGenerator::ApplyYourself(const G4HadProjectile & thePrimary, G4Nucleus &theNucleus)
{
  // init particle change
  theParticleChange->Clear();
  theParticleChange->SetStatusChange(stopAndKill);
  
	SetP0AntiDeuteron(thePrimary);
  // check if models have been registered, and use default, in case this is not true @@
  
  const G4DynamicParticle aPart(thePrimary.GetDefinition(),thePrimary.Get4Momentum().vect());

  if ( theQuasielastic ) {
  
     if ( theQuasielastic->GetFraction(theNucleus, aPart) > G4UniformRand() )
     {
       G4cout<<"___G4TheoFSGenerator: before Scatter (1) QE=" << theQuasielastic<<G4endl;
       G4KineticTrackVector *result= theQuasielastic->Scatter(theNucleus, aPart);
       G4cout << "^^G4TheoFSGenerator: after Scatter (1) " << G4endl;
       if (result)
       {
	    for(unsigned int  i=0; i<result->size(); i++)
	    {
	      G4DynamicParticle * aNew = 
		 new G4DynamicParticle(result->operator[](i)->GetDefinition(),
                        	       result->operator[](i)->Get4Momentum().e(),
                        	       result->operator[](i)->Get4Momentum().vect());
	      theParticleChange->AddSecondary(aNew);
	      delete result->operator[](i);
	    }
	    delete result;
	   
       } else 
       {
	    theParticleChange->SetStatusChange(isAlive);
	    theParticleChange->SetEnergyChange(thePrimary.GetKineticEnergy());
	    theParticleChange->SetMomentumChange(thePrimary.Get4Momentum().vect().unit());
 
       }
	return theParticleChange;
     } 
  }

 // get result from high energy model

  G4KineticTrackVector * theInitialResult =
               theHighEnergyGenerator->Scatter(theNucleus, aPart);

#define DEBUG_initial_result
  #ifdef DEBUG_initial_result
  	  G4double E_out(0);
  	  G4IonTable * ionTable=G4ParticleTable::GetParticleTable()->GetIonTable();
  	  std::vector<G4KineticTrack *>::iterator ir_iter;
  	  for(ir_iter=theInitialResult->begin(); ir_iter!=theInitialResult->end(); ir_iter++)
  	  {
  		  //G4cout << "TheoFS secondary, mom " << (*ir_iter)->GetDefinition()->GetParticleName() << " " << (*ir_iter)->Get4Momentum() << G4endl;
  		  E_out += (*ir_iter)->Get4Momentum().e();
  	  }
  	  G4double init_mass= ionTable->GetIonMass(theNucleus.GetZ_asInt(),theNucleus.GetA_asInt());
          G4double init_E=aPart.Get4Momentum().e();
  	  // residual nucleus

  	  const std::vector<G4Nucleon> & thy = theHighEnergyGenerator->GetWoundedNucleus()->GetNucleons();

  	  G4int resZ(0),resA(0);
	  G4double delta_m(0);
  	  for(size_t them=0; them<thy.size(); them++)
  	  {
   	     if(thy[them].AreYouHit()) {
  	       ++resA;
  	       if ( thy[them].GetDefinition() == G4Proton::Proton() ) {
	          ++resZ;
		  delta_m +=G4Proton::Proton()->GetPDGMass();
	       } else {
	          delta_m +=G4Neutron::Neutron()->GetPDGMass();
	       }  
  	     }
	  }

  	  G4double final_mass(0);
	  if ( theNucleus.GetA_asInt() ) {
	   final_mass=ionTable->GetIonMass(theNucleus.GetZ_asInt()-resZ,theNucleus.GetA_asInt()- resA);
  	  }
	  G4double E_excit=init_mass + init_E - final_mass - E_out;
	  G4cout << " Corrected delta mass " << init_mass - final_mass - delta_m << G4endl;
  	  G4cout << "initial E, mass = " << init_E << ", " << init_mass << G4endl;
  	  G4cout << "  final E, mass = " << E_out <<", " << final_mass << "  excitation_E " << E_excit << G4endl;
  #endif

  G4ReactionProductVector * theTransportResult = NULL;
  G4ReactionProductVector * theTransportResultnew = NULL;

// Uzhi Nov. 2012
  G4V3DNucleus* theProjectileNucleus = theHighEnergyGenerator->GetProjectileNucleus(); 
if(theProjectileNucleus == 0)                                       // Uzhi Nov. 2012
{                                                                   // Uzhi Nov. 2012
	G4cout<<"theProjectileNucleus == 0"<<G4endl;
  G4int hitCount = 0;
  const std::vector<G4Nucleon>& they = theHighEnergyGenerator->GetWoundedNucleus()->GetNucleons();
  for(size_t them=0; them<they.size(); them++)
  {
    if(they[them].AreYouHit()) hitCount ++;
  }
	G4cout<<"hitCount:"<<hitCount<<" "<<theHighEnergyGenerator->GetWoundedNucleus()->GetMassNumber()<<G4endl;
  if(hitCount != theHighEnergyGenerator->GetWoundedNucleus()->GetMassNumber() )
  {
    theTransport->SetPrimaryProjectile(thePrimary);	// For Bertini Cascade
    theTransportResult = 
               theTransport->Propagate(theInitialResult, theHighEnergyGenerator->GetWoundedNucleus());
    if ( !theTransportResult ) {
       G4cout << "G4TheoFSGenerator: null ptr from transport propagate " << G4endl;
       throw G4HadronicException(__FILE__, __LINE__, "Null ptr from transport propagate");

    }
			G4cout<<"hitCount!=GetMassNumber()"<<G4endl; 
  }
  else
  {
    theTransportResult = theDecay.Propagate(theInitialResult, theHighEnergyGenerator->GetWoundedNucleus());
    if ( !theTransportResult ) {
       G4cout << "G4TheoFSGenerator: null ptr from decay propagate " << G4endl;
       throw G4HadronicException(__FILE__, __LINE__, "Null ptr from decay propagate");
    }

			G4cout<<"hitCount==GetMassNumber()"<<G4endl;    
  }

} else                                                              // Uzhi Nov. 2012
{                                                                   // Uzhi Nov. 2012
    theTransport->SetPrimaryProjectile(thePrimary);
    theTransportResult = 
    theTransport->PropagateNuclNucl(theInitialResult, 
                            theHighEnergyGenerator->GetWoundedNucleus(),
                            theHighEnergyGenerator->GetProjectileNucleus());
    if ( !theTransportResult ) {
       G4cout << "G4TheoFSGenerator: null ptr from transport propagate " << G4endl;
       throw G4HadronicException(__FILE__, __LINE__, "Null ptr from transport propagate");
    } 
}                                                                   // Uzhi Nov. 2012


    //G4DynamicParticle * newdeut = new G4DynamicParticle(deuteron,0.,0.);
	//G4ReactionProduct * finaldeut = new G4ReactionProduct();


/*
	G4String name = theTransportResult->operator[](i)->GetDefinition()->GetParticleName();
	G4cout<<name<<"  "<<theTransportResult->operator[](i)->GetTotalEnergy()<<G4endl;

	if (theTransportResult->operator[](i)->GetDefinition()->GetPDGEncoding() == 2212) {

		G4ParticleDefinition* deuteron = G4ParticleTable::GetParticleTable()->FindParticle("deuteron");		
		G4ReactionProduct * finaldeut = new G4ReactionProduct();//theTransportResult->operator[](i);
		finaldeut->SetDefinition(deuteron);
		finaldeut->SetMomentum(theTransportResult->operator[](i)->GetMomentum());
		finaldeut->SetTotalEnergy(theTransportResult->operator[](i)->GetTotalEnergy());
		finaldeut->SetMass(deuteron->GetPDGMass());
		theTransportResult->push_back(finaldeut);
		theTransportResult->erase(theTransportResult->begin()+i);
	}	
*/

	theTransportResultnew = GenerateDeuterons(theTransportResult);
  // Fill particle change
  unsigned int i;
  for(i=0; i<theTransportResultnew->size(); i++)
  {
	G4String name = theTransportResultnew->operator[](i)->GetDefinition()->GetParticleName();
	G4cout<<name<<G4endl;

    G4DynamicParticle * aNew = 
       new G4DynamicParticle(theTransportResultnew->operator[](i)->GetDefinition(),
                             theTransportResultnew->operator[](i)->GetTotalEnergy(),
                             theTransportResultnew->operator[](i)->GetMomentum());

    // @@@ - overkill? G4double newTime = theParticleChange->GetGlobalTime(theTransportResult->operator[](i)->GetFormationTime());
    theParticleChange->AddSecondary(aNew);
    delete theTransportResultnew->operator[](i);
  }
	//G4cout<<GetP0AntiDeuteron()<<G4endl;  
  // some garbage collection
  delete theTransportResultnew;
  
  // Done
  return theParticleChange;
}

std::pair<G4double, G4double> G4TheoFSGenerator::GetEnergyMomentumCheckLevels() const
{
  if ( theHighEnergyGenerator ) {
	 return theHighEnergyGenerator->GetEnergyMomentumCheckLevels();
  } else {
	 return std::pair<G4double, G4double>(DBL_MAX, DBL_MAX);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4ReactionProductVector * G4TheoFSGenerator::GenerateDeuterons(G4ReactionProductVector * result) const
{
//
// Clusters are made with the first nucleon pair that fulfill
// the coalescence conditions, starting with the protons
//
//
// a deuteron is a pair (i,j) where i is the proton and j the neutron in current event
// with the relative momentum less than p0 (within a sphere of radius p0)
//
//

	//std::vector<std::pair<int, int>> moves
	//std::pair <G4int, G4ThreeVector> foo;
	vector<std::pair<G4int, G4ThreeVector>> * bar= new vector<std::pair<G4int, G4ThreeVector>>;
	vector<G4int> * proton = new vector<G4int>;
	vector<G4ThreeVector> * protonMom = new vector<G4ThreeVector>;
	vector<G4int> * antiproton = new vector<G4int>;
	vector<G4ThreeVector> * antiprotonMom = new vector<G4ThreeVector>;
	vector<G4int> * neutron = new vector<G4int>;
	vector<G4ThreeVector> * neutronMom = new vector<G4ThreeVector>;
	vector<G4int> * antineutron = new vector<G4int>;
	vector<G4ThreeVector> * antineutronMom = new vector<G4ThreeVector>;

	vector<G4int> * protom_del = new vector<G4int>;
	vector<G4int> * neutron_del = new vector<G4int>;
	vector<G4int> * antiprotom_del = new vector<G4int>;
	vector<G4int> * antineutron_del = new vector<G4int>;

   //i;
	//G4LorentzVector v1;
  for(unsigned int i=0; i<result->size(); i++)
  {

	G4int pdgid = result->operator[](i)->GetDefinition()->GetPDGEncoding();

	if (pdgid == 2212){
		//proton->push_back(i);
		protonMom->push_back(result->operator[](i)->GetMomentum());
		//G4double mproton = result->operator[](i)->GetDefinition()->GetPDGMass();
		//G4cout<<"proton encontrado"<<G4endl;
		//vec.SetVect(theTransportResult->operator[](i)->GetMomentum());
		//vec.SetE(theTransportResult->operator[](i)->GetTotalEnergy());
		delete result->operator[](i);
	}

	if (pdgid == -2212){
		//antiproton->push_back(i);
		antiprotonMom->push_back(result->operator[](i)->GetMomentum());
		//G4cout<<"antiproton encontrado"<<G4endl;
		delete result->operator[](i);
	}

	if (pdgid == 2112){
		//neutron->push_back(i);
		neutronMom->push_back(result->operator[](i)->GetMomentum());
		//G4double mneutron = result->operator[](i)->GetDefinition()->GetPDGMass();
		//G4cout<<"neutron encontrado"<<G4endl;
		delete result->operator[](i);
	}

	if (pdgid == -2112){
		//antineutron->push_back(i);
		antineutronMom->push_back(result->operator[](i)->GetMomentum());
		//G4cout<<"antineutron encontrado"<<G4endl;
		delete result->operator[](i);
	}

  }

	for(unsigned int i=0; i<proton->size(); ++i) // center of the sphere
	{
		if(proton->at(i)==-1) continue;  // with next proton
		
		G4ThreeVector p1 = protonMom->at(i);

		int partner1 = this->FindPartner(p1, G4Proton::Proton()->GetPDGMass(), neutron, neutronMom, G4Neutron::Neutron()->GetPDGMass());
		
		if(partner1 == -1) continue; // with next proton
		
		G4ThreeVector p2 = neutronMom->at(partner1);

		this->PushDeuteron(proton->at(i), neutron->at(partner1), p1, p2, result, 1);

		// tag the bound neutron
		neutron->at(partner1) = -1;
		
		//++npart;
	}

	for(unsigned int i=0; i<antiproton->size(); ++i) // center of the sphere
	{
		if(antiproton->at(i)==-1) continue;  // with next proton
		
		G4ThreeVector p1 = antiprotonMom->at(i);

		int partner1 = this->FindPartner(p1, G4Proton::Proton()->GetPDGMass(), antineutron, antineutronMom, G4Neutron::Neutron()->GetPDGMass());
		
		if(partner1 == -1) continue; // with next proton
		
		G4ThreeVector p2 = antineutronMom->at(partner1);

		this->PushDeuteron(antiproton->at(i), antineutron->at(partner1), p1, p2, result, -1);
		
		// tag the bound neutron
		antineutron->at(partner1) = -1;
		
		//++npart;
	}


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4int G4TheoFSGenerator::FindPartner(const G4ThreeVector& p1, G4double m1, vector<G4int> * Neutron, vector<G4ThreeVector> * NeutronMom, G4double m2) const
{
//
// find a nucleon partner within a sphere of radius p0 center at the proton
// and exclude nucleon nx
//

		for(unsigned int j=0; j<Neutron->size(); ++j)
		{
			if(Neutron->at(j) == -1) continue;  // with next nucleon
		
			G4ThreeVector p2 = NeutronMom->at(j);
		
			if(!this->Coalescence(p1,m1,p2,m2)) continue; // with next nucleon
		
			return j;
		}
	
		return -1; // no partner
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool G4TheoFSGenerator::Coalescence( G4double p1x, G4double p1y, G4double p1z, G4double m1
                                 , G4double p2x, G4double p2y, G4double p2z, G4double mass2) const
{
//
// returns true if the nucleons are inside of an sphere of radius p0
// (assume the nucleons are in the same place e.g. PYTHIA, PHOJET,...)
//

	G4double deltaP = this->GetPcm( p1x, p1y, p1z, m1, p2x, p2y, p2z, mass2);
	//G4cout<<fP0<<G4endl;
	return (deltaP < fP0);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4bool G4TheoFSGenerator::Coalescence(const G4ThreeVector& p1, G4double m1, const G4ThreeVector& p2, G4double mass2) const
{
//
// returns true if the nucleons are inside of an sphere of radius p0
// (assume the nucleons are in the same place e.g. PYTHIA, PHOJET,...)
//
	return this->Coalescence(  p1.x(), p1.y(), p1.z(), m1
	                         , p2.x(), p2.y(), p2.z(), mass2);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4TheoFSGenerator::PushDeuteron(G4int i, G4int j, const G4ThreeVector& p1, const G4ThreeVector& p2, G4ReactionProductVector * result, G4int charge) const
{

	if(charge > 0) {
		G4ParticleDefinition* deuteron = G4ParticleTable::GetParticleTable()->FindParticle("deuteron");		
		G4ReactionProduct * finaldeut = new G4ReactionProduct();//theTransportResult->operator[](i);
		finaldeut->SetDefinition(deuteron);

		G4ThreeVector pT = p1+p2;
		G4double massd = deuteron->GetPDGMass();
		G4double E = std::sqrt(pT.mag()*pT.mag()+massd*massd);	
		finaldeut->SetMomentum(pT);
		finaldeut->SetTotalEnergy(E);
		finaldeut->SetMass(massd);
		//result->erase(result->begin()+i);
		//result->erase(result->begin()+j);
		result->push_back(finaldeut);
	} else {
		G4ParticleDefinition* antideuteron = G4ParticleTable::GetParticleTable()->FindAntiParticle("deuteron");	
		G4ReactionProduct * finalantideut = new G4ReactionProduct();//theTransportResult->operator[](i);
		finalantideut->SetDefinition(antideuteron);

		G4ThreeVector pT = p1+p2;
		G4double massd = antideuteron->GetPDGMass();
		G4double E = std::sqrt(pT.mag()*pT.mag()+massd*massd);	
		finalantideut->SetMomentum(pT);
		finalantideut->SetTotalEnergy(E);
		finalantideut->SetMass(massd);
		//result->erase(result->begin()+i);
		//result->erase(result->begin()+j);
		result->push_back(finalantideut);
	}

}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
void G4TheoFSGenerator::SetP0AntiDeuteron(const G4HadProjectile & thePrimary)
{
	//Colaescence condition just available for proton-A, antiproton-A collisions.

	G4cout<<"Projectile: "<<thePrimary.GetDefinition()->GetParticleName()<<G4endl;

	if(std::abs(thePrimary.GetDefinition()->GetPDGEncoding())==2212) {

		//G4int n = target->GetN();
		//G4String targetName = target->GetName();
		G4double mproj = thePrimary.GetDefinition()->GetPDGMass();
		G4double pz = thePrimary.Get4Momentum().z();
		G4double rs = std::sqrt(this->GetS(0.0, 0.0, pz, mproj, 0.0, 0.0, 0.0, mproj));		
		G4cout<<pz<<" "<<rs<<G4endl;
		//G4double A = 2.92*n+188.27;
		//fP0 = A*std::tanh(0.0000226*rs);
		if(thePrimary.GetDefinition()->GetPDGEncoding()==2212){

			if(rs>6.){
				fP0 = 188.5/(1+std::exp(6.09-std::log(0.001*rs)/0.48));
			} else fP0 = 0.0;

		} else {

			if(rs>4.){
            	fP0 = 188.5/(1+std::exp(6.09-std::log(0.001*rs+2.)/0.48));
			} else fP0 = 0.0;
		}

		if(fP0<0.0) fP0 = 0.0;

	} else fP0 = 0.0;

	G4cout<<"Coalescence parameter p0: "<<fP0<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4TheoFSGenerator::GetS(G4double p1x, G4double p1y, G4double p1z, G4double m1, G4double p2x, G4double p2y, G4double p2z, G4double mass2) const
{
//
// square of center of mass energy of two particles from LAB values
//
	G4double E1 = std::sqrt( p1x*p1x + p1y*p1y + p1z*p1z + m1*m1);
	G4double E2 = std::sqrt( p2x*p2x + p2y*p2y + p2z*p2z + mass2*mass2);	

	return (E1+E2)*(E1+E2) - ((p1x+p2x)*(p1x+p2x) + (p1y+p2y)*(p1y+p2y) + (p1z+p2z)*(p1z+p2z));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4TheoFSGenerator::GetPcm(G4double p1x, G4double p1y, G4double p1z, G4double m1, G4double p2x, G4double p2y, G4double p2z, G4double mass2) const
{
//
// momentum in the CM frame of two particles from LAB values
//
	G4double scm = this->GetS(p1x, p1y, p1z, m1, p2x, p2y, p2z, mass2);
	
	return std::sqrt( (scm-(m1-mass2)*(m1-mass2))*(scm-(m1+mass2)*(m1+mass2)) )/(2.*std::sqrt(scm));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
G4double G4TheoFSGenerator::GetPcm(const G4ThreeVector& p1, G4double m1, const G4ThreeVector& p2, G4double mass2) const
{
//
// momentum in the CM frame of two particles from LAB values
//
	return this->GetPcm(p1.x(),p1.y(),p1.z(),m1,p2.x(),p2.y(),p2.z(),mass2);
}


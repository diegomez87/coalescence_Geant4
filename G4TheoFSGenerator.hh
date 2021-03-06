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
//
// $Id$
//

// Class Description
// Final state production model for theoretical models of hadron inelastic
// scattering in geant4;
// To be used in your physics list in case you need this physics.
// In this case you want to register an object of this class with 
// the corresponding process.
// Note: This class is part of an implementation framework. You need to
// register corresponding high energy generators and transport codes to 
// fill it with life; decay of strong resonances is done directly,
// in case there is no residual nucleus. 
// Class Description - End

#ifndef G4TheoFSGenerator_h
#define G4TheoFSGenerator_h 1

#include "G4VIntraNuclearTransportModel.hh"
#include "G4QuasiElasticChannel.hh"
#include "G4HadronicInteraction.hh"
#include "G4VHighEnergyGenerator.hh"
#include "G4DecayStrongResonances.hh"
#include "G4HadFinalState.hh"
#include "G4QuasiElasticChannel.hh"

using namespace std;

class G4TheoFSGenerator : public G4HadronicInteraction 

{
  public:
      G4TheoFSGenerator(const G4String& name = "TheoFSGenerator");
      ~G4TheoFSGenerator();

  private:
      G4TheoFSGenerator(const G4TheoFSGenerator &right);
      const G4TheoFSGenerator & operator=(const G4TheoFSGenerator &right);
      int operator==(const G4TheoFSGenerator &right) const;
      int operator!=(const G4TheoFSGenerator &right) const;

  public:
      G4HadFinalState * ApplyYourself(const G4HadProjectile & thePrimary, G4Nucleus & theNucleus);
      void SetTransport(G4VIntraNuclearTransportModel *const  value);
      void SetHighEnergyGenerator(G4VHighEnergyGenerator *const  value);
      void SetQuasiElasticChannel(G4QuasiElasticChannel *const value);
      virtual std::pair<G4double, G4double> GetEnergyMomentumCheckLevels() const;
      void ModelDescription(std::ostream& outFile) const;

	  void SetP0AntiDeuteron(const G4HadProjectile & thePrimary);
      G4double GetS(G4double p1x, G4double p1y, G4double p1z, G4double m1, G4double p2x, G4double p2y, G4double p2z, G4double mass2) const;
	  G4double GetPcm(G4double p1x, G4double p1y, G4double p1z, G4double m1, G4double p2x, G4double p2y, G4double p2z, G4double mass2) const;
	  G4double GetPcm(const G4ThreeVector& p1, G4double m1, const G4ThreeVector& p2, G4double mass2) const;
	  G4ReactionProductVector * GenerateDeuterons(G4ReactionProductVector * result) const;
	  G4int FindPartner(const G4ThreeVector& p1, G4double m1, vector<G4int> * Neutron, vector<G4ThreeVector> * NeutronMom, G4double m2) const;
	  G4bool Coalescence( G4double p1x, G4double p1y, G4double p1z, G4double m1, G4double p2x, G4double p2y, G4double p2z, G4double mass2) const;
	  G4bool Coalescence(const G4ThreeVector& p1, G4double m1, const G4ThreeVector& p2, G4double mass2) const;
	  void PushDeuteron(G4int i, G4int j, const G4ThreeVector& p1, const G4ThreeVector& p2, G4ReactionProductVector * result, G4int charge) const;

  private:
      const G4VIntraNuclearTransportModel * GetTransport() const;
      const G4VHighEnergyGenerator * GetHighEnergyGenerator() const;
      const G4HadFinalState * GetFinalState() const;

  private: 
      G4VIntraNuclearTransportModel * theTransport;
      G4VHighEnergyGenerator * theHighEnergyGenerator;
      G4DecayStrongResonances theDecay;
      G4HadFinalState * theParticleChange;
      G4QuasiElasticChannel * theQuasielastic;

	  G4double fP0;
};

inline const G4VIntraNuclearTransportModel * G4TheoFSGenerator::GetTransport() const
{
  return theTransport;
}

inline void G4TheoFSGenerator::SetTransport(G4VIntraNuclearTransportModel *const  value)
{
  theTransport = value;
}

inline const G4VHighEnergyGenerator * G4TheoFSGenerator::GetHighEnergyGenerator() const
{
  return theHighEnergyGenerator;
}

inline void G4TheoFSGenerator::SetHighEnergyGenerator(G4VHighEnergyGenerator *const  value)
{
  theHighEnergyGenerator= value;
}

inline void G4TheoFSGenerator::SetQuasiElasticChannel(G4QuasiElasticChannel *const value)
{
  theQuasielastic = value;
}

inline const G4HadFinalState * G4TheoFSGenerator::GetFinalState() const
{
  return theParticleChange;
}

//inline const G4double G4TheoFSGenerator::GetP0AntiDeuteron() const
//{
//  return fP0;
//}

#endif



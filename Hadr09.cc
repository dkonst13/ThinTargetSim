
/// \file Hadr09.cc
/// \brief Main program of the hadronic/Hadr09 example
//
//------------------------------------------------------------------------
// This program shows how to use the class Hadronic Generator.
// The class HadronicGenerator is a kind of "hadronic generator", i.e.
// provides Geant4 final states (i.e. secondary particles) produced by
// hadron-nuclear inelastic collisions.
// Please see the class itself for more information.
//
// The use of the class Hadronic Generator is very simple:
// the constructor needs to be invoked only once - specifying the name
// of the Geant4 "physics case" to consider ("FTFP_BERT" will be
// considered as default is the name is not specified) - and then one
// method needs to be called at each collision, specifying the type of
// collision (hadron, energy, direction, material) to be simulated.
// The class HadronicGenerator is expected to work also in a
// multi-threaded environment with "external" threads (i.e. threads
// that are not necessarily managed by Geant4 run-manager):
// each thread should have its own instance of the class.
//
// See the string "***LOOKHERE***" below for the setting of parameters
// of this example: the "physics case", the set of possibilities from
// which to sample the projectile (i.e. whether the projectile is a
// hadron or an ion - in the case of hadron projectile, a list of hadrons
// is possible from which to sample at each collision; in the case of
// ion projectile, only one type of ion needs to be specified),
// the kinetic energy of the projectile (which can be sampled within
// an interval), whether the direction of the projectile is fixed or
// sampled at each collision, the target material (a list of materials
// is possible, from which the target material can be sampled at each
// collision, and then from this target material, the target nucleus
// will be chosen randomly by Geant4 itself), and whether to print out
// some information or not and how frequently.
// Once a well-defined type of hadron-nucleus or nucleus-nucleus
// inelastic collision has been chosen, the method
//   HadronicGenerator::GenerateInteraction
// returns the secondaries produced by that interaction (in the form
// of a G4VParticleChange object).
// Some information about this final-state is printed out as an example.
//
// Usage:  Hadr09
//------------------------------------------------------------------------

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "CLHEP/Random/Randomize.h"
#include "CLHEP/Random/Ranlux64Engine.h"
#include "HadronicGenerator.hh"

#include "G4GenericIon.hh"
#include "G4HadronicParameters.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4ProcessManager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4VParticleChange.hh"
#include "G4ios.hh"
#include "globals.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int, char**)
{
  G4cout << "=== Test of the HadronicGenerator ===" << G4endl;

  // Enable light hypernuclei and anti-hypernuclei
  G4HadronicParameters::Instance()->SetEnableHyperNuclei(true);

  // See the HadronicGenerator class for the possibilities and meaning of the "physics cases".
  // ( In short, it is the name of the Geant4 hadronic model used for the simulation of
  //   the collision, with the possibility of having a transition between two models in
  //   a given energy interval, as in physics lists. )
  const G4String namePhysics = "FTFP_BERT";  

  // The kinetic energy of the projectile will be sampled randomly, with flat probability
  // in the interval [minEnergy, maxEnergy].
  G4double minEnergy = 1.0 * CLHEP::GeV;  //***LOOKHERE***  HADRON PROJECTILE MIN Ekin
  G4double maxEnergy = 30.0 * CLHEP::GeV;  //***LOOKHERE***  HADRON PROJECTILE MAX Ekin

  const G4int numCollisions = 1000;  //***LOOKHERE***  NUMBER OF COLLISIONS

  // Enable or disable the print out of this program: if enabled, the number of secondaries
  // produced in each collisions is printed out; moreover, once every "printingGap"
  // collisions, the list of secondaries is printed out.
  const G4bool isPrintingEnabled = true;  //***LOOKHERE***  PRINT OUT ON/OFF
  const G4int printingGap = 100;  //***LOOKHERE***  GAP IN PRINTING

  // Vector of Geant4 names of hadron projectiles: one of this will be sampled randomly
  // (with uniform probability) for each collision, when the projectile is not a generic ion
  // (note that the 6 light hypernuclei and anti-hypernuclei are treated here as for the
  // other hadrons, not as generic ions).
  // Note: comment out the corresponding line in order to exclude a particle.
  // Use only proton as projectile for testing
  std::vector<G4String> vecProjectiles;
  vecProjectiles.push_back("proton");

  G4ParticleDefinition* projectileNucleus = nullptr;
  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager(new G4ProcessManager(gion));
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  G4IonTable* ions = partTable->GetIonTable();
  partTable->SetReadiness();
  ions->CreateAllIon();
  ions->CreateAllIsomer();

  const G4bool isProjectileIon = false;  //***LOOKHERE***  HADRON (false) OR ION (true) PROJECTILE?
  if (isProjectileIon) {
    minEnergy = 40.0 * 13.0 * CLHEP::GeV;  //***LOOKHERE***  ION PROJECTILE MIN Ekin
    maxEnergy = 40.0 * 13.0 * CLHEP::GeV;  //***LOOKHERE***  ION PROJECTILE MAX Ekin
    G4int ionZ = 18, ionA = 40;  //***LOOKHERE***  ION PROJECTILE (Z, A)
    projectileNucleus = partTable->GetIonTable()->GetIon(ionZ, ionA, 0.0);
  }

  // Vector of Geant4 NIST names of materials: one of this will be sampled randomly
  // (with uniform probability) for each collision and used as target material.
  // Note: comment out the corresponding line in order to exclude a material;
  //       or, vice versa, add a new line to extend the list with another material.
  G4String selectedMaterial = "G4_Fe"; // default
  std::vector<G4String> vecMaterials;
  vecMaterials.push_back(selectedMaterial);

  const G4int numProjectiles = vecProjectiles.size();
  const G4int numMaterials = vecMaterials.size();

  G4cout << G4endl << "=================  Configuration ==================" << G4endl
         << "Model: " << namePhysics << G4endl << "Ekin: [ " << minEnergy / CLHEP::GeV << " , "
         << maxEnergy / CLHEP::GeV << " ] GeV" << G4endl
         << "Number of collisions:  " << numCollisions << G4endl
         << "Number of hadron projectiles: " << numProjectiles << G4endl
         << "Number of materials:   " << numMaterials << G4endl
         << "IsIonProjectile: " << (projectileNucleus != nullptr ? "true \t" : "false")
         << (projectileNucleus != nullptr ? projectileNucleus->GetParticleName() : G4String(""))
         << G4endl << "===================================================" << G4endl << G4endl;

  CLHEP::Ranlux64Engine defaultEngine(1234567, 4);
  CLHEP::HepRandom::setTheEngine(&defaultEngine);
  G4int seed = time(NULL);
  CLHEP::HepRandom::setTheSeed(seed);
  G4cout << G4endl << " Initial seed = " << seed << G4endl << G4endl;

  // Instanciate the HadronicGenerator providing the name of the "physics case"
  HadronicGenerator* theHadronicGenerator = new HadronicGenerator(namePhysics);
  //****************************************************************************

  if (theHadronicGenerator == nullptr) {
    G4cerr << "ERROR: theHadronicGenerator is NULL !" << G4endl;
    return 1;
  }
  else if (!theHadronicGenerator->IsPhysicsCaseSupported()) {
    G4cerr << "ERROR: this physics case is NOT supported !" << G4endl;
    return 2;
  }

  // Loop over the collisions
  G4double rnd1, rnd2, rnd3, rnd4, rnd5, rnd6, normalization, projectileEnergy;
  G4VParticleChange* aChange = nullptr;
  for (G4int i = 0; i < numCollisions; ++i) {
    // Draw some random numbers to select the hadron-nucleus interaction:
    // projectile hadron, projectile kinetic energy, projectile direction, and target material.
    rnd1 = CLHEP::HepRandom::getTheEngine()->flat();
    rnd2 = CLHEP::HepRandom::getTheEngine()->flat();
    rnd3 = CLHEP::HepRandom::getTheEngine()->flat();
    rnd4 = CLHEP::HepRandom::getTheEngine()->flat();
    rnd5 = CLHEP::HepRandom::getTheEngine()->flat();
    rnd6 = CLHEP::HepRandom::getTheEngine()->flat();
    // Sample the projectile kinetic energy
    projectileEnergy = minEnergy + rnd1 * (maxEnergy - minEnergy);
    if (projectileEnergy <= 0.0) projectileEnergy = minEnergy;
    // Sample the projectile direction
    normalization = 1.0 / std::sqrt(rnd2 * rnd2 + rnd3 * rnd3 + rnd4 * rnd4);
    const G4bool isOnSmearingDirection = true;  //***LOOKHERE***
    G4ThreeVector aDirection = G4ThreeVector(0.0, 0.0, 1.0);  //***LOOKHERE***
    if (isOnSmearingDirection) {
      aDirection = G4ThreeVector(normalization * rnd2, normalization * rnd3, normalization * rnd4);
    }
    // Sample the projectile hadron from the vector vecProjectiles
    G4int index_projectile = std::trunc(rnd5 * numProjectiles);
    G4String nameProjectile = vecProjectiles[index_projectile];
    G4ParticleDefinition* projectile = partTable->FindParticle(nameProjectile);
    if (projectileNucleus) {
      nameProjectile = projectileNucleus->GetParticleName();
      projectile = projectileNucleus;
    }
    // Sample the target material from the vector vecMaterials
    // (Note: the target nucleus will be sampled by Geant4)
    G4int index_material = std::trunc(rnd6 * numMaterials);
    G4String nameMaterial = vecMaterials[index_material];
    G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(nameMaterial);
    if (material == nullptr) {
      G4cerr << "ERROR: Material " << nameMaterial << " is not found !" << G4endl;
      return 3;
    }
    if (isPrintingEnabled) {
      G4cout << "\t Collision " << i << " ; projectile=" << nameProjectile;
      if (projectileNucleus) {
        G4cout << " ; Ekin[MeV]/nucleon="
               << projectileEnergy
                    / static_cast<G4double>(std::abs(projectileNucleus->GetBaryonNumber()));
      }
      else {
        G4cout << " ; Ekin[MeV]=" << projectileEnergy;
      }
      G4cout << " ; direction=" << aDirection << " ; material=" << nameMaterial;
    }

    // Call here the "hadronic generator" to get the secondaries produced by the hadronic collision
    aChange = theHadronicGenerator->GenerateInteraction(
      projectile, projectileEnergy,
      /* ********************************************** */ aDirection, material);

    G4int nsec = aChange ? aChange->GetNumberOfSecondaries() : 0;
    G4bool isPrintingOfSecondariesEnabled = false;
    if (isPrintingEnabled) {
      G4cout << G4endl << "\t --> #secondaries=" << nsec
             << " ; impactParameter[fm]=" << theHadronicGenerator->GetImpactParameter() / fermi
             << " ; #projectileSpectatorNucleons="
             << theHadronicGenerator->GetNumberOfProjectileSpectatorNucleons()
             << " ; #targetSpectatorNucleons="
             << theHadronicGenerator->GetNumberOfTargetSpectatorNucleons()
             << " ; #NNcollisions=" << theHadronicGenerator->GetNumberOfNNcollisions() << G4endl;
      if (i % printingGap == 0) {
        isPrintingOfSecondariesEnabled = true;
        G4cout << "\t \t List of produced secondaries: " << G4endl;
      }
    }
    // Loop over produced secondaries and eventually print out some information.
    for (G4int j = 0; j < nsec; ++j) {
      const G4DynamicParticle* sec = aChange->GetSecondary(j)->GetDynamicParticle();
      if (isPrintingOfSecondariesEnabled) {
        G4cout << "\t \t \t j=" << j << "\t" << sec->GetDefinition()->GetParticleName()
               << "\t p=" << sec->Get4Momentum() << " MeV" << G4endl;
      }
      delete aChange->GetSecondary(j);
    }
    if (aChange) aChange->Clear();
  }

  G4cout << G4endl << " Final random number = " << CLHEP::HepRandom::getTheEngine()->flat()
         << G4endl << "=== End of test ===" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

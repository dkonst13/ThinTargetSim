#include "CLHEP/Units/SystemOfUnits.h"
#include "HadronicGenerator.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
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

int main(int, char**) {
  G4cout << "=== Test of the HadronicGenerator ===" << G4endl;

  // Enable hypernuclei (not used here, but harmless to keep)
  G4HadronicParameters::Instance()->SetEnableHyperNuclei(true);

  const G4String namePhysics = "FTFP";
  const G4int numCollisions = 1000;

  // Fixed configuration
  const G4String nameProjectile = "proton";
  const G4String nameMaterial = "G4_C";
  const G4double projectileEnergy = 5.0 * CLHEP::GeV;
  const G4ThreeVector aDirection(0.0, 0.0, 1.0);

  // Construct standard Geant4 particles
  G4LeptonConstructor().ConstructParticle();
  G4MesonConstructor().ConstructParticle();
  G4BaryonConstructor().ConstructParticle();
  G4IonConstructor().ConstructParticle();
  G4BosonConstructor().ConstructParticle();
  G4ShortLivedConstructor().ConstructParticle();

  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager(new G4ProcessManager(gion));
  partTable->GetIonTable()->CreateAllIon();
  partTable->GetIonTable()->CreateAllIsomer();

  G4ParticleDefinition* projectile = partTable->FindParticle(nameProjectile);
  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(nameMaterial);

  if (!projectile) {
    G4cerr << "ERROR: Could not find projectile particle: " << nameProjectile << G4endl;
    return 1;
  }
  if (!material) {
    G4cerr << "ERROR: Could not find material: " << nameMaterial << G4endl;
    return 2;
  }

  G4cout << "\n=================  Configuration ==================" << G4endl
         << "Model: " << namePhysics << G4endl
         << "Projectile: " << nameProjectile << G4endl
         << "Target material: " << nameMaterial << G4endl
         << "Projectile Ekin: " << projectileEnergy / CLHEP::GeV << " GeV" << G4endl
         << "Direction: " << aDirection << G4endl
         << "Number of collisions: " << numCollisions << G4endl
         << "===================================================" << G4endl << G4endl;

  // Initialize the hadronic generator
  HadronicGenerator* theHadronicGenerator = new HadronicGenerator(namePhysics);
  if (!theHadronicGenerator) {
    G4cerr << "ERROR: theHadronicGenerator is NULL !" << G4endl;
    return 3;
  }
  if (!theHadronicGenerator->IsPhysicsCaseSupported()) {
    G4cerr << "ERROR: this physics case is NOT supported !" << G4endl;
    return 4;
  }

  // Enable or disable all printing
  const G4bool isPrintingEnabled = true;


G4ThreeVector mom(0, 0, 1*GeV);  // z-direction, 31 GeV/c
  // Loop over collisions
  for (G4int i = 0; i < numCollisions; ++i) {
    //G4VParticleChange* aChange = theHadronicGenerator->GenerateInteraction(
    //      projectile, projectileEnergy, aDirection, material);
auto aChange = theHadronicGenerator->GenerateInteraction(projectile, mom, material);

    G4int nsec = aChange ? aChange->GetNumberOfSecondaries() : 0;

    if (isPrintingEnabled) {

       G4cout << "\n================================================================"
              << "\n Collision #" << i
              << "\n ----------------------------------------------------------------"
              << "\n   Projectile:       " << nameProjectile
              << "\n   Kinetic Energy:   " << std::fixed << std::setprecision(2)
              << projectileEnergy / CLHEP::GeV << " GeV"
              << "\n   Direction:        " << aDirection
              << "\n   Target Material:  " << nameMaterial

              << "\n\n   Impact Parameter:           "
              << theHadronicGenerator->GetImpactParameter() / fermi << " fm"
              << "\n   Projectile Spectators:       "
              << theHadronicGenerator->GetNumberOfProjectileSpectatorNucleons()
              << "\n   Target Spectators:           "
              << theHadronicGenerator->GetNumberOfTargetSpectatorNucleons()
              << "\n   Number of NN Collisions:     "
              << theHadronicGenerator->GetNumberOfNNcollisions()
              << "\n   Number of Secondaries:       " << nsec
              << "\n ----------------------------------------------------------------\n";


      G4cout << "\n"
             << "    j    Name         px [MeV]     py [MeV]     pz [MeV]     E [MeV]\n"
             << "-----------------------------------------------------------------------"
             << G4endl;

      for (G4int j = 0; j < nsec; ++j) {
        const G4DynamicParticle* sec = aChange->GetSecondary(j)->GetDynamicParticle();
        const G4String pname = sec->GetDefinition()->GetParticleName();
        const G4LorentzVector p4 = sec->Get4Momentum();

        G4cout << std::setw(5) << j << "  "
               << std::setw(10) << pname << "  "
               << std::fixed << std::setprecision(2)
               << std::setw(11) << p4.px() << "  "
               << std::setw(11) << p4.py() << "  "
               << std::setw(11) << p4.pz() << "  "
               << std::setw(9)  << p4.e()  << G4endl;
      }
    }

    if (aChange) aChange->Clear();
  }

  G4cout << "\n=== End of test ===" << G4endl;
  return 0;
}

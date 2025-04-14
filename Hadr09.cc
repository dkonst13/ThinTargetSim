#include "CLHEP/Units/SystemOfUnits.h"
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

int main(int, char**)
{
  G4cout << "=== Test of the HadronicGenerator ===" << G4endl;

  // Enable hypernuclei (not used here, but harmless to keep)
  G4HadronicParameters::Instance()->SetEnableHyperNuclei(true);

  const G4String namePhysics = "FTFP_BERT";
  const G4int numCollisions = 1000;

  // Fixed configuration
  const G4String nameProjectile = "proton";
  const G4String nameMaterial = "G4_Fe";
  const G4double projectileEnergy = 13.0 * CLHEP::GeV;
  const G4ThreeVector aDirection(0.0, 0.0, 1.0);

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

  // Enable output control
  const G4bool isPrintingEnabled = true;
  const G4int printingGap = 100;

  // Loop over collisions
  for (G4int i = 0; i < numCollisions; ++i) {
    if (isPrintingEnabled) {
      G4cout << "\t Collision " << i
             << " ; projectile=" << nameProjectile
             << " ; Ekin[MeV]=" << projectileEnergy
             << " ; direction=" << aDirection
             << " ; material=" << nameMaterial;
    }

    G4VParticleChange* aChange = theHadronicGenerator->GenerateInteraction(
      projectile, projectileEnergy, aDirection, material);

    G4int nsec = aChange ? aChange->GetNumberOfSecondaries() : 0;
    G4bool isPrintingOfSecondariesEnabled = false;

    if (isPrintingEnabled) {
      G4cout << G4endl
             << "\t --> #secondaries=" << nsec
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

    for (G4int j = 0; j < nsec; ++j) {
      const G4DynamicParticle* sec = aChange->GetSecondary(j)->GetDynamicParticle();
      if (isPrintingOfSecondariesEnabled) {
        G4cout << "\t \t \t j=" << j << "\t"
               << sec->GetDefinition()->GetParticleName()
               << "\t p=" << sec->Get4Momentum() << " MeV" << G4endl;
      }
      delete aChange->GetSecondary(j);
    }

    if (aChange) aChange->Clear();
  }

  G4cout << "\n=== End of test ===" << G4endl;
  return 0;
}

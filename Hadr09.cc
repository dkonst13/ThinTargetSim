// STL
#include <iomanip>

// CLHEP
#include "CLHEP/Units/SystemOfUnits.h"

#include "G4Nucleus.hh"
#include "G4NucleiProperties.hh"

// Geant4 Core
#include "globals.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"

// Geant4 Materials and Particles
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4GenericIon.hh"
#include "G4ProcessManager.hh"

// Geant4 Constructors
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4ShortLivedConstructor.hh"

// Geant4 Hadronics
#include "G4VParticleChange.hh"
#include "G4HadronicParameters.hh"

// Local
#include "HadronicGenerator.hh"


int main(int, char**) {

  // Compute energy and kinetic energy from momentum and mass
  G4ParticleTable* partTable = G4ParticleTable::GetParticleTable();
  partTable->SetReadiness();

  // Construct standard Geant4 particles
  G4LeptonConstructor().ConstructParticle();
  G4MesonConstructor().ConstructParticle();
  G4BaryonConstructor().ConstructParticle();
  G4IonConstructor().ConstructParticle();
  G4BosonConstructor().ConstructParticle();
  G4ShortLivedConstructor().ConstructParticle();

  G4GenericIon* gion = G4GenericIon::GenericIon();
  gion->SetProcessManager(new G4ProcessManager(gion));
  partTable->GetIonTable()->CreateAllIon();
  partTable->GetIonTable()->CreateAllIsomer();



  G4cout << "=== Test of the HadronicGenerator ===" << G4endl;

  // Enable hypernuclei (not used here, but harmless to keep)
  G4HadronicParameters::Instance()->SetEnableHyperNuclei(true);

  const G4String namePhysics = "FTFP";
  const G4int numCollisions = 1000;

  // === Fixed configuration parameters ===
  const G4String nameProjectile = "proton";
  const G4String nameMaterial   = "G4_C";

  // Define projectile momentum in z-direction (magnitude in GeV/c)
  const G4double projectileMomentumZ = 5.0 * CLHEP::GeV;
  const G4ThreeVector projectileMomentum(0.0, 0.0, projectileMomentumZ);
  const G4ThreeVector projectileDirection = projectileMomentum.unit();  // Normalize for direction

  G4ParticleDefinition* projectile = partTable->FindParticle(nameProjectile);
  if (!projectile) {
    G4cerr << "ERROR: Could not find projectile particle: " << nameProjectile << G4endl;
    return 1;
  }

  const G4double mass = projectile->GetPDGMass();
  const G4double totalEnergy = std::sqrt(projectileMomentum.mag2() + mass * mass);
  const G4double projectileEnergy = totalEnergy - mass;  // kinetic energy

  G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(nameMaterial);
  if (!material) {
    G4cerr << "ERROR: Could not find material: " << nameMaterial << G4endl;
    return 2;
  }

  // --- Extract effective (Z, A) from material's primary element ---
  const G4Element* element = material->GetElement(0);
  const G4int Z = static_cast<G4int>(element->GetZ());  // Atomic number
  G4int A = 0;  // Mass number

  // Try to get A from isotope data if available
  if (element->GetNumberOfIsotopes() > 0) {
    A = element->GetIsotope(0)->GetN();
  } else {
    // Estimate A from molar mass (in g/mole)
    A = G4lrint(element->GetA() / (CLHEP::g / CLHEP::mole));
  }

  // Get nucleus mass from Geant4 nuclear database (includes binding energy)
  const G4double nucleusMass = G4NucleiProperties::GetNuclearMass(A, Z);

  G4cout << "\n--- Target Nucleus Info ----------------------------------\n"
         << "  Element name:      " << element->GetName() << "\n"
         << "  Atomic number Z:   " << Z << "\n"
         << "  Mass number A:     " << A << "\n"
         << "  Nucleus mass:      " << G4BestUnit(nucleusMass, "Energy") << "\n";

  // Approximate mass from material A (not subtracting binding energy)
  const G4double Af = material->GetA();
  const G4double massFromA = Af / CLHEP::Avogadro * CLHEP::c_squared;

  G4cout << "  Material mass (A): " << G4BestUnit(massFromA, "Energy")
         << "\n----------------------------------------------------------\n";

  G4cout << "\n=================  Configuration ==================" << G4endl
         << "Model:              " << namePhysics << G4endl
         << "Projectile:         " << nameProjectile << G4endl
         << "Target material:    " << nameMaterial << G4endl
         << "Projectile Ekin:    " << projectileEnergy / CLHEP::GeV << " GeV" << G4endl
         << "Momentum:           " << G4BestUnit(projectileMomentum.mag(), "Energy") << G4endl
         << "Direction:          " << projectileDirection << G4endl
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


  // Loop over collisions
  for (G4int i = 0; i < numCollisions; ++i) {

    auto aChange = theHadronicGenerator->GenerateInteraction(projectile, projectileMomentum, material);

    G4int nsec = aChange ? aChange->GetNumberOfSecondaries() : 0;

    if (isPrintingEnabled) {

       G4cout << "\n================================================================"
              << "\n Collision #" << i
              << "\n ----------------------------------------------------------------"
              << "\n   Projectile:       " << nameProjectile
              << "\n   Kinetic Energy:   " << std::fixed << std::setprecision(2)
              << projectileEnergy / CLHEP::GeV << " GeV"
              << "\n   Direction:        " << projectileDirection
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

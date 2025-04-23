// main.cc (refactored from your original main)
#include "HadronicAnalysis.hh"
#include "Observables.hh"
#include "HadronicAnalysisLoader.hh"
#include "HadronicGenerator.hh"
#include "G4HadronicParameters.hh"

#include <G4ParticleTable.hh>
#include <G4SystemOfUnits.hh>
#include <G4NistManager.hh>
#include <G4Element.hh>
#include <G4Material.hh>
#include <G4NucleiProperties.hh>
#include <G4LeptonConstructor.hh>
#include <G4MesonConstructor.hh>
#include <G4BaryonConstructor.hh>
#include <G4IonConstructor.hh>
#include <G4BosonConstructor.hh>
#include <G4ShortLivedConstructor.hh>
#include <getopt.h>
#include <iostream>

#include "YODA/WriterYODA.h"

#include <dlfcn.h>

int main(int argc, char** argv) {
    std::string analysisName;
    G4int numCollisions = 1000000;

    int opt;
    while ((opt = getopt(argc, argv, "a:n:")) != -1) {
        if (opt == 'a') analysisName = optarg;
        else if (opt == 'n') numCollisions = std::stoi(optarg);
    }

    if (analysisName.empty()) {
        std::cerr << "Usage: " << argv[0] << " -a <AnalysisName> [-n Ncoll]" << std::endl;
        return 1;
    }

    void* handle = nullptr;
    std::string libPath;

    // Try LD_LIBRARY_PATH first
    libPath = "lib" + analysisName + ".so";
    handle = dlopen(libPath.c_str(), RTLD_LAZY);

    if (!handle) {
       // Fallback: try local ./plugins directory
       libPath = "./plugins/lib" + analysisName + ".so";
       handle = dlopen(libPath.c_str(), RTLD_LAZY);
   }


    HadronicAnalysis* analysis = LoadAnalysis(libPath, &handle);
    if (!analysis) return 2;

    analysis->Initialize(numCollisions);

    // Standard Geant4 init
    G4ParticleTable::GetParticleTable()->SetReadiness();
    G4LeptonConstructor().ConstructParticle();
    G4MesonConstructor().ConstructParticle();
    G4BaryonConstructor().ConstructParticle();
    G4IonConstructor().ConstructParticle();
    G4BosonConstructor().ConstructParticle();
    G4ShortLivedConstructor().ConstructParticle();

    G4HadronicParameters::Instance()->SetEnableHyperNuclei(true);

    G4String namePhysics = "QGSP";
    G4String nameProjectile = "proton";
    G4String nameMaterial = "G4_C";

    G4ParticleDefinition* projectile = G4ParticleTable::GetParticleTable()->FindParticle(nameProjectile);
    if (!projectile) return 1;

    G4double projectileMomentumZ = 31.0 * CLHEP::GeV;
    G4ThreeVector projectileMomentum(0., 0., projectileMomentumZ);
    G4double mass = projectile->GetPDGMass();
    G4double totalEnergy = std::sqrt(projectileMomentum.mag2() + mass * mass);

    G4Material* material = G4NistManager::Instance()->FindOrBuildMaterial(nameMaterial);
    const G4Element* element = material->GetElement(0);
    G4int Z = static_cast<G4int>(element->GetZ());
    G4int A = (element->GetNumberOfIsotopes() > 0) ? element->GetIsotope(0)->GetN() : G4lrint(element->GetA() / (CLHEP::g / CLHEP::mole));
    G4double nucleusMass = G4NucleiProperties::GetNuclearMass(A, Z);

    G4LorentzVector proj4mom_lab(projectileMomentum, totalEnergy);
    G4LorentzVector targ4mom_lab(G4ThreeVector(0., 0., 0.), nucleusMass);
    G4LorentzVector labv = proj4mom_lab + targ4mom_lab;
    G4ThreeVector cmsBoost = labv.boostVector();

    HadronicGenerator* theHadronicGenerator = new HadronicGenerator(namePhysics);
    if (!theHadronicGenerator->IsPhysicsCaseSupported()) return 3;

    for (G4int i = 0; i < numCollisions; ++i) {
        auto aChange = theHadronicGenerator->GenerateInteraction(projectile, projectileMomentum, material);
        G4int nsec = aChange ? aChange->GetNumberOfSecondaries() : 0;
        for (G4int j = 0; j < nsec; ++j) {
            const auto* sec = aChange->GetSecondary(j)->GetDynamicParticle();
            const auto* pd = sec->GetDefinition();
            const G4LorentzVector p4 = sec->Get4Momentum();
            Observables obs = computeObservables(p4, pd, cmsBoost, labv.mag());
            analysis->Fill(obs, pd);
        }
        if (aChange) aChange->Clear();
    }

    analysis->Finalize();
    UnloadAnalysis(analysis, handle);
    return 0;
}

#ifdef WITH_YODA
#include "YODA/WriterYODA.h"
// Force the linker to keep libYODA.so by actually using a symbol
static auto __force_yoda_link = YODA::WriterYODA::create();
#endif


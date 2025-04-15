#ifndef OBSERVABLES_HH
#define OBSERVABLES_HH

#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"
#include "G4ParticleDefinition.hh"

struct Observables {
    G4LorentzVector p4_lab;
    G4ThreeVector   p_lab;
    G4double        E_lab;
    G4double        T_lab;
    G4double        pt_lab;
    G4double        theta_lab;

    G4LorentzVector p4_cms;
    G4ThreeVector   p_cms;
    G4double        E_cms;
    G4double        T_cms;
    G4double        pt_cms;
    G4double        theta_cms;
    G4double        y_cms;
    G4double        xF;
};

Observables computeObservables(const G4LorentzVector& p4_lab,
                               const G4ParticleDefinition* pd,
                               const G4ThreeVector& boostToCMS,
                               G4double sqrtS);

#endif

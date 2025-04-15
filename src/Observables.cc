#include "Observables.hh"

Observables computeObservables(const G4LorentzVector& p4_lab,
                               const G4ParticleDefinition* pd,
                               const G4ThreeVector& boostToCMS,
                               G4double sqrtS)
{
    Observables obs;
    obs.p4_lab     = p4_lab;
    obs.p_lab      = p4_lab.vect();
    obs.E_lab      = p4_lab.e();
    obs.T_lab      = obs.E_lab - pd->GetPDGMass();
    obs.pt_lab     = obs.p_lab.perp();
    obs.theta_lab  = obs.p_lab.theta();

    obs.p4_cms     = p4_lab;
    obs.p4_cms.boost(-boostToCMS);
    obs.p_cms      = obs.p4_cms.vect();
    obs.E_cms      = obs.p4_cms.e();
    obs.T_cms      = obs.E_cms - pd->GetPDGMass();
    obs.pt_cms     = obs.p_cms.perp();
    obs.theta_cms  = obs.p_cms.theta();
    obs.y_cms      = obs.p4_cms.rapidity();
    obs.xF         = 2. * obs.p4_cms.z() / sqrtS;

    return obs;
}

#ifndef HADRONIC_ANALYSIS_HH
#define HADRONIC_ANALYSIS_HH

#include "G4VParticleChange.hh"
#include "G4ParticleDefinition.hh"
#include "G4LorentzVector.hh"

#include "Observables.hh"

#include <string>

// Abstract base class for analyses (like Rivet::Analysis)
class HadronicAnalysis {
public:
    virtual ~HadronicAnalysis() {}

    /// Called once at the beginning of the run
    virtual void Initialize(G4int numCollisions) = 0;

    /// Called for every secondary particle
    virtual void Fill(const Observables& obs, const G4ParticleDefinition* pd) = 0;

    /// Called once at the end of the run
    virtual void Finalize() = 0;

    /// Return name of the analysis
    virtual std::string GetName() const = 0;
};

// Factory function signature used by plugins
extern "C" HadronicAnalysis* CreateAnalysis();

#endif

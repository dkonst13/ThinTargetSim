#include "HadronicAnalysis.hh"
#include <vector>
#include <iostream>
#include <iomanip>
#include <memory>
#include "YODA/Histo.h"
#include "YODA/WriterYODA.h"
#include "YODA/ReaderYODA.h"
#include "HistogramUtils.hh"

class NA61_2009_I151002703 : public HadronicAnalysis {
public:
    void Initialize(G4int numCollisions) override {
        std::cout << GetName() << std::endl;
        std::vector<YODA::AnalysisObject*> aovec;
        YODA::ReaderYODA::create().read(GetName() + ".yoda", aovec);

        for (auto* ao : aovec) {
            const auto* est = dynamic_cast<const YODA::Estimate1D*>(ao);
            if (!est) {
                std::cerr << "Warning: Skipping non-Estimate1D object." << std::endl;
                continue;
            }
            if (!ao->hasAnnotation("Theta range") || !ao->hasAnnotation("Secondary")) continue;

            auto* hist = new YODA::Histo1D(est->xEdges());
            std::cout << "Looping over annotations for histogram..." << std::endl;
            for (const auto& key : ao->annotations()) {
                hist->addAnnotation(key, ao->annotation(key));
                std::cout << key << ": " << ao->annotation(key) << std::endl;
            }

            std::cout << "Histo added" << std::endl;
            _histos.push_back(hist);
        }

        _nCollisions = numCollisions;
    }

    void Fill(const Observables& obs, const G4ParticleDefinition* pd) override {
        if (pd->GetPDGEncoding() != -211 && pd->GetPDGEncoding() != 211) return;

        const G4double theta_mrad = obs.theta_lab * 1000.0;
        const G4double p = obs.p_lab.mag() / CLHEP::GeV;

        for (auto* hist : _histos) {
            auto [theta_min, theta_max] = parseThetaRange(hist->annotation("Theta range"));
            int pdg = pdgFromName(hist->annotation("Secondary"));
            if (pd->GetPDGEncoding() != pdg) continue;
            if (theta_mrad >= theta_min && theta_mrad < theta_max) {
                hist->fill(p, 1.0);
                break;
            }
        }
    }

    void Finalize() override {
        const G4double dtheta_rad = 60.0 * CLHEP::milliradian;
        std::vector<YODA::AnalysisObject*> out;
        for (auto* hist : _histos) {
            out.push_back(hist);
        }
        YODA::WriterYODA::create().write(GetName() + ".yoda", out);
        std::cout << "Saved YODA histogram" << std::endl;
    }

    std::string GetName() const override {
        return "NA61_2009_I151002703";
    }

private:
    G4int _nCollisions = 0;
    std::vector<YODA::Histo1D*> _histos;
};

extern "C" HadronicAnalysis* CreateAnalysis() {
    return new NA61_2009_I151002703();
}

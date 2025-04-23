#include "HadronicAnalysis.hh"
#include <vector>
#include <iostream>
#include <iomanip>
#include <memory>
#include "YODA/Histo.h"
#include "YODA/WriterYODA.h"

class NA61_2024 : public HadronicAnalysis {
public:
    void Initialize(G4int numCollisions) override {
        _nCollisions = numCollisions;

        // Extract edges from bin pairs
        std::vector<double> edges;
        for (const auto& b : _binEdges) edges.push_back(b.first);
        edges.push_back(_binEdges.back().second); // add last upper edge

        _histo = std::make_shared<YODA::Histo1D>(edges);
        _histo->setPath("/NA61_2024/piMinus_180_240_mrad");
        _histo->setTitle("π⁻ d²n/dpdθ in 180–240 mrad");
    }

    void Fill(const Observables& obs, const G4ParticleDefinition* pd) override {
        if (pd->GetPDGEncoding() != -211) return;  // pi⁻ only

        const G4double theta_mrad = obs.theta_lab * 1000.0;
        if (theta_mrad < 180.0 || theta_mrad >= 240.0) return;

        const G4double p = obs.p_lab.mag() / CLHEP::GeV;
        _histo->fill(p, 1.0);
    }

    void Finalize() override {
        const G4double dtheta_rad = 60.0 * CLHEP::milliradian;

/*
        // Normalize the histogram manually
        for (auto& b : _histo->bins()) {
            double dp = b.xMax() - b.xMin();
            double scale = 1.0 / (_nCollisions * dp * dtheta_rad);
            b.scaleW(scale);
        }
*/
        // Save YODA file
        std::vector<YODA::AnalysisObject*> out = { _histo.get() };
        YODA::WriterYODA::create().write("NA61_2024.yoda", out);

        std::cout << "Saved YODA histogram: NA61_2024.yoda" << std::endl;
    }

    std::string GetName() const override {
        return "NA61_2024";
    }

private:
    G4int _nCollisions = 0;
    std::shared_ptr<YODA::Histo1D> _histo;

    const std::vector<std::pair<G4double, G4double>> _binEdges = {
        {0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8},
        {0.8, 0.9}, {0.9, 1.0}, {1.0, 1.2}, {1.2, 1.4}, {1.4, 1.6}, {1.6, 1.8},
        {1.8, 2.0}, {2.0, 2.2}, {2.2, 2.4}, {2.4, 2.6}, {2.6, 2.8}, {2.8, 3.0},
        {3.0, 3.2}, {3.2, 3.4}, {3.4, 3.6}, {3.6, 3.8}, {3.8, 4.0}, {4.0, 4.2},
        {4.2, 4.4}, {4.4, 4.6}, {4.6, 4.8}, {4.8, 5.0}
    };
};

extern "C" HadronicAnalysis* CreateAnalysis() {
    return new NA61_2024();
}


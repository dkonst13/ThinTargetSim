#include "HadronicAnalysis.hh"
#include <vector>
#include <iostream>
#include <iomanip>

class NA61_2024 : public HadronicAnalysis {
public:
    void Initialize(G4int numCollisions) override {
        _nCollisions = numCollisions;
        _bins.resize(28, 0);
    }

    void Fill(const Observables& obs, const G4ParticleDefinition* pd) override {
        if (pd->GetPDGEncoding() != -211) return;  // pi-

        const G4double theta_mrad = obs.theta_lab * 1000.0;
        if (theta_mrad < 180.0 || theta_mrad >= 240.0) return;

        const G4double p = obs.p_lab.mag() / CLHEP::GeV;

        for (size_t i = 0; i < _binEdges.size(); ++i) {
            if (p >= _binEdges[i].first && p < _binEdges[i].second) {
                _bins[i]++;
                break;
            }
        }
    }

    void Finalize() override {
        const G4double dtheta = 60.0 * CLHEP::milliradian;  // 60 mrad

        std::cout << "\n=== NA61_2024 pi- Histogram in 180-240 mrad ===" << std::endl;
        for (size_t i = 0; i < _binEdges.size(); ++i) {
            G4double dp = _binEdges[i].second - _binEdges[i].first;
            double norm = _bins[i] / (_nCollisions * dp * dtheta);
            std::cout << std::fixed << std::setprecision(2)
                      << _binEdges[i].first << " - " << _binEdges[i].second << " GeV/c : "
                      << "count = " << std::setw(5) << _bins[i]
                      << ", d²n/dpdθ = " << std::setprecision(4) << norm
                      << " (rad·GeV/c)^-1" << std::endl;
        }
    }

    std::string GetName() const override {
        return "NA61_2024";
    }

private:
    G4int _nCollisions = 0;
    std::vector<int> _bins;

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

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
        
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_0_10mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_0_10mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_10_20mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_10_20mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_20_40mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_20_40mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_40_60mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_40_60mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_60_100mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_60_100mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_100_140mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_100_140mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_140_180mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_140_180mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_180_240mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_180_240mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_240_300mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_240_300mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_300_360mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_300_360mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_360_420mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_360_420mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_0_10mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_0_10mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_10_20mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_10_20mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_20_40mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_20_40mrad);
        histoInit(aovec, "h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_40_60mrad", h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_40_60mrad);

std::string thetaRange = h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_0_10mrad->annotation("Theta range");
std::string secondary = h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_0_10mrad->annotation("Secondary");
std::cout << "DIMA Annotation Name: " << thetaRange << "  " <<secondary  <<std::endl;
auto [theta_min, theta_max] = parseThetaRange(thetaRange);
std::cout << theta_min << "  " << theta_max << std::endl;        
        for (const auto* ao : aovec) {
            if (ao->path().find("piPlus") == std::string::npos) continue;
            
            
            for (const auto& key : ao->annotations()) {
                std::cout << "Annotation: " << key << " = " << ao->annotation(key) << std::endl;
            }
            
            
            const auto* est1d = dynamic_cast<const YODA::Estimate1D*>(ao);
            if (!est1d) continue;
            
            std::cout << "Reference histogram: " << ao->path() << std::endl;
            for (const auto& bin : est1d->bins()) {
                auto val = bin.val();
                auto stat = bin.err("stat");
                auto syst = bin.err("syst");
                
                std::cout << "  Bin " << bin.xMin() << "-" << bin.xMax()
                << " : val = " << val
                << ", stat = " << stat.first << "±" << stat.second
                << ", syst = " << syst.first << "±" << syst.second
                << std::endl;
            }
        }
        
        
        
        
        _nCollisions = numCollisions;
        
        // Extract edges from bin pairs
        std::vector<double> edges;
        for (const auto& b : _binEdges) edges.push_back(b.first);
        edges.push_back(_binEdges.back().second); // add last upper edge
        
        _histo = std::make_shared<YODA::Histo1D>(edges);
        _histo->setPath("/NA61_2009_I151002703/piMinus_180_240_mrad");
        _histo->setTitle("π⁻ d²n/dpdθ in 180–240 mrad");
        _histo->addAnnotation("Name", "Name");
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
        std::vector<YODA::AnalysisObject*> out = {
            h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_0_10mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_10_20mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_20_40mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_40_60mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_60_100mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_100_140mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_140_180mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_180_240mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_240_300mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_300_360mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_360_420mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_0_10mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_10_20mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_20_40mrad.get(),
            h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_40_60mrad.get()
        };
        
        YODA::WriterYODA::create().write("NA61_2024.yoda", out);
        
        std::cout << "Saved YODA histogram: NA61_2024.yoda" << std::endl;
    }
    
    std::string GetName() const override {
        return "NA61_2009_I151002703";
    }
    
private:
    G4int _nCollisions = 0;
    std::shared_ptr<YODA::Histo1D> _histo;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_0_10mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_10_20mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_20_40mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_40_60mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_60_100mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_100_140mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_140_180mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_180_240mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_240_300mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_300_360mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piPlus_theta_360_420mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_0_10mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_10_20mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_20_40mrad;
    std::shared_ptr<YODA::Histo1D> h_d2sigma_dpdtheta_beam30GeV_targetC_piMinus_theta_40_60mrad;
    const std::vector<std::pair<G4double, G4double>> _binEdges = {
        {0.2, 0.3}, {0.3, 0.4}, {0.4, 0.5}, {0.5, 0.6}, {0.6, 0.7}, {0.7, 0.8},
        {0.8, 0.9}, {0.9, 1.0}, {1.0, 1.2}, {1.2, 1.4}, {1.4, 1.6}, {1.6, 1.8},
        {1.8, 2.0}, {2.0, 2.2}, {2.2, 2.4}, {2.4, 2.6}, {2.6, 2.8}, {2.8, 3.0},
        {3.0, 3.2}, {3.2, 3.4}, {3.4, 3.6}, {3.6, 3.8}, {3.8, 4.0}, {4.0, 4.2},
        {4.2, 4.4}, {4.4, 4.6}, {4.6, 4.8}, {4.8, 5.0}
    };
};

extern "C" HadronicAnalysis* CreateAnalysis() {
    return new NA61_2009_I151002703();
}



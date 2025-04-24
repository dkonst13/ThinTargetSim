// HistogramUtils.hh
#pragma once
#include <YODA/Estimate.h>
#include <YODA/Histo.h>
#include <memory>
#include <vector>
#include <string>
#include <regex>
#include <stdexcept>
#include <utility>

inline void histoInit(const std::vector<YODA::AnalysisObject*>& aovec,
                      const std::string& name,
                      std::shared_ptr<YODA::Histo1D>& histo) {
    for (const auto* ao : aovec) {
        if (!ao->hasAnnotation("Name")) continue;
        if (ao->annotation("Name") != name) continue;
        
        const auto* est = dynamic_cast<const YODA::Estimate1D*>(ao);
        if (!est) continue;
        
        // Create histogram with bin edges from Estimate1D
        histo = std::make_shared<YODA::Histo1D>(est->xEdges());
        
        // Copy path, title, and all annotations
        histo->setPath("/TMP/" + name);
        histo->setTitle(name);
        
        for (const auto& key : ao->annotations()) {
            histo->addAnnotation(key, ao->annotation(key));
        }
        
        return;
    }
    
    throw std::runtime_error("No matching AO found for histogram: " + name);
}

inline std::pair<double, double> parseThetaRange(const std::string& s) {
    std::regex re(R"((\d+(?:\.\d+)?)[\s\-]+(\d+(?:\.\d+)?)\s*mrad)");
    std::smatch match;
    if (std::regex_search(s, match, re)) {
        double theta_min = std::stod(match[1]);
        double theta_max = std::stod(match[2]);
        return {theta_min, theta_max};
    } else {
        throw std::invalid_argument("Could not parse theta range from string: " + s);
    }
}

inline int pdgFromName(const std::string& name) {
    if (name == "pi+") return 211;
    if (name == "pi-") return -211;
    if (name == "proton") return 2212;
    if (name == "kaon_") return 321;
    if (name == "kaon-") return -321;
    throw std::runtime_error("Unknown secondary particle name: " + name);
}

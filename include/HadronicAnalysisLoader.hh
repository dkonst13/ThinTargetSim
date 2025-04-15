#ifndef HADRONIC_ANALYSIS_LOADER_HH
#define HADRONIC_ANALYSIS_LOADER_HH

#include <string>
#include "HadronicAnalysis.hh"

HadronicAnalysis* LoadAnalysis(const std::string& libPath, void** handleOut);
void UnloadAnalysis(HadronicAnalysis* analysis, void* handle);

#endif

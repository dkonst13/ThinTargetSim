#include "HadronicAnalysisLoader.hh"
#include <dlfcn.h>
#include <iostream>

typedef HadronicAnalysis* (*CreateFunc)();

HadronicAnalysis* LoadAnalysis(const std::string& libPath, void** handleOut) {
    void* handle = dlopen(libPath.c_str(), RTLD_LAZY);
    if (!handle) {
        std::cerr << "ERROR: Failed to load " << libPath << ": " << dlerror() << std::endl;
        return nullptr;
    }

    CreateFunc create = (CreateFunc) dlsym(handle, "CreateAnalysis");
    if (!create) {
        std::cerr << "ERROR: Cannot find CreateAnalysis in " << libPath << std::endl;
        dlclose(handle);
        return nullptr;
    }

    *handleOut = handle;
    return create();
}

void UnloadAnalysis(HadronicAnalysis* analysis, void* handle) {
    delete analysis;
    if (handle) dlclose(handle);
}

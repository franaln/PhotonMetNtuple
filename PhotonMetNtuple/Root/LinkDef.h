#ifndef PhotonMetNtuple_LINKDEF_H
#define PhotonMetNtuple_LINKDEF_H

#include <vector>
#include <string>
#include <map>
#include <utility>

#include <PhotonMetNtuple/xAODAnalysis.h>
#include <PhotonMetNtuple/xAODTruthAnalysis.h>
#include <PhotonMetNtuple/xAODCountEwkProcesses.h>
#include <PhotonMetNtuple/MiniTree2.h>
#include <PhotonMetNtuple/TruthTree.h>

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

#endif

#ifdef __CINT__
#pragma link C++ class xAODAnalysis+;
#pragma link C++ class xAODTruthAnalysis+;
#pragma link C++ class xAODCountEwkProcesses+;
#pragma link C++ class MiniTree2+;
#pragma link C++ class TruthTree+;
#pragma link C++ class pair<std::string,TTree >+;
#endif

#endif

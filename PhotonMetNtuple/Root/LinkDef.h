#ifndef PhotonMetNtuple_LINKDEF_H
#define PhotonMetNtuple_LINKDEF_H

#include <vector>
#include <string>
#include <map>
#include <utility>

#include <PhotonMetNtuple/xAODAnalysis.h>
#include <PhotonMetNtuple/OutTree.h>
#include <PhotonMetNtuple/Utils.h>

#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclass;

#endif

#ifdef __CINT__
#pragma link C++ class xAODAnalysis+;
#pragma link C++ class OutTree+;
#pragma link C++ class Utils+;
#pragma link C++ class pair<string,TTree >+;
#endif

#endif

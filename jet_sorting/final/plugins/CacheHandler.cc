/*
Rewritten by chatGPT so that Rory could use it with CMSSW 13.3.1
*/

#include "CacheHandler.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"

CacheHandler::CacheHandler(const edm::FileInPath GraphPath) {
  // Load the graph definition from file
  std::string FullPath = GraphPath.fullPath();
  graph_.reset(tensorflow::loadGraphDef(FullPath));

  // Create session with default options
  tensorflow::Options options;  // this is the correct struct in CMSSW_13
  session_ = tensorflow::createSession(graph_.get(), options);
}

CacheHandler::~CacheHandler() {
  tensorflow::closeSession(session_);
}


#ifndef CACHE_HANDLER_H
#define CACHE_HANDLER_H

#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "tensorflow/core/util/memmapped_file_system.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"

class CacheHandler{
 public:
  //  CacheHandler(std::string GraphName, std::string GraphPath);
  CacheHandler(edm::FileInPath GraphPath);

  ~CacheHandler();
  
  tensorflow::Session& getSession() const { return *session_; }
  const tensorflow::GraphDef& getGraph() const { return *graph_; }

 private:
  std::shared_ptr<tensorflow::GraphDef> graph_;
  tensorflow::Session* session_;
  std::unique_ptr<tensorflow::MemmappedEnv> memmappedEnv_;

};
#endif


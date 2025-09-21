#include "BESTEvaluation.h"
#include <fstream>
#include <sstream>

//========================================================================================
// Loads the BEST Neural Network using the tensforflow interface -------------------------
//----------------------------------------------------------------------------------------

void BESTEvaluation::configure(const edm::ParameterSet& iConfig){
  if(isConfigured_) return;
  

  BESTname_   = iConfig.getParameter<std::string>("BESTname");  // Used for debug messages
  BESTpath_   = iConfig.getParameter<edm::FileInPath>("BESTpath"); // Path to BEST network
  BESTscale_  = iConfig.getParameter<edm::FileInPath>("BESTscale"); // Path to BEST Scalar Parameters file
  runYear_    = iConfig.getParameter<std::string>("year");

  std::string fullParamPath = BESTscale_.fullPath();
  std::ifstream inputParamsFile(fullParamPath);
  std::string line, word, besVar, scaler;
  float param1, param2;

  while (std::getline(inputParamsFile, line)) {
    std::istringstream ss(line);
    std::getline(ss,word,',');
    besVar = word;
    std::getline(ss,word,',');
    scaler = word;
    std::getline(ss,word,',');
    param1 = std::stof(word);
    std::getline(ss,word,',');
    param2 = std::stof(word);
    listOfBESTVars.push_back(besVar); // ordered list of bes vars
    paramDict_[besVar] = std::make_tuple(scaler, param1, param2); // dictionary containing the scaler and parameters for each bes var
  }
  inputParamsFile.close();

  // The input name needs to match the name of the input layer in the data/constantgraph.pb file
  
  // old method for this
  /*if      (runYear_ == "2015")  inputNames_.push_back("input_1");
  else if (runYear_ == "2016")  inputNames_.push_back("input_2");
  else if (runYear_ == "2017")  inputNames_.push_back("input_3");
  else if (runYear_ == "2018")  inputNames_.push_back("input_4");

  // The output name needs to match the name of the output layer in the data/constantgraph.pb file
  if      (runYear_ == "2015")  outputName_ = "dense_9/Softmax";
  else if (runYear_ == "2016")  outputName_ = "dense_18/Softmax";
  else if (runYear_ == "2017")  outputName_ = "dense_27/Softmax";
  else if (runYear_ == "2018")  outputName_ = "dense_36/Softmax"; */


  // single NN for every year
  inputNames_.push_back("input_1");
  outputName_ = "dense_9/Sigmoid";



  inputShapes_.push_back(tensorflow::TensorShape{1, NumBESTInputs_});
  inputTensors_.resize(inputShapes_.size()); // Set the internal tensor list to the size of your inputs

  // Now make each element have the correct name and shape
  for (size_t i=0; i<inputShapes_.size(); i++) 
  {
    inputTensors_[i] = tensorflow::NamedTensor(inputNames_[i], tensorflow::Tensor(tensorflow::DT_FLOAT, inputShapes_.at(i)));
  }

  // This class must be constructed after the global cache was created and passed to it
  auto graph = cache_->getGraph(); // Get the graph from the cache

  // The rest of this code is a sanity check
  for (size_t i=0; i<inputShapes_.size(); i++)
  {
    // abbott: this looks like it grabs the input and output names from the BEST model .pb file
    //    update later to just use these so we don't need to manually update every time the model changes
    const auto& name = graph.node(i).name();
    auto it = std::find(inputNames_.begin(), inputNames_.end(), name);

    if (it == inputNames_.end()) 
    { // Check if input layer name is in the graph
      throw cms::Exception("BESTEvaluation") 
      << "Processing graph: " << BESTname_ << ".\nUnknown input name: " << name;  
    }

    const auto& shape = graph.node(i).attr().at("shape").shape();
    int j = std::distance(inputNames_.begin(),it); // Do inputs in correct order, even if they weren't declared in correct order
    
    for (int d=1; d<inputShapes_.at(j).dims(); d++) 
    {
      if (shape.dim(d).size() != inputShapes_.at(j).dim_size(d)) 
      {
        throw cms::Exception("BESTEvaluation")
        << "Number of inputs in graph does not match those expected for " << BESTname_ << ".\n"
        << "Expected input " << j << " dim " << d << " = " << inputShapes_.at(j).dim_size(d) << ".\n"
        << "Found " << shape.dim(d).size() << ".";
      }
    }
  }

  const auto& outName = graph.node(graph.node_size() - 1).name();
  if (outName != outputName_) {
    throw cms::Exception("BESTEvaluation")
    << "Processing graph " << BESTname_ << ".\n"
    << "Unexpected output name. Expected " << outputName_ << " found " << outName << ".";
  }

  isConfigured_ = true;
}




//========================================================================================
// Get a prediction using the BEST Neural Network ----------------------------------------
//----------------------------------------------------------------------------------------
// std::vector<float> BESTEvaluation::getPrediction(const std::vector<float> &BESTInputs){
std::vector<float> BESTEvaluation::getPrediction(std::map<std::string, float> &BESTInputs){
  std::vector<tensorflow::Tensor> pred_vector; // vector of predictions to allow for evaluation multiple jets, but only using one at the moment
  tensorflow::Tensor prediction;
  std::vector<float> BESTScores;

  inputTensors_.at(kBEST_).second.flat<float>().setZero();


  //std::cout << "------------------------ New Superjet --------------------------" << std::endl;
  //std::cout << "Event number " << std::fixed << int(BESTInputs["eventNumber"]) << std::endl;

  // Scale the BESTInputs and load them into the input tensor in the correct order
  for (int n=0; n < NumBESTInputs_; n++){
    float value = BESTInputs[listOfBESTVars[n]];
    auto [scalerString, param1, param2] = paramDict_[listOfBESTVars[n]];

    //std::cout  <<  listOfBESTVars[n] << " " << getScaledValue(value, scalerStringToInt_[scalerString], param1, param2) << std::endl; 
    inputTensors_.at(kBEST_).second.matrix<float>()(0, n) = getScaledValue(value, scalerStringToInt_[scalerString], param1, param2); 
  }

  // Evaluate network
  tensorflow::run(&(cache_->getSession()), inputTensors_, {outputName_}, &pred_vector);

   // Assuming the model outputs a single sigmoid value
  float pred = pred_vector[0].flat<float>()(0);
  if (!(pred >= 0 && pred <= 1)) {
    throw cms::Exception("BESTEvaluation") << "invalid prediction = " << pred;
  }


  // abbott: Why fill prediction.matrix instead of BESTScores directly? 

  BESTScores.push_back(pred);

  return BESTScores;
}






/* OLD VERSION
// std::vector<float> BESTEvaluation::getPrediction(const std::vector<float> &BESTInputs){
std::vector<float> BESTEvaluation::getPrediction(std::map<std::string, float> &BESTInputs){
  std::vector<tensorflow::Tensor> pred_vector; // vector of predictions to allow for evaluation multiple jets, but only using one at the moment
  tensorflow::Tensor prediction;
  std::vector<float> BESTScores;

  inputTensors_.at(kBEST_).second.flat<float>().setZero();


  //std::cout << "------------------------ New Superjet --------------------------" << std::endl;
  //std::cout << "Event number " << std::fixed << int(BESTInputs["eventNumber"]) << std::endl;

  // Scale the BESTInputs and load them into the input tensor in the correct order
  for (int n=0; n < NumBESTInputs_; n++)
  {
    float value = BESTInputs[listOfBESTVars[n]];
    auto [scalerString, param1, param2] = paramDict_[listOfBESTVars[n]];

    //std::cout  <<  listOfBESTVars[n] << " " << getScaledValue(value, scalerStringToInt_[scalerString], param1, param2) << std::endl; 
    inputTensors_.at(kBEST_).second.matrix<float>()(0, n) = getScaledValue(value, scalerStringToInt_[scalerString], param1, param2); 
  }

  // Evaluate network
  tensorflow::run(&(cache_->getSession()), inputTensors_, {outputName_}, &pred_vector);

  prediction = tensorflow::Tensor(tensorflow::DT_FLOAT, {1, nCategories}); 

  for (int k = 0; k < nCategories; ++k) 
  {
    const float pred = pred_vector[0].flat<float>()(k);
    if (!(pred >= 0 && pred <= 1)) 
    {
      throw cms::Exception("BESTEvaluation")  << "invalid prediction = " << pred << " for pred_index = " << k;
    }
    prediction.matrix<float>()(0, k) = pred;
  }

  // abbott: Why fill prediction.matrix instead of BESTScores directly? 

  for (int i = 0; i < nCategories; i++)
  {
    BESTScores.push_back(prediction.matrix<float>()(0,i)); // 0,5 determined by number of outputs
  }

  return BESTScores;
} */

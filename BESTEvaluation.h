/**
 * @file BESTEvaluation.h
 * @brief Header file containing the declaration of the BESTEvaluation class.
 * 
 * This file defines the BESTEvaluation class, which is responsible for evaluating BEST (Boosted Event Shape Tagger) predictions.
 * It provides methods for configuring the evaluation, obtaining predictions, and storing related information.
 * This class is used by the BESTAnalyzer class to evaluate BEST predictions for each event.
 */
#ifndef BEST_EVALUATION_H
#define BEST_EVALUATION_H

#include <Math/VectorUtil.h>
//#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
//#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
//#include "DataFormats/Common/interface/AssociationMap.h"
#include "CacheHandler.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
//#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

/**
 * @class BESTEvaluation
 * @brief Class for loading and evaluating BEST (Boosted Event Shape Tagger) predictions.
 * 
 * The BESTEvaluation class provides methods for configuring the BEST model,
 * obtaining predictions based on input variables, and managing the model parameters.
 * It also includes member variables for storing the list of BEST input features (BESvars),
 * the configuration status, and the cache handler.
 */
class BESTEvaluation {
  public:
    BESTEvaluation(const CacheHandler* cache)
    : isConfigured_(false),
    cache_(cache)
    {}
    void configure(const edm::ParameterSet&);
    //std::vector<float> getPrediction(const std::vector<float> &BESTInputs);
    std::vector<float> getPrediction(std::map<std::string, float> &BESTInputs);

    ~BESTEvaluation() {}
    std::vector<std::string> listOfBESTVars;

  private:

    // This function scales the input variables to the BEST Neural Network to be used in the BESTAnalyzer.
    // The ScalerParameters.txt (loaded in BEST::configure) contains each BESvar's name, scaler, and scaler parameters.
    // More info on the scalers used from sci-kit learn: https://scikit-learn.org/stable/modules/preprocessing.html
    static float getScaledValue(float value, int scaler, float param1, float param2){
      // std::unordered_map<std::string,int> const scalerInt = { {"NoScale",0}, {"MinMax",1}, {"Standard",2}, {"MaxAbs",3} };
      
      float result;
      switch(scaler){

        case 0: // NoScale: result = value
          result = value; break;

        case 1: // MinMax: result = ( value - var_min[param1] ) / ( var_max[param2] - var_min[param1] )
          result = (value-param1)/(param2-param1); break;

        case 2: // Standard: result = (value - mean[param2]) / std_dev[param1]
          result = (value-param2)/param1; break;
        
        case 3: // MaxAbs: result = value / max_value[param1]
          result = value/param1; break;

        default: // If one of the four scalers is not found, something bad is happening
          throw cms::Exception("BESTEvaluation") << "invalid scaler = " << scaler; break;
      }
      return result;
    }

    bool isConfigured_;
    const CacheHandler* cache_;

    std::string BESTname_;
    edm::FileInPath BESTpath_;
    edm::FileInPath BESTscale_;
    std::string runYear_;

    std::map< std::string, std::tuple<std::string,float,float> > paramDict_;
    std::unordered_map<std::string,int> scalerStringToInt_ = { {"NoScale",0}, {"MinMax",1}, {"Standard",2}, {"MaxAbs",3} };

    size_t kBEST_ = 0;  // kBEST is the number of different input branches
    std::vector<tensorflow::TensorShape> inputShapes_;
    tensorflow::NamedTensorList inputTensors_;
    
    int nCategories = 2;

    int NumBESTInputs_ = 110; 
    std::vector<std::string> inputNames_; 
    std::string outputName_;
};
#endif

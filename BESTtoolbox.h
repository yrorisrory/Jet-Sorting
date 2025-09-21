//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// BESTtoolbox.h --------------------------------------------------------------
//=================================================================================
// Header file containing functions for use with CMS EDAnalyzer and EDProducer ----
///////////////////////////////////////////////////////////////////////////////////

// make sure the functions are not declared more than once
#ifndef BESTtoolbox_H
#define BESTtoolbox_H

// include files
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "PhysicsTools/CandUtils/interface/EventShapeVariables.h"
#include "PhysicsTools/CandUtils/interface/Thrust.h"
#include "TMath.h"
#include "TLorentzVector.h"

// Fast Jet Include files
#include <fastjet/JetDefinition.hh>
#include <fastjet/PseudoJet.hh>
#include "fastjet/tools/Filter.hh"
#include <fastjet/ClusterSequence.hh>
#include <fastjet/ClusterSequenceArea.hh>

///////////////////////////////////////////////////////////////////////////////////
// Functions ----------------------------------------------------------------------
///////////////////////////////////////////////////////////////////////////////////

// calculate Legendre Polynomials
float LegendreP(float x, int order);

// calculate Fox Wolfram moments
int FWMoments(std::vector<TLorentzVector> particles, double (&outputs)[5] );

#endif


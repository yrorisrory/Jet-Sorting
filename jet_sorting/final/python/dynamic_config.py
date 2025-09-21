import FWCore.ParameterSet.VarParsing as VarParsing

options = VarParsing.VarParsing('analysis')
options.parseArguments()

from PhysicsTools.PatAlgos.tools.helpers import getPatAlgosToolsTask
import FWCore.ParameterSet.Config as cms

process = cms.Process("analysis")
process.patAlgosToolsTask = getPatAlgosToolsTask(process)

process.load('Configuration.StandardSequences.Services_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('JetMETCorrections.Configuration.JetCorrectors_cff')
process.load('JetMETCorrections.Configuration.CorrectedJetProducers_cff')
process.load('JetMETCorrections.Configuration.CorrectedJetProducersDefault_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '106X_mcRun2_asymptotic_v17', '')

from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_v17'
process.options = cms.untracked.PSet( allowUnscheduled = cms.untracked.bool(True) )
from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
################# JEC ################
corrLabels = ['L2Relative', 'L3Absolute']
from PhysicsTools.PatAlgos.tools.jetTools import *
from RecoBTag.ONNXRuntime.pfDeepBoostedJet_cff import *
updateJetCollection(
 process,
 jetSource = cms.InputTag('slimmedJetsAK8'),
 labelName = 'AK8',
 jetCorrections = ('AK8PFPuppi', cms.vstring(corrLabels), 'None'), #previous corrections: 'L2Relative', 'L3Absolute', 'L2L3Residual'
 postfix = 'UpdatedJEC',
 printWarning = False
)
corrLabels_AK4 = ['L1FastJet', 'L2Relative', 'L3Absolute']
updateJetCollection(
 process,
 jetSource = cms.InputTag('slimmedJets'),
 labelName = 'AK4',
 jetCorrections = ('AK4PFchs', cms.vstring(corrLabels_AK4), 'None'),
 postfix = 'UpdatedJEC',
 printWarning = False
)  
################# Jet PU ID ################
from RecoJets.JetProducers.PileupJetID_cfi import pileupJetId
process.load("RecoJets.JetProducers.PileupJetID_cfi")
process.pileupJetIdUpdated = process.pileupJetId.clone( 
jets=cms.InputTag("updatedPatJetsAK4UpdatedJEC"),  #should be the name of the post-JEC jet collection
inputIsCorrected=True,
applyJec=False,
vertexes=cms.InputTag("offlineSlimmedPrimaryVertices"),
)
process.patAlgosToolsTask.add(process.pileupJetIdUpdated)
updateJetCollection(    # running in unscheduled mode, need to manually update collection
 process,
 labelName = "PileupJetID",
 jetSource = cms.InputTag("updatedPatJetsAK4UpdatedJEC"),
)
process.updatedPatJetsPileupJetID.userData.userInts.src = ["pileupJetIdUpdated:fullId"]
process.content = cms.EDAnalyzer("EventContentAnalyzer")
##############################################################################
process.leptonVeto = cms.EDFilter("leptonVeto",
 muonCollection= cms.InputTag("slimmedMuons"),
 electronCollection = cms.InputTag("slimmedElectrons"),
 metCollection = cms.InputTag("slimmedMETs"),
 tauCollection = cms.InputTag("slimmedTaus")
)
from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
process.prefiringweight = l1PrefiringWeightProducer.clone(
TheJets = cms.InputTag('updatedPatJetsPileupJetID'), #updatedPatJetsAK4UpdatedJEC 
DataEraECAL = cms.string('None'), #Use 2016BtoH for 2016
DataEraMuon = cms.string('20172018'), #Use 2016 for 2016
UseJetEMPt = cms.bool(False),
PrefiringRateSystematicUnctyECAL = cms.double(0.2),
PrefiringRateSystematicUnctyMuon = cms.double(0.2)
)
process.clusteringAnalyzerAll_nom = cms.EDAnalyzer("matchingTest",
 runType = cms.string("Suu4_chi1_HTHT"), #types: QCDMC1000to1500, QCDMC1500to2000, QCDMC2000toInf, TTbarMC, DataA, etc. , Suu8_chi3, etc.
 genPartCollection = cms.InputTag("prunedGenParticles"),
 genParticles = cms.InputTag("genParticles"),
# Rory had to change this input tag to match his gen particle tags
 packedGenParticles = cms.InputTag("packedGenParticles", "", "PAT"),
 BESTname = cms.string('BESTGraph'),  BESTpath = cms.FileInPath('data/BEST_models/constantgraph_combine.pb'),
 BESTscale = cms.FileInPath('data/BESTScalerParameters_all_mass_combine.txt'),
 PUfile_path = cms.FileInPath('data/POG/LUM/2018_UL/puWeights.json'),
 bTagEff_path = cms.FileInPath('data/btaggingEffMaps/btag_efficiency_map_SuuToChiChi_combined_2018.root'),
 bTagSF_path = cms.FileInPath('data/bTaggingSFs/2018_UL/btagging.json'),
 JECUncert_AK8_path = cms.FileInPath("data/JEC/2018_UL/MC/Summer19UL18_V5_MC_Uncertainty_AK8PFPuppi.txt"),
 JECUncert_AK4_path = cms.FileInPath("data/JEC/2018_UL/MC/Summer19UL18_V5_MC_Uncertainty_AK4PFchs.txt"),
 fatJetCollection = cms.InputTag("updatedPatJetsAK8UpdatedJEC"),
 jetCollection = cms.InputTag("updatedPatJetsPileupJetID"),
 muonCollection= cms.InputTag("slimmedMuons"),
 electronCollection = cms.InputTag("slimmedElectrons"),
 metCollection = cms.InputTag("slimmedMETs"),
 tauCollection = cms.InputTag("slimmedTaus"),
 pileupCollection = cms.InputTag("slimmedAddPileupInfo"),
 systematicType = cms.string("nom"),
 jetVetoMapFile = cms.FileInPath("data/jetVetoMaps/hotjets-UL18.root"),
 jetVetoMapName = cms.string("h2hot_ul18_plus_hem1516_and_hbp2m1"),
 includeAllBranches = cms.bool(False),
 slimmedSelection   = cms.bool(True),
 verbose            = cms.bool(False),
 runSideband            = cms.bool(False),
 year = cms.string("2018"), #types: 2015,2016,2017,2018
 genEventInfoTag=cms.InputTag("generator"),
 lheEventInfoTag=cms.InputTag("externalLHEProducer"),
 bits = cms.InputTag("TriggerResults", "", "HLT"),
 triggers = cms.string("HLT_PFHT900_v"),
 doPUID = cms.bool(True),
 doPDF = cms.bool(False)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(options.inputFiles)
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(options.outputFile)
)

process.options = cms.untracked.PSet(
 wantSummary = cms.untracked.bool(True),
)
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(400))
process.p = cms.Path(  process.pileupJetIdUpdated * process.leptonVeto * process.clusteringAnalyzerAll_nom #* process.prefiringweight
)
process.pathRunPatAlgos = cms.Path(process.patAlgosToolsTask)

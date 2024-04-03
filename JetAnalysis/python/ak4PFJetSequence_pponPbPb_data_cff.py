import FWCore.ParameterSet.Config as cms

from HeavyIonsAnalysis.JetAnalysis.inclusiveJetAnalyzer_cff import *

ak4PFJetAnalyzer = inclusiveJetAnalyzer.clone(
    jetTag = cms.InputTag("slimmedJets"),
    rParam = 0.4,
    fillGenJets = False,
    isMC = False,
    jetName = cms.untracked.string("ak4PF"),
    hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
    )

# ak2PFJetAnalyzer = inclusiveJetAnalyzer.clone(
#     jetTag = cms.InputTag("slimmedJets"),
#     rParam = 0.2,
#     fillGenJets = False,
#     isMC = False,
#     jetName = cms.untracked.string("ak2PF"),
#     hltTrgResults = cms.untracked.string('TriggerResults::'+'HISIGNAL'),
#     )

akFlowPu4PFJetAnalyzer = ak4PFJetAnalyzer.clone()

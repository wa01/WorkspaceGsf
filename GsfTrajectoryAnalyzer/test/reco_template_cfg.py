# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RECO -s RAW2DIGI,L1Reco,RECO,ALCA:EcalUncalWElectron+EcalUncalZElectron+EcalESAlign+HcalCalIterativePhiSym+HcalCalIsoTrkFilter,EI,PAT,DQM:@allForPrompt --runUnscheduled --nThreads 4 --data --era Run2_2016 --scenario pp --conditions 80X_dataRun2_2016LegacyRepro_v3 --eventcontent AOD,MINIAOD,DQM --datatier AOD,MINIAOD,DQMIO --customise Configuration/DataProcessing/RecoTLR.customisePostEra_Run2_2016,RecoTracker/Configuration/customizeMinPtForHitRecoveryInGluedDet.customizeHitRecoveryInGluedDetTkSeedsOnly --filein file:pippo.root -n 100 --python_filename=recoskim_Run2016B-v2_SingleElectron.py --no_exec
import FWCore.ParameterSet.Config as cms

from Configuration.StandardSequences.Eras import eras

process = cms.Process('RECO',eras.Run2_2016)

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.L1Reco_cff')
process.load('Configuration.StandardSequences.Reconstruction_Data_cff')
process.load('Configuration.StandardSequences.AlCaRecoStreams_cff')
process.load('CommonTools.ParticleFlow.EITopPAG_cff')
process.load('PhysicsTools.PatAlgos.slimming.metFilterPaths_cff')
process.load('DQMOffline.Configuration.DQMOffline_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3)
)
#process.MessageLogger.FrameworkJobReport.default.limit = 1000

# Input source
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/t/tpook/public/dEta_03Feb_raw.root'),
    fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016B-03Feb2017_Rec_1Ele+X_DEta_110_274999_1970753919.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016B-03Feb2017_Rec_1Ele+X_DEta_133_274200_720352580.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016B-03Feb2017_Rec_1Ele+X_DEta_14_275125_927024589.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016B-03Feb2017_Rec_1Ele+X_DEta_17_274244_596153446.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016B-03Feb2017_Rec_1Ele+X_DEta_186_275310_349349861.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016B-03Feb2017_Rec_1Ele+X_DEta_227_275125_188399335.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016B-03Feb2017_Rec_1Ele+X_DEta_38_274999_1460665511.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016B-03Feb2017_Rec_1Ele+X_DEta_79_274388_2781933717.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016B-03Feb2017_Rec_1Ele+X_DEta_95_274388_1422660397.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016C-03Feb2017-v1_Rec_1Ele+X_DEta_104_275911_327701147.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016C-03Feb2017-v1_Rec_1Ele+X_DEta_160_276242_2737025591.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016C-03Feb2017-v1_Rec_1Ele+X_DEta_19_276242_1731111950.root',
                                      'file:/afs/cern.ch/user/t/tpook/public/03Feb_raw/pick_DoubleEGRun2016C-03Feb2017-v1_Rec_1Ele+X_DEta_216_276097_186490714.root'
                                      ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(
    allowUnscheduled = cms.untracked.bool(True)
)

## Production Info
#process.configurationMetadata = cms.untracked.PSet(
#    annotation = cms.untracked.string('RECO nevts:100'),
#    name = cms.untracked.string('Applications'),
#    version = cms.untracked.string('$Revision: 1.19 $')
#)

process.TFileService = cms.Service("TFileService", 
                                   fileName = cms.string("histo.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )

process.contentAnalyzer = cms.EDAnalyzer("EventContentAnalyzer")
process.trajectoryAnalyzer = cms.EDAnalyzer("GsfTrajectoryAnalyzer",
                                            trajectories = cms.InputTag("electronGsfTracks")
)
process.contentAnalyzer_step = cms.Path(#process.contentAnalyzer+
                                        process.trajectoryAnalyzer)

# Output definition

process.RECOoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('RECO'),
        filterName = cms.untracked.string('')
    ),
#    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('RECO_RAW2DIGI_L1Reco_RECO.root'),
    outputCommands = process.RECOEventContent.outputCommands
)

process.AODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('AOD'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fileName = cms.untracked.string('RECO_RAW2DIGI_L1Reco_RECO_ALCA_EI_PAT_DQM.root'),
    outputCommands = process.AODEventContent.outputCommands
)

process.MINIAODoutput = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAOD'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('RECO_RAW2DIGI_L1Reco_RECO_ALCA_EI_PAT_DQM_inMINIAOD.root'),
    outputCommands = process.MINIAODEventContent.outputCommands,
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)

process.DQMoutput = cms.OutputModule("DQMRootOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('DQMIO'),
        filterName = cms.untracked.string('')
    ),
    fileName = cms.untracked.string('RECO_RAW2DIGI_L1Reco_RECO_ALCA_EI_PAT_DQM_inDQM.root'),
    outputCommands = process.DQMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition
process.ALCARECOStreamEcalESAlign = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOEcalESAlign')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('ALCARECO'),
        filterName = cms.untracked.string('EcalESAlign')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('EcalESAlign.root'),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep ESDCCHeaderBlocksSorted_ecalPreshowerDigis_*_*', 
        'keep ESDigiCollection_ecalPreshowerDigis_*_*', 
        'keep ESKCHIPBlocksSorted_ecalPreshowerDigis_*_*', 
        'keep SiPixelClusteredmNewDetSetVector_ecalAlCaESAlignTrackReducer_*_*', 
        'keep SiStripClusteredmNewDetSetVector_ecalAlCaESAlignTrackReducer_*_*', 
        'keep TrackingRecHitsOwned_ecalAlCaESAlignTrackReducer_*_*', 
        'keep recoTrackExtras_ecalAlCaESAlignTrackReducer_*_*', 
        'keep recoTracks_ecalAlCaESAlignTrackReducer_*_*', 
        'keep recoBeamSpot_offlineBeamSpot_*_*')
)
process.ALCARECOStreamEcalUncalWElectron = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOEcalUncalWElectron')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('ALCARECO'),
        filterName = cms.untracked.string('EcalUncalWElectron')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('EcalUncalWElectron.root'),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_pfMet_*_*', 
        'keep *_kt6PFJetsForRhoCorrection_rho_*', 
        'keep *_kt6PFJets_rho_*', 
        'keep recoVertexs_offlinePrimaryVertices_*_*', 
        'keep recoVertexs_offlinePrimaryVerticesWithBS_*_*', 
        'keep *BeamSpot_offlineBeamSpot_*_*', 
        'keep *_allConversions_*_*', 
        'keep *_conversions_*_*', 
        'keep *GsfTrack*_electronGsfTracks_*_*', 
        'keep *GsfTrack*_uncleanedOnlyElectronGsfTracks_*_*', 
        'keep *_generator_*_*', 
        'keep *_addPileupInfo_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGsfElectron*_gsfElectron*_*_*', 
        'keep recoGsfElectron*_gedGsfElectrons_*_*', 
        'keep recoGsfElectron*_gedGsfElectronsTmp_*_*', 
        'keep recoGsfElectron*_gedGsfElectronCores_*_*', 
        'keep recoPhoton*_gedPhoton_*_*', 
        'keep recoCaloClusters_hfEMClusters_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_alCaIsolatedElectrons_*_*', 
        'keep recoCaloClusters_cleanedHybridSuperClusters_*_*', 
        'keep recoCaloClusters_hybridSuperClusters_*_*', 
        'keep recoCaloClusters_uncleanedHybridSuperClusters_*_*', 
        'keep recoCaloClusters_multi5x5BasicClustersCleaned_*_*', 
        'keep recoCaloClusters_multi5x5BasicClustersUncleaned_*_*', 
        'keep recoCaloClusters_multi5x5SuperClusters_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoSuperClusters_SCselector_*_*', 
        'keep recoSuperClusters_cleanedHybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_hfEMClusters_*_*', 
        'keep recoSuperClusters_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_mergedSuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_uncleanedHybridSuperClusters_*_*', 
        'keep recoSuperClusters_uncleanedOnlyCorrectedHybridSuperClusters_*_*', 
        'keep recoSuperClusters_uncleanedOnlyCorrectedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_uncleanedOnlyMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersCleaned_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersUncleaned_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoPreshowerCluster*_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerCluster*_uncleanedOnlyMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerCluster*_multi5x5PreshowerClusterShape_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_alcaElectronTracksReducer_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep edmConditionsInEventBlock_conditionsInEdm_*_*', 
        'keep edmConditionsInLumiBlock_conditionsInEdm_*_*', 
        'keep edmConditionsInRunBlock_conditionsInEdm_*_*', 
        'keep *_TriggerResults_*_*', 
        'keep *_hltTriggerSummaryAOD_*_HLT', 
        'keep *EcalRecHit*_alCaIsolatedElectrons_*_*', 
        'keep *EcalRecHit*_reducedEcalRecHitsES_alCaRecHitsES_*', 
        'keep *_ecalDigis_*_*', 
        'keep *EcalTriggerPrimitiveDigi*_ecalDigis_*_*', 
        'keep *_ecalPreshowerDigis_*_*', 
        'keep *_ecalDetIdToBeRecovered_*_*', 
        'keep reco*Clusters_pfElectronTranslator_*_*', 
        'drop recoCaloClusters_*_*_*', 
        'drop recoSuperClusters_*_*_*', 
        'drop recoPreshowerCluster*_*_*_*', 
        'drop *EcalRecHit*_reducedEcalRecHitsES*_*_*', 
        'keep reco*Clusters_pfElectronTranslator_*_*')
)
process.ALCARECOStreamEcalUncalZElectron = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOEcalUncalZElectron', 
            'pathALCARECOEcalUncalZSCElectron')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('ALCARECO'),
        filterName = cms.untracked.string('EcalUncalZElectron')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('EcalUncalZElectron.root'),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_pfMet_*_*', 
        'keep *_kt6PFJetsForRhoCorrection_rho_*', 
        'keep *_kt6PFJets_rho_*', 
        'keep recoVertexs_offlinePrimaryVertices_*_*', 
        'keep recoVertexs_offlinePrimaryVerticesWithBS_*_*', 
        'keep *BeamSpot_offlineBeamSpot_*_*', 
        'keep *_allConversions_*_*', 
        'keep *_conversions_*_*', 
        'keep *GsfTrack*_electronGsfTracks_*_*', 
        'keep *GsfTrack*_uncleanedOnlyElectronGsfTracks_*_*', 
        'keep *_generator_*_*', 
        'keep *_addPileupInfo_*_*', 
        'keep *_genParticles_*_*', 
        'keep recoGsfElectron*_gsfElectron*_*_*', 
        'keep recoGsfElectron*_gedGsfElectrons_*_*', 
        'keep recoGsfElectron*_gedGsfElectronsTmp_*_*', 
        'keep recoGsfElectron*_gedGsfElectronCores_*_*', 
        'keep recoPhoton*_gedPhoton_*_*', 
        'keep recoCaloClusters_hfEMClusters_*_*', 
        'keep recoCaloClusters_particleFlowEGamma_*_*', 
        'keep recoCaloClusters_alCaIsolatedElectrons_*_*', 
        'keep recoCaloClusters_cleanedHybridSuperClusters_*_*', 
        'keep recoCaloClusters_hybridSuperClusters_*_*', 
        'keep recoCaloClusters_uncleanedHybridSuperClusters_*_*', 
        'keep recoCaloClusters_multi5x5BasicClustersCleaned_*_*', 
        'keep recoCaloClusters_multi5x5BasicClustersUncleaned_*_*', 
        'keep recoCaloClusters_multi5x5SuperClusters_*_*', 
        'keep recoCaloClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoSuperClusters_SCselector_*_*', 
        'keep recoSuperClusters_cleanedHybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedHybridSuperClusters_*_*', 
        'keep recoSuperClusters_correctedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_hfEMClusters_*_*', 
        'keep recoSuperClusters_hybridSuperClusters_*_*', 
        'keep recoSuperClusters_mergedSuperClusters_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_particleFlowEGamma_*_*', 
        'keep recoSuperClusters_uncleanedHybridSuperClusters_*_*', 
        'keep recoSuperClusters_uncleanedOnlyCorrectedHybridSuperClusters_*_*', 
        'keep recoSuperClusters_uncleanedOnlyCorrectedMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_uncleanedOnlyMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersCleaned_*_*', 
        'keep recoSuperClusters_multi5x5SuperClustersUncleaned_*_*', 
        'keep recoSuperClusters_multi5x5SuperClusters_*_*', 
        'keep recoSuperClusters_particleFlowSuperClusterECAL_*_*', 
        'keep recoPreshowerCluster*_multi5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerCluster*_uncleanedOnlyMulti5x5SuperClustersWithPreshower_*_*', 
        'keep recoPreshowerCluster*_multi5x5PreshowerClusterShape_*_*', 
        'keep *_pfElectronTranslator_*_*', 
        'keep *_alcaElectronTracksReducer_*_*', 
        'keep *_l1extraParticles_*_*', 
        'keep L1GlobalTriggerReadoutRecord_gtDigis_*_*', 
        'keep *_l1L1GtObjectMap_*_*', 
        'keep edmConditionsInEventBlock_conditionsInEdm_*_*', 
        'keep edmConditionsInLumiBlock_conditionsInEdm_*_*', 
        'keep edmConditionsInRunBlock_conditionsInEdm_*_*', 
        'keep *_TriggerResults_*_*', 
        'keep *_hltTriggerSummaryAOD_*_HLT', 
        'keep *EcalRecHit*_alCaIsolatedElectrons_*_*', 
        'keep *EcalRecHit*_reducedEcalRecHitsES_alCaRecHitsES_*', 
        'keep *_ecalDigis_*_*', 
        'keep *EcalTriggerPrimitiveDigi*_ecalDigis_*_*', 
        'keep *_ecalPreshowerDigis_*_*', 
        'keep *_ecalDetIdToBeRecovered_*_*', 
        'keep reco*Clusters_pfElectronTranslator_*_*', 
        'drop recoCaloClusters_*_*_*', 
        'drop recoSuperClusters_*_*_*', 
        'drop recoPreshowerCluster*_*_*_*', 
        'drop *EcalRecHit*_reducedEcalRecHitsES*_*_*', 
        'keep reco*Clusters_pfElectronTranslator_*_*')
)
process.ALCARECOStreamHcalCalIsoTrkFilter = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOHcalCalIsoTrkFilter')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('ALCARECO'),
        filterName = cms.untracked.string('HcalCalIsoTrkFilter')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('HcalCalIsoTrkFilter.root'),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_hbhereco_*_*', 
        'keep *_ecalRecHit_*_*', 
        'keep *_offlineBeamSpot_*_*', 
        'keep *_hltTriggerSummaryAOD_*_*', 
        'keep *_TriggerResults_*_*', 
        'keep *_generalTracks_*_*', 
        'keep *_generalTracksExtra_*_*', 
        'keep *_offlinePrimaryVertices_*_*', 
        'keep *_TkAlIsoProdFilter_*_*')
)
process.ALCARECOStreamHcalCalIterativePhiSym = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('pathALCARECOHcalCalIterativePhiSym')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('ALCARECO'),
        filterName = cms.untracked.string('ALCARECOHcalCalIterativePhiSym')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('ALCARECOHcalCalIterativePhiSym.root'),
    outputCommands = cms.untracked.vstring('drop *', 
        'keep *_horeco_*_*', 
        'keep *_hfreco_*_*', 
        'keep *_hbhereco_*_*', 
        'keep *_offlinePrimaryVertices_*_*', 
        'keep edmTriggerResults_*_*_HLT')
)

# Other statements
process.ALCARECOEventContent.outputCommands.extend(process.OutALCARECOEcalESAlign_noDrop.outputCommands)
process.ALCARECOEventContent.outputCommands.extend(process.OutALCARECOEcalUncalWElectron_noDrop.outputCommands)
process.ALCARECOEventContent.outputCommands.extend(process.OutALCARECOHcalCalIsoTrkFilter_noDrop.outputCommands)
process.ALCARECOEventContent.outputCommands.extend(process.OutALCARECOHcalCalIterativePhiSym_noDrop.outputCommands)
process.ALCARECOEventContent.outputCommands.extend(process.OutALCARECOEcalUncalZElectron_noDrop.outputCommands)
from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016LegacyRepro_v3', '')
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')

# Path and EndPath definitions
process.raw2digi_step = cms.Path(process.RawToDigi)
process.L1Reco_step = cms.Path(process.L1Reco)
process.reconstruction_step = cms.Path(process.reconstruction)
process.eventinterpretaion_step = cms.Path(process.EIsequence)
process.Flag_trackingFailureFilter = cms.Path(process.goodVertices+process.trackingFailureFilter)
process.Flag_goodVertices = cms.Path(process.primaryVertexFilter)
process.Flag_CSCTightHaloFilter = cms.Path(process.CSCTightHaloFilter)
process.Flag_trkPOGFilters = cms.Path(process.trkPOGFilters)
process.Flag_HcalStripHaloFilter = cms.Path(process.HcalStripHaloFilter)
process.Flag_trkPOG_logErrorTooManyClusters = cms.Path(~process.logErrorTooManyClusters)
process.Flag_EcalDeadCellTriggerPrimitiveFilter = cms.Path(process.EcalDeadCellTriggerPrimitiveFilter)
process.Flag_ecalLaserCorrFilter = cms.Path(process.ecalLaserCorrFilter)
process.Flag_globalSuperTightHalo2016Filter = cms.Path(process.globalSuperTightHalo2016Filter)
process.Flag_eeBadScFilter = cms.Path(process.eeBadScFilter)
process.Flag_METFilters = cms.Path(process.metFilters)
process.Flag_chargedHadronTrackResolutionFilter = cms.Path(process.chargedHadronTrackResolutionFilter)
process.Flag_globalTightHalo2016Filter = cms.Path(process.globalTightHalo2016Filter)
process.Flag_CSCTightHaloTrkMuUnvetoFilter = cms.Path(process.CSCTightHaloTrkMuUnvetoFilter)
process.Flag_HBHENoiseIsoFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseIsoFilter)
process.Flag_BadChargedCandidateSummer16Filter = cms.Path(process.BadChargedCandidateSummer16Filter)
process.Flag_hcalLaserEventFilter = cms.Path(process.hcalLaserEventFilter)
process.Flag_BadPFMuonFilter = cms.Path(process.BadPFMuonFilter)
process.Flag_HBHENoiseFilter = cms.Path(process.HBHENoiseFilterResultProducer+process.HBHENoiseFilter)
process.Flag_trkPOG_toomanystripclus53X = cms.Path(~process.toomanystripclus53X)
process.Flag_EcalDeadCellBoundaryEnergyFilter = cms.Path(process.EcalDeadCellBoundaryEnergyFilter)
process.Flag_BadChargedCandidateFilter = cms.Path(process.BadChargedCandidateFilter)
process.Flag_trkPOG_manystripclus53X = cms.Path(~process.manystripclus53X)
process.Flag_BadPFMuonSummer16Filter = cms.Path(process.BadPFMuonSummer16Filter)
process.Flag_muonBadTrackFilter = cms.Path(process.muonBadTrackFilter)
process.Flag_CSCTightHalo2015Filter = cms.Path(process.CSCTightHalo2015Filter)
process.dqmoffline_step = cms.Path(process.DQMOfflineCommon)
process.dqmoffline_1_step = cms.Path(process.DQMOfflineMuon)
process.dqmoffline_2_step = cms.Path(process.DQMOfflineHcal)
process.dqmoffline_3_step = cms.Path(process.DQMOfflineJetMET)
process.dqmoffline_4_step = cms.Path(process.DQMOfflineEcal)
process.dqmofflineOnPAT_step = cms.EndPath(process.PostDQMOffline)
process.RECOoutput_step = cms.EndPath(process.RECOoutput)
#process.AODoutput_step = cms.EndPath(process.AODoutput)
#process.MINIAODoutput_step = cms.EndPath(process.MINIAODoutput)
#process.DQMoutput_step = cms.EndPath(process.DQMoutput)
process.ALCARECOStreamEcalESAlignOutPath = cms.EndPath(process.ALCARECOStreamEcalESAlign)
process.ALCARECOStreamEcalUncalWElectronOutPath = cms.EndPath(process.ALCARECOStreamEcalUncalWElectron)
process.ALCARECOStreamEcalUncalZElectronOutPath = cms.EndPath(process.ALCARECOStreamEcalUncalZElectron)
process.ALCARECOStreamHcalCalIsoTrkFilterOutPath = cms.EndPath(process.ALCARECOStreamHcalCalIsoTrkFilter)
process.ALCARECOStreamHcalCalIterativePhiSymOutPath = cms.EndPath(process.ALCARECOStreamHcalCalIterativePhiSym)

# Schedule definition
process.schedule = cms.Schedule(process.raw2digi_step,process.L1Reco_step,process.reconstruction_step,process.contentAnalyzer_step,process.RECOoutput_step
#,process.pathALCARECOEcalESAlign,process.pathALCARECOEcalUncalWElectron,process.pathALCARECOHcalCalIsoTrkFilter,process.pathALCARECOHcalCalIterativePhiSym,process.pathALCARECOEcalUncalZElectron,process.pathALCARECOEcalUncalZSCElectron,process.eventinterpretaion_step,process.Flag_HBHENoiseFilter,process.Flag_HBHENoiseIsoFilter,process.Flag_CSCTightHaloFilter,process.Flag_CSCTightHaloTrkMuUnvetoFilter,process.Flag_CSCTightHalo2015Filter,process.Flag_globalTightHalo2016Filter,process.Flag_globalSuperTightHalo2016Filter,process.Flag_HcalStripHaloFilter,process.Flag_hcalLaserEventFilter,process.Flag_EcalDeadCellTriggerPrimitiveFilter,process.Flag_EcalDeadCellBoundaryEnergyFilter,process.Flag_goodVertices,process.Flag_eeBadScFilter,process.Flag_ecalLaserCorrFilter,process.Flag_trkPOGFilters,process.Flag_chargedHadronTrackResolutionFilter,process.Flag_muonBadTrackFilter,process.Flag_BadChargedCandidateFilter,process.Flag_BadPFMuonFilter,process.Flag_BadChargedCandidateSummer16Filter,process.Flag_BadPFMuonSummer16Filter,process.Flag_trkPOG_manystripclus53X,process.Flag_trkPOG_toomanystripclus53X,process.Flag_trkPOG_logErrorTooManyClusters,process.Flag_METFilters,process.dqmoffline_step,process.dqmoffline_1_step,process.dqmoffline_2_step,process.dqmoffline_3_step,process.dqmoffline_4_step,process.dqmofflineOnPAT_step,process.AODoutput_step,process.MINIAODoutput_step,process.DQMoutput_step,process.ALCARECOStreamEcalESAlignOutPath,process.ALCARECOStreamEcalUncalWElectronOutPath,process.ALCARECOStreamEcalUncalZElectronOutPath,process.ALCARECOStreamHcalCalIsoTrkFilterOutPath,process.ALCARECOStreamHcalCalIterativePhiSymOutPath
)

#Setup FWK for multithreaded
process.options.numberOfThreads=cms.untracked.uint32(2)
process.options.numberOfStreams=cms.untracked.uint32(0)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.DataProcessing.RecoTLR
from Configuration.DataProcessing.RecoTLR import customisePostEra_Run2_2016 

#call to customisation function customisePostEra_Run2_2016 imported from Configuration.DataProcessing.RecoTLR
process = customisePostEra_Run2_2016(process)

# Automatic addition of the customisation function from RecoTracker.Configuration.customizeMinPtForHitRecoveryInGluedDet
from RecoTracker.Configuration.customizeMinPtForHitRecoveryInGluedDet import customizeHitRecoveryInGluedDetTkSeedsOnly 

#call to customisation function customizeHitRecoveryInGluedDetTkSeedsOnly imported from RecoTracker.Configuration.customizeMinPtForHitRecoveryInGluedDet
process = customizeHitRecoveryInGluedDetTkSeedsOnly(process)

# End of customisation functions
#do not add changes to your config after this point (unless you know what you are doing)
from FWCore.ParameterSet.Utilities import convertToUnscheduled
process=convertToUnscheduled(process)
process.load('Configuration.StandardSequences.PAT_cff')
from FWCore.ParameterSet.Utilities import cleanUnscheduled
process=cleanUnscheduled(process)

# customisation of the process.

# Automatic addition of the customisation function from PhysicsTools.PatAlgos.slimming.miniAOD_tools
from PhysicsTools.PatAlgos.slimming.miniAOD_tools import miniAOD_customizeAllData 

#call to customisation function miniAOD_customizeAllData imported from PhysicsTools.PatAlgos.slimming.miniAOD_tools
process = miniAOD_customizeAllData(process)

# End of customisation functions

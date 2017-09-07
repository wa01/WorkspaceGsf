from reco_template_cfg import process

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016LegacyRepro_v3', '')
process.RECOoutput.fileName = cms.untracked.string('reco1_RECO.root')
process.TFileService.fileName = cms.string("reco1_histo.root")

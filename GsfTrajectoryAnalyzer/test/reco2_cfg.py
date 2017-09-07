from reco_template_cfg import process

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '80X_dataRun2_2016SeptRepro_v7', '')
process.RECOoutput.fileName = cms.untracked.string('reco2_RECO.root')
process.TFileService.fileName = cms.string("reco2_histo.root")

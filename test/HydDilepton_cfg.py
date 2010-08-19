import FWCore.ParameterSet.Config as cms

process = cms.Process("GenAna")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use
    duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    fileNames = cms.untracked.vstring(
#    'dcache:/pnfs/cmsaf.mit.edu/t2bat/cms/store/user/davidlw/Hydjet_Quenched_MinBias_2800GeV_GEN_341_020510/Hydjet_Quenched_MinBias_2800GeV_GEN_341_020510/2cd3c91d56ca6f732b8741d2d845caeb/Hydjet_Quenched_MinBias_2800GeV_cfi_py_GEN_1.root'
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0002/B69C6E44-B027-DF11-9D6E-00248C0BE005.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0002/AEC8C081-B127-DF11-A364-00304867C16A.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0002/9295A0A7-3B28-DF11-AAD2-0018F3D09688.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0002/84331DB1-BC27-DF11-8FA5-00304867C0F6.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0002/8275181C-AD27-DF11-ABC7-002618FDA265.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0002/5ADD7B74-A727-DF11-92CB-002618943927.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0002/0E80FA8C-AA27-DF11-BEE5-001A92971B8A.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0001/F6D35C5B-A227-DF11-BE5B-0030486790A0.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0001/EADA70E6-9E27-DF11-AA7E-003048D25B68.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0001/D6D2A1B1-9B27-DF11-9DF2-002354EF3BE1.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0001/9650F38C-9527-DF11-9A86-00304866C398.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0001/84CEDEE3-A027-DF11-AE40-002618943943.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0001/78409A6D-A527-DF11-BD07-002618943962.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0001/5458719F-9E27-DF11-977D-003048D3FC94.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0001/3896BD83-9727-DF11-AAD6-00304867915A.root',
  '/store/relval/CMSSW_3_5_3/RelValHydjetQ_MinBias_4TeV/GEN-SIM-RAW/MC_3XY_V24-v1/0001/085B27DD-A327-DF11-B7F0-003048679182.root'
	 )
)

process.GenAna = cms.EDAnalyzer("HydDileptonAnalyzer",
	#genSource = cms.untracked.InputTag("genParticles"),    
	genSource = cms.untracked.InputTag("hiGenParticles"),    
	hOutputFile = cms.untracked.string("test.root") 
)
process.p = cms.Path(process.GenAna)

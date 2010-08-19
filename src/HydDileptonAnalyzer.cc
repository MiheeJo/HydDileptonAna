// -*- C++ -*-
//
// Package:    HydDileptonAnalyzer
// Class:      HydDileptonAnalyzer
// 
/**\class HydDileptonAnalyzer HydDileptonAnalyzer.cc UserCode/DMoonAna/src/HydDileptonAnalyzer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Dong Ho Moon
//         Created:  Mon Oct  6 14:12:25 CEST 2008
// $Id: HydDileptonAnalyzer.cc,v 1.5 2010/03/01 17:33:52 dmoon Exp $
//
//


// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"


#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "PhysicsTools/UtilAlgos/interface/TFileService.h"

#include <vector>
#include <iostream>


#include "DataFormats/Candidate/interface/Particle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"


#include "FWCore/ParameterSet/interface/InputTag.h"
//#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Ref.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"


//
// class decleration
//

#include <TH1.h>
#include <TFile.h>
#include <TString.h>
#include <stdio.h>
#include <TH2F.h>
#include <TNtuple.h>

#include <TLorentzVector.h>

using namespace edm;
using namespace reco;
using namespace std;

class HydDileptonAnalyzer : public edm::EDAnalyzer {
public:
  explicit HydDileptonAnalyzer(const edm::ParameterSet&);
//  explicit HydDileptonAnalyzer(const edm::EventSetup&);
  ~HydDileptonAnalyzer();
  
  const static int xNum=1000;
  const static int yNum=3;
  
  
private:
//  virtual void beginJob(const edm::EventSetup&) ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
//  virtual void endJob() ;

  std::string fOutputFile_;
  edm::InputTag genSource_; 
  Int_t m_run, m_evt, nmom, nPar, ndau, m_nMom[xNum], m_nDau[xNum];
  Int_t m_pId[xNum], momId[xNum][yNum], dauId[xNum][yNum];

  double p_pt[xNum], p_eta[xNum], p_phi[xNum], p_st[xNum], p_mass[xNum], p_mtum[xNum], 
         p_vx[xNum], p_vy[xNum], p_vz[xNum], p_y[xNum];
  double m_pt[xNum][yNum], m_eta[xNum][yNum], m_phi[xNum][yNum], m_st[xNum][yNum], 
         m_mass[xNum][yNum], m_mtum[xNum][yNum], m_vx[xNum][yNum], m_vy[xNum][yNum], 
         m_vz[xNum][yNum], m_y[xNum][yNum];
  double d_pt[xNum][yNum], d_eta[xNum][yNum], d_phi[xNum][yNum], d_st[xNum][yNum], 
         d_mass[xNum][yNum], d_mtum[xNum][yNum], d_vx[xNum][yNum], d_vy[xNum][yNum], 
         d_vz[xNum][yNum], d_y[xNum][yNum]; 
  
  TFile *hOutputFile;
  TTree *rtree2;
  
  
  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//

HydDileptonAnalyzer::~HydDileptonAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)
    
    hOutputFile->SetCompressionLevel(2);
    hOutputFile->cd();

    hOutputFile->Write();
    hOutputFile->Close();

}

/*HydDileptonAnalyzer::HydDileptonAnalyzer(const edm::ParameterSet& iConfig):
{ }
*/

//HydDileptonAnalyzer::HydDileptonAnalyzer(const edm::EventSetup&)
HydDileptonAnalyzer::HydDileptonAnalyzer(const edm::ParameterSet& iConfig)
//HydDileptonAnalyzer::beginJob(const edm::EventSetup&)
{
    fOutputFile_ = iConfig.getUntrackedParameter<std::string>("hOutputFile") ;
    genSource_ = iConfig.getUntrackedParameter<edm::InputTag>("genSource") ;
    hOutputFile   = new TFile( fOutputFile_.c_str(), "RECREATE" ) ;
//    hOutputFile = new TFile( "hOutputfile", "RECREATE");

    rtree2 = new TTree("ana","HydDilepton analysis ntuple");

    rtree2->Branch("run",        &m_run,       "run/I");
    rtree2->Branch("evt",        &m_evt,       "evt/I");
    rtree2->Branch("nPar",       &nPar,        "nPar/I");

    rtree2->Branch("pId",        &m_pId,       "pId/I");
    rtree2->Branch("p_pt",       &p_pt,        "p_pt[nPar]/D");
    rtree2->Branch("p_y",        &p_y,         "p_y[nPar]/D");
    rtree2->Branch("p_eta",      &p_eta,       "p_eta[nPar]/D");
    rtree2->Branch("p_phi",      &p_phi,       "p_phi[nPar]/D");
    rtree2->Branch("p_st",       &p_st,        "p_st[nPar]/D");
    rtree2->Branch("p_mass",     &p_mass,      "p_mass[nPar]/D");
    //rtree2->Branch("p_vx",       &p_vx,        "p_vx[nPar]/D");
    //rtree2->Branch("p_vy",       &p_vy,        "p_vy[nPar]/D");
    rtree2->Branch("p_vz",       &p_vz,        "p_vz[nPar]/D");

    rtree2->Branch("nMom",       &m_nMom,      "nMom[nPar]/I");
    rtree2->Branch("momId",      &momId,       "momId[nPar][3]/I");
    rtree2->Branch("m_pt",       &m_pt,        "m_pt[nPar][3]/D");
    rtree2->Branch("m_y",        &m_y,        "m_y[nPar][3]/D");
    rtree2->Branch("m_eta",      &m_eta,       "m_eta[nPar][3]/D");
    rtree2->Branch("m_phi",      &m_phi,       "m_phi[nPar][3]/D");
    rtree2->Branch("m_st",       &m_st,        "m_st[nPar][3]/D");
    rtree2->Branch("m_mass",     &m_mass,      "m_mass[nPar][3]/D");
    //rtree2->Branch("m_vx",       &m_vx,        "m_vx[nPar][3]/D");
    //rtree2->Branch("m_vy",       &m_vy,        "m_vy[nPar][3]/D");
    rtree2->Branch("m_vz",       &m_vz,        "m_vz[nPar][3]/D");

    rtree2->Branch("nDau",       &m_nDau,      "nDau[nPar]/I");
    rtree2->Branch("dauId",      &dauId,       "dauId[nPar][3]/I");
    rtree2->Branch("d_pt",       &d_pt,        "d_pt[nPar][3]/D");
    rtree2->Branch("d_y",        &d_y,         "d_y[nPar][3]/D");
    rtree2->Branch("d_eta",      &d_eta,       "d_eta[nPar][3]/D");
    rtree2->Branch("d_phi",      &d_phi,       "d_phi[nPar][3]/D");
    rtree2->Branch("d_st",       &d_st,        "d_st[nPar][3]/D");
    rtree2->Branch("d_mass",     &d_mass,      "d_mass[nPar][3]/D");
    //rtree2->Branch("d_vx",       &p_vx,        "d_vx[nPar][3]/D");
    //rtree2->Branch("d_vy",       &p_vy,        "d_vy[nPar][3]/D");
    rtree2->Branch("d_vz",       &p_vz,        "d_vz[nPar][3]/D");

}

//
// member functions
//

// ------------ method called to for each event  ------------
    void
HydDileptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

#ifdef THIS_IS_AN_EVENT_EXAMPLE
    Handle<ExampleData> pIn;
    iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
    ESHandle<SetupData> pSetup;
    iSetup.get<SetupRecord>().get(pSetup);
#endif

    Handle<GenParticleCollection> genParticles;
    iEvent.getByLabel(genSource_, genParticles);

    m_evt = iEvent.id().event();
    m_run = iEvent.id().run();

    for (Int_t i=0; i<xNum; ++i) {
        m_pId[i]=0;
        p_pt[i]=0;
        p_y[i]=0;
        p_eta[i]=0;
        p_phi[i]=0;
        p_st[i]=0;
        p_mass[i]=0;
        //p_vx[i]=0;
        //p_vy[i]=0;
        p_vz[i]=0;
        m_nMom[i]=0;

        for(Int_t j=0; j<yNum; ++j) {
            momId[i][j]=0;
            m_pt[i][j]=0;
            m_y[i][j]=0;
            m_eta[i][j]=0;
            m_phi[i][j]=0;
            m_st[i][j]=0;
            m_mass[i][j]=0;
            //m_vx[i][j]=0;
            //m_vy[i][j]=0;
            m_vz[i][j]=0;

            dauId[i][j]=0;
            d_pt[i][j]=0;
            d_y[i][j]=0;
            d_eta[i][j]=0;
            d_phi[i][j]=0;
            d_st[i][j]=0;
            d_mass[i][j]=0;
            //d_vx[i][j]=0;
            //d_vy[i][j]=0;
            d_vz[i][j]=0;
        }
    }


    nPar = 0;

    std::vector<const Candidate*> sons;
    for(size_t i = 0; i < (genParticles->size()); ++ i) {

      const Candidate & p = (*genParticles)[i];
      if (! &p ) continue;
      int id = p.pdgId();

      //if(abs(id) == 443 || abs(id) == 553 || abs(id) == 23 || abs(id) ==13){
      if(abs(id) == 13){
            m_pId[nPar] = id;
            p_pt[nPar] = p.pt();
            p_y[nPar] = p.rapidity();
            p_eta[nPar]= p.eta();
            p_phi[nPar] = p.phi();
            p_mass[nPar] = p.mass();
            p_st[nPar] = p.status();
            //p_vx[nPar] = p.vx();
            //p_vy[nPar] = p.vy();
            p_vz[nPar] = p.vz();
            nmom = p.numberOfMothers();
            m_nMom[nPar] = nmom;

            // for(int j=0; j<m_nMom[nPar];j++){
            //		if (nmom>=2) continue;
            nmom = p.numberOfMothers();

            for(int j=0; j<nmom;++j){
                const Candidate * mom = p.mother(j);
                if (!mom) continue;
                momId[nPar][j] = mom->pdgId();
                //if (mom->pdgId() == 443 && id != 443)
                //sons.push_back(&p);
                m_pt[nPar][j] = mom->pt();
                m_y[nPar][j] = mom->rapidity();
                m_eta[nPar][j] = mom->eta();
                m_phi[nPar][j] = mom->phi();
                m_st[nPar][j] = mom->status();
                m_mass[nPar][j] = mom->mass();
                //m_vx[nPar][j] = mom->vx();
                //m_vy[nPar][j] = mom->vy();
                m_vz[nPar][j] = mom->vz();
                //++nMom;

            }

            ndau = p.numberOfDaughters();

            m_nDau[nPar] = ndau;
            for(int l=0; l<ndau;l++){
                const Candidate * dau = p.daughter(l);
                dauId[nPar][l] = dau->pdgId();
                d_pt[nPar][l] = dau->pt();
                d_y[nPar][l] = dau->rapidity();
                d_eta[nPar][l]=dau->eta();
                d_phi[nPar][l] = dau->phi();
                d_st[nPar][l] = dau->status();
                //d_vx[nPar][l] = dau->vx();
                //d_vy[nPar][l] = dau->vy();
                d_vz[nPar][l] = dau->vz();
                //++nDau;
            }

            ++nPar;
      }


    }
        //cout << "WP :: " << myMass(sons) << endl;


        rtree2->Fill();    

}

//define this as a plug-in
DEFINE_FWK_MODULE(HydDileptonAnalyzer);

// -*- C++ -*-
//
// Package:    Workspace/GsfTrajectoryAnalyzer
// Class:      GsfTrajectoryAnalyzer
// 
/**\class GsfTrajectoryAnalyzer GsfTrajectoryAnalyzer.cc Workspace/GsfTrajectoryAnalyzer/plugins/GsfTrajectoryAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Wolfgang Adam
//         Created:  Wed, 09 Aug 2017 13:54:40 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "TTree.h"
#include "TVector3.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

using std::cout;
using std::endl;

class GsfTrajectoryAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit GsfTrajectoryAnalyzer(const edm::ParameterSet&);
      ~GsfTrajectoryAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  void fillGlobal(const TrajectoryStateOnSurface& tsos, TVector3& position, TVector3& momentum);
  void fillParameters(const TrajectoryStateOnSurface& tsos, 
		      std::vector<float>& parameters, std::vector<float>& errors);
  void fillComponent(const std::vector<TrajectoryStateOnSurface>& components, size_t ic,
		     float& weight, std::vector<float>& parameters, std::vector<float>& errors);

      // ----------member data ---------------------------
private:
  edm::InputTag trajectoryTag_;

  TTree* tree_;
  int run_;
  int lumi_;
  long long evt_;
  int itraj_;
  int itm_[2];
  int ncs_[3];
  int ic_;
  float wgts_[3];
  TVector3 fwPredGPos_;
  TVector3 fwPredGMom_;
  std::vector<float> fwPredLPar_;
  std::vector<float> fwPredLErr_;
  TVector3 bwPredGPos_;
  TVector3 bwPredGMom_;
  std::vector<float> bwPredLPar_;
  std::vector<float> bwPredLErr_;
  TVector3 updGPos_;
  TVector3 updGMom_;
  std::vector<float> updLPar_;
  std::vector<float> updLErr_;
  TVector3 hitGPos_;
  std::vector<float> hitPar_;
  std::vector<float> hitErr_;
  // std::vector<float> projMat_;

  AlgebraicVector5 aVec5_;
  AlgebraicSymMatrix55 aSymMat55_;

  // typedef typename AlgebraicROOTObject<1,5>::Matrix Mat15;
  // typedef typename AlgebraicROOTObject<5,1>::Matrix Mat51;
  typedef typename AlgebraicROOTObject<1,1>::SymMatrix SMat11;
  typedef typename AlgebraicROOTObject<1>::Vector Vec1;
  ProjectMatrix<double,5,1>  projMat51_;
  Vec1 resVec1_, resMeas1_; 
  SMat11 matV11_, matVMeas11_;
  KfComponentsHolder holder1_;

  // typedef typename AlgebraicROOTObject<2,5>::Matrix Mat25;
  // typedef typename AlgebraicROOTObject<5,2>::Matrix Mat52;
  typedef typename AlgebraicROOTObject<2,2>::SymMatrix SMat22;
  typedef typename AlgebraicROOTObject<2>::Vector Vec2;
  ProjectMatrix<double,5,2>  projMat52_;
  Vec2 resVec2_, resMeas2_; 
  SMat22 matV22_, matVMeas22_;
  KfComponentsHolder holder2_;

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
GsfTrajectoryAnalyzer::GsfTrajectoryAnalyzer(const edm::ParameterSet& iConfig) :
  trajectoryTag_(iConfig.getParameter<edm::InputTag>("trajectories")),
  fwPredLPar_(5,0.),fwPredLErr_(5,0.),
  bwPredLPar_(5,0.),bwPredLErr_(5,0.),
  updLPar_(5,0.),updLErr_(5,0.),
  hitPar_(2,0.),hitErr_(2,0.),
  matV11_(ROOT::Math::SMatrixNoInit{}),matVMeas11_(ROOT::Math::SMatrixNoInit{}),
  matV22_(ROOT::Math::SMatrixNoInit{}),matVMeas22_(ROOT::Math::SMatrixNoInit{})
{
   // //now do what ever initialization is needed
  usesResource("TFileService");

  consumes< std::vector<Trajectory> >(trajectoryTag_);

  using ROOT::Math::SMatrixNoInit;
 
 
  holder1_.template setup<1>(&resVec1_, &matV11_, &projMat51_, &resMeas1_, &matVMeas11_, aVec5_, aSymMat55_);
  holder2_.template setup<2>(&resVec2_, &matV22_, &projMat52_, &resMeas2_, &matVMeas22_, aVec5_, aSymMat55_);

}


GsfTrajectoryAnalyzer::~GsfTrajectoryAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

void
GsfTrajectoryAnalyzer::fillGlobal(const TrajectoryStateOnSurface& tsos, TVector3& position, TVector3& momentum)
{
  GlobalPoint gPos = tsos.globalPosition();
  GlobalVector gMom = tsos.globalMomentum();
  position.SetXYZ(gPos.x(),gPos.y(),gPos.z());
  momentum.SetXYZ(gMom.x(),gMom.y(),gMom.z());
}

void
GsfTrajectoryAnalyzer::fillParameters(const TrajectoryStateOnSurface& tsos, 
				      std::vector<float>& parameters, std::vector<float>& errors)
{
       const LocalTrajectoryParameters& localUpdPars = tsos.localParameters();
       const LocalTrajectoryError& localUpdErrs = tsos.localError();
       for ( size_t i=0; i<5; ++i ) {
	 parameters[i] = localUpdPars.vector()[i];
	 errors[i] = sqrt(localUpdErrs.matrix()[i][i]);
       }
}

void
GsfTrajectoryAnalyzer::fillComponent(const std::vector<TrajectoryStateOnSurface>& components, size_t ic,
				     float& weight, std::vector<float>& parameters, std::vector<float>& errors)
{
  if ( ic<components.size() ) {
    const TrajectoryStateOnSurface& tsos = components[ic];
    weight = tsos.weight();
    fillParameters(tsos,parameters,errors);
  }
  else {
    weight = -1.;
    for ( size_t i=0; i<5; ++i ) {
      parameters[i] = 0.;
      errors[i] = 0.;
    }
  }
}				     

// ------------ method called for each event  ------------
void
GsfTrajectoryAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   run_ = iEvent.eventAuxiliary().run();
   lumi_ = iEvent.eventAuxiliary().luminosityBlock();
   evt_ = iEvent.eventAuxiliary().event();

   Handle< std::vector<Trajectory> > trajectoryHandle;
   iEvent.getByLabel(trajectoryTag_,trajectoryHandle);

   itraj_ = 0;
   // const std::vector<Trajectory>& trajectories = *trajectoryHandle.product();
   for ( std::vector<Trajectory>::const_iterator traj=trajectoryHandle->begin();
	 traj!=trajectoryHandle->end(); ++traj,++itraj_ ) {
     // cout << endl << "New trajectory" << endl;
     GlobalPoint gPos;
     GlobalVector gMom;
     
     size_t ntm = traj->measurements().size();
     itm_[0] = 0;
     itm_[1] = -ntm+1;
     for ( Trajectory::DataContainer::const_iterator itm=traj->measurements().begin();
	   itm!=traj->measurements().end(); ++itm,++itm_[0],++itm_[1] ) {
       ic_ = -1;
       for ( size_t i=0; i<3; ++i )  wgts_[i] = -1.;

       fillGlobal(itm->forwardPredictedState(),fwPredGPos_,fwPredGMom_);
       fillParameters(itm->forwardPredictedState(),fwPredLPar_,fwPredLErr_);
       ncs_[0] = itm->forwardPredictedState().components().size();

       fillGlobal(itm->backwardPredictedState(),bwPredGPos_,bwPredGMom_);
       fillParameters(itm->backwardPredictedState(),bwPredLPar_,bwPredLErr_);
       ncs_[1] = itm->backwardPredictedState().components().size();
       // cout << "  Position (pred bwd) " << gPos.perp() << " " << gPos.z() << endl;
       // cout << "  Momentum (pred bwd) " << gMom.perp() << " " << gMom.eta() << " " << gMom.phi() << endl;

       fillGlobal(itm->updatedState(),updGPos_,updGMom_);
       fillParameters(itm->updatedState(),updLPar_,updLErr_);
       ncs_[2] = itm->updatedState().components().size();

       const TrackingRecHit& recHit = *itm->recHit();
       for ( int i=0; i<2; ++i ) {
	 hitPar_[i] = 0.;
	 hitErr_[i] = 0.;
       }
       // projMat_.clear();
       if ( recHit.isValid() && (recHit.dimension()==1 || recHit.dimension()==2) && ic_<0 ) {
	 gPos = recHit.globalPosition();
	 hitGPos_.SetXYZ(gPos.x(),gPos.y(),gPos.z());
	 switch ( recHit.dimension() ) {
	 case 1:
	   holder1_.template setup<1>(&resVec1_, &matV11_, &projMat51_, &resMeas1_, &matVMeas11_, 
				      aVec5_, aSymMat55_);
	   recHit.getKfComponents(holder1_);
	   hitPar_[0] = holder1_.params<1>()[0];
	   hitErr_[0] = sqrt(holder1_.errors<1>()[0][0]);
	   // projMat_.resize(5,0.);
	   // for ( size_t i=0; i<5; ++i )  projMat_[i] = projMat51_.matrix()[0][i];
	   break;
	 case 2:
	   holder2_.template setup<2>(&resVec2_, &matV22_, &projMat52_, &resMeas2_, &matVMeas22_, 
				      aVec5_, aSymMat55_);
	   cout << "Global position " << bwPredGPos_.Pt() << " " << bwPredGPos_.z() << endl;
	   recHit.getKfComponents(holder2_);
	   for ( int i=0; i<2; ++i ) {
	     hitPar_[i] = holder2_.params<2>()[i];
	     hitErr_[i] = sqrt(holder2_.errors<2>()[i][i]);
	   }
	   double hitRho = holder2_.errors<2>()[0][1]/hitErr_[0]/hitErr_[1];
	   cout << "hitRho(1) = " << hitRho << endl;
	   if ( fabs(hitRho)>1.e-5 ) {
	     cout << "hitRho = " << hitRho << endl;
	     std::vector<float> par(hitPar_);
	     std::vector<float> err(hitErr_);
	     cout << "Old parameters / errors :" 
		  << " " << hitPar_[0] << " +- " << hitErr_[0] << "  ;  "
		  << " " << hitPar_[1] << " +- " << hitErr_[1] << endl;
	     for ( int i=0; i<2; ++i ) {
	       float pref;
	       if ( fwPredLErr_[5-i]<bwPredLErr_[5-i] ) {
		 pref = fwPredLPar_[5-i];
	       }
	       else {
		 pref = bwPredLPar_[5-i];
	       }
	       cout << "RefPar for " << 1-i << " is " << pref << endl;
	       hitPar_[i] += err[i]*hitRho/err[1-i]*(pref-par[1-i]);
	       hitErr_[i] *= sqrt(1-hitRho*hitRho);
	       cout << "New parameters / errors for " << i << " :" 
		    << " " << hitPar_[i] << " +- " << hitErr_[i] << endl;
	     }
	   }
	   // projMat_.resize(10,0.);
	   // for ( size_t j=0; j<2; ++j ) {
	   //   for ( size_t i=0; i<5; ++i )  projMat_[i+5*j] = projMat52_.matrix()[j][i];
	   // }
	   break;
	 }
	 // cout << "  Position (hit)      " << gPos.perp() << " " << gPos.z() << endl;
       }
       else {
	 hitGPos_.SetXYZ(0.,0.,0.);
	 // cout << "  Position (hit)      invalid hit" << endl;
       }
       // cout << "  Position (updated)  " << gPos.perp() << " " << gPos.z() << endl;
       // cout << "  Momentum (updated)  " << gMom.perp() << " " << gMom.eta() << " " << gMom.phi() << endl;
       // cout << "  ---" << endl;
       tree_->Fill();

       int ncmax = ncs_[0];
       if ( ncmax<ncs_[1] )  ncmax = ncs_[1];
       if ( ncmax<ncs_[2] )  ncmax = ncs_[2];

       std::vector<TrajectoryStateOnSurface> fwPredComps = itm->forwardPredictedState().components();
       std::vector<TrajectoryStateOnSurface> bwPredComps = itm->backwardPredictedState().components();
       std::vector<TrajectoryStateOnSurface> updComps = itm->updatedState().components();
       for ( size_t i=0; i<(size_t)ncmax; ++i ) {
	 ic_ = i;
	 fillComponent(fwPredComps,i,wgts_[0],fwPredLPar_,fwPredLErr_);
	 fillComponent(bwPredComps,i,wgts_[1],bwPredLPar_,bwPredLErr_);
	 fillComponent(updComps,i,wgts_[2],updLPar_,updLErr_);
	 tree_->Fill();
       }

     }
   }
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
void 
GsfTrajectoryAnalyzer::beginJob()
{
  edm::Service<TFileService> fs;

  tree_ = fs->make<TTree>("GsfTree","GsfTree");
  tree_->Branch("run",&run_,"run/I");
  tree_->Branch("lumi",&lumi_,"lumi/I");
  tree_->Branch("evt",&evt_,"evt/L");
  tree_->Branch("itraj",&itraj_,"itraj/I");
  tree_->Branch("itm",itm_,"itmf/I:itmr/I");
  tree_->Branch("ncs",ncs_,"ncFwPred/I:ncBwPred/I:ncUpd/I");
  tree_->Branch("ic",&ic_,"ic/I");
  tree_->Branch("wgts",wgts_,"wgtFwPred/F:wgtBwPred/F:wgtUpd/F");
  tree_->Branch("fwPredGPos",&fwPredGPos_);
  tree_->Branch("fwPredGMom",&fwPredGMom_);
  tree_->Branch("fwPredLPar",&fwPredLPar_);
  tree_->Branch("fwPredLErr",&fwPredLErr_);
  tree_->Branch("bwPredGPos",&bwPredGPos_);
  tree_->Branch("bwPredGMom",&bwPredGMom_);
  tree_->Branch("bwPredLPar",&bwPredLPar_);
  tree_->Branch("bwPredLErr",&bwPredLErr_);
  tree_->Branch("updGPos",&updGPos_);
  tree_->Branch("updGMom",&updGMom_);
  tree_->Branch("updLPar",&updLPar_);
  tree_->Branch("updLErr",&updLErr_);
  tree_->Branch("hitGPos",&hitGPos_);
  tree_->Branch("hitPar",&hitPar_);
  tree_->Branch("hitErr",&hitErr_);
  // tree_->Branch("projMat",&projMat_);

}

// ------------ method called once each job just after ending the event loop  ------------
void 
GsfTrajectoryAnalyzer::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
GsfTrajectoryAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(GsfTrajectoryAnalyzer);

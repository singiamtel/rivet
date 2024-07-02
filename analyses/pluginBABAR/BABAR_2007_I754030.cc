// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  B -> omega K
  class BABAR_2007_I754030 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2007_I754030);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable(310);
      BP.addStable(223);
      declare(BP, "BP");
      // histos
      book(_h,1,1,1);
    }

    void findChildren(const Particle & p, Particles & pim, Particles & pip,
		      Particles & pi0, unsigned int &ncount) {
      for( const Particle &child : p.children()) {
	if(child.pid()==PID::PIPLUS) {
	  pip.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::PIMINUS) {
	  pim.push_back(child);
	  ncount+=1;
	}
	else if(child.pid()==PID::PI0) {
	  pi0.push_back(child);
	  ncount+=1;
	}
	else if(child.children().empty()) {
	  ncount+=1;
	}
    	else
    	  findChildren(child,pim,pip,pi0,ncount);
      }
    }

    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 223,1}, { 321,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 223,1}, {-321,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
	if (BP.modeMatches(ix,2,mode) || BP.modeMatches(ix,2,modeCC)) {
	  const Particle & omega = BP.decayProducts()[ix].at(223)[0];
	  // children of the omega
	  unsigned int ncount=0;
	  Particles pip,pim,pi0;
	  findChildren(omega,pim,pip,pi0,ncount);
	  if( ncount!=3 || !(pim.size()==1 && pip.size()==1 && pi0.size()==1)) continue;
	  LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix].momentum().betaVec());
	  FourMomentum pomega = boost1.transform(omega.momentum());
	  FourMomentum pPip   = boost1.transform(pip[0].momentum());
	  FourMomentum pPim   = boost1.transform(pim[0].momentum());
	  // boost to omega frame
	  LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(pomega.betaVec());
	  pPip   = boost3.transform(pPip);
	  pPim   = boost3.transform(pPim);
	  Vector3 axisOmega = pPip.p3().cross(pPim.p3()).unit();
	  // helicity angle
	  _h->fill(abs(pomega.p3().unit().dot(axisOmega)));
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h,1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h;
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2007_I754030);

}

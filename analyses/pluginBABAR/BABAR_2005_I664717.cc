// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B+ -> omega rho+ 
  class BABAR_2005_I664717 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2005_I664717);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable( 213);
      BP.addStable(-213);
      BP.addStable( 223);
      declare(BP, "BP");
      // histos
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],1,1,1+ix);
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
      static const map<PdgId,unsigned int> & mode   = { { 223,1}, { 213,1}};
      static const map<PdgId,unsigned int> & modeCC = { { 223,1}, {-213,1}};
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
	  // rho decay
	  int sign = BP.decaying()[ix].pid()/BP.decaying()[ix].abspid();
	  const Particle & rho = BP.decayProducts()[ix].at(sign*213)[0];
	  if(rho.children().size()!=2) continue;
	  Particle pi;
	  if(rho.children()[0].pid()==sign*211 &&
	     rho.children()[1].pid()==111)
	    pi = rho.children()[0];
	  else if (rho.children()[1].pid()==sign*211 &&
		   rho.children()[0].pid()==111)
	    pi = rho.children()[1];
	  else continue;
	  // boosts
	  LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix].momentum().betaVec());
	  FourMomentum pomega = boost1.transform(omega .momentum());
	  FourMomentum prho   = boost1.transform(  rho .momentum());
	  FourMomentum pPip   = boost1.transform(pip[0].momentum());
	  FourMomentum pPim   = boost1.transform(pim[0].momentum());
	  FourMomentum pPi    = boost1.transform(pi    .momentum());
	  // boost to omega frame
	  LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pomega.betaVec());
	  pPip   = boost2.transform(pPip);
	  pPim   = boost2.transform(pPim);
	  Vector3 axisOmega = pPip.p3().cross(pPim.p3()).unit();
	  // helicity angle
	  _h[0]->fill(abs(pomega.p3().unit().dot(axisOmega)));
	  // boost to rho rest frame
	  LorentzTransform boost3 = LorentzTransform::mkFrameTransformFromBeta(prho.betaVec());
	  pPi = boost3.transform(pPi);
	  _h[1]->fill(prho.p3().unit().dot(pPi.p3().unit()));
      	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	normalize(_h[ix],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2005_I664717);

}

// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief  B0 -> K*0 eta' 
  class BABAR_2007_I722710 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2007_I722710);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511);
      declare(ufs, "UFS");
      DecayedParticles B0(ufs);
      B0.addStable(PID::ETAPRIME);
      declare(B0, "B0");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	book(_h[ix],ix+1,1,1);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1},{-211,1}, { 331,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1},{ 211,1}, { 331,1}};
      DecayedParticles B0 = apply<DecayedParticles>(event, "B0");
      // loop over particles
      for(unsigned int ix=0;ix<B0.decaying().size();++ix) {
      	int sign = 1;
      	if (B0.decaying()[ix].pid()>0 && B0.modeMatches(ix,3,mode)) {
      	  sign=1;
      	}
      	else if  (B0.decaying()[ix].pid()<0 && B0.modeMatches(ix,3,modeCC)) {
      	  sign=-1;
      	}
      	else
      	  continue;
	const Particle & Kp  = B0.decayProducts()[ix].at( sign*321)[0];
	const Particle & pim = B0.decayProducts()[ix].at(-sign*211)[0];
	double mKpi  = (Kp.momentum()+pim.momentum()).mass();
	if(mKpi<.755 || mKpi>1.305) continue;
	_h[0]->fill(mKpi);
	// boost to b rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(B0.decaying()[ix].momentum().betaVec());
	FourMomentum pK  = boost1.transform(Kp.momentum());
	FourMomentum ppi = boost1.transform(pim.momentum());
	FourMomentum pKS = pK+ppi;
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKS.betaVec());
	Vector3 axis = pKS.p3().unit();
	_h[1]->fill(axis.dot(boost2.transform(pK).p3().unit()));
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


  RIVET_DECLARE_PLUGIN(BABAR_2007_I722710);

}

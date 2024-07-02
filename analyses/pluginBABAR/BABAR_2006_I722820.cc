// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> eta K*
  class BABAR_2006_I722820 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2006_I722820);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // projection
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==511||
						Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BB(ufs);
      BB.addStable(PID::ETA);
      BB.addStable(PID::K0S);
      BB.addStable(PID::PI0);
      declare(BB, "BB");
      // histograms
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  book(_h[ix][iy],1,1+ix,1+iy);
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode1   = { { 321,1},{-211,1}, { 221,1}};
      static const map<PdgId,unsigned int> & mode1CC = { {-321,1},{ 211,1}, { 221,1}};
      static const map<PdgId,unsigned int> & mode2   = { { 321,1},{ 111,1}, { 221,1}};
      static const map<PdgId,unsigned int> & mode2CC = { {-321,1},{ 111,1}, { 221,1}};
      static const map<PdgId,unsigned int> & mode3   = { { 310,1},{ 211,1}, { 221,1}};
      static const map<PdgId,unsigned int> & mode3CC = { { 310,1},{-211,1}, { 221,1}};
      DecayedParticles BB = apply<DecayedParticles>(event, "BB");
      // loop over particles
      for(unsigned int ix=0;ix<BB.decaying().size();++ix) {
      	int imode = 1, iK=321, iPi=-211;
       	if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode1)) {
	  imode=0;
	  iK  = 321;
	  iPi =-211;
      	}
      	else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode1CC)) {
	  imode=0;
	  iK  =-321;
	  iPi = 211;
      	}
	else if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode2)) {
	  iK  = 321;
	  iPi = 111;
      	}
      	else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode2CC)) {
	  iK  =-321;
	  iPi = 111;
      	}
	else if (BB.decaying()[ix].pid()>0 && BB.modeMatches(ix,3,mode3)) {
	  iK  = 310;
	  iPi = 211;
      	}
      	else if  (BB.decaying()[ix].pid()<0 && BB.modeMatches(ix,3,mode3CC)) {
	  iK  = 310;
	  iPi =-211;
      	}
      	else
      	  continue;
	const Particle & KK = BB.decayProducts()[ix].at(iK )[0];
	const Particle & pi = BB.decayProducts()[ix].at(iPi)[0];
	double mKpi  = (KK.momentum()+pi.momentum()).mass();
	_h[0][imode]->fill(mKpi);
	// boost to b rest frame
	LorentzTransform boost1 = LorentzTransform::mkFrameTransformFromBeta(BB.decaying()[ix].momentum().betaVec());
	FourMomentum pK  = boost1.transform(KK.momentum());
	FourMomentum ppi = boost1.transform(pi.momentum());
	FourMomentum pKS = pK+ppi;
	LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKS.betaVec());
	Vector3 axis = pKS.p3().unit();
	_h[1][imode]->fill(axis.dot(boost2.transform(pK).p3().unit()));
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0;ix<2;++ix)
	for(unsigned int iy=0;iy<2;++iy)
	  normalize(_h[ix][iy],1.,false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[2][2];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2006_I722820);

}

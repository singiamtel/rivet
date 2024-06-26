// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> K pi pi
  class BABAR_2004_I626730 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2004_I626730);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==521);
      declare(ufs, "UFS");
      DecayedParticles BP(ufs);
      BP.addStable(PID::K0S);
      declare(BP, "BP");
      for(unsigned int ix=0; ix<3; ++ix) {
	book(_h_mass2[ix],2,1,1+ix);
	book(_h_angle[ix],2,2,1+ix);
	if(ix==2) continue;
	book(_h_mass [ix],1,1,1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      static const map<PdgId,unsigned int> & mode   = { { 321,1}, { 211,1}, {-211,1}};
      static const map<PdgId,unsigned int> & modeCC = { {-321,1}, { 211,1}, {-211,1}};
      DecayedParticles BP = apply<DecayedParticles>(event, "BP");
      // loop over particles
      for(unsigned int ix=0;ix<BP.decaying().size();++ix) {
      	int sign = BP.decaying()[ix].pid()>0 ? 1 : -1;
       	if((sign>0 and BP.modeMatches(ix,3,mode)) ||
	   (sign<0 and BP.modeMatches(ix,3,mode))) {
	  // boost to B rest frame
	  LorentzTransform boost =
	    LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix]. momentum().betaVec());
	  // momenta
	  FourMomentum pip  = boost.transform(BP.decayProducts()[ix].at( 211*sign)[0].momentum());
	  FourMomentum pim  = boost.transform(BP.decayProducts()[ix].at(-211*sign)[0].momentum());
	  FourMomentum Kp   = boost.transform(BP.decayProducts()[ix].at( 321*sign)[0].momentum());
	  // mass distributions
	  double mKpi = (Kp+pim).mass();
	  _h_mass[0]->fill(mKpi);
	  double mpipi = (pip+pim).mass();
	  _h_mass[1]->fill(mpipi);
	  // regions for the helicity analysis
	  unsigned int imode=0;
	  if(mKpi>0.816 && mKpi<0.976 && mpipi>1.5)   imode=0;
	  else if(mpipi>0.9 && mpipi<1.1)             imode=1;
	  else if(mKpi>0.976&& mKpi<1.8 && mpipi>1.5) imode=2;
	  else continue;
	  // pi+pi- resonance
	  if(imode==1) {
	    _h_mass2[imode]->fill(mpipi);
	    FourMomentum ppipi = pim+pip;
	    LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(ppipi.betaVec());
	    double cTheta = boost2.transform(pim).p3().unit().dot(Kp.p3().unit());
	    _h_angle[imode]->fill(cTheta);
	  }
	  // K pi- resonance
	  else {
	    _h_mass2[imode]->fill(mKpi);
	    FourMomentum pKpim = Kp+pim;
	    LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(pKpim.betaVec());
	    double cTheta = boost2.transform(Kp).p3().unit().dot(pip.p3().unit());
	    _h_angle[imode]->fill(cTheta);
	  }
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int ix=0; ix<3; ++ix) {
	normalize(_h_mass2[ix],1.,false);
	normalize(_h_angle[ix],1.,false);
	if(ix==2) continue;
	normalize(_h_mass[ix],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[2],_h_mass2[3],_h_angle[3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2004_I626730);

}

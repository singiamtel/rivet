// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief B -> K pi pi
  class BABAR_2008_I782405 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(BABAR_2008_I782405);


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
      // histos
      for(unsigned int iy=0;iy<3;++iy) {
	book(_h_mass[iy],1,1,1+iy);
	for(unsigned int ix=0;ix<2;++ix)
	  book(_h_mass2[ix][iy],2,1+iy,1+ix);
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
       	if((sign>0 and BP.modeMatches(ix,3,mode  )) ||
	   (sign<0 and BP.modeMatches(ix,3,modeCC))) {
	  LorentzTransform boost =
	    LorentzTransform::mkFrameTransformFromBeta(BP.decaying()[ix]. momentum().betaVec());
	  // momenta
	  FourMomentum pip  = boost.transform(BP.decayProducts()[ix].at( 211*sign)[0].momentum());
	  FourMomentum pim  = boost.transform(BP.decayProducts()[ix].at(-211*sign)[0].momentum());
	  FourMomentum Kp   = boost.transform(BP.decayProducts()[ix].at( 321*sign)[0].momentum());
	  // mass distributions
	  double mKpi = (Kp+pim).mass();
	  double mpipi = (pip+pim).mass();
	  // vetos
	  if(mKpi>1.756 && mKpi<1.931) continue;
	  if(mpipi>1.660 && mpipi<1.800) continue;
	  if(mpipi>3.019 && mpipi<3.179) continue;
	  if(mpipi>3.627 && mpipi<3.747) continue;
	  // mass distributions
	  if(mpipi>2.) _h_mass[0]->fill(mKpi);
	  if(mKpi>2.) {
	    _h_mass[1]->fill(mpipi);
	    _h_mass[2]->fill(mpipi);
	  }
	  FourMomentum ppipi = pim+pip;
	  LorentzTransform boost2 = LorentzTransform::mkFrameTransformFromBeta(ppipi.betaVec());
	  double cTheta = boost2.transform(pim).p3().unit().dot(Kp.p3().unit());
	  _h_mass2[(sign+1)/2][0]->fill(mpipi);
	  if(cTheta>0)
	    _h_mass2[(sign+1)/2][1]->fill(mpipi);
	  else
	    _h_mass2[(sign+1)/2][2]->fill(mpipi);
	}
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for(unsigned int iy=0;iy<3;++iy) {
	normalize(_h_mass[iy],1.,false);
	for(unsigned int ix=0;ix<2;++ix)
	  normalize(_h_mass2[ix][iy],1.,false);
      }
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h_mass[3],_h_mass2[2][3];
    /// @}


  };


  RIVET_DECLARE_PLUGIN(BABAR_2008_I782405);

}

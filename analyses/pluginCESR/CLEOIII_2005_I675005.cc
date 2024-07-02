// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/UnstableParticles.hh"
#include "Rivet/Projections/DecayedParticles.hh"

namespace Rivet {


  /// @brief Hadronic mass spectra in four hadron tau- decays with K+/-
  class CLEOIII_2005_I675005 : public Analysis {
  public:

    /// Constructor
    RIVET_DEFAULT_ANALYSIS_CTOR(CLEOIII_2005_I675005);


    /// @name Analysis methods
    /// @{

    /// Book histograms and initialise projections before the run
    void init() {
      // Initialise and register projections
      UnstableParticles ufs = UnstableParticles(Cuts::abspid==15);
      declare(ufs, "UFS");
      DecayedParticles TAU(ufs);
      TAU.addStable(310);
      TAU.addStable(223);
      TAU.addStable(111);
      declare(TAU, "TAU");
      // histos
      for (unsigned int ix=0; ix<3; ++ix) {
	      book(_h[ix], 1, 1, 1+ix);
      }
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      DecayedParticles TAU = apply<DecayedParticles>(event, "TAU");
      // loop over particles
      for (unsigned int ix=0;ix<TAU.decaying().size();++ix) {
      	int sign = 1, imode=0;
      	if (TAU.decaying()[ix].pid()>0 &&
      	  TAU.modeMatches(ix,5,mode1)) {
          sign=1;
          imode=0;
        }
      	else if (TAU.decaying()[ix].pid()<0 &&
      		TAU.modeMatches(ix,5,mode1CC)) {
          sign=-1;
          imode=0;
        }
        else if(TAU.decaying()[ix].pid()>0 &&
      	  TAU.modeMatches(ix,3,mode2)) {
          sign=1;
          imode=1;
        }
      	else if (TAU.decaying()[ix].pid()<0 &&
      		TAU.modeMatches(ix,3,mode2CC)) {
          sign=-1;
          imode=1;
        }
        else if(TAU.decaying()[ix].pid()>0 &&
                 TAU.modeMatches(ix,5,mode3)) {
          sign=1;
          imode=2;
        }
      	else if (TAU.decaying()[ix].pid()<0 &&
      		TAU.modeMatches(ix,5,mode3CC)) {
          sign=-1;
          imode=2;
        }
      	else {
      	  continue;
        }
      	const Particle& nu = TAU.decayProducts()[ix].at(16*sign)[0];
        const double mass = (TAU.decaying()[ix].mom()-nu.mom()).mass();
        _h[imode]->fill(mass);
      }
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      normalize(_h, 1.0, false);
    }

    /// @}


    /// @name Histograms
    /// @{
    Histo1DPtr _h[3];
    const map<PdgId,unsigned int> mode1   = { {-321,1}, { 211,1}, {-211,1}, {111,1}, { 16,1}};
    const map<PdgId,unsigned int> mode1CC = { { 321,1}, { 211,1}, {-211,1}, {111,1}, {-16,1}};
    const map<PdgId,unsigned int> mode2   = { {-321,1}, { 223,1}, { 16,1}};
    const map<PdgId,unsigned int> mode2CC = { { 321,1}, { 223,1}, {-16,1}};
    const map<PdgId,unsigned int> mode3   = { { 321,1}, {-321,1}, {-211,1}, {111,1}, { 16,1}};
    const map<PdgId,unsigned int> mode3CC = { { 321,1}, {-321,1}, { 211,1}, {111,1}, {-16,1}};
    /// @}


  };


  RIVET_DECLARE_PLUGIN(CLEOIII_2005_I675005);

}
